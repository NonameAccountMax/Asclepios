#!/usr/bin/env python3
"""Run the validation decks with asclepios and compare with the Abaqus reference results.

For every case of cases/manifest.json: runs bin/asclepios (output in validation/out/), extracts the same quantities
as abaqus/extract_odb.py (displacements and rotations of the monitored nodes, total reaction force on the FIX node
set, energy histories), and compares them with reference/<case>.json when it exists.

Errors are relative to the peak reference value of each quantity over the whole history:
  U, UR : max over frames and monitored nodes of |x - x_ref| / max |x_ref|
  RF    : same on the total reaction force of the FIX node set
  E     : max over time of |E - E_ref| for ALLIE, ALLKE, ALLWK (or the energies listed by the case), divided by the
          largest of these energies
Writes RESULTS.md and returns a non-zero status if a case fails (or, with --strict, has no reference).

usage: compare.py [case ...] [--no-run] [--jobs N] [--strict] [--update-readme]
"""
import base64
import bisect
import concurrent.futures
import csv
import json
import math
import os
import re
import struct
import subprocess
import sys
import time

import analytic

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, 'bin', 'asclepios')
CASES = os.path.join(HERE, 'cases')
REF = os.path.join(HERE, 'reference')
OUT = os.path.join(HERE, 'out')
ENERGIES = ('ALLIE', 'ALLKE', 'ALLWK')


def vtu_arrays(path):
    s = open(path).read()
    res = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="(\w+)"[^>]*format="binary"[^>]*>\s*([^<]*)<', s):
        typ, name, data = m.groups()
        b = base64.b64decode(data.strip())
        n = struct.unpack('<I', b[:4])[0]
        code = {'Float32': 'f', 'Float64': 'd', 'Int32': 'i'}.get(typ)
        if code:
            res[name] = struct.unpack('<%d%s' % (n // struct.calcsize(code), code), b[4:4 + n])
    return res


def mises(s):
    if len(s) == 3:  # plane stress 11, 22, 12
        s = [s[0], s[1], 0.0, s[2], 0.0, 0.0]
    a, b, c, d, e, f = s
    return (0.5 * ((a - b) ** 2 + (b - c) ** 2 + (c - a) ** 2) + 3 * (d * d + e * e + f * f)) ** 0.5


def element_values(a, case):
    """Same element quantities as the Abaqus extraction, from the cell data of an asclepios frame."""
    mone, etype = set(case.get('mone') or []), (case.get('mone_type') or '')
    out = {}
    if not mone or 'ElementLabel' not in a:
        return out
    for i, lab in enumerate(a['ElementLabel']):
        if lab not in mone:
            continue
        s6 = [a[k][i] for k in ('S11', 'S22', 'S33', 'S12', 'S13', 'S23')]
        e = {}
        if etype.startswith('B3'):
            e['sf1'] = a['SF'][i]
        elif etype.startswith(('S4', 'S3')):
            e['smem'] = mises([s6[0], s6[1], s6[3]])
            e['mises'] = a['S_Mises'][i]
        else:
            e['mises'] = a['S_Mises'][i]
            e['smem'] = mises(s6)
            e['savg'] = s6
        e['peeq'] = a['PEEQ'][i]
        out[str(lab)] = e
    return out


def asclepios_results(case):
    """Same layout as the Abaqus extraction: frames (t, U, UR, RF, RM, E), histories, number of increments."""
    job = os.path.join(OUT, case['name'])
    pvd = open(job + '.pvd').read()
    frames = []
    mon, fix = set(case['mon']), set(case['fix'])
    for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', pvd):
        a = vtu_arrays(os.path.join(OUT, f))
        labels = a['NodeLabel']
        fr = {'t': float(t), 'U': {}, 'UR': {}, 'RF': [0.0, 0.0, 0.0], 'RM': [0.0, 0.0, 0.0]}
        for i, lab in enumerate(labels):
            if lab in mon:
                fr['U'][str(lab)] = list(a['U'][3 * i:3 * i + 3])
                if 'UR' in a:
                    fr['UR'][str(lab)] = list(a['UR'][3 * i:3 * i + 3])
            if lab in fix:
                for k in range(3):
                    fr['RF'][k] += a['RF'][3 * i + k]
                    if 'RM' in a:
                        fr['RM'][k] += a['RM'][3 * i + k]
        fr['E'] = element_values(a, case)
        if 'NT11' in a:
            fr['NT'] = {str(lab): [a['NT11'][i]] for i, lab in enumerate(labels) if lab in mon}
        frames.append(fr)
    hist = {k: [] for k in ENERGIES + ('ALLPD', 'ALLVD', 'ETOTAL')}
    with open(job + '_history.csv') as f:
        for r in csv.DictReader(f):
            t = float(r['total_time'])
            for k in hist:
                if k in r:
                    hist[k].append([t, float(r[k])])
    mon_hist = None
    if os.path.exists(job + '_mon.csv'):
        mon_hist = {}
        with open(job + '_mon.csv') as f:
            for r in csv.DictReader(f):
                mon_hist.setdefault(r['node'], []).append(
                    [float(r['total_time'])] + [float(r[k]) for k in ('U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3')])
    rf_hist = None
    if os.path.exists(job + '_rf.csv'):
        with open(job + '_rf.csv') as f:
            rf_hist = [[float(r['total_time']), float(r['RF1']), float(r['RF2']), float(r['RF3'])] for r in csv.DictReader(f)]
    inc = None
    m = re.search(r'COMPLETED SUCCESSFULLY\s+\((\d+) increments', open(job + '.dat').read())
    if m:
        inc = int(m.group(1))
    return {'frames': frames, 'history': hist, 'increments': inc, 'rf_hist': rf_hist, 'mon_hist': mon_hist}


def interp(series, t):
    """Linear interpolation in a [[t, value], ...] list (values may be vectors)."""
    ts = [p[0] for p in series]
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return series[0][1]
    if i >= len(ts):
        return series[-1][1]
    t0, t1 = ts[i - 1], ts[i]
    w = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
    a, b = series[i - 1][1], series[i][1]
    if isinstance(a, list):
        return [x + w * (y - x) for x, y in zip(a, b)]
    return a + w * (b - a)


def norm(v):
    return sum(x * x for x in v) ** 0.5


def field_peak(ref, key):
    return max((norm(v) for fr in ref['frames'] for v in fr.get(key, {}).values()), default=0.0)


def force_scale(ref):
    """Characteristic force of the problem: largest energy / largest displacement of the reference."""
    e = max((abs(v) for k in ENERGIES for _, v in ref.get('history', {}).get(k, [])), default=0.0)
    u = field_peak(ref, 'U')
    return e / u if u > 0 else 0.0


def canonical_rotation(v):
    """Rotation vector written as Abaqus reports it: half angle brought into (-pi, pi] (angle in (-2 pi, 2 pi]), so
    that a total rotation of -7.53 rad about y reads +5.03 rad about y."""
    phi = norm(v)
    if phi < 1e-12:
        return list(v)
    half = 0.5 * phi
    half -= 2 * math.pi * round(half / (2 * math.pi))
    return [x / phi * 2 * half for x in v]


def history_error(mine, ref, key):
    """Explicit cases with the motion of the monitored nodes at every increment on both sides: peak-normalised error
    after a moving average over RF_WINDOW reference increments (same filter as the reactions). None when the
    histories are missing or the rotations are large (compared at the frames, in the Abaqus representation)."""
    rh, mh = ref.get('mon_hist'), mine.get('mon_hist')
    if not rh or not mh or key not in ('U', 'UR'):
        return None
    c0 = 1 if key == 'U' else 4
    out_err, peak = 0.0, 0.0
    for lab, rows in rh.items():
        if lab not in mh or len(rows) < 3 or len(rows[0]) < c0 + 3:
            return None
        rser = [[r[0]] + r[c0:c0 + 3] for r in rows]
        mser = [[r[0]] + r[c0:c0 + 3] for r in mh[lab]]
        if key == 'UR' and max(norm(r[1:]) for r in rser) > 1.0:
            return None
        if key == 'UR' and max(norm(r[1:]) for r in rser) < 1e-3:  # numerical zero at this node
            continue
        dts = sorted(b[0] - a[0] for a, b in zip(rser[:-1], rser[1:]) if b[0] > a[0])
        w = RF_WINDOW * dts[len(dts) // 2]
        times = window_times(rser, mser, w)
        fr, fm = window_average(rser, times, w), window_average(mser, times, w)
        peak = max(peak, max(norm(v) for v in fr))
        out_err = max(out_err, max(norm([a - b for a, b in zip(x, y)]) for x, y in zip(fm, fr)))
    return out_err / peak if peak > 0 else None


def field_error(mine, ref, key):
    """Peak-normalised error of a nodal vector field over all reference frames and monitored nodes. Rotations that
    stay below 1e-3 rad in the reference (numerical zero, e.g. a flat impact), for the case or at a monitored node,
    are not compared."""
    peak, err = field_peak(ref, key), 0.0
    if peak == 0.0 or (key == 'UR' and peak < 1e-3):
        return None
    h = history_error(mine, ref, key)
    if h is not None:
        return h
    small = set()  # monitored nodes whose rotation stays a numerical zero (below 1e-3 rad)
    if key == 'UR':
        labs = {lab for fr in ref['frames'] for lab in fr.get(key, {})}
        small = {lab for lab in labs if max(norm(fr[key][lab]) for fr in ref['frames'] if lab in fr.get(key, {})) < 1e-3}
    for fr in ref['frames']:
        for lab, v in fr.get(key, {}).items():
            if lab in small:
                continue
            series = [[m['t'], m[key][lab]] for m in mine['frames'] if lab in m.get(key, {})]
            if not series:
                return float('inf')
            got = interp(series, fr['t'])
            if key == 'UR':  # same representation of large rotations on both sides
                got, v = canonical_rotation(got), canonical_rotation(v)
            err = max(err, norm([a - b for a, b in zip(got, v)]))
    return err / peak


RF_WINDOW = 30  # moving average of explicit reaction histories, in reference increments


def window_average(series, times, w):
    """Moving average over [t - w/2, t + w/2] of a piecewise linear history [[t, v1, v2, v3], ...] (exact integral
    of the linear interpolant; the history is held constant outside its time range)."""
    ts = [r[0] for r in series]
    cum = [[0.0, 0.0, 0.0]]
    for a, b in zip(series[:-1], series[1:]):
        h = b[0] - a[0]
        cum.append([c + 0.5 * h * (a[k + 1] + b[k + 1]) for k, c in enumerate(cum[-1])])
    def integral(t):
        if t <= ts[0]:
            return [(t - ts[0]) * v for v in series[0][1:]]
        if t >= ts[-1]:
            return [c + (t - ts[-1]) * v for c, v in zip(cum[-1], series[-1][1:])]
        i = bisect.bisect_right(ts, t) - 1
        a, b = series[i], series[i + 1]
        h, x = b[0] - a[0], t - a[0]
        return [cum[i][k] + x * a[k + 1] + 0.5 * x * x * (b[k + 1] - a[k + 1]) / h for k in range(3)]
    out = []
    for t in times:
        lo, hi = integral(t - 0.5 * w), integral(t + 0.5 * w)
        out.append([(h - l) / w for l, h in zip(lo, hi)])
    return out


def window_times(ref, mine, w):
    """Reference times at which the moving averages are compared: none past half a window before the end of the
    histories (the last value would be held over half a window and the ringing at that instant would not be
    averaged); before the start, the initial value is the true one."""
    tend = min(ref[-1][0], mine[-1][0]) - 0.5 * w
    return [r[0] for r in ref if r[0] <= tend] or [ref[-1][0]]


def rf_error(mine, ref):
    """Total reaction force; not compared when it is a numerical zero (pure moment, self-equilibrated loads):
    below 1e-3 of the characteristic force of the problem. Explicit cases with reaction histories at every
    increment (both sides): compared after a moving average over RF_WINDOW reference increments, which removes the
    ringing near the highest frequencies of the mesh (sampled at the frames with an arbitrary phase); see
    window_times."""
    rh, mh = ref.get('rf_hist'), mine.get('rf_hist')
    if rh and mh and len(rh) > 2 and len(mh) > 2:
        dts = sorted(b[0] - a[0] for a, b in zip(rh[:-1], rh[1:]) if b[0] > a[0])
        w = RF_WINDOW * dts[len(dts) // 2]
        times = window_times(rh, mh, w)
        fr_ref, fr_mine = window_average(rh, times, w), window_average(mh, times, w)
        peak = max(norm(v) for v in fr_ref)
        if peak == 0.0 or peak < 1e-3 * force_scale(ref):
            return None
        return max(norm([a - b for a, b in zip(x, y)]) for x, y in zip(fr_mine, fr_ref)) / peak
    if not all('RF' in fr for fr in ref['frames']):
        return None
    peak = max((norm(fr['RF']) for fr in ref['frames']), default=0.0)
    if peak == 0.0 or peak < 1e-3 * force_scale(ref):
        return None
    series = [[m['t'], m['RF']] for m in mine['frames']]
    return max(norm([a - b for a, b in zip(interp(series, fr['t']), fr['RF'])]) for fr in ref['frames']) / peak


def rm_error(mine, ref):
    """Total reaction moment of the FIX nodes (shells, beams); a numerical zero is not compared."""
    if not all('RM' in fr for fr in ref['frames']) or field_peak(ref, 'UR') < 1e-4:
        return None  # no significant rotation (membrane problem): the moments are numerical noise
    peak = max((norm(fr['RM']) for fr in ref['frames']), default=0.0)
    emax = max((abs(v) for k in ENERGIES for _, v in ref.get('history', {}).get(k, [])), default=0.0)
    if peak == 0.0 or peak < 1e-3 * emax:
        return None
    series = [[m['t'], m['RM']] for m in mine['frames']]
    return max(norm([a - b for a, b in zip(interp(series, fr['t']), fr['RM'])]) for fr in ref['frames']) / peak


def element_error(mine, ref, shell=False):
    """Element quantities of the monitored elements (largest von Mises stress, membrane or mean-stress von Mises,
    PEEQ, axial force of beams): largest error relative to the peak reference value of each quantity."""
    worst = None
    peaks = {q: max((abs(e[q]) for fr in ref['frames'] for e in fr.get('E', {}).values() if e.get(q) is not None),
                    default=0.0) for q in ('mises', 'smem', 'peeq', 'sf1')}
    smax = max(peaks['mises'], peaks['smem'])
    for q in ('mises', 'smem', 'peeq', 'sf1'):
        peak = peaks[q]
        # numerical zeros are not compared: no plastic strain, membrane stress of a plate in pure bending
        # the membrane stress of a bending-dominated shell (below 5 % of the peak stress) is not compared
        if peak < (1e-4 if q == 'peeq' else 1e-9) or (q in ('mises', 'smem') and peak < 5e-2 * smax):
            continue
        err = 0.0
        for fr in ref['frames']:
            for lab, e in fr.get('E', {}).items():
                if e.get(q) is None:
                    continue
                # the reference membrane stress of a shell is the mean of its bottom and top section points: it is
                # the membrane stress only while the stress varies linearly through the thickness (no plasticity)
                if shell and q == 'smem' and (e.get('peeq') or 0.0) > 0.0:
                    continue
                series = [[m['t'], m['E'][lab][q]] for m in mine['frames'] if lab in m.get('E', {}) and
                          m['E'][lab].get(q) is not None]
                if not series:
                    break
                err = max(err, abs(interp(series, fr['t']) - e[q]))
        worst = max(worst or 0.0, err / peak)
    return worst


def energy_error(mine, ref, energies=None):
    keys = [k for k in (energies or ENERGIES) if ref['history'].get(k) and mine['history'].get(k)]
    scale = max((abs(v) for k in keys for _, v in ref['history'][k]), default=0.0)
    if scale == 0.0:
        return None
    err = 0.0
    for k in keys:  # at the times of the coarser history, the finer one interpolated (no extrapolation)
        rs, ms = ref['history'][k], mine['history'][k]
        inside = [(t, v) for t, v in ms if rs[0][0] - 1e-12 <= t <= rs[-1][0] + 1e-12]
        if len(inside) >= len(rs):
            pts = [(abs(interp(ms, t) - v)) for t, v in rs if ms[0][0] - 1e-12 <= t <= ms[-1][0] + 1e-12]
        else:
            pts = [abs(v - interp(rs, t)) for t, v in inside]
        err = max([err] + pts)
    return err / scale


def run_case(case):
    t0 = time.time()
    # energy histories at every increment: compared with the reference at its own output times (an Abaqus history
    # point may fall one increment after the requested time, in the middle of an impact)
    env = dict(os.environ, ABQ_HISTORY_EVERY_INCREMENT='1')
    r = subprocess.run([BIN, os.path.join(CASES, case['name'] + '.inp'), '-q', '-j', '2'], cwd=OUT,
                       capture_output=True, text=True, env=env)
    ok = r.returncode == 0 and os.path.exists(os.path.join(OUT, case['name'] + '.pvd'))
    with open(os.path.join(OUT, case['name'] + '.status'), 'w') as f:
        f.write('ok\n' if ok else 'failed\n')
    return ok, time.time() - t0, (r.stderr.strip().splitlines() or [''])[-1]


def evaluate(case, ran):
    res = {'name': case['name'], 'status': 'RUN FAILED', 'U': None, 'UR': None, 'RF': None, 'RM': None, 'S': None,
           'E': None, 'inc': None}
    if not case.get('supported', True):
        res['status'] = 'NOT SUPPORTED'
        return res
    if not ran:
        return res
    try:
        mine = asclepios_results(case)
    except (OSError, KeyError, ValueError) as e:
        res['message'] = str(e)
        return res
    if case.get('reference') == 'literature':
        return literature_check(case, mine, res)
    if case.get('reference') == 'analytic':
        ref = analytic.reference(case)
    else:
        ref_path = os.path.join(REF, case['name'] + '.json')
        if not os.path.exists(ref_path):
            res['status'] = 'NO REFERENCE'
            return res
        ref = json.load(open(ref_path))
    res['U'] = field_error(mine, ref, 'U')
    res['UR'] = field_error(mine, ref, 'UR') if case['rot'] else None
    res['RF'] = rf_error(mine, ref) if case['fix'] else None
    res['RM'] = rm_error(mine, ref) if case['fix'] and case['rot'] else None
    nt = field_error(mine, ref, 'NT') if any('NT' in fr for fr in ref['frames']) else None
    if nt is not None:
        res['U'] = max(res['U'] or 0.0, nt)
    res['S'] = element_error(mine, ref, shell=(case.get('mone_type') or '').startswith(('S3', 'S4')))
    res['E'] = energy_error(mine, ref, case.get('energies'))
    if mine.get('increments') and ref.get('increments'):
        res['inc'] = mine['increments'] / float(ref['increments'])
    tol = dict(case['tol'])
    tol.setdefault('S', 2 * tol['U'])
    tol['UR'], tol['RM'] = tol['U'], tol['RF']
    skip = set(case.get('known_quantities') or [])
    checks = [(res[q], tol[q]) for q in ('U', 'UR', 'RF', 'RM', 'S', 'E') if q not in skip]
    res['status'] = 'PASS' if all(v is None or v <= t for v, t in checks) else 'FAIL'
    return res


def last_frame_arrays(case):
    pvd = open(os.path.join(OUT, case['name'] + '.pvd')).read()
    files = re.findall(r'file="([^"]+)"', pvd)
    return vtu_arrays(os.path.join(OUT, files[-1]))


def extrapolate(dist, vals, deg):
    """Value at distance 0 of the least-squares polynomial of degree deg through (dist, vals) (interpolating when
    deg = len(dist) - 1)."""
    n = deg + 1
    a = [[sum(d ** (i + j) for d in dist) for j in range(n)] + [sum(v * d ** i for d, v in zip(dist, vals))]
         for i in range(n)]
    for k in range(n):  # Gauss-Jordan elimination with partial pivoting
        p = max(range(k, n), key=lambda r: abs(a[r][k]))
        a[k], a[p] = a[p], a[k]
        for r in range(n):
            if r != k:
                f = a[r][k] / a[k][k]
                a[r] = [x - f * y for x, y in zip(a[r], a[k])]
    return a[0][n] / a[0][0]


def literature_check(case, mine, res):
    """Published benchmarks: every check of the manifest entry at the last frame, error relative to the published
    value; U collects the displacement checks, S the stress checks."""
    lit = case['literature']
    last = mine['frames'][-1]
    arrays = None
    res['checks'] = []
    ok = True
    for c in lit['checks']:
        if c['kind'] == 'U':
            v = c['sign'] * last['U'][str(c['node'])][c['dof'] - 1]
        else:
            if arrays is None:
                arrays = last_frame_arrays(case)
            idx = {lab: i for i, lab in enumerate(arrays['ElementLabel'])}
            vals = [sum(arrays[c['comp']][idx[e]] for e in g) / len(g) for g in c['groups']]
            v = extrapolate(c['dist'], vals, c.get('deg', len(vals) - 1))
        err = abs(v - c['value']) / abs(c['value'])
        res['checks'].append((c['label'], v, c['value'], err, c['tol']))
        key = 'U' if c['kind'] == 'U' else 'S'
        res[key] = max(res[key] or 0.0, err)
        ok = ok and err <= c['tol']
    res['status'] = 'PASS' if ok else 'FAIL'
    return res


def with_known(case, res):
    """Known limitations (manifest 'known') are reported as KNOWN FAIL and do not fail the run; application
    validation cases (manifest 'deferred') are reported as DEFERRED and do not count either."""
    if case.get('deferred') and res['status'] not in ('NOT SUPPORTED',):
        res['status'] = 'DEFERRED (%s)' % res['status']
        return res
    if case.get('known') and res['status'] in ('FAIL', 'RUN FAILED'):
        res['status'] = 'KNOWN FAIL'
    elif case.get('known') and res['status'] == 'PASS':
        res['status'] = 'PASS (known limitation fixed: remove it from KNOWN in gen_cases.py)'
    return res


def pct(v):
    return '' if v is None else ('inf' if v == float('inf') else '%.2f %%' % (100 * v))


def write_results(cases, results, elapsed):
    by = {r['name']: r for r in results}
    counts = {}
    for r in results:
        counts[r['status']] = counts.get(r['status'], 0) + 1
    lines = ['# Validation results', '',
             'Generated by `validation/compare.py` (`make validate`). Reference: Abaqus results of the same decks '
             '(`validation/reference/`). Errors are relative to the peak reference value over the whole history '
             '(see `validation/README.md`).', '',
             'Summary: ' + ', '.join('%d %s' % (n, s) for s, n in sorted(counts.items())) +
             ' (%d cases, %.0f s).' % (len(results), elapsed), '',
             ]
    titles = {'element': 'Elements and deformation modes', 'patch': 'Patch tests, distorted meshes, robustness',
              'benchmark': 'Standard benchmarks', 'literature': 'Published benchmarks (refined meshes)', 'input': 'Input file features', 'contact': 'Contact',
              'constraint': 'Constraints and connectors', 'load': 'Loads, boundary conditions, steps',
              'material': 'Materials', 'explicit': 'Options specific to the explicit procedure',
              'standard': 'Procedures specific to the implicit solver', 'procedure': 'Other procedures',
              'application': 'Application models', 'parallel': 'Parallel runs and performance'}
    order = ['element', 'patch', 'literature', 'benchmark', 'input', 'contact', 'constraint', 'load', 'material', 'explicit',
             'standard', 'procedure', 'application', 'parallel']
    def section(c):
        return 'element' if c.get('category', 'element') == 'element' and 'contact' not in c['mode'] else (
            'contact' if 'contact' in c['mode'] else c.get('group', 'element'))
    cur = None
    for c in sorted(cases, key=lambda c: (order.index(section(c)) if section(c) in order else 99)):
        r = by.get(c['name'])
        if not r:
            continue
        if section(c) != cur:
            cur = section(c)
            lines += ['', '## ' + titles.get(cur, cur), '',
                      '| Case | Element | Mode | Procedure | Reference | U | UR | RF | RM | Stress | Energy | '
                      'Increments asclepios / ref | Tolerance (U / RF / S / E) | Status |',
                      '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
        t = c['tol']
        lines.append('| `%s` | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %g / %g / %g / %g %% | %s |' % (
            c['name'], c['element'], c['mode'], c['procedure'], {'analytic': 'closed form', 'literature': 'published'}
            .get(c.get('reference'), 'Abaqus'), pct(r['U']), pct(r['UR']), pct(r['RF']), pct(r.get('RM')), pct(r.get('S')), pct(r['E']),
            '%.2f' % r['inc'] if r.get('inc') else '', 100 * t['U'], 100 * t['RF'], 100 * t.get('S', 2 * t['U']),
            100 * t['E'], r['status']))
    lits = [c for c in cases if c.get('reference') == 'literature' and by.get(c['name'], {}).get('checks')]
    if lits:
        lines += ['', '## Published benchmarks: values', '', '| Case | Quantity | asclepios | Published | Error | '
                  'Tolerance | Source |', '|---|---|---|---|---|---|---|']
        for c in lits:
            for label, v, ref, err, tol in by[c['name']]['checks']:
                lines.append('| `%s` | %s | %.5g | %.5g | %s | %g %% | %s |' % (c['name'], label, v, ref, pct(err),
                                                                           100 * tol, c['literature']['source']))
    known = [c for c in cases if c.get('known') and c['name'] in by]
    partial = [c for c in cases if c.get('known_quantities') and c['name'] in by]
    if known or partial:
        lines += ['', 'Known limitations:', '']
        lines += ['- `%s`: %s.' % (c['name'], c['known']) for c in known]
        lines += ['- `%s` (%s not checked): %s.' % (c['name'], ', '.join(c['known_quantities']), c['known_quantities_note'])
                  for c in partial]
    with open(os.path.join(HERE, 'RESULTS.md'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return counts


def update_readme(counts, n):
    path = os.path.join(ROOT, 'README.md')
    s = open(path).read()
    merged = {}
    for k, v in counts.items():
        key = 'DEFERRED' if k.startswith('DEFERRED') else k
        merged[key] = merged.get(key, 0) + v
    counts = merged
    words = {'PASS': 'pass', 'FAIL': 'fail', 'RUN FAILED': 'did not run to completion',
             'DEFERRED': 'application validation cases (deferred)',
             'NO REFERENCE': 'awaiting Abaqus reference results', 'KNOWN FAIL': 'known limitations',
             'NOT SUPPORTED': 'not supported by asclepios yet (run in Abaqus only)'}
    summary = ', '.join('%d %s' % (k, words.get(v, v.lower())) for v, k in sorted(counts.items(), key=lambda kv: kv[0]))
    block = ('<!-- validation:start -->\nLast run of the %d cases: %s. Details: '
             '[validation/RESULTS.md](validation/RESULTS.md).\n<!-- validation:end -->' % (n, summary))
    s2 = re.sub(r'<!-- validation:start -->.*?<!-- validation:end -->', block, s, flags=re.S)
    if s2 != s:
        open(path, 'w').write(s2)


def main():
    args = sys.argv[1:]
    jobs = int(args[args.index('--jobs') + 1]) if '--jobs' in args else max(1, (os.cpu_count() or 4) // 2)
    names = [a for i, a in enumerate(args) if not a.startswith('--') and (i == 0 or args[i - 1] != '--jobs')]
    cases = json.load(open(os.path.join(CASES, 'manifest.json')))
    if names:
        cases = [c for c in cases if c['name'] in names]
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    ran = {}
    todo = [c for c in cases if c.get('supported', True)]
    if '--no-run' in args:
        def status(c):
            p = os.path.join(OUT, c['name'] + '.status')
            return os.path.exists(p) and open(p).read().strip() == 'ok'
        ran = {c['name']: status(c) for c in cases}
    else:
        with concurrent.futures.ThreadPoolExecutor(jobs) as ex:
            for c, (ok, dt, msg) in zip(todo, ex.map(run_case, todo)):
                ran[c['name']] = ok
                print('%-28s %s %6.1f s %s' % (c['name'], 'ran   ' if ok else 'FAILED', dt, '' if ok else msg))
    results = [with_known(c, evaluate(c, ran.get(c['name'], False))) for c in cases]
    for r in results:
        print('%-30s U %-9s UR %-9s RF %-9s RM %-9s S %-9s E %-9s %s' % (
            r['name'], pct(r['U']), pct(r['UR']), pct(r['RF']), pct(r.get('RM')), pct(r.get('S')), pct(r['E']),
            r['status']))
    counts = write_results(cases, results, time.time() - t0)
    if '--update-readme' in args and not names:
        update_readme(counts, len(results))
    print('summary: ' + ', '.join('%d %s' % (n, s) for s, n in sorted(counts.items())))
    bad = counts.get('FAIL', 0) + counts.get('RUN FAILED', 0) + (counts.get('NO REFERENCE', 0) if '--strict' in args else 0)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
