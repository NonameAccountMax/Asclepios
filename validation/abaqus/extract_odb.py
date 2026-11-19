# -*- coding: utf-8 -*-
"""Extract the validation quantities from an Abaqus output database (runs with `abaqus python`, Python 2.7 or 3).

usage: abaqus python extract_odb.py job.odb reference.json manifest.json case_name

Writes, for every field output frame of every step (total time t):
  U, UR : vectors of the monitored nodes (node set MON of the deck),
  RF, RM: total reaction force and moment of the node set FIX,
  E     : per monitored element (element set MONE): largest von Mises stress over its output points, von Mises
          stress of the mean in-plane stress (shells: membrane), mean stress tensor (solids), largest PEEQ,
          section force SF1 (axial force of beams);
the whole-model histories (ALLIE, ALLKE, ALLWK, ALLPD, ALLVD, ALLAE, ETOTAL, DT, DMASS) against total time, and the
number of increments and iterations read from the status file (<job>.sta).
"""
from __future__ import print_function

import json
import math
import os
import sys

from odbAccess import openOdb

HISTORY = ('ALLIE', 'ALLKE', 'ALLWK', 'ALLPD', 'ALLVD', 'ALLAE', 'ETOTAL', 'DT', 'DMASS')


def vec(v):
    return [float(x) for x in v]


def mises(s):
    """von Mises stress of a tensor given as (11, 22, 33, 12, 13, 23), or (11, 22, 12) for plane stress."""
    if len(s) == 3:
        s = [s[0], s[1], 0.0, s[2], 0.0, 0.0]
    elif len(s) == 4:
        s = [s[0], s[1], s[2], s[3], 0.0, 0.0]
    a, b, c, d, e, f = (list(s) + [0.0] * 6)[:6]
    return math.sqrt(0.5 * ((a - b) ** 2 + (b - c) ** 2 + (c - a) ** 2) + 3 * (d * d + e * e + f * f))


def element_values(fo, mone):
    out = {}
    if 'S' in fo.keys():
        acc = {}
        for v in fo['S'].values:
            if v.elementLabel in mone:
                acc.setdefault(v.elementLabel, []).append(vec(v.data))
        for lab, vals in acc.items():
            n = len(vals[0])
            mean = [sum(x[i] for x in vals) / len(vals) for i in range(n)]
            if n == 3:  # plane stress shell output: 11, 22, 12
                smem = mises(mean)
            elif n == 4:  # 11, 22, 33, 12
                smem = mises(mean)
            else:
                smem = mises(mean)
            out.setdefault(str(lab), {}).update({'mises': max(mises(x) for x in vals), 'smem': smem,
                                                 'savg': mean if n == 6 else None})
    if 'PEEQ' in fo.keys():
        acc = {}
        for v in fo['PEEQ'].values:
            if v.elementLabel in mone:
                acc[v.elementLabel] = max(acc.get(v.elementLabel, 0.0), float(v.data))
        for lab, p in acc.items():
            out.setdefault(str(lab), {})['peeq'] = p
    if 'SF' in fo.keys():
        acc = {}
        for v in fo['SF'].values:
            if v.elementLabel in mone:
                acc.setdefault(v.elementLabel, []).append(float(v.data[0]))
        for lab, vals in acc.items():
            out.setdefault(str(lab), {})['sf1'] = sum(vals) / len(vals)
    return out


def status_counts(sta):
    """(increments, iterations) from an Abaqus status file; iterations only for Abaqus/Standard."""
    if not os.path.exists(sta):
        return None, None
    lines = open(sta).read().splitlines()
    standard = any('STEP' in l and 'INC' in l and 'ATT' in l for l in lines)
    inc, it, last, total = 0, 0, 0, 0
    for l in lines:
        p = l.split()
        if not p or not p[0].isdigit():
            continue
        if standard and len(p) >= 9:
            if not p[2].endswith('U'):  # converged attempts only (a 'U' marks a cutback)
                inc += 1
                try:
                    it += int(p[5])
                except ValueError:
                    pass
        elif not standard and len(p) >= 6 and 'E' in p[1]:
            n = int(p[0])
            if n < last:
                total += last
            last = n
    if not standard:
        inc, it = total + last, None
    return inc, it


def main():
    odb_path, out_path, manifest, name = sys.argv[1:5]
    case = [c for c in json.load(open(manifest)) if c['name'] == name][0]
    mon, fix, mone = set(case['mon']), set(case['fix']), set(case.get('mone') or [])
    odb = openOdb(odb_path, readOnly=True)
    steps = sorted(odb.steps.values(), key=lambda s: s.totalTime)
    frames = []
    history = dict((k, []) for k in HISTORY)
    for step in steps:
        for frame in step.frames:
            fo = frame.fieldOutputs
            fr = {'t': step.totalTime + frame.frameValue, 'U': {}, 'UR': {}, 'RF': [0.0, 0.0, 0.0],
                  'RM': [0.0, 0.0, 0.0]}
            for key in ('U', 'UR'):
                if key in fo.keys():
                    for v in fo[key].values:
                        if v.nodeLabel in mon:
                            fr[key][str(v.nodeLabel)] = vec(v.data)
            if 'NT11' in fo.keys():
                fr['NT'] = dict((str(v.nodeLabel), [float(v.data)]) for v in fo['NT11'].values if v.nodeLabel in mon)
            for key in ('RF', 'RM'):
                if key in fo.keys():
                    for v in fo[key].values:
                        if v.nodeLabel in fix:
                            for k in range(3):
                                fr[key][k] += float(v.data[k])
            if 'RM' not in fo.keys():
                del fr['RM']
            if mone:
                fr['E'] = element_values(fo, mone)
            frames.append(fr)
        for rname in step.historyRegions.keys():
            if not rname.startswith('Assembly'):
                continue
            region = step.historyRegions[rname]
            for key in HISTORY:
                if key in region.historyOutputs.keys():
                    data = region.historyOutputs[key].data
                    if not data:
                        continue
                    # history times are step times in the database: shift them to total time
                    shift = step.totalTime if data[0][0] < step.totalTime - 1e-12 or step.totalTime == 0.0 else 0.0
                    history[key] += [[float(t) + shift, float(v)] for t, v in data]
    info = {'version': str(odb.jobData.version), 'precision': str(odb.jobData.precision)}
    odb.close()
    inc, it = status_counts(os.path.splitext(odb_path)[0] + '.sta')
    with open(out_path, 'w') as f:
        json.dump({'case': name, 'abaqus': info, 'increments': inc, 'iterations': it, 'frames': frames,
                   'history': history}, f)
    print('%s: %d frames written to %s' % (name, len(frames), out_path))


if __name__ == '__main__':
    main()
