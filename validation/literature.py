#!/usr/bin/env python3
"""Published benchmarks: mesh-convergence series and figures (make literature).

For every benchmark of gen_literature.py, runs asclepios on a series of meshes and compares the result with the
published value; writes one figure per benchmark (deformed shape of the validation mesh and convergence of the
normalised result) to validation/literature/ and the page validation/LITERATURE.md.

usage: literature.py [--no-run]      (--no-run: reuse the results already in validation/out/literature)
"""
import json
import math
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection   # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import deck                                           # noqa: E402
import compare                                        # noqa: E402
import gen_literature as gl                           # noqa: E402

ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, 'bin', 'asclepios')
WORK = os.path.join(HERE, 'out', 'literature')
FIG = os.path.join(HERE, 'literature')
CMAP = 'viridis'

# benchmark: (title, [(series label, generator, [argument tuples]), ...], index of the validation mesh in each series)
BENCHMARKS = [
    ('scordelis_lo', 'Scordelis-Lo roof', [('S4R', gl.scordelis_lo, [('S4R', n) for n in (4, 8, 16, 32)]),
                                           ('S3R', gl.scordelis_lo, [('S3R', n) for n in (8, 16, 32, 48)])]),
    ('pinched_cylinder', 'Pinched cylinder with rigid diaphragms',
     [('S4R', gl.pinched_cylinder, [(n,) for n in (8, 16, 32, 48)])]),
    ('pinched_hemisphere', 'Pinched hemisphere with an 18 degree hole',
     [('S4R', gl.pinched_hemisphere, [(n,) for n in (4, 8, 16, 32)])]),
    ('twisted_beam', 'Twisted beam', [('S4R, load along the clamped width', gl.twisted_beam,
                                       [('y', 6 * m, m) for m in (2, 4, 8, 16)]),
                                      ('S4R, load normal to the clamped end', gl.twisted_beam,
                                       [('z', 6 * m, m) for m in (2, 4, 8, 16)])]),
    ('cook_membrane', "Cook's membrane", [('S4R', gl.cook_membrane, [('S4R', n) for n in (4, 8, 16, 32, 48)]),
                                          ('S3R', gl.cook_membrane, [('S3R', n) for n in (4, 8, 16, 32, 48)])]),
    ('nafems_le1', 'NAFEMS LE1 elliptic membrane', [('S4R', gl.nafems_le1, [('S4R', 2 * m, m) for m in (6, 12, 24, 32)]),
                                                    ('S3R', gl.nafems_le1, [('S3R', 2 * m, m) for m in (6, 12, 24, 32)])]),
    ('clamped_plate', 'Clamped square plate under pressure', [('S4R', gl.clamped_plate, [('S4R', n) for n in (4, 8, 16, 32)]),
                                                             ('S3R', gl.clamped_plate, [('S3R', n) for n in (4, 8, 16, 32)])]),
    ('lame_cylinder', 'Thick cylinder under internal pressure (Lame)',
     [('C3D8R', gl.lame_cylinder, [('C3D8R', m, 2 * m) for m in (2, 4, 8, 16)]),
      ('C3D4', gl.lame_cylinder, [('C3D4', m, 2 * m) for m in (2, 4, 8, 16)])]),
    ('elastica', 'Elastica: cantilever under a large end load', [('B31', gl.elastica, [(5.0, n) for n in (5, 10, 20, 40)])]),
]


def tag(args):
    return '_'.join(str(a).replace('.', 'p') for a in args)


def build_and_run(gen, args, run):
    """Writes the deck of one mesh under a unique name, runs it, returns (case entry, results)."""
    deck.OUT = WORK
    deck.CASES.clear()
    gen(*args)
    case = dict(deck.CASES[-1])
    name = '%s__%s' % (case['name'], tag(args))
    os.replace(os.path.join(WORK, case['name'] + '.inp'), os.path.join(WORK, name + '.inp'))
    case['name'] = name
    if run:
        r = subprocess.run([BIN, name + '.inp', '-q', '-j', '4'], cwd=WORK, capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError('%s: %s' % (name, r.stderr.strip()[-300:]))
    compare.OUT = WORK
    res = compare.literature_check(case, compare.asclepios_results(case), {'U': None, 'S': None})
    return case, res


def elements_per_side(args):
    nums = [a for a in args if isinstance(a, int)]
    return nums[0] if nums else 0


def read_mesh(path):
    """Nodes and element connectivities (by type) of a deck."""
    nodes, elems, mode, etype = {}, [], None, None
    for line in open(path):
        if line.startswith('*'):
            low = line.lower()
            mode = None
            if low.startswith('*node') and not low.startswith('*node output'):
                mode = 'node'
            elif low.startswith('*element') and not low.startswith('*element output'):
                mode = 'elem'
                etype = re.search(r'type=(\w+)', line, re.I).group(1).upper()
            continue
        if not line.strip():
            continue
        p = [v.strip() for v in line.split(',') if v.strip()]
        if mode == 'node':
            nodes[int(p[0])] = tuple(float(v) for v in p[1:4])
        elif mode == 'elem':
            elems.append((etype, int(p[0]), [int(v) for v in p[1:]]))
    return nodes, elems


def last_fields(case):
    a = compare.last_frame_arrays(case)
    lab = a['NodeLabel']
    u = {l: a['U'][3 * i:3 * i + 3] for i, l in enumerate(lab)}
    cells = {}
    if 'ElementLabel' in a:
        for i, l in enumerate(a['ElementLabel']):
            cells[l] = {k: a[k][i] for k in ('S11', 'S22', 'S12', 'S_Mises') if k in a}
    return u, cells


def deformed(nodes, u, scale):
    return {n: tuple(x + scale * d for x, d in zip(p, u[n])) for n, p in nodes.items()}


def deformation_factor(nodes, u, frac=0.12):
    size = max(max(p[k] for p in nodes.values()) - min(p[k] for p in nodes.values()) for k in range(3))
    umax = max(math.sqrt(sum(c * c for c in v)) for v in u.values())
    return frac * size / umax if umax > 0 else 1.0


def faces(elems, nodes):
    """Surface polygons: shells as they are, top faces (largest z) of solids."""
    out = []
    for etype, lab, c in elems:
        if etype.startswith(('S4', 'S3')):
            out.append((lab, c[:4] if etype.startswith('S4') else c[:3]))
        elif etype == 'C3D8R':
            top = c[4:8] if sum(nodes[n][2] for n in c[4:8]) > sum(nodes[n][2] for n in c[:4]) else c[:4]
            out.append((lab, top))
    return out


def draw_surface(ax, case, inp, value='umag', title='', view=(25, -60), two_d=False, mark=None):
    nodes, elems = read_mesh(inp)
    u, cells = last_fields(case)
    sc = deformation_factor(nodes, u) if value == 'umag' else 0.0  # stresses on the undeformed shape
    xd = deformed(nodes, u, sc)
    polys = faces(elems, nodes)
    if value == 'umag':
        vals = [sum(math.sqrt(sum(c * c for c in u[n])) for n in f) / len(f) for _, f in polys]
        label = '|U|'
    else:
        vals = [cells[lab][value] for lab, _ in polys]
        label = {'S22': 'sigma_yy (MPa)'}.get(value, value)
    vmin, vmax = min(vals), max(vals)
    cmap = plt.get_cmap(CMAP)
    colors = [cmap((v - vmin) / (vmax - vmin) if vmax > vmin else 0.5) for v in vals]
    if two_d:
        ax.add_collection(PolyCollection([[xd[n][:2] for n in f] for _, f in polys], facecolors=colors,
                                         edgecolors='k', linewidths=0.15))
        if sc:
            ax.add_collection(PolyCollection([[nodes[n][:2] for n in f] for _, f in polys], facecolors='none',
                                             edgecolors='0.6', linewidths=0.3, linestyles=':'))
        ax.autoscale()
        ax.set_aspect('equal')
        if mark:
            ax.plot(*mark, 'o', color='crimson', ms=5)
    else:
        ax.add_collection3d(Poly3DCollection([[nodes[n] for n in f] for _, f in polys], facecolors=(0, 0, 0, 0),
                                             edgecolors='0.7', linewidths=0.2))
        ax.add_collection3d(Poly3DCollection([[xd[n] for n in f] for _, f in polys], facecolors=colors,
                                             edgecolors='k', linewidths=0.15))
        pts = list(nodes.values()) + list(xd.values())
        lo = [min(p[k] for p in pts) for k in range(3)]
        hi = [max(p[k] for p in pts) for k in range(3)]
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect([max(h - l, 1e-9) for l, h in zip(lo, hi)])
        ax.view_init(*view)
        ax.set_axis_off()
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
    plt.colorbar(sm, ax=ax, shrink=0.7, label=label)
    ax.set_title('%s\n(displacements x %.3g)' % (title, sc) if sc else title, fontsize=9)


def elastica_shape(k, n=4000):
    """Exact deformed shape of the inextensible cantilever under a dead end load P L^2 / EI = k (shooting)."""
    def shoot(c0, steps):
        h = 1.0 / steps
        y = [0.0, c0, 0.0, 0.0]
        path = [(0.0, 0.0)]
        f = lambda y: [y[1], -k * math.cos(y[0]), math.cos(y[0]), math.sin(y[0])]
        for _ in range(steps):
            k1 = f(y); k2 = f([a + h / 2 * b for a, b in zip(y, k1)])
            k3 = f([a + h / 2 * b for a, b in zip(y, k2)]); k4 = f([a + h * b for a, b in zip(y, k3)])
            y = [a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4)]
            path.append((y[2], y[3]))
        return y, path
    a, b = 0.0, 1.5 * k
    for _ in range(60):
        m = 0.5 * (a + b)
        if shoot(m, 2000)[0][1] > 0:
            b = m
        else:
            a = m
    return shoot(0.5 * (a + b), n)[1]


def draw_elastica(ax, case, inp):
    nodes, elems = read_mesh(inp)
    u, _ = last_fields(case)
    L = max(p[0] for p in nodes.values())
    ex = elastica_shape(5.0)
    ax.plot([L * x for x, _ in ex], [L * y for _, y in ex], '-', color='0.3', lw=2.5, label='exact elastica')
    order = sorted(nodes)
    ax.plot([nodes[n][0] + u[n][0] for n in order], [nodes[n][1] + u[n][1] for n in order], 'o-', color='crimson',
            ms=3, lw=1, label='asclepios B31')
    ax.plot([0, L], [0, 0], ':', color='0.6', label='undeformed')
    ax.set_aspect('equal')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title('Deformed shape (true scale)', fontsize=9)


def figure(key, title, series, runs, val_case):
    fig = plt.figure(figsize=(10.5, 4.2))
    two_d = key in ('cook_membrane', 'nafems_le1', 'lame_cylinder', 'elastica')
    ax1 = fig.add_subplot(1, 2, 1) if two_d else fig.add_subplot(1, 2, 1, projection='3d')
    case, inp = val_case
    if key == 'elastica':
        draw_elastica(ax1, case, inp)
    elif key == 'nafems_le1':
        draw_surface(ax1, case, inp, value='S22', title='sigma_yy, validation mesh (point D in red)', two_d=True,
                     mark=(2.0, 0.0))
    else:
        view = {'scordelis_lo': (25, -130), 'pinched_cylinder': (25, -35), 'pinched_hemisphere': (20, 30),
                'twisted_beam': (20, -70)}.get(key, (25, -60))
        draw_surface(ax1, case, inp, title='Deformed shape, validation mesh', view=view, two_d=two_d)
    ax2 = fig.add_subplot(1, 2, 2)
    markers = 'osd^v'
    tol = None
    for k, (label, pts) in enumerate(runs):
        nchk = len(pts[0][1])
        for j in range(nchk):
            xs = [p[0] for p in pts]
            ys = [p[1][j][1] / p[1][j][2] for p in pts]
            lab = label if nchk == 1 else '%s, %s' % (label, pts[0][1][j][0])
            ax2.plot(xs, ys, markers[(k * nchk + j) % 5] + '-', label=lab, ms=5)
            tol = pts[0][1][j][4]
    ax2.axhline(1.0, color='k', lw=1)
    if tol:
        ax2.axhspan(1 - tol, 1 + tol, color='0.85', zorder=0, label='tolerance +/- %g %%' % (100 * tol))
    ax2.set_xscale('log', base=2)
    ax2.set_xlabel('elements along the beam' if key == 'elastica' else 'elements per side')
    ax2.set_ylabel('asclepios / published value')
    ax2.set_title('Convergence', fontsize=9)
    ax2.legend(fontsize=7)
    ax2.grid(True, which='both', alpha=0.3)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    path = os.path.join(FIG, key + '.png')
    fig.savefig(path, dpi=90, metadata={'Software': None})
    plt.close(fig)
    return os.path.relpath(path, HERE)


def main():
    run = '--no-run' not in sys.argv
    os.makedirs(WORK, exist_ok=True)
    validation = {c['name']: c for c in json.load(open(os.path.join(HERE, 'cases', 'manifest.json')))
                  if c.get('reference') == 'literature'}
    md = ['# Published benchmarks', '',
          'Generated by `validation/literature.py` (`make literature`). Classic benchmarks of the finite element '
          'literature, run with asclepios on series of meshes and compared with the converged value of the cited '
          'publication. Left: deformed shape (or stress) of the mesh used in the validation suite (`make validate`, '
          'cases `lit_*`); right: result divided by the published value against the mesh density, with the '
          'tolerance of the validation case. Decks: `validation/gen_literature.py`.', '']
    for key, title, series in BENCHMARKS:
        runs, rows, val_case, source = [], [], None, None
        for label, gen, arglist in series:
            pts = []
            for args in arglist:
                case, res = build_and_run(gen, args, run)
                base = case['name'].split('__')[0]
                source = case['literature']['source']
                inp = os.path.join(WORK, case['name'] + '.inp')
                if base in validation and val_case is None and \
                        len(read_mesh(inp)[0]) == len(read_mesh(os.path.join(HERE, 'cases', base + '.inp'))[0]):
                    val_case = (case, inp)  # same mesh as the validation case
                pts.append((elements_per_side(args), res['checks']))
                for lbl, v, ref, err, tol in res['checks']:
                    rows.append((label, '%s' % ' x '.join(str(a) for a in args if isinstance(a, int)), lbl, v, ref, err))
            runs.append((label, pts))
        if val_case is None:  # first mesh of the first series that matches the validation mesh size
            label, gen, arglist = series[0]
            case, _ = build_and_run(gen, arglist[-2], False)
            val_case = (case, os.path.join(WORK, case['name'] + '.inp'))
        img = figure(key, title, series, runs, val_case)
        md += ['## ' + title, '', 'Source: ' + source + '.', '', '![%s](%s)' % (title, img), '',
               '| Element | Mesh | Quantity | asclepios | Published | Error |', '|---|---|---|---|---|---|']
        md += ['| %s | %s | %s | %.5g | %.5g | %.2f %% |' % (a, m, q, v, r, 100 * e) for a, m, q, v, r, e in rows]
        md.append('')
        print('%-20s done' % key)
    with open(os.path.join(HERE, 'LITERATURE.md'), 'w') as f:
        f.write('\n'.join(md) + '\n')


if __name__ == '__main__':
    main()
