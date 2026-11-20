"""Benchmarks with published reference values: classic problems of the finite element literature, run on refined
meshes and compared with the converged value given in the cited publication (reference = 'literature' in the
manifest). Each case lists its checks: a nodal displacement component, or a stress extrapolated to a point from the
element values next to it. Tolerances cover the discretisation error left on the refined mesh and the spread of the
published values.
"""
import math

from deck import Deck, fmt, near, quad_plate, register, static_step, tris_from_quads, faces_on
from gen_more import beam, cylinder_nodes, hex_ok, nodal_area_loads, shell, solid


def check_u(node, dof, value, tol, sign=1.0, label=''):
    return {'kind': 'U', 'node': node, 'dof': dof, 'sign': sign, 'value': value, 'tol': tol, 'label': label}


def check_s(groups, dist, comp, value, tol, deg, label=''):
    """Stress at a point, extrapolated (least-squares polynomial of degree deg through the group means) from groups
    of elements whose centroids lie at the given distances from the point."""
    return {'kind': 'S', 'groups': groups, 'dist': dist, 'comp': comp, 'value': value, 'tol': tol, 'deg': deg,
            'label': label}


def lit(d, element, mode, description, source, checks, fix, rot=True, extra_mon=()):
    mon = [c['node'] for c in checks if c['kind'] == 'U'] + list(extra_mon)
    register(d, 'static', element, mode, description, mon, fix, rot, 'published', 'linear', reference='literature',
             category='feature', group='literature', twin=False,
             literature={'source': source, 'checks': checks})


def quads_of(ids, ni, nj):
    return [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(ni) for j in range(nj)]


def scordelis_lo(etype, n):
    d = Deck('lit_scordelis_lo_%s' % etype.lower(), '%s Scordelis-Lo roof, %d x %d quarter model' % (etype, n, n))
    ids = cylinder_nodes(d, n, n, 25.0, 25.0, 0.0, math.radians(40.0))
    quads = quads_of(ids, n, n)
    shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 0.25, 'ROOF')
    d.nset('DIAPH', [ids[n, j] for j in range(n + 1)])
    d.nset('XS', [ids[0, j] for j in range(n + 1)])
    d.nset('YS', [ids[i, 0] for i in range(n + 1)])
    d.steps += static_step('Load', ['*Boundary', 'DIAPH, 2, 3', 'DIAPH, 4, 4', 'XS, XSYMM', 'YS, YSYMM', '*Cload'] +
                           nodal_area_loads(d, quads, -90.0), False, rot=True, nframes=1)
    lit(d, etype, 'Scordelis-Lo roof', 'Cylindrical roof on rigid diaphragms under self-weight (R = 25, L = 50, '
        't = 0.25, E = 4.32e8, nu = 0, 90 per unit area): vertical displacement at the middle of the free edge.',
        'MacNeal and Harder (1985), Finite Elem. Anal. Des. 1:3-20: 0.3024',
        [check_u(ids[0, n], 3, 0.3024, 0.02, sign=-1.0, label='w at the middle of the free edge')],
        [ids[n, j] for j in range(n + 1)])


def pinched_cylinder(n):
    d = Deck('lit_pinched_cylinder_s4r', 'S4R pinched cylinder with rigid end diaphragms, %d x %d octant model' % (n, n))
    ids = cylinder_nodes(d, n, n, 300.0, 300.0, 0.0, math.pi / 2)
    shell(d, 'S4R', quads_of(ids, n, n), 3.0, 'CYL')
    d.nset('DIAPH', [ids[n, j] for j in range(n + 1)])
    d.nset('XS', [ids[0, j] for j in range(n + 1)])
    d.nset('YS', [ids[i, 0] for i in range(n + 1)])
    d.nset('ZS', [ids[i, n] for i in range(n + 1)])
    d.steps += static_step('Load', ['*Boundary', 'DIAPH, 2, 3', 'DIAPH, 4, 4', 'XS, XSYMM', 'YS, YSYMM', 'ZS, ZSYMM',
                                    '*Cload', '%d, 3, -0.25' % ids[0, 0]], False, rot=True, nframes=1)
    lit(d, 'S4R', 'pinched cylinder', 'Cylinder (R = 300, L = 600, t = 3, E = 3e6, nu = 0.3) pinched by two opposite '
        'unit loads at mid-length: radial displacement under the load.',
        'Belytschko et al. (1985), Comput. Methods Appl. Mech. Engrg. 51:221-258: 1.8248e-5',
        [check_u(ids[0, 0], 3, 1.8248e-5, 0.02, sign=-1.0, label='displacement under the load')],
        [ids[n, j] for j in range(n + 1)])


def pinched_hemisphere(n):
    d = Deck('lit_pinched_hemisphere_s4r', 'S4R pinched hemisphere with an 18 degree hole, %d x %d quarter model' % (n, n))
    ids = {}
    for i in range(n + 1):
        phi = math.radians(18.0 + 72.0 * i / n)
        for j in range(n + 1):
            th = math.pi / 2 * j / n
            ids[i, j] = d.node(10.0 * math.sin(phi) * math.cos(th), 10.0 * math.sin(phi) * math.sin(th),
                               10.0 * math.cos(phi))
    quads = [[ids[i, j], ids[i, j + 1], ids[i + 1, j + 1], ids[i + 1, j]] for i in range(n) for j in range(n)]
    shell(d, 'S4R', quads, 0.04, 'HEMI')
    d.nset('YS', [ids[i, 0] for i in range(n + 1)])
    d.nset('XS', [ids[i, n] for i in range(n + 1)])
    d.steps += static_step('Load', ['*Boundary', 'YS, YSYMM', 'XS, XSYMM', '%d, 3, 3' % ids[n, 0], '*Cload',
                                    '%d, 1, 1.' % ids[n, 0], '%d, 2, -1.' % ids[n, n]], False, rot=True, nframes=1)
    lit(d, 'S4R', 'pinched hemisphere', 'Hemispherical shell (R = 10, t = 0.04, E = 6.825e7, nu = 0.3) with an 18 '
        'degree hole, loaded by alternating radial forces of 2 at the equator: displacement under the load.',
        'MacNeal and Harder (1985), Finite Elem. Anal. Des. 1:3-20: 0.094',
        [check_u(ids[n, 0], 1, 0.094, 0.02, label='outward displacement under the load'),
         check_u(ids[n, n], 2, 0.094, 0.02, sign=-1.0, label='inward displacement under the load')],
        [ids[n, 0]])


def twisted_beam(direction, nx, ny):
    d = Deck('lit_twisted_beam_%s_s4r' % direction, 'S4R twisted beam (90 degree pretwist), tip load along %s, '
             '%d x %d mesh' % (direction, nx, ny))
    ids = {}
    for i in range(nx + 1):
        a = math.radians(90.0 * i / nx)
        for j in range(ny + 1):
            s = -0.55 + 1.1 * j / ny
            ids[i, j] = d.node(12.0 * i / nx, s * math.cos(a), s * math.sin(a))
    shell(d, 'S4R', quads_of(ids, nx, ny), 0.32, 'TWIST')
    fix = [ids[0, j] for j in range(ny + 1)]
    tip = [ids[nx, j] for j in range(ny + 1)]
    dof = 2 if direction == 'y' else 3
    w = [0.5 if j in (0, ny) else 1.0 for j in range(ny + 1)]
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, %d, %s' % (nd, dof, fmt(wi / sum(w))) for nd, wi in zip(tip, w)], False, rot=True,
                           nframes=1)
    value = 1.754e-3 if direction == 'y' else 5.424e-3
    lit(d, 'S4R', 'twisted beam', 'Cantilever (L = 12, width 1.1, t = 0.32, E = 29e6, nu = 0.22) twisted by 90 degrees '
        'over its length, unit tip load %s: tip displacement along the load.' % (
            'in the width direction of the clamped end' if direction == 'y' else 'normal to the clamped end'),
        'MacNeal and Harder (1985), Finite Elem. Anal. Des. 1:3-20: %g' % value,
        [check_u(ids[nx, ny // 2], dof, value, 0.02, label='tip displacement')], fix)


def cook_membrane(etype, n):
    d = Deck('lit_cook_membrane_%s' % etype.lower(), "%s Cook's membrane, %d x %d mesh" % (etype, n, n))
    ids = {}
    for i in range(n + 1):
        x = 48.0 * i / n
        y0, y1 = 44.0 * x / 48.0, 44.0 + 16.0 * x / 48.0
        for j in range(n + 1):
            ids[i, j] = d.node(x, y0 + (y1 - y0) * j / n)
    quads = quads_of(ids, n, n)
    shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 1.0, 'COOK')
    d.nset('NALL', list(d.nodes))
    fix = [ids[0, j] for j in range(n + 1)]
    right = [ids[n, j] for j in range(n + 1)]
    d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 2', '*Cload'] +
                           ['%d, 2, %s' % (nd, fmt((0.5 if k in (0, n) else 1.0) / n)) for k, nd in enumerate(right)],
                           False, rot=True, nframes=1)
    lit(d, etype, "Cook's membrane", 'Tapered panel (E = 1, nu = 1/3, plane stress, thickness 1) clamped on one edge, '
        'unit shear on the other: vertical displacement of the middle of the loaded edge (48, 52).',
        'Cook (1974), J. Struct. Div. ASCE 100:1851-1863; converged value 23.96 at the middle of the loaded edge',
        [check_u(ids[n, n // 2], 2, 23.96, 0.02, label='vertical displacement of the middle of the loaded edge')],
        fix, extra_mon=[ids[n, n]])


def nafems_le1(etype, nt, ns):
    d = Deck('lit_nafems_le1_%s' % etype.lower(), '%s NAFEMS LE1 elliptic membrane, %d x %d mesh' % (etype, nt, ns))
    ids = {}
    for i in range(nt + 1):
        th = math.pi / 2 * i / nt
        for j in range(ns + 1):
            s = j / ns
            ids[i, j] = d.node((2.0 + 1.25 * s) * math.cos(th), (1.0 + 1.75 * s) * math.sin(th), 0.0)
    quads = [[ids[i, j], ids[i, j + 1], ids[i + 1, j + 1], ids[i + 1, j]] for i in range(nt) for j in range(ns)]
    conns = tris_from_quads(quads) if etype == 'S3R' else quads
    labels = d.elements(etype, 'EALL', conns)
    d.materials.add('STEEL')
    d.model += ['*Shell Section, elset=EALL, material=STEEL', '0.1, 5']
    per = 2 if etype == 'S3R' else 1
    groups, dist = [], []
    npts, deg = 3, 2  # quadratic through the first three elements (pairs of triangles) of the strip along x
    for j in range(npts):  # first elements (pairs of triangles) of the strip next to the x axis, from point D
        k = j * per  # quads are ordered with j fastest for i = 0
        groups.append(labels[k:k + per])
        q = quads[j]
        dist.append(sum(d.nodes[nd][0] for nd in q) / 4 - 2.0)
    d.nset('NALL', list(d.nodes))
    d.nset('XS', [ids[nt, j] for j in range(ns + 1)])
    d.nset('YS', [ids[0, j] for j in range(ns + 1)])
    loads, acc = [], {}
    for i in range(nt):  # 10 MPa normal tension on the outer edge, consistent nodal forces
        a, b = d.nodes[ids[i, ns]], d.nodes[ids[i + 1, ns]]
        nx_, ny_ = (b[1] - a[1]), -(b[0] - a[0])
        for nd in (ids[i, ns], ids[i + 1, ns]):
            f = acc.setdefault(nd, [0.0, 0.0])
            f[0] += 0.5 * 10.0 * 0.1 * nx_
            f[1] += 0.5 * 10.0 * 0.1 * ny_
    for nd, (fx, fy) in sorted(acc.items()):  # (nx_, ny_) is the outward normal times the edge length
        loads += ['%d, 1, %s' % (nd, fmt(fx)), '%d, 2, %s' % (nd, fmt(fy))]
    d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'XS, XSYMM', 'YS, YSYMM', '*Cload'] + loads, False,
                           rot=True, nframes=1)
    lit(d, etype, 'NAFEMS LE1', 'Elliptic membrane (E = 210000, nu = 0.3, thickness 0.1) under 10 MPa outer normal '
        'tension: hoop stress sigma_yy at the inner point D on the x axis, extrapolated from the elements next to D.',
        'NAFEMS (1990), The Standard NAFEMS Benchmarks, test LE1: 92.7 MPa',
        [check_s(groups, dist, 'S22', 92.7, 0.03, deg, label='sigma_yy at D')],
        [ids[nt, j] for j in range(ns + 1)] + [ids[0, j] for j in range(ns + 1)], extra_mon=[ids[0, 0]])


def clamped_plate(etype, n):
    d = Deck('lit_clamped_plate_%s' % etype.lower(), '%s clamped square plate under pressure (a / t = 100), %d x %d '
             'mesh' % (etype, n, n))
    quads, ids = quad_plate(d, n, n, 20.0, 20.0)
    shell(d, etype, tris_from_quads(quads) if etype.startswith('S3') else quads, 0.2, 'STEEL')
    d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
    edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
    q = 0.01
    dflex = 210000.0 * 0.2 ** 3 / (12 * (1 - 0.3 ** 2))
    w = 0.00126532 * q * 20.0 ** 4 / dflex
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, %s' % fmt(q)], False, rot=True,
                           nframes=1)
    lit(d, etype, 'clamped plate', 'Clamped square plate (a = 20, t = 0.2, E = 210000, nu = 0.3) under uniform '
        'pressure 0.01: deflection of the centre (thin plate theory, w = 0.00126532 q a^4 / D).',
        'Timoshenko and Woinowsky-Krieger (1959), Theory of Plates and Shells, sec. 44; series value from Taylor and '
        'Govindjee (2004), Commun. Numer. Meth. Engng 20:757-765: 0.00126532 q a^4 / D',
        [check_u(ids[n // 2, n // 2], 3, w, 0.02, sign=-1.0, label='centre deflection')], edge)


def thick_plate(etype, n):
    d = Deck('lit_thick_plate_%s' % etype.lower(), '%s simply supported thick square plate under pressure (a / t = 10), '
             '%d x %d mesh' % (etype, n, n))
    a, t, q, nu = 10.0, 1.0, 0.1, 0.3
    quads, ids = quad_plate(d, n, n, a, a)
    shell(d, etype, tris_from_quads(quads) if etype.startswith('S3') else quads, t, 'STEEL')
    d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
    # hard simple support: w and the rotation about the edge normal held on every edge
    d.nset('EX', d.select(lambda x, y, z: near(x, 0) or near(x, a)))
    d.nset('EY', d.select(lambda x, y, z: near(y, 0) or near(y, a)))
    d.nset('NALL', list(d.nodes))
    edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, a) or near(y, a))
    dflex = 210000.0 * t ** 3 / (12 * (1 - nu ** 2))
    kgt = 5.0 / 6.0 * 210000.0 / (2 * (1 + nu)) * t
    w = 0.00406235 * q * a ** 4 / dflex + 0.0736713 * q * a * a / kgt
    d.steps += static_step('Load', ['*Boundary', 'NALL, 1, 2', 'EX, 3, 4', 'EY, 3, 3', 'EY, 5, 5',
                                    '*Dsload', 'STOP, P, %s' % fmt(q)], False, rot=True, nframes=1)
    lit(d, etype, 'thick plate', 'Simply supported (hard) square plate (a = 10, t = 1, E = 210000, nu = 0.3, shear '
        'factor 5/6) under uniform pressure 0.1: deflection of the centre, Mindlin plate theory (shear part 5 %).',
        'Wang, Reddy and Lee (2000), Shear Deformable Beams and Plates, Elsevier: w = w_K + M_K / (k G t) for simply '
        'supported polygonal plates; w_K = 0.00406235 q a^4 / D (Timoshenko and Woinowsky-Krieger 1959, sec. 30), '
        'Marcus moment at the centre M_K = 0.0736713 q a^2; w = 0.4273 q a^4 / (100 D)',
        [check_u(ids[n // 2, n // 2], 3, w, 0.02, sign=-1.0, label='centre deflection')], edge)


def lame_cylinder(etype, nr, nt):
    d = Deck('lit_lame_cylinder_%s' % etype.lower(), '%s thick-walled cylinder under internal pressure (plane strain, '
             'quarter model), %d x %d mesh' % (etype, nr, nt))
    ri, ro, p, hz = 10.0, 20.0, 100.0, 1.0
    layers = []
    for k in range(2):
        ids = {}
        for i in range(nr + 1):
            r = ri + (ro - ri) * i / nr
            for j in range(nt + 1):
                th = math.pi / 2 * j / nt
                ids[i, j] = d.node(r * math.cos(th), r * math.sin(th), hz * k)
        layers.append(ids)
    lo, hi = layers
    hexes = [hex_ok(d, [lo[i, j], lo[i + 1, j], lo[i + 1, j + 1], lo[i, j + 1],
                        hi[i, j], hi[i + 1, j], hi[i + 1, j + 1], hi[i, j + 1]]) for i in range(nr) for j in range(nt)]
    from deck import tets_from_hexes
    conns = tets_from_hexes(d, hexes) if etype == 'C3D4' else hexes
    el = solid(d, etype, conns, 'STEEL')
    inner = faces_on(d, etype, el, conns, lambda x, y, z: near(math.hypot(x, y), ri))
    d.model += ['*Surface, type=ELEMENT, name=SIN'] + ['%d, S%d' % f for f in inner]
    d.nset('NALL', list(d.nodes))
    d.nset('XS', d.select(lambda x, y, z: near(x, 0)))
    d.nset('YS', d.select(lambda x, y, z: near(y, 0)))
    d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 3', 'XS, 1, 1', 'YS, 2, 2', '*Dsload', 'SIN, P, %s' % fmt(p)],
                           False, nframes=1)
    E, nu = 210000.0, 0.3
    c = (1 + nu) * p * ri ** 2 / (E * (ro ** 2 - ri ** 2))
    u = lambda r: c * ((1 - 2 * nu) * r + ro ** 2 / r)
    lit(d, etype, 'Lame cylinder', 'Thick-walled cylinder (ri = 10, ro = 20, E = 210000, nu = 0.3) under 100 MPa '
        'internal pressure, plane strain: radial displacement of the inner and outer surfaces.',
        'Lame solution, e.g. Timoshenko and Goodier (1970), Theory of Elasticity, sec. 28: u(r) = (1 + nu) p ri^2 / '
        '(E (ro^2 - ri^2)) ((1 - 2 nu) r + ro^2 / r)',
        [check_u(lo[0, 0], 1, u(ri), 0.01, label='inner radial displacement'),
         check_u(lo[nr, 0], 1, u(ro), 0.01, label='outer radial displacement')],
        d.nsets['YS'], rot=False)


def elastica(k, n):
    d = Deck('lit_elastica_b31', 'B31 cantilever under a large end load (P L^2 / EI = %g), %d elements' % (k, n))
    L, EI = 100.0, 210000.0 / 12.0
    nodes = [d.node(L * i / n, 0.0, 0.0) for i in range(n + 1)]
    beam(d, [[nodes[i], nodes[i + 1]] for i in range(n)], 'RECT', '1., 1.')
    P = k * EI / L ** 2
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % nodes[0], '*Cload', '%d, 2, %s' % (nodes[-1], fmt(P))],
                           True, rot=True, nframes=10)
    lit(d, 'B31', 'elastica', 'Slender cantilever (L = 100, square section 1 x 1, E = 210000) under a dead end load '
        'P L^2 / EI = %g: tip displacements along and across the load (inextensible elastica).' % k,
        'Elliptic-integral solution of the elastica, tables in Mattiasson (1981), Int. J. Numer. Meth. Engng '
        '17:145-153: 0.71379 L and 0.38763 L',
        [check_u(nodes[-1], 2, 0.71379 * L, 0.01, label='tip displacement along the load'),
         check_u(nodes[-1], 1, 0.38763 * L, 0.01, sign=-1.0, label='tip shortening')], [nodes[0]])


def generate():
    scordelis_lo('S4R', 16)
    scordelis_lo('S3R', 48)
    pinched_cylinder(32)
    pinched_hemisphere(16)
    twisted_beam('y', 24, 4)
    twisted_beam('z', 24, 4)
    for etype in ('S4R', 'S3R'):
        cook_membrane(etype, 48)
        nafems_le1(etype, 48, 24)
    clamped_plate('S4R', 16)
    clamped_plate('S3R', 32)
    clamped_plate('S3RS', 16)
    thick_plate('S3RS', 16)
    lame_cylinder('C3D8R', 8, 16)
    lame_cylinder('C3D4', 16, 32)
    elastica(5.0, 40)
