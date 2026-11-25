"""More validation cases: patch tests and distorted meshes, standard benchmarks, input-file features, Abaqus features
not covered elsewhere (elements, materials, contact options, connectors, amplitudes, procedures), application models
and performance runs. Closed-form references where the nodal solution is known by hand (see analytic.py), Abaqus
results otherwise.
"""
import math

from deck import (CURRENT, MATERIALS, MODE, QS_DURATION, Deck, beam_line, edge_shares, explicit_step, faces_on, fmt,
                  hex_block, near, output_block, quad_plate, register, static_step, tris_from_quads)

NS = 'not supported by asclepios yet: '

MATERIALS.update({
    'ROOF': '*Material, name=ROOF\n*Density\n1e-09\n*Elastic\n4.32e+08, 0.\n',
    'CYL': '*Material, name=CYL\n*Density\n1e-09\n*Elastic\n3e+06, 0.3\n',
    'HEMI': '*Material, name=HEMI\n*Density\n1e-09\n*Elastic\n6.825e+07, 0.3\n',
    'COOK': '*Material, name=COOK\n*Density\n1e-09\n*Elastic\n1., 0.333333333333\n',
    'TWIST': '*Material, name=TWIST\n*Density\n1e-09\n*Elastic\n2.9e+07, 0.22\n',
    'COPPER': ('*Material, name=COPPER\n*Density\n8.93e-09\n*Elastic\n117000., 0.35\n'
               '*Plastic\n400., 0.\n500., 1.\n'),
    'RUBBER_INC': '*Material, name=RUBBER_INC\n*Density\n1.1e-09\n*Hyperelastic, neo hooke\n0.5, 0.0001\n',
    'RUBBER_495': '*Material, name=RUBBER_495\n*Density\n1.1e-09\n*Hyperelastic, neo hooke\n0.5, 0.0202\n',
    'OGDEN': ('*Material, name=OGDEN\n*Density\n1.1e-09\n*Hyperelastic, ogden, n=2\n'
              '0.6, 2.5, -0.1, -2., 0.05, 0.\n'),
    'YEOH': '*Material, name=YEOH\n*Density\n1.1e-09\n*Hyperelastic, yeoh\n0.5, -0.02, 0.003, 0.04, 0., 0.\n',
    'ORTHO': ('*Material, name=ORTHO\n*Density\n1.6e-09\n*Elastic, type=ENGINEERING CONSTANTS\n'
              '150000., 10000., 10000., 0.3, 0.3, 0.4, 5000., 5000.\n3500.\n'),
    'HGO': ('*Material, name=HGO\n*Density\n1.1e-09\n*Anisotropic Hyperelastic, holzapfel, local directions=2\n'
            '0.05, 0.02, 2., 10., 0.1\n'),
    'ABQ_SUPER_ELASTIC_N3D': ('*Material, name=ABQ_SUPER_ELASTIC_N3D\n*Density\n6.5e-09\n*Depvar\n31\n'
                               '*User Material, constants=15\n'
                               '60000., 0.3, 30000., 0.3, 0.05, 6.7, 520., 600.\n'
                               '37., 6.7, 300., 200., 600., 0.05, 0.\n'),
    'STEEL_RATE': ('*Material, name=STEEL_RATE\n*Density\n7.85e-09\n*Elastic\n210000., 0.3\n'
                   '*Plastic\n250., 0.\n350., 0.2\n*Rate Dependent, type=POWER LAW\n40.4, 5.\n'),
    'STEEL_DMG': ('*Material, name=STEEL_DMG\n*Density\n7.85e-09\n*Elastic\n210000., 0.3\n'
                  '*Plastic\n250., 0.\n300., 0.2\n*Damage Initiation, criterion=DUCTILE\n0.1, 0.333, 0.\n'
                  '*Damage Evolution, type=DISPLACEMENT\n0.05\n'),
    'COND': ('*Material, name=COND\n*Density\n7.85e-09\n*Conductivity\n50.\n*Specific Heat\n4.6e+08\n'
             '*Elastic\n210000., 0.3\n*Expansion\n1.2e-05\n'),
})


def shell(d, etype, conns, t, mat, elset='EALL', extra=''):
    d.elements(etype, elset, conns)
    d.materials.add(mat)
    d.model += ['*Shell Section, elset=%s, material=%s%s' % (elset, mat, extra), '%s, 5' % fmt(t)]


def solid(d, etype, conns, mat, elset='EALL'):
    el = d.elements(etype, elset, conns)
    d.materials.add(mat)
    d.model += ['*Solid Section, elset=%s, material=%s' % (elset, mat), ',']
    return el


def beam(d, conns, section, dims, mat='STEEL', elset='EALL', n1='0., 0., -1.'):
    d.elements('B31', elset, conns)
    d.materials.add(mat)
    d.model += ['*Beam Section, elset=%s, material=%s, section=%s' % (elset, mat, section), dims, n1]


def general_contact(d, friction=None, extra_model=None):
    d.model += ['*Surface Interaction, name=IP'] + (extra_model or []) + (['*Friction', fmt(friction)] if friction else [])
    return ['*Contact', '*Contact Inclusions, ALL EXTERIOR', '*Contact Property Assignment', ' , , IP']


def hex_ok(d, h):
    """Orients an 8-node hexahedron so that its volume is positive (nodes 1-4 one face, 5-8 the opposite one)."""
    p = [d.nodes[n] for n in h]
    c = [sum(q[i] for q in p) / 8 for i in range(3)]
    b = [sum(p[k][i] for k in range(4)) / 4 for i in range(3)]
    t = [sum(p[k][i] for k in range(4, 8)) / 4 for i in range(3)]
    u = [p[1][i] - p[0][i] for i in range(3)]
    v = [p[3][i] - p[0][i] for i in range(3)]
    nrm = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    if sum(nrm[i] * (t[i] - b[i]) for i in range(3)) < 0:
        return [h[0], h[3], h[2], h[1], h[4], h[7], h[6], h[5]]
    return h


# --------------------------------------------------------------------------------- patch tests, robustness
MH_INNER = [(0.04, 0.02), (0.18, 0.03), (0.16, 0.08), (0.08, 0.08)]   # MacNeal-Harder membrane patch (x 100)


def membrane_patch_mesh(d):
    c = [d.node(0., 0.), d.node(24., 0.), d.node(24., 12.), d.node(0., 12.)]
    i = [d.node(100 * x, 100 * y) for x, y in MH_INNER]
    quads = [[c[0], c[1], i[1], i[0]], [c[1], c[2], i[2], i[1]], [c[2], c[3], i[3], i[2]], [c[3], c[0], i[0], i[3]],
             [i[0], i[1], i[2], i[3]]]
    return c, i, quads


def patch_features():
    for etype in ('S4R', 'S4', 'S3R'):
        d = Deck('patch_membrane_%s' % etype.lower(), '%s membrane patch test (distorted elements)' % etype)
        c, i, quads = membrane_patch_mesh(d)
        shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 0.1, 'STEEL')
        d.nset('NALL', list(d.nodes))
        body = ['*Boundary', 'NALL, 3, 5', '*Boundary']
        for n in c:
            x, y, z = d.nodes[n]
            body += ['%d, 1, 1, %s' % (n, fmt(1e-3 * (x + y / 2))), '%d, 2, 2, %s' % (n, fmt(1e-3 * (y + x / 2)))]
        d.steps += static_step('Patch', body, False, rot=True, nframes=1)
        register(d, 'static', etype, 'membrane patch test', 'Constant strain imposed on the corners of a patch of '
                 'distorted elements: interior nodes and stresses must be exact.', i, [], True, 'exact', 'linear',
                 reference='analytic', category='feature', group='patch')

    for etype in ('S4R', 'S4'):
        d = Deck('patch_plate_%s' % etype.lower(), '%s plate bending patch test (constant curvature, distorted elements)'
                 % etype)
        c, i, quads = membrane_patch_mesh(d)
        shell(d, etype, quads, 0.5, 'STEEL')
        d.nset('NALL', list(d.nodes))
        body = ['*Boundary', 'NALL, 1, 2', '*Boundary']
        for n in c:
            x, y, z = d.nodes[n]
            body += ['%d, 3, 3, %s' % (n, fmt(1e-3 * (x * x + x * y + y * y) / 2)),
                     '%d, 4, 4, %s' % (n, fmt(1e-3 * (x / 2 + y))), '%d, 5, 5, %s' % (n, fmt(-1e-3 * (x + y / 2)))]
        d.steps += static_step('Patch', body, False, rot=True, nframes=1)
        register(d, 'static', etype, 'plate bending patch test', 'Constant curvature w = (x^2 + xy + y^2) / 2000 imposed '
                 'on the corners: interior deflections, rotations and bending stresses must be exact.', i, [], True,
                 'exact', 'linear', reference='analytic', category='feature', group='patch')

    inner = [(0.249, 0.342, 0.192), (0.826, 0.288, 0.288), (0.850, 0.649, 0.263), (0.273, 0.750, 0.230),
             (0.320, 0.186, 0.643), (0.677, 0.305, 0.683), (0.788, 0.693, 0.644), (0.165, 0.745, 0.702)]
    corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for etype in ('C3D8R', 'C3D8', 'C3D8I', 'C3D4'):
        d = Deck('patch_solid_%s' % etype.lower(), '%s solid patch test (distorted elements)' % etype)
        c = [d.node(*[10.0 * v for v in p]) for p in corners]
        if etype == 'C3D4':
            ii = [d.node(4.0, 5.5, 4.5)]
            faces = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
            tets = []
            for f in faces:
                for tri in ((f[0], f[1], f[2]), (f[0], f[2], f[3])):
                    t = [c[tri[0]], c[tri[1]], c[tri[2]], ii[0]]
                    from deck import vol6
                    if vol6(d, t) < 0:
                        t[1], t[2] = t[2], t[1]
                    tets.append(t)
            solid(d, etype, tets, 'STEEL')
        else:
            ii = [d.node(*[10.0 * v for v in p]) for p in inner]
            pairs = [((0, 3, 2, 1), (0, 3, 2, 1)), ((4, 5, 6, 7), (4, 5, 6, 7)), ((0, 1, 5, 4), (0, 1, 5, 4)),
                     ((1, 2, 6, 5), (1, 2, 6, 5)), ((2, 3, 7, 6), (2, 3, 7, 6)), ((3, 0, 4, 7), (3, 0, 4, 7))]
            hexes = [hex_ok(d, [ii[k] for k in range(4)] + [ii[k] for k in range(4, 8)])]
            for outer, inn in pairs:
                hexes.append(hex_ok(d, [c[k] for k in outer] + [ii[k] for k in inn]))
            solid(d, etype, hexes, 'STEEL')
        body = ['*Boundary']
        for n in c:
            x, y, z = d.nodes[n]
            body += ['%d, 1, 1, %s' % (n, fmt(1e-3 * (2 * x + y + z) / 2)),
                     '%d, 2, 2, %s' % (n, fmt(1e-3 * (x + 2 * y + z) / 2)),
                     '%d, 3, 3, %s' % (n, fmt(1e-3 * (x + y + 2 * z) / 2))]
        d.steps += static_step('Patch', body, False, nframes=1)
        register(d, 'static', etype, 'solid patch test', 'Constant strain imposed on the corners of a cube meshed '
                 'with distorted elements: interior nodes and stresses must be exact.', ii, [], False, 'exact',
                 'linear', reference='analytic', category='feature', group='patch')

    for etype in ('S4R', 'S3R'):
        d = Deck('thin_shell_moment_%s' % etype.lower(), '%s very thin strip (t = 0.01 mm) under an end moment' % etype)
        quads, ids = quad_plate(d, 20, 2, 40.0, 4.0)
        shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 0.01, 'STEEL0')
        fix = [ids[0, j] for j in range(3)]
        tip = [ids[20, j] for j in range(3)]
        ei = 210000.0 * 4.0 * 0.01 ** 3 / 12
        m = 0.01 / 40.0 * ei
        d.steps += static_step('Moment', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 5, %s' % (n, fmt(-m * w)) for n, w in edge_shares(d, tip)], False, rot=True,
                               nframes=1)
        register(d, 'static', etype, 'thin shell, pure bending', 'Strip 4000 times thinner than long under a constant '
                 'moment: exact arc expected (no shear locking).', tip + [ids[10, 1]], fix, True, 'exact', 'linear',
                 reference='analytic', category='feature', group='patch')

    def cantilever(name, title, skew, mode, desc, enhanced=False):
        d = Deck(name, title)
        ids = {}
        for i, x in enumerate([0., 10., 20., 30., 40., 50., 60.]):
            ids[i, 0] = d.node(x - skew[i] / 2, 0.)
            ids[i, 1] = d.node(x + skew[i] / 2, 2.)
        quads = [[ids[i, 0], ids[i + 1, 0], ids[i + 1, 1], ids[i, 1]] for i in range(6)]
        shell(d, 'S4R', quads, 0.1, 'STEEL', extra=', controls=HGE' if enhanced else '')
        if enhanced:
            d.model += ['*Section Controls, name=HGE, hourglass=ENHANCED']
        d.nset('NALL', list(d.nodes))
        fix, tip = [ids[0, 0], ids[0, 1]], [ids[6, 0], ids[6, 1]]
        d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 2', '*Cload', '%d, 2, 0.5' % tip[0],
                                        '%d, 2, 0.5' % tip[1]], False, rot=True, nframes=1)
        register(d, 'static', 'S4R', mode, desc, tip, fix, True, 'linear', 'linear', category='feature',
                 group='patch', twin=False)

    trapezoids = [0., 4., -4., 4., -4., 4., 0.]
    cantilever('distorted_cantilever_s4r', 'S4R in-plane cantilever with trapezoidal elements (MacNeal-Harder beam)',
               trapezoids, 'distorted elements', 'In-plane cantilever meshed with trapezoidal elements, tip shear '
               '(MacNeal-Harder straight beam).')
    cantilever('straight_cantilever_s4r', 'S4R in-plane cantilever with rectangular elements (MacNeal-Harder beam)',
               [0.] * 7, 'in-plane bending, membrane hourglass', 'In-plane cantilever, one row of six rectangular '
               'elements, tip shear: in-plane bending carried by the membrane hourglass stiffness.')
    cantilever('hourglass_enhanced_s4r', 'S4R in-plane cantilever, rectangular elements, enhanced hourglass control',
               [0.] * 7, '*SECTION CONTROLS, HOURGLASS=ENHANCED', 'straight_cantilever_s4r with the enhanced '
               'hourglass control: hourglass stiffness of a physical element (beam theory expected).', enhanced=True)
    cantilever('hourglass_enhanced_distorted_s4r', 'S4R in-plane cantilever, trapezoidal elements, enhanced hourglass '
               'control', trapezoids, '*SECTION CONTROLS, HOURGLASS=ENHANCED', 'distorted_cantilever_s4r with the '
               'enhanced hourglass control.', enhanced=True)

    d = Deck('twisted_beam_s4r', 'S4R twisted beam (90 degree pretwist, MacNeal-Harder)')
    nx, ny = 12, 2
    ids = {}
    for i in range(nx + 1):
        a = math.radians(90.0 * i / nx)
        for j in range(ny + 1):
            s = -5.5 + 11.0 * j / ny
            ids[i, j] = d.node(120.0 * i / nx, s * math.cos(a), s * math.sin(a))
    quads = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(nx) for j in range(ny)]
    shell(d, 'S4R', quads, 3.2, 'TWIST')
    fix = [ids[0, j] for j in range(ny + 1)]
    tip = [ids[nx, j] for j in range(ny + 1)]
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 2, %s' % (n, fmt(1.0 * w)) for n, w in zip(tip, (0.25, 0.5, 0.25))], False, rot=True,
                           nframes=1)
    register(d, 'static', 'S4R', 'warped elements', 'Twisted beam: every element is warped; unit tip load in the '
             'direction of the width at the tip.', tip + [ids[6, 1]], fix, True, 'linear', 'linear',
             category='feature', group='patch', twin=False)

    d = Deck('spin_two_turns_s4r', 'S4R free plate spinning two full turns (explicit)')
    quads, ids = quad_plate(d, 6, 6, 20.0, 20.0, x0=(-10.0, -10.0, 0.0))
    shell(d, 'S4R', quads, 0.5, 'STEEL')
    w = 500.0
    d.initial += ['*Initial Conditions, type=VELOCITY']
    for n, (x, y, z) in sorted(d.nodes.items()):
        d.initial += ['%d, 1, %s' % (n, fmt(-w * y)), '%d, 2, %s' % (n, fmt(w * x)), '%d, 6, %s' % (n, fmt(w))]
    d.steps += explicit_step('Spin', 4 * math.pi / w, [], rot=True, nframes=24)
    register(d, 'explicit', 'S4R', 'rotation beyond 360 degrees', 'Free plate spinning at 500 rad/s for two full turns: '
             'positions must follow the rigid rotation exactly.', [ids[6, 6], ids[6, 3]], [], True, 'exact',
             reference='analytic', category='feature', group='patch')

    for etype in ('C3D8H', 'C3D4H', 'C3D8R'):
        d = Deck('incompressible_%s' % etype.lower(), '%s nearly incompressible rubber cantilever (K / mu = 2e4)' % etype)
        hexes = hex_block(d, 10, 2, 2, 20.0, 4.0, 4.0)
        from deck import tets_from_hexes
        solid(d, etype, tets_from_hexes(d, hexes) if etype.startswith('C3D4') else hexes,
              'RUBBER_INC' if etype.endswith('H') else 'RUBBER_495')
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, 20.0))
        QS_DURATION[d.name] = 0.05
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 3, %s' % (n, fmt(-1.0 / len(tip))) for n in tip], True)
        register(d, 'static', etype, 'near incompressibility', 'Rubber cantilever with a bulk modulus 20000 times the '
                 'shear modulus, tip load (volumetric locking).', tip[:3], fix, False, 'nonlinear', category='feature',
                 group='patch', supported=etype == 'C3D8R', reason=None if etype == 'C3D8R' else
                 NS + 'hybrid elements (%s)' % etype, twin=etype == 'C3D8R')

    for n in (4, 16):
        d = Deck('s4r_pressure_n%d' % n, 'S4R clamped square plate under pressure, %d x %d mesh (convergence)' % (n, n))
        quads, ids = quad_plate(d, n, n, 20.0, 20.0)
        shell(d, 'S4R', quads, 0.2, 'STEEL')
        d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
        edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
        QS_DURATION[d.name] = 0.002
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.5'], True, rot=True)
        register(d, 'static', 'S4R', 'mesh convergence', 'Same plate as s4r_pressure with a %d x %d mesh.' % (n, n),
                 [ids[n // 2, n // 2], ids[n // 4, n // 4]], edge, True, 'nonlinear', category='feature',
                 group='patch')


# ----------------------------------------------------------------------------------------------- benchmarks
def cylinder_nodes(d, nx, nt, length, radius, theta0, theta1):
    ids = {}
    for i in range(nx + 1):
        for j in range(nt + 1):
            th = theta0 + (theta1 - theta0) * j / nt
            ids[i, j] = d.node(length * i / nx, radius * math.sin(th), radius * math.cos(th))
    return ids


def nodal_area_loads(d, quads, fz):
    """Consistent nodal shares (quarter of each element area) of a uniform load per unit area along z."""
    acc = {}
    for q in quads:
        p = [d.nodes[n] for n in q]
        a = [p[2][i] - p[0][i] for i in range(3)]
        b = [p[3][i] - p[1][i] for i in range(3)]
        cr = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        area = 0.5 * math.sqrt(sum(c * c for c in cr))
        for n in q:
            acc[n] = acc.get(n, 0.0) + area / 4
    return ['%d, 3, %s' % (n, fmt(fz * a)) for n, a in sorted(acc.items())]


def benchmark_features():
    for etype in ('S4R', 'S3R'):
        d = Deck('scordelis_lo_%s' % etype.lower(), '%s Scordelis-Lo roof (quarter model, self-weight)' % etype)
        n = 8
        ids = cylinder_nodes(d, n, n, 25.0, 25.0, 0.0, math.radians(40.0))
        quads = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(n) for j in range(n)]
        shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 0.25, 'ROOF')
        d.nset('DIAPH', [ids[n, j] for j in range(n + 1)])
        d.nset('XS', [ids[0, j] for j in range(n + 1)])
        d.nset('YS', [ids[i, 0] for i in range(n + 1)])
        d.steps += static_step('Load', ['*Boundary', 'DIAPH, 2, 3', 'DIAPH, 4, 4', 'XS, XSYMM', 'YS, YSYMM', '*Cload'] +
                               nodal_area_loads(d, quads, -90.0), False, rot=True, nframes=1)
        register(d, 'static', etype, 'Scordelis-Lo roof', 'Cylindrical roof on rigid diaphragms under self-weight '
                 '(published converged value: 0.3024 at the middle of the free edge).',
                 [ids[0, n], ids[n // 2, n], ids[0, n // 2]], [ids[n, j] for j in range(n + 1)], True, 'linear',
                 'linear', category='feature', group='benchmark', twin=False)

    d = Deck('pinched_cylinder_s4r', 'S4R pinched cylinder with rigid end diaphragms (octant model)')
    n = 8
    ids = cylinder_nodes(d, n, n, 300.0, 300.0, 0.0, math.pi / 2)
    quads = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(n) for j in range(n)]
    shell(d, 'S4R', quads, 3.0, 'CYL')
    d.nset('DIAPH', [ids[n, j] for j in range(n + 1)])
    d.nset('XS', [ids[0, j] for j in range(n + 1)])
    d.nset('YS', [ids[i, 0] for i in range(n + 1)])
    d.nset('ZS', [ids[i, n] for i in range(n + 1)])
    d.steps += static_step('Load', ['*Boundary', 'DIAPH, 2, 3', 'DIAPH, 4, 4', 'XS, XSYMM', 'YS, YSYMM', 'ZS, ZSYMM',
                                    '*Cload', '%d, 3, -0.25' % ids[0, 0]], False, rot=True, nframes=1)
    register(d, 'static', 'S4R', 'pinched cylinder', 'Cylinder pinched by two opposite unit loads (published converged '
             'value: 1.8248e-5 under the load).', [ids[0, 0], ids[0, n], ids[n // 2, 0]],
             [ids[n, j] for j in range(n + 1)], True, 'linear', 'linear', category='feature', group='benchmark',
             twin=False)

    d = Deck('pinched_hemisphere_s4r', 'S4R pinched hemisphere with an 18 degree hole (quarter model)')
    n = 8
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
    register(d, 'static', 'S4R', 'pinched hemisphere', 'Hemispherical shell loaded by opposite radial forces at the '
             'equator (published converged value: 0.093 in the load direction).', [ids[n, 0], ids[n, n], ids[0, 0]],
             [ids[n, 0]], True, 'linear', 'linear', category='feature', group='benchmark', twin=False)

    for etype in ('S4R', 'S3R'):
        d = Deck('cook_membrane_%s' % etype.lower(), "%s Cook's membrane (tapered panel, shear load)" % etype)
        n = 8
        ids = {}
        for i in range(n + 1):
            x = 48.0 * i / n
            y0, y1 = 44.0 * x / 48.0, 44.0 + 16.0 * x / 48.0
            for j in range(n + 1):
                ids[i, j] = d.node(x, y0 + (y1 - y0) * j / n)
        quads = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(n) for j in range(n)]
        shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 1.0, 'COOK')
        d.nset('NALL', list(d.nodes))
        fix = [ids[0, j] for j in range(n + 1)]
        right = [ids[n, j] for j in range(n + 1)]
        d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 2', '*Cload'] +
                               ['%d, 2, %s' % (nd, fmt(1.0 * (0.5 if k in (0, n) else 1.0) / n))
                                for k, nd in enumerate(right)], False, rot=True, nframes=1)
        register(d, 'static', etype, "Cook's membrane", 'Tapered membrane clamped on one edge, unit shear on the '
                 'other (published converged value: 23.96 at the middle of the loaded edge).', [ids[n, n], ids[n, 0]], fix, True,
                 'linear', 'linear', category='feature', group='benchmark', twin=False)

    def ellipse_mesh(d, nt, ns, z=0.0):
        ids = {}
        for i in range(nt + 1):
            th = math.pi / 2 * i / nt
            for j in range(ns + 1):
                s = j / ns
                ids[i, j] = d.node((2.0 + 1.25 * s) * math.cos(th), (1.0 + 1.75 * s) * math.sin(th), z)
        return ids

    for etype in ('S4R', 'S3R'):
        d = Deck('nafems_le1_%s' % etype.lower(), '%s NAFEMS LE1 elliptic membrane under outer tension' % etype)
        nt, ns = 12, 4
        ids = ellipse_mesh(d, nt, ns)
        quads = [[ids[i, j], ids[i, j + 1], ids[i + 1, j + 1], ids[i + 1, j]] for i in range(nt) for j in range(ns)]
        shell(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, 0.1, 'STEEL')
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
        for nd, (fx, fy) in sorted(acc.items()):
            loads += ['%d, 1, %s' % (nd, fmt(-fx)), '%d, 2, %s' % (nd, fmt(-fy))]
        d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'XS, XSYMM', 'YS, YSYMM', '*Cload'] + loads, False,
                               rot=True, nframes=1)
        register(d, 'static', etype, 'NAFEMS LE1', 'Elliptic membrane under a uniform normal load of 10 MPa on its '
                 'outer edge (geometry of NAFEMS LE1, load reversed; see lit_nafems_le1_* for the published value).', [ids[0, 0], ids[0, ns], ids[nt, 0]],
                 [ids[nt, j] for j in range(ns + 1)] + [ids[0, j] for j in range(ns + 1)], True, 'linear', 'linear',
                 category='feature', group='benchmark', twin=False)

    for etype in ('C3D8R', 'C3D4'):
        d = Deck('nafems_le10_%s' % etype.lower(), '%s NAFEMS LE10 thick elliptic plate under pressure' % etype)
        nt, ns, nz = 8, 4, 4
        layers = [ellipse_mesh(d, nt, ns, -0.3 + 0.6 * k / nz) for k in range(nz + 1)]
        hexes = []
        for k in range(nz):
            for i in range(nt):
                for j in range(ns):
                    lo, hi = layers[k], layers[k + 1]
                    hexes.append(hex_ok(d, [lo[i, j], lo[i, j + 1], lo[i + 1, j + 1], lo[i + 1, j],
                                            hi[i, j], hi[i, j + 1], hi[i + 1, j + 1], hi[i + 1, j]]))
        from deck import tets_from_hexes
        conns = tets_from_hexes(d, hexes) if etype == 'C3D4' else hexes
        el = solid(d, etype, conns, 'STEEL')
        top = faces_on(d, etype, el, conns, lambda x, y, z: near(z, 0.3))
        d.model += ['*Surface, type=ELEMENT, name=STOP'] + ['%d, S%d' % f for f in top]
        d.nset('XS', d.select(lambda x, y, z: near(x, 0)))
        d.nset('YS', d.select(lambda x, y, z: near(y, 0)))
        outer = [l[i, ns] for l in layers for i in range(nt + 1)]
        d.nset('OUTER', outer)
        d.nset('OUTMID', [n for n in outer if near(d.nodes[n][2], 0.0)])
        d.steps += static_step('Load', ['*Boundary', 'XS, 1, 1', 'YS, 2, 2', 'OUTER, 1, 2', 'OUTMID, 3, 3',
                                        '*Dsload', 'STOP, P, 1.'], False, nframes=1)
        register(d, 'static', etype, 'NAFEMS LE10', 'Thick elliptic plate under 1 MPa pressure (published value: '
                 'sigma_yy = -5.38 MPa at point D, top face, inner edge on the x axis).',
                 [layers[nz][0, 0], layers[nz][nt, 0]], outer, False, 'linear', 'linear', category='feature',
                 group='benchmark', twin=False)

    d = Deck('taylor_impact_c3d8r', 'Taylor-type impact of a copper bar on a rigid wall (explicit, quarter model)')
    hexes = hex_block(d, 3, 3, 20, 3.2, 3.2, 32.4)
    solid(d, 'C3D8R', hexes, 'COPPER')
    d.nset('WALL', d.select(lambda x, y, z: near(z, 0)))
    d.nset('XS', d.select(lambda x, y, z: near(x, 0)))
    d.nset('YS', d.select(lambda x, y, z: near(y, 0)))
    moving = d.select(lambda x, y, z: z > 1e-9)
    d.nset('MOVING', moving)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'MOVING, 3, -227000.']
    d.steps += explicit_step('Impact', 8e-05, ['*Boundary', 'WALL, 3, 3', 'XS, XSYMM', 'YS, YSYMM'], time_marks=True)
    register(d, 'explicit', 'C3D8R', 'Taylor impact', 'Copper bar hitting a rigid frictionless wall at 227 m/s: '
             'mushrooming, large plastic strains (explicit dynamics).',
             d.select(lambda x, y, z: near(z, 0) and near(y, 0) and near(x, 3.2)) +
             d.select(lambda x, y, z: near(z, 32.4) and near(x, 0) and near(y, 0)), d.nsets['WALL'], False, 'dynamic',
             category='feature', group='benchmark', twin=False)


# ---------------------------------------------------------------------------------------- input-file features
class AssemblyDeck(Deck):
    """Two parts, each instanced once (translation, rotation), sets defined in the assembly."""

    def text(self):
        out = ['*Heading', self.title]
        for part, (nodes, elems, mat) in self.parts.items():
            out += ['*Part, name=%s' % part, '*Node']
            out += ['%d, %s, %s, %s' % (i, fmt(p[0]), fmt(p[1]), fmt(p[2])) for i, p in sorted(nodes.items())]
            out += ['*Element, type=C3D8R, elset=EALL'] + ['%d, %s' % (e, ', '.join(map(str, c))) for e, c in elems]
            out += ['*Solid Section, elset=EALL, material=%s' % mat, ',', '*End Part']
        out += ['*Assembly, name=Assembly']
        for inst, (part, trans, rot) in self.instances.items():
            out += ['*Instance, name=%s, part=%s' % (inst, part), '%s, %s, %s' % tuple(fmt(v) for v in trans)]
            if rot:
                out.append(', '.join(fmt(v) for v in rot))
            out.append('*End Instance')
        for name, (inst, ids) in self.asets.items():
            out.append('*Nset, nset=%s, instance=%s' % (name, inst))
            out += [', '.join(str(i) for i in ids[k:k + 16]) for k in range(0, len(ids), 16)]
        out.append('*End Assembly')
        for m in sorted(self.materials):
            out.append(MATERIALS[m].rstrip('\n'))
        out += self.steps
        return '\n'.join(out) + '\n'


def input_features():
    d = AssemblyDeck('input_parts_instances', 'Two parts instanced with a translation and a rotation (tension)')
    d.parts, d.instances, d.asets = {}, {}, {}
    mon, fix = [], []
    for part, base, inst, trans, rot in (('PA', 0, 'IA', (0., 10., 0.), None),
                                         ('PB', 1000, 'IB', (20., 0., 0.), (0., 0., 0., 0., 0., 1., 90.))):
        nodes, elems, idx = {}, [], {}
        for k in range(2):
            for j in range(2):
                for i in range(5):
                    lab = base + 1 + i + 5 * (j + 2 * k)
                    nodes[lab] = (10.0 * i, 4.0 * j, 4.0 * k)
                    idx[i, j, k] = lab
        for i in range(4):
            elems.append((base + 1 + i, [idx[i, 0, 0], idx[i + 1, 0, 0], idx[i + 1, 1, 0], idx[i, 1, 0],
                                         idx[i, 0, 1], idx[i + 1, 0, 1], idx[i + 1, 1, 1], idx[i, 1, 1]]))
        d.parts[part] = (nodes, elems, 'STEEL')
        d.instances[inst] = (part, trans, rot)
        d.asets[inst + 'X0'] = (inst, [idx[0, j, k] for j in range(2) for k in range(2)])
        d.asets[inst + 'Y0'] = (inst, [idx[i, 0, k] for i in range(5) for k in range(2)])
        d.asets[inst + 'Z0'] = (inst, [idx[i, j, 0] for i in range(5) for j in range(2)])
        d.asets[inst + 'TIP'] = (inst, [idx[4, j, k] for j in range(2) for k in range(2)])
        mon += [idx[4, 1, 1], idx[2, 1, 1]]
        fix += [idx[0, j, k] for j in range(2) for k in range(2)]
    d.materials.add('STEEL')
    d.steps += static_step('Load', ['*Boundary', 'IAX0, 1, 1', 'IAY0, 2, 2', 'IAZ0, 3, 3', 'IBX0, 2, 2', 'IBY0, 1, 1',
                                    'IBZ0, 3, 3', '*Boundary', 'IATIP, 1, 1, 0.04', 'IBTIP, 2, 2, 0.04'], False,
                           nframes=1)
    register(d, 'static', 'C3D8R', 'parts and instances', 'Two parts, one instance translated, the other rotated by 90 '
             'degrees; homogeneous tension of each bar along its own axis, boundary conditions on assembly sets.', mon,
             fix, False, 'exact', 'linear', reference='analytic', category='feature', group='input')

    d = Deck('input_include_parameter', 'Bar in tension: mesh in an included file, values from *PARAMETER')
    hexes = hex_block(d, 4, 1, 1, 40.0, 4.0, 4.0)
    d.elements('C3D8R', 'EALL', hexes)
    fx = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 40.0))
    d.include_mesh = True
    d.model += ['*Solid Section, elset=EALL, material=PARAM_STEEL', ',',
                '*Material, name=PARAM_STEEL', '*Density', '7.85e-09', '*Elastic', '<YOUNG>, <NU>']
    d.steps += static_step('Load', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary', 'TIP, 1, 1, <TIPU>'],
                           False, nframes=1)
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.nset('TIP', tip)
    register(d, 'static', 'C3D8R', '*INCLUDE, *PARAMETER', 'Tension bar whose nodes, elements and sets are read from an '
             'included file, with modulus and displacement given by parameters.', tip, fx, False, 'exact', 'linear',
             reference='analytic', category='feature', group='input')

    d = Deck('contact_selected_surfaces', 'S4R sheet dropped on a plate, general contact restricted to two surfaces')
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    shell(d, 'S4R', quads, 0.5, 'STEEL', 'PLATE')
    squads, sids = quad_plate(d, 3, 3, 6.0, 6.0, x0=(7.0, 7.0, 1.0))
    shell(d, 'S4R', squads, 0.5, 'ALU', 'SHEET')
    edge = [ids[i, j] for (i, j) in ids if i in (0, 8) or j in (0, 8)]
    d.model += ['*Surface, type=ELEMENT, name=SPLATE', 'PLATE, SPOS', 'PLATE, SNEG',
                '*Surface, type=ELEMENT, name=SSHEET', 'SHEET, SPOS', 'SHEET, SNEG', '*Surface Interaction, name=IP']
    d.nset('SHEETN', list(sids.values()))
    d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE', '*Contact',
                                              '*Contact Inclusions', 'SSHEET, SPLATE',
                                              '*Contact Property Assignment', ' , , IP'], rot=True, nframes=200,
                             time_marks=True)
    register(d, 'explicit', 'S4R', 'contact on selected surfaces', 'Sheet dropped on a plate; general contact limited '
             'to the pair of surfaces instead of all exterior faces.', [ids[4, 4], sids[1, 1]], edge, True, 'contact',
             category='feature', group='input', twin=False)

    d = Deck('section_controls_hourglass', 'S4R square plate twisted by a corner load, hourglass stiffness x 0.1')
    quads, ids = quad_plate(d, 6, 6, 20.0, 20.0)
    shell(d, 'S4R', quads, 0.5, 'STEEL', extra=', controls=HG')
    d.model += ['*Section Controls, name=HG, hourglass=STIFFNESS', '0.1, 0.1, 0.1']
    c00, c10, c01, c11 = ids[0, 0], ids[6, 0], ids[0, 6], ids[6, 6]
    d.steps += static_step('Load', ['*Boundary', '%d, 1, 3' % c00, '%d, 2, 3' % c10, '%d, 3, 3' % c01,
                                    '*Cload', '%d, 3, -2.' % c11], False, rot=True, nframes=1)
    register(d, 'static', 'S4R', '*SECTION CONTROLS', 'Twisting test of s4r_twist with the hourglass stiffness scaled '
             'by 0.1.', [c11, ids[3, 3]], [c00, c10, c01], True, 'linear', 'linear', category='feature', group='input',
             twin=False)

    # *NSET, ELSET=: the nodes of an element set (here the first column of elements, clamped)
    d = Deck('nset_from_elset', 'S4R cantilever plate clamped through a node set built from an element set')
    quads, ids = quad_plate(d, 8, 2, 20.0, 5.0)
    shell(d, 'S4R', quads, 0.5, 'STEEL')
    root = [k + 1 for k, q in enumerate(quads) if min(d.nodes[n][0] for n in q) < 1e-9]
    d.model += ['*Elset, elset=ROOT'] + [', '.join(map(str, root))] + ['*Nset, nset=CLAMP, elset=ROOT']
    fix = sorted({n for k in root for n in quads[k - 1]})
    tip = [ids[8, j] for j in range(3)]
    d.steps += static_step('Load', ['*Boundary', 'CLAMP, ENCASTRE', '*Cload'] +
                           ['%d, 3, %s' % (n, fmt(-1.0 * w)) for n, w in zip(tip, (0.25, 0.5, 0.25))], False,
                           rot=True, nframes=1)
    register(d, 'static', 'S4R', '*NSET, ELSET=', 'Cantilever plate clamped on the nodes of its first column of '
             'elements, given as *NSET, ELSET= (node set built from an element set), tip load.', tip, fix, True,
             'linear', 'linear', category='feature', group='input', twin=False)

    d = Deck('transverse_shear_thick', 'Thick S4R cantilever with user transverse shear stiffness')
    quads, ids = quad_plate(d, 10, 2, 20.0, 4.0)
    shell(d, 'S4R', quads, 4.0, 'STEEL')
    d.model += ['*Transverse Shear', '20000., 20000., 0.']
    fix = [ids[0, j] for j in range(3)]
    tip = [ids[10, j] for j in range(3)]
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 3, %s' % (n, fmt(-100.0 * w)) for n, w in zip(tip, (0.25, 0.5, 0.25))], False,
                           rot=True, nframes=1)
    register(d, 'static', 'S4R', '*TRANSVERSE SHEAR', 'Short thick cantilever whose deflection is dominated by the '
             'user-defined transverse shear stiffness.', tip + [ids[5, 1]], fix, True, 'linear', 'linear',
             category='feature', group='input', twin=False)

    for npt in (3, 9):
        d = Deck('section_points_%d' % npt, 'S4R plastic bending, %d section points (explicit)' % npt)
        quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
        d.elements('S4R', 'EALL', quads)
        d.materials.add('STEEL_PL')
        d.model += ['*Shell Section, elset=EALL, material=STEEL_PL', '1., %d' % npt]
        fix = [ids[0, j] for j in range(3)]
        tip = [ids[10, j] for j in range(3)]
        d.nset('TIP', tip)
        d.model += ['*Amplitude, name=BEND, definition=SMOOTH STEP', '0., 0., 0.02, 1.']
        d.steps += explicit_step('Bend', 0.025, ['*Boundary', 'FIX, ENCASTRE', '*Boundary, amplitude=BEND',
                                                 'TIP, 3, 3, -8.'], rot=True)
        register(d, 'explicit', 'S4R', 'section points (%d)' % npt, 'Plastic bending of a plate strip with %d '
                 'integration points through the thickness.' % npt, tip[:2] + [ids[5, 1]], fix, True, 'dynamic',
                 category='feature', group='input', twin=False)


# ----------------------------------------------------------------------------------------- more features
def cantilever_solid(name, etype, title, mat='STEEL', conns_fn=None):
    d = Deck(name, title)
    hexes = hex_block(d, 10, 2, 2, 40.0, 4.0, 4.0)
    solid(d, etype, conns_fn(d, hexes) if conns_fn else hexes, mat)
    fix = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 40.0) and (near(y, 0) or near(y, 4.0)) and (near(z, 0) or near(z, 4.0)))
    return d, fix, tip


def uniaxial_block(name, mat, title, tip_u, steps=None, nlgeom=True, nframes=10):
    d = Deck(name, title)
    hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
    solid(d, 'C3D8R', hexes, mat)
    fx = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 4.0))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.nset('TIP', tip)
    d.steps += static_step('Load', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary',
                                    'TIP, 1, 1, %s' % fmt(tip_u)], nlgeom, nframes=nframes,
                           inc=0.02 if nlgeom else 0.1, maxinc=0.02 if nlgeom else None)
    return d, fx, tip


def more_features():
    from gen_features import tet10  # noqa
    for etype in ('C3D8I', 'C3D10'):
        if etype == 'C3D10':
            from deck import tets_from_hexes
            d, fix, tip = cantilever_solid('elem_c3d10', 'C3D10', 'C3D10 cantilever, tip load',
                                           conns_fn=lambda d, h: tet10(d, tets_from_hexes(d, h)))
        else:
            d, fix, tip = cantilever_solid('elem_c3d8i', 'C3D8I', 'C3D8I (incompatible modes) cantilever, tip load')
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 3, -50.' % n for n in tip], False, nframes=1)
        register(d, 'static', etype, 'element %s' % etype, 'Cantilever under a 200 N tip load.', tip, fix, False,
                 'linear', 'linear', category='feature', group='element', supported=etype == 'C3D8I',
                 reason=None if etype == 'C3D8I' else NS + '%s elements' % etype, twin=False)

    for sec, dims, title in (('BOX', '4., 2., 0.2, 0.2, 0.2, 0.2', 'box'), ('I', '1.5, 3., 2., 2., 0.2, 0.2, 0.2', 'I')):
        d = Deck('beam_section_%s' % title.lower(), 'B31 cantilever with a %s section, tip loads' % title)
        segs, ids = beam_line(d, 20, 100.0)
        beam(d, segs, sec, dims)
        d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload', '%d, 2, 1.' % ids[-1],
                                        '%d, 3, -2.' % ids[-1], '%d, 4, 50.' % ids[-1]], False, rot=True, nframes=1)
        register(d, 'static', 'B31', '%s section' % title, 'Thin-walled %s section cantilever under bending in both '
                 'planes and torsion.' % title, [ids[-1], ids[10]], [ids[0]], True, 'linear', 'linear',
                 category='feature', group='element', supported=False, reason=NS + '%s beam sections' % sec,
                 twin=False)

    # homogeneous material tests (closed form)
    for name, mat, title, u in (('mat_ogden', 'OGDEN', 'Ogden (N = 2) block, uniaxial stretch 1.5', 2.0),
                                ('mat_yeoh', 'YEOH', 'Yeoh block, uniaxial stretch 1.5', 2.0)):
        d, fx, tip = uniaxial_block(name, mat, title, u)
        register(d, 'static', 'C3D8R', title.split(' block')[0], 'Homogeneous uniaxial stretch to 1.5.', tip[:3], fx,
                 False, 'exact', reference='analytic', category='feature', group='material', supported=True)
    d, fx, tip = uniaxial_block('mat_orthotropic', 'ORTHO', 'Orthotropic block (engineering constants), tension 0.1 %',
                                0.004, nlgeom=False, nframes=1)
    d.model = ['*Orientation, name=OR', '1., 0., 0., 0., 1., 0.'] + [
        l + ', orientation=OR' if l.startswith('*Solid Section') else l for l in d.model]
    register(d, 'static', 'C3D8R', 'orthotropic elasticity', 'Homogeneous tension along the stiff direction of an '
             'orthotropic material.', tip[:3], fx, False, 'exact', 'linear', reference='analytic', category='feature',
             group='material', supported=True)

    d = Deck('mat_hgo', 'Holzapfel-Gasser-Ogden block with two fibre families at +/- 30 degrees, uniaxial stretch')
    hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
    d.elements('C3D8R', 'EALL', hexes)
    d.materials.add('HGO')
    d.model += ['*Orientation, name=ORI, local directions=2', '1., 0., 0., 0., 1., 0.', '3, 0.',
                '0.866025404, 0.5, 0.', '0.866025404, -0.5, 0.',
                '*Solid Section, elset=EALL, material=HGO, orientation=ORI', ',']
    fx = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 4.0))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.nset('TIP', tip)
    QS_DURATION['mat_hgo'] = 0.01
    d.steps += static_step('Load', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary', 'TIP, 1, 1, 1.2'],
                           True)
    register(d, 'static', 'C3D8R', 'Holzapfel-Gasser-Ogden', 'Fibre-reinforced (arterial) hyperelastic block, two '
             'fibre families, uniaxial stretch to 1.3.', tip[:3], fx, False, 'nonlinear', category='feature',
             group='material', supported=True)

    d = Deck('mat_superelastic', 'Superelastic NiTi bar (ABQ_SUPER_ELASTIC), loading and unloading through transformation')
    hexes = hex_block(d, 4, 1, 1, 20.0, 2.0, 2.0)
    solid(d, 'C3D8', hexes, 'ABQ_SUPER_ELASTIC_N3D')
    fx = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 20.0))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.nset('TIP', tip); d.nset('NALL', list(d.nodes))
    d.initial += ['*Initial Conditions, type=TEMPERATURE', 'NALL, 37.']
    QS_DURATION['mat_superelastic'] = 0.002
    d.steps += static_step('Load', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary', 'TIP, 1, 1, 1.4'],
                           True, inc=0.02, maxinc=0.02)
    d.steps += static_step('Unload', ['*Boundary', 'TIP, 1, 1, 0.'], True, inc=0.02, maxinc=0.02)
    register(d, 'static', 'C3D8', 'superelastic nitinol', 'NiTi bar stretched to 7 % (austenite to martensite '
             'transformation) and unloaded (reverse transformation, hysteresis).', tip[:3], fx, False, 'nonlinear',
             category='feature', group='material', supported=True)

    for name, mat, title, reason in (('mat_rate_dependent', 'STEEL_RATE', 'rate-dependent plasticity (power law)',
                                      'rate-dependent plasticity'),
                                     ('mat_damage_deletion', 'STEEL_DMG', 'ductile damage and element deletion',
                                      'damage and element deletion')):
        d = Deck(name, 'Steel bar pulled at high speed, %s (explicit)' % title)
        hexes = hex_block(d, 10, 2, 2, 20.0, 2.0, 2.0)
        solid(d, 'C3D8R', hexes, mat)
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, 20.0))
        d.nset('TIP', tip)
        d.model += ['*Amplitude, name=PULL, definition=SMOOTH STEP', '0., 0., 0.0002, 1.']
        d.steps += explicit_step('Pull', 0.0004, ['*Boundary', 'FIX, ENCASTRE', '*Boundary, type=VELOCITY, '
                                                  'amplitude=PULL', 'TIP, 1, 1, 10000.', 'TIP, 2, 3, 0.'])
        register(d, 'explicit', 'C3D8R', title, 'Clamped bar pulled at 10 m/s, %s.' % title, tip[:2] +
                 d.select(lambda x, y, z: near(x, 10.0) and near(y, 0) and near(z, 0)), fix, False, 'dynamic',
                 category='feature', group='material', supported=False, reason=NS + reason, twin=False)

    # constraints, connectors and options
    d = Deck('connector_nonlinear', 'CONN3D2 AXIAL spring with a nonlinear force-displacement law')
    a, b = d.node(0, 0, 0), d.node(10, 0, 0)
    d.elements('CONN3D2', 'CONN', [[a, b]])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL', '*Connector Behavior, name=SPRING',
                '*Connector Elasticity, nonlinear, component=1', '-100., -1.', '0., 0.', '100., 0.5', '300., 1.']
    d.steps += static_step('Load', ['*Boundary', '%d, 1, 3' % a, '%d, 2, 3' % b, '*Cload', '%d, 1, 200.' % b], False,
                           nframes=4, inc=0.05, maxinc=0.05)
    register(d, 'static', 'CONN3D2', 'nonlinear connector', 'Axial spring with a piecewise linear force-displacement '
             'law, loaded beyond its kink.', [b], [a], False, 'exact', reference='analytic', category='feature',
             group='constraint', twin=False)

    for name, kw, data, title, reason in (
            ('connector_damping', '*Connector Damping, component=1', '0.05', 'damping', 'connector damping'),
            ('connector_stop', '*Connector Stop, component=1', '9.95, 10.05', 'stop', 'connector stops')):
        d = Deck(name, 'Mass on an axial connector with %s (explicit)' % title)
        hexes = hex_block(d, 1, 1, 1, 2.0, 2.0, 2.0, x0=(10.0, -1.0, -1.0))
        solid(d, 'C3D8R', hexes, 'HEAVY', 'MASS')
        face = d.select(lambda x, y, z: near(x, 10.0))
        ground = [d.node(0.0, d.nodes[n][1], d.nodes[n][2]) for n in face]
        d.elements('CONN3D2', 'CONN', [[g, n] for g, n in zip(ground, face)])
        d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL', '*Connector Behavior, name=SPRING',
                    '*Connector Elasticity, component=1', '25.', kw, data]
        cube = [n for n in d.nodes if n not in set(ground)]
        d.nset('CUBE', cube); d.nset('GROUND', ground)
        d.initial += ['*Initial Conditions, type=VELOCITY', 'CUBE, 1, 100.']
        d.steps += explicit_step('Oscillate', 0.01, ['*Boundary', 'GROUND, 1, 3', 'CUBE, 2, 3'])
        register(d, 'explicit', 'CONN3D2', 'connector %s' % title, 'Block on four axial connectors with %s, released '
                 'with 100 mm/s.' % title, cube[:2], ground, False, 'dynamic', category='feature', group='constraint',
                 supported=name != 'connector_stop', reason=None if name != 'connector_stop' else NS + reason, twin=False)

    for name, ctype, title in (('connector_beam', 'BEAM', 'rigid BEAM connector'),
                               ('connector_bushing', 'BUSHING', 'BUSHING connector')):
        d = Deck(name, 'Two B31 cantilevers joined by a %s' % title)
        s1, b1 = beam_line(d, 10, 20.0)
        s2, b2 = beam_line(d, 10, 20.0, x0=(21.0, 0.0, 0.0))
        beam(d, s1 + s2, 'RECT', '2., 1.')
        d.elements('CONN3D2', 'JOINT', [[b1[-1], b2[0]]])
        if ctype == 'BEAM':
            d.model += ['*Connector Section, elset=JOINT', 'BEAM']
        else:
            d.model += ['*Orientation, name=BO', '1., 0., 0., 0., 1., 0.', '*Connector Section, elset=JOINT, '
                        'behavior=BUSH', 'BUSHING', 'BO', '*Connector Behavior, name=BUSH']
            for comp, k in ((1, 1000.), (2, 500.), (3, 500.), (4, 1e5), (5, 1e5), (6, 1e5)):
                d.model += ['*Connector Elasticity, component=%d' % comp, fmt(k)]
        d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % b1[0], '*Cload', '%d, 3, -1.' % b2[-1],
                                        '%d, 2, 0.5' % b2[-1]], False, rot=True, nframes=1)
        register(d, 'static', 'B31 + CONN3D2', title, 'Two cantilever segments joined across a 1 mm gap by a %s.'
                 % title, [b2[-1], b1[-1], b2[0]], [b1[0]], True, 'linear', 'linear', category='feature',
                 group='constraint', supported=True, twin=False)

    d = Deck('spring_nonlinear', 'Point mass on a nonlinear SPRINGA spring (explicit)')
    g, m = d.node(0.0, 0.0, 0.0), d.node(10.0, 0.0, 0.0)
    d.elements('SPRINGA', 'ESP', [[g, m]])
    d.elements('MASS', 'EM', [[m], [g]])
    d.model += ['*Spring, elset=ESP, nonlinear', '', '-200., -1.', '0., 0.', '100., 1.', '400., 2.', '*Mass, elset=EM',
                '1e-04']
    d.initial += ['*Initial Conditions, type=VELOCITY', '%d, 1, 150.' % m]
    d.steps += explicit_step('Oscillate', 0.03, ['*Boundary', '%d, 1, 3' % g, '%d, 2, 3' % m], fixed_dt=1e-05)
    register(d, 'explicit', 'SPRINGA + MASS', 'nonlinear spring', 'Point mass on a spring stiffer in tension than in '
             'compression, released with 150 mm/s.', [m], [g], False, 'dynamic', category='feature', group='element',
             supported=True, twin=False)

    # contact options
    def sheet_on_plate(name, title, step_extra=None, model_extra=None, gap=1.0, sheet_offset=None):
        d = Deck(name, title)
        quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
        shell(d, 'S4R', quads, 0.5, 'STEEL', 'PLATE')
        squads, sids = quad_plate(d, 3, 3, 6.0, 6.0, x0=(7.0, 7.0, gap))
        shell(d, 'S4R', squads, 0.5, 'ALU', 'SHEET')
        edge = [ids[i, j] for (i, j) in ids if i in (0, 8) or j in (0, 8)]
        d.model += ['*Surface, type=ELEMENT, name=SPLATE', 'PLATE, SPOS', '*Surface, type=ELEMENT, name=SSHEET',
                    'SHEET, SNEG']
        d.nset('SHEETN', list(sids.values()))
        d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -2000.']
        body = ['*Boundary', 'FIX, ENCASTRE'] + general_contact(d, extra_model=model_extra) + (step_extra or [])
        d.steps += explicit_step('Impact', 1e-3, body, rot=True, nframes=200, time_marks=True)
        return d, ids, sids, edge
    d, ids, sids, edge = sheet_on_plate('contact_initial_overclosure', 'S4R sheet initially overclosing a plate by '
                                        '0.1 mm (default treatment of initial overclosures)', gap=0.4)
    register(d, 'explicit', 'S4R', 'initial overclosure', 'Sheet initially penetrating the plate by 0.1 mm (thickness '
             'offsets), then dropped: default treatment of initial overclosures (Abaqus: strain-free adjustment).',
             [ids[4, 4], sids[1, 1]], edge, True, 'contact', category='feature', group='contact', twin=False)
    d, ids, sids, edge = sheet_on_plate('contact_exclusions', 'S4R sheet falling on a plate, pair excluded from '
                                        'general contact', step_extra=['*Contact Exclusions', 'SSHEET, SPLATE'])
    register(d, 'explicit', 'S4R', 'contact exclusions', 'Sheet dropped on a plate with the sheet-plate pair excluded '
             'from general contact: the sheet passes through.', [ids[4, 4], sids[1, 1]], edge, True, 'contact',
             category='feature', group='contact', supported=True, twin=False)
    d, ids, sids, edge = sheet_on_plate('contact_softened', 'S4R sheet dropped on a plate, exponential '
                                        'pressure-overclosure law', model_extra=['*Surface Behavior, '
                                                                                 'pressure-overclosure=TABULAR',
                                                                                 '0., 0.', '5., 0.1', '500., 0.2'])
    register(d, 'explicit', 'S4R', 'softened contact', 'Sheet dropped on a plate with a tabular softened '
             'pressure-overclosure law (pressure grows from zero over the first 0.1 mm of overclosure).', [ids[4, 4], sids[1, 1]],
             edge, True, 'contact', category='feature', group='contact', supported=True, twin=False)

    d = Deck('contact_softened_single', 'S4R sheet dropped on one large clamped S4R element, linear softened '
             'pressure-overclosure law (no reverse contact)')
    pq, pids = quad_plate(d, 1, 1, 1.0, 1.0)
    shell(d, 'S4R', pq, 0.2, 'STEEL', 'PLATE')
    sq, sids = quad_plate(d, 1, 1, 0.4, 0.4, x0=(0.31, 0.33, 0.25))
    shell(d, 'S4R', sq, 0.2, 'ALU', 'SHEET')
    d.nset('SHEETN', list(sids.values()))
    d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -1000.']
    d.steps += explicit_step('Drop', 1e-4, ['*Boundary', 'FIX, ENCASTRE'] +
                             general_contact(d, extra_model=['*Surface Behavior, pressure-overclosure=LINEAR', '1000.']),
                             rot=True, nframes=200, time_marks=True)
    register(d, 'explicit', 'S4R', 'softened contact, single pass', 'Small sheet (one element) dropped at 1 m/s on '
             'one large clamped element with a linear softened law: no plate node lies under the sheet, so only the '
             'sheet nodes are slaves and carry the full pressure of their tributary area.', [sids[0, 0], sids[1, 1]],
             list(pids.values()), True, 'contact', category='feature', group='contact', twin=False)

    d = Deck('contact_analytical_cylinder', 'S4R strip pressed onto an analytical rigid cylinder (explicit)')
    quads, ids = quad_plate(d, 16, 2, 32.0, 4.0, x0=(-16.0, -2.0, 5.5))
    shell(d, 'S4R', quads, 0.5, 'STEEL')
    rp = d.node(0.0, 0.0, 0.0)
    # generator along Y (third point): the profile plane is (X, Z), local y = x cross generator = Z; the profile runs
    # clockwise, so that its normal (to the left) points outwards, towards the strip
    d.model += ['*Surface, type=CYLINDER, name=RCYL', '0., 0., 0., 1., 0., 0.', '0., 1., 0.', 'START, -5., 0.',
                'CIRCL, 0., 5., 0., 0.', 'CIRCL, 5., 0., 0., 0.', '*Rigid Body, ref node=%d, analytical surface=RCYL' % rp]
    ends = [ids[0, j] for j in range(3)] + [ids[16, j] for j in range(3)]
    d.nset('ENDS', ends)
    d.model += ['*Amplitude, name=PRESS, definition=SMOOTH STEP', '0., 0., 0.004, 1.']
    d.steps += explicit_step('Press', 0.005, ['*Boundary', '%d, ENCASTRE' % rp, '*Boundary, amplitude=PRESS',
                                              'ENDS, 3, 3, -3.', 'ENDS, 2, 2, 0.'] + general_contact(d), rot=True)
    register(d, 'explicit', 'S4R + analytical rigid surface', 'analytical rigid surface', 'Strip whose ends are pulled '
             'down, wrapping it around a rigid cylinder of radius 5.', [ids[8, 1], ids[4, 1]], [rp], True, 'contact',
             category='feature', group='contact', supported=True, twin=False)

    def tube(d, nt=16, nz=6, r=5.0, length=10.0, elset='TUBE', mat='STEEL_PL', t=0.2):
        ids = {}
        for i in range(nt):
            a = 2 * math.pi * i / nt
            for k in range(nz + 1):
                ids[i, k] = d.node(r * math.cos(a), r * math.sin(a), length * k / nz)
        quads = [[ids[i, k], ids[(i + 1) % nt, k], ids[(i + 1) % nt, k + 1], ids[i, k + 1]]
                 for i in range(nt) for k in range(nz)]
        shell(d, 'S4R', quads, t, mat, elset)
        return ids

    d = Deck('crimp_sfm3d4', 'S4R tube crimped by a shrinking cylinder of SFM3D4 surface elements (explicit)')
    # thick tube (t/r = 0.16): a thin one (t = 0.2) buckles in plastic hoop compression between the contacts, a
    # bifurcation that differs from run to run (Abaqus and asclepios alike)
    tids = tube(d, t=0.8)
    cids = {}
    for i in range(24):
        a = 2 * math.pi * i / 24
        for k in range(4):
            cids[i, k] = d.node(6.0 * math.cos(a), 6.0 * math.sin(a), -1.0 + 4.0 * k)
    d.elements('SFM3D4R', 'CRIMPER', [[cids[i, k], cids[(i + 1) % 24, k], cids[(i + 1) % 24, k + 1], cids[i, k + 1]]
                                     for i in range(24) for k in range(3)])
    # mass per unit area of a 0.2 mm steel sheet: without it the penalty of the massless crimper nodes is zero and
    # the crimper passes through the tube (Abaqus/Explicit, first reference of this case)
    d.model += ['*Surface Section, elset=CRIMPER, density=1.5e-09', '*Amplitude, name=CRIMP, definition=SMOOTH STEP',
                '0., 0., 0.002, 1.']
    body = ['*Boundary, amplitude=CRIMP']
    for (i, k), n in sorted(cids.items()):
        x, y, z = d.nodes[n]
        body += ['%d, 1, 1, %s' % (n, fmt(-x / 4)), '%d, 2, 2, %s' % (n, fmt(-y / 4)), '%d, 3, 3, 0.' % n]
    d.nset('TUBEZ0', [tids[i, 0] for i in range(16)])
    d.steps += explicit_step('Crimp', 0.0025, body + ['*Boundary', 'TUBEZ0, 3, 3'] + general_contact(d, 0.1), rot=True)
    register(d, 'explicit', 'S4R + SFM3D4', 'crimping with surface elements', 'Tube of radius 5 crimped by a cylinder '
             'of surface elements (mass per unit area 1.5e-9) shrinking from radius 6 to 4.5 (thick tube, t = 0.8).',
             [tids[0, 3], tids[4, 3]],
             [], True, 'contact', category='feature', group='contact', supported=True, twin=False)

    d = Deck('crimp_rigid_plates', 'S4R tube crimped by four rigid R3D4 plates moving radially (explicit)')
    tids = tube(d)
    rps = []
    # four 5 mm plates moving from the radius 6 to 4.5: the tube is pressed into four lobes that bulge freely between
    # the plates (eight plates leave narrow slots at the corners of the octagon, through which the tube escapes in a
    # way that depends on round-off)
    for p in range(4):
        a = 2 * math.pi * p / 4
        cx, cy = math.cos(a), math.sin(a)
        tx, ty = -cy, cx
        corners = [d.node(6.0 * cx + s * 2.5 * tx, 6.0 * cy + s * 2.5 * ty, z) for s, z in
                   ((-1, -1.0), (1, -1.0), (1, 11.0), (-1, 11.0))]
        d.elements('R3D4', 'PLATE%d' % p, [corners[::-1]])
        rp = d.node(6.0 * cx, 6.0 * cy, 5.0)
        d.model += ['*Rigid Body, ref node=%d, elset=PLATE%d' % (rp, p)]
        rps.append((rp, cx, cy))
    d.model += ['*Amplitude, name=CRIMP, definition=SMOOTH STEP', '0., 0., 0.002, 1.']
    body = ['*Boundary, amplitude=CRIMP']
    for rp, cx, cy in rps:
        body += ['%d, 1, 1, %s' % (rp, fmt(-1.5 * cx)), '%d, 2, 2, %s' % (rp, fmt(-1.5 * cy)), '%d, 3, 6, 0.' % rp]
    d.nset('TUBEZ0', [tids[i, 0] for i in range(16)])
    d.steps += explicit_step('Crimp', 0.0025, body + ['*Boundary', 'TUBEZ0, 3, 3'] + general_contact(d, 0.1), rot=True)
    register(d, 'explicit', 'S4R + R3D4', 'crimping with rigid plates', 'Tube of radius 5 crimped by four rigid '
             '5 mm plates moving radially from radius 6 to 4.5 (crimping force of one plate).',
             [tids[0, 3], tids[1, 3], tids[2, 3]], [rps[0][0]], True, 'contact', category='feature', group='contact',
             supported=True, twin=False)

    # shell section options
    d = Deck('shell_offset', 'S4R cantilever strip with the reference surface on its bottom face (OFFSET=SPOS)')
    quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
    shell(d, 'S4R', quads, 1.0, 'STEEL', extra=', offset=SPOS')
    fix = [ids[0, j] for j in range(3)]
    tip = [ids[10, j] for j in range(3)]
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 1, %s' % (n, fmt(100.0 * w)) for n, w in zip(tip, (0.25, 0.5, 0.25))], False,
                           rot=True, nframes=1)
    register(d, 'static', 'S4R', 'shell offset', 'Axial pull applied on the bottom (reference) surface of an offset '
             'strip: extension and bending.', tip + [ids[5, 1]], fix, True, 'linear', 'linear', category='feature',
             group='element', supported=False, reason=NS + 'shell offsets', twin=False)

    d = Deck('shell_composite', 'S4R cantilever strip with a two-layer composite section')
    quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
    d.elements('S4R', 'EALL', quads)
    d.materials.update(['STEEL', 'ALU'])
    d.model += ['*Shell Section, elset=EALL, composite', '0.3, 3, STEEL', '0.7, 3, ALU']
    fix = [ids[0, j] for j in range(3)]
    tip = [ids[10, j] for j in range(3)]
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 3, %s' % (n, fmt(-1.0 * w)) for n, w in zip(tip, (0.25, 0.5, 0.25))], False,
                           rot=True, nframes=1)
    register(d, 'static', 'S4R', 'composite section', 'Cantilever strip made of a steel and an aluminium layer, tip '
             'load (bending-extension coupling).', tip + [ids[5, 1]], fix, True, 'linear', 'linear',
             category='feature', group='element', supported=False, reason=NS + 'composite shell sections', twin=False)

    # amplitudes and acceleration
    for name, amp, title in (('amp_decay', ['*Amplitude, name=AMP, definition=DECAY', '1., -1., 0., 0.3'], 'DECAY'),
                             ('amp_equally_spaced', ['*Amplitude, name=AMP, definition=EQUALLY SPACED, '
                                                     'fixed interval=0.25', '0., 0.5, 1., 0.8, 0.6'],
                              'EQUALLY SPACED')):
        d = Deck(name, 'S4R strip, tip displacement driven by a %s amplitude' % title)
        quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
        shell(d, 'S4R', quads, 0.5, 'STEEL')
        fix = [ids[0, j] for j in range(3)]
        tip = [ids[10, j] for j in range(3)]
        d.nset('TIP', tip); d.nset('NALL', list(d.nodes))
        d.model += amp
        d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 1', '%d, 2, 2' % ids[0, 0],
                                        '*Boundary, amplitude=AMP', 'TIP, 1, 1, 0.04'], False, rot=True, inc=0.05,
                               nframes=8)
        register(d, 'static', 'S4R', '%s amplitude' % title.lower(), 'Homogeneous membrane strain following a %s '
                 'amplitude.' % title, tip, fix, True, 'exact', 'linear', reference='analytic', category='feature',
                 group='load', supported=True, twin=False)

    d = Deck('bc_acceleration', 'Free C3D8R block driven by a prescribed acceleration (explicit)')
    hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
    solid(d, 'C3D8R', hexes, 'STEEL')
    d.nset('NALL', list(d.nodes))
    d.steps += explicit_step('Accelerate', 0.001, ['*Boundary, type=ACCELERATION', 'NALL, 1, 1, 100000.'])
    register(d, 'explicit', 'C3D8R', 'prescribed acceleration', 'Block uniformly accelerated at 100 m/s^2 along x: '
             'u = a t^2 / 2, reaction m a.', list(d.nodes)[:2], list(d.nodes), False, 'exact', reference='analytic',
             category='feature', group='load', supported=True, twin=False)

    # procedures
    d, fix, tip = cantilever_solid('static_stabilize', 'C3D8R', 'C3D8R cantilever, static step with automatic '
                                   'stabilization')
    d.steps += ['*Step, name=Load, nlgeom=YES, inc=1000', '*Static, stabilize=0.0002', '0.1, 1., 1e-08, 0.1',
                '*Boundary', 'FIX, ENCASTRE', '*Cload'] + ['%d, 3, -50.' % n for n in tip] + \
        output_block(True, False, 10) + ['*End Step']
    register(d, 'static', 'C3D8R', 'stabilization', 'Cantilever loaded in a static step with automatic stabilization '
             '(dissipated energy fraction 2e-4).', tip, fix, False, 'nonlinear', category='feature',
             group='procedure', twin=False)

    d, fix, tip = cantilever_solid('static_perturbation', 'C3D8R', 'C3D8R cantilever, linear perturbation step')
    d.steps += ['*Step, name=Perturb, perturbation', '*Static', '*Boundary', 'FIX, ENCASTRE', '*Cload'] + \
        ['%d, 3, -50.' % n for n in tip] + ['*Output, field', '*Node Output', 'U, RF', '*End Step']
    register(d, 'static', 'C3D8R', 'linear perturbation', 'Cantilever load case solved as a linear perturbation step.',
             tip, fix, False, 'linear', 'linear', category='feature', group='procedure', supported=False,
             reason=NS + 'linear perturbation steps', twin=False)

    d = Deck('visco_relaxation', 'Viscoelastic block, step strain then relaxation (*VISCO)')
    hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
    solid(d, 'C3D8R', hexes, 'VISCO')
    fx = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 4.0))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.nset('TIP', tip)
    d.steps += ['*Step, name=Stretch, nlgeom=NO', '*Static', '1., 1.', '*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3',
                'TIP, 1, 1, 0.2', '*Output, field', '*Node Output', 'U, RF', '*End Step',
                '*Step, name=Relax, nlgeom=NO, inc=1000', '*Visco, cetol=0.0001', '1e-05, 0.01, 1e-09, 0.0005',
                '*Output, field, number interval=10', '*Node Output', 'U, RF', '*End Step']
    register(d, 'visco', 'C3D8R', 'viscoelastic relaxation', 'Block stretched by 5 % then held: stress relaxation in '
             'a *VISCO step.', tip[:3], fx, False, 'nonlinear', category='feature', group='procedure',
             supported=False, reason=NS + '*VISCO steps', twin=False)

    from gen_features import cantilever_beam
    d, ids = cantilever_beam('modal_dynamic', 'B31 cantilever, frequency extraction then modal dynamics')
    d.model += ['*Amplitude, name=PULSE', '0., 0., 0.001, 1., 0.002, 0., 1., 0.']
    d.steps += ['*Step, name=Modes, perturbation', '*Frequency, eigensolver=Lanczos', '10, , , , ,',
                '*Boundary', '%d, ENCASTRE' % ids[0], '*End Step',
                '*Step, name=Response', '*Modal Dynamic', '0.0001, 0.05', '*Modal Damping', '1, 10, 0.02',
                '*Cload, amplitude=PULSE', '%d, 3, -1.' % ids[-1], '*Output, field, number interval=20',
                '*Node Output', 'U', '*End Step']
    register(d, 'modal dynamic', 'B31', 'modal dynamics', 'Response of a cantilever to a short tip pulse by mode '
             'superposition (10 modes, 2 % damping).', [ids[-1], ids[5]], [ids[0]], True, 'dynamic',
             category='feature', group='procedure', supported=False, reason=NS + '*MODAL DYNAMIC', twin=False)

    d, ids = cantilever_beam('steady_state_dynamics', 'B31 cantilever, harmonic response (*STEADY STATE DYNAMICS)')
    d.steps += ['*Step, name=Modes, perturbation', '*Frequency, eigensolver=Lanczos', '10, , , , ,',
                '*Boundary', '%d, ENCASTRE' % ids[0], '*End Step',
                '*Step, name=Harmonic, perturbation', '*Steady State Dynamics', '10., 1000., 40, 1.',
                '*Modal Damping', '1, 10, 0.02', '*Cload', '%d, 3, -1.' % ids[-1], '*Output, field',
                '*Node Output', 'U', '*End Step']
    register(d, 'steady state dynamics', 'B31', 'harmonic response', 'Frequency response of a cantilever between 10 '
             'and 1000 Hz under a unit harmonic tip load.', [ids[-1], ids[5]], [ids[0]], True, 'dynamic',
             category='feature', group='procedure', supported=False, reason=NS + '*STEADY STATE DYNAMICS', twin=False)

    d = Deck('heat_transfer_steady', 'DC3D8 bar, steady heat conduction between two imposed temperatures')
    hexes = hex_block(d, 10, 1, 1, 40.0, 4.0, 4.0)
    d.elements('DC3D8', 'EALL', hexes)
    d.materials.add('COND')
    d.model += ['*Solid Section, elset=EALL, material=COND', ',']
    d.nset('HOT', d.select(lambda x, y, z: near(x, 40.0)))
    cold = d.select(lambda x, y, z: near(x, 0))
    d.steps += ['*Step, name=Conduct', '*Heat Transfer, steady state', '1., 1.', '*Boundary', 'FIX, 11, 11, 20.',
                'HOT, 11, 11, 100.', '*Output, field', '*Node Output', 'NT', '*End Step']
    register(d, 'heat transfer', 'DC3D8', 'heat transfer', 'Steady conduction along a bar between 20 and 100 degrees: '
             'linear temperature profile.', d.select(lambda x, y, z: near(y, 0) and near(z, 0) and x in (10.0, 20.0,
                                                                                                              30.0)),
             cold, False, 'exact', 'linear', reference='analytic', category='feature', group='procedure',
             supported=False, reason=NS + 'heat transfer', twin=False)

    d = Deck('coupled_temp_displacement', 'C3D8T bar heated at one end, free expansion (coupled analysis)')
    hexes = hex_block(d, 10, 1, 1, 40.0, 4.0, 4.0)
    d.elements('C3D8T', 'EALL', hexes)
    d.materials.add('COND')
    d.model += ['*Solid Section, elset=EALL, material=COND', ',']
    fx = d.select(lambda x, y, z: near(x, 0))
    d.nset('HOT', d.select(lambda x, y, z: near(x, 40.0)))
    d.nset('NALL', list(d.nodes))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.initial += ['*Initial Conditions, type=TEMPERATURE', 'NALL, 20.']
    d.steps += ['*Step, name=Heat, inc=1000', '*Coupled Temperature-Displacement, deltmx=10.', '0.1, 10., 1e-06, 1.',
                '*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', 'FIX, 11, 11, 20.', 'HOT, 11, 11, 100.',
                '*Output, field, number interval=10', '*Node Output', 'U, RF, NT', '*End Step']
    register(d, 'coupled temperature-displacement', 'C3D8T', 'thermomechanical coupling', 'Transient conduction '
             'from a hot end and the resulting expansion.', d.select(lambda x, y, z: near(x, 40.0))[:2], fx, False,
             'nonlinear', category='feature', group='procedure', supported=False,
             reason=NS + 'coupled temperature-displacement analysis', twin=False)

    d = Deck('ale_adaptive_mesh', 'Taylor-type impact with adaptive meshing (ALE, explicit)')
    hexes = hex_block(d, 3, 3, 20, 3.2, 3.2, 32.4)
    solid(d, 'C3D8R', hexes, 'COPPER')
    d.nset('WALL', d.select(lambda x, y, z: near(z, 0)))
    d.nset('XS', d.select(lambda x, y, z: near(x, 0))); d.nset('YS', d.select(lambda x, y, z: near(y, 0)))
    moving = d.select(lambda x, y, z: z > 1e-9)
    d.nset('MOVING', moving)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'MOVING, 3, -227000.']
    d.steps += explicit_step('Impact', 8e-05, ['*Boundary', 'WALL, 3, 3', 'XS, XSYMM', 'YS, YSYMM'],
                             extra=['*Adaptive Mesh, elset=EALL, frequency=10'])
    register(d, 'explicit', 'C3D8R', 'adaptive meshing', 'Taylor impact of taylor_impact_c3d8r with adaptive '
             'meshing every 10 increments.', d.select(lambda x, y, z: near(z, 32.4) and near(x, 0) and near(y, 0)),
             d.nsets['WALL'], False, 'dynamic', category='feature', group='procedure', supported=False,
             reason=NS + 'adaptive meshing (ALE)', twin=False)

    # transfer between procedures
    for stage in ('a', 'b'):
        d = Deck('import_%s_%s' % (stage, 'static' if stage == 'a' else 'explicit'),
                 'Strip bent statically (Abaqus/Standard), then released in Abaqus/Explicit (*IMPORT)')
        quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
        fix = [ids[0, j] for j in range(3)]
        tip = [ids[10, j] for j in range(3)]
        if stage == 'a':
            shell(d, 'S4R', quads, 0.5, 'STEEL')
            d.nset('TIP', tip)
            d.steps += static_step('Bend', ['*Boundary', 'FIX, ENCASTRE', '*Boundary', 'TIP, 3, 3, -4.',
                                            '*Restart, write, number interval=1'], True, rot=True)
            register(d, 'static', 'S4R', 'transfer Standard to Explicit (1/2)', 'Strip bent by a 4 mm tip '
                     'displacement; state imported by import_b_explicit.', tip, fix, True, 'nonlinear',
                     category='feature', group='procedure', twin=False)
        else:
            d.title += '\n** validation: oldjob=import_a_static'
            d.import_only = True
            d.steps += ['*Import, step=1, increment=10, state=yes, update=no', 'EALL',
                        '*Step, name=Release, nlgeom=YES', '*Dynamic, Explicit', ', 0.002',
                        '*Boundary', 'FIX, ENCASTRE'] + output_block(False, True, 20, 0.002) + ['*End Step']
            register(d, 'explicit', 'S4R', 'transfer Standard to Explicit (2/2)', 'Bent strip imported from '
                     'import_a_static (stresses and deformed shape) and released: free vibration.', tip, fix, True,
                     'dynamic', category='feature', group='procedure', supported=False,
                     reason=NS + '*IMPORT (transfer between procedures)', twin=False)


# --------------------------------------------------------------------------------------- application models
def application_features():
    from deck import MATERIALS as M  # noqa
    d = Deck('app_mini_valve', 'Two Fung S4R leaflets closing under pressure against each other (explicit)')
    c, s = 1.0, 0.0
    lids = []
    for k, x0 in enumerate((0.0, 4.0)):
        ids = {}
        for i in range(11):
            for j in range(9):
                ids[i, j] = d.node(x0, 10.0 * i / 10, 8.0 * j / 8)
        quads = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for i in range(10) for j in range(8)]
        if k == 1:
            quads = [q[::-1] for q in quads]
        d.elements('S4R', 'LEAF%d' % (k + 1), quads)
        lids.append(ids)
    d.materials.add('FUNG')
    d.model += ['*Orientation, name=ORI, local directions=2', '0., 1., 0., 0., 0., 1.',
                '*Shell Section, elset=LEAF1, material=FUNG, orientation=ORI', '0.4, 5',
                '*Shell Section, elset=LEAF2, material=FUNG, orientation=ORI', '0.4, 5']
    ring = []
    for k in range(2):
        segs = [[lids[k][i, 0], lids[k][i + 1, 0]] for i in range(10)]
        d.elements('B31', 'RING%d' % (k + 1), segs)
        ring += [lids[k][i, 0] for i in range(11)]
    d.materials.add('STEEL')
    d.model += ['*Beam Section, elset=RING1, material=STEEL, section=CIRC', '0.5', '1., 0., 0.',
                '*Beam Section, elset=RING2, material=STEEL, section=CIRC', '0.5', '1., 0., 0.']
    d.model += ['*Surface, type=ELEMENT, name=P1', 'LEAF1, SNEG', '*Surface, type=ELEMENT, name=P2', 'LEAF2, SNEG',
                '*Amplitude, name=CLOSE, definition=SMOOTH STEP', '0., 0., 0.004, 1.']
    ends = [lids[k][0, 0] for k in range(2)] + [lids[k][10, 0] for k in range(2)]
    d.nset('ENDS', ends)
    d.steps += explicit_step('Close', 0.005, ['*Boundary', 'ENDS, ENCASTRE', '*Dsload, amplitude=CLOSE',
                                              'P1, P, 0.002', 'P2, P, 0.002'] + general_contact(d, 0.1), rot=True)
    register(d, 'explicit', 'S4R Fung + B31', 'valve closure', 'Two anisotropic hyperelastic leaflets on flexible '
             'annulus segments, closed by pressure until their free edges meet (contact with friction).',
             [lids[0][5, 8], lids[1][5, 8], lids[0][5, 4]], ends, True, 'contact', category='feature',
             group='application', twin=False)

    d = Deck('app_stent_crimp', 'B31 zig-zag stent ring crimped by a shrinking S4R tube (explicit)')
    nz = 12
    pts = []
    for i in range(2 * nz):
        a = 2 * math.pi * i / (2 * nz)
        pts.append(d.node(4.0 * math.cos(a), 4.0 * math.sin(a), 1.0 if i % 2 else -1.0))
    d.elements('B31', 'STENT', [[pts[i], pts[(i + 1) % (2 * nz)]] for i in range(2 * nz)])
    d.materials.update(['STEEL_PL', 'STEEL'])
    d.model += ['*Beam Section, elset=STENT, material=STEEL_PL, section=CIRC', '0.1', '0., 0., 1.']
    tids = {}
    for i in range(24):
        a = 2 * math.pi * i / 24
        for k in range(5):
            tids[i, k] = d.node(5.0 * math.cos(a), 5.0 * math.sin(a), -2.0 + k)
    d.elements('S4R', 'TUBE', [[tids[i, k], tids[(i + 1) % 24, k], tids[(i + 1) % 24, k + 1], tids[i, k + 1]]
                               for i in range(24) for k in range(4)])
    d.model += ['*Shell Section, elset=TUBE, material=STEEL', '0.1, 3',
                '*Amplitude, name=CRIMP, definition=SMOOTH STEP', '0., 0., 0.004, 1.']
    body = ['*Boundary, amplitude=CRIMP']
    for (i, k), n in sorted(tids.items()):
        x, y, z = d.nodes[n]
        body += ['%d, 1, 1, %s' % (n, fmt(-0.6 * x)), '%d, 2, 2, %s' % (n, fmt(-0.6 * y)), '%d, 3, 6, 0.' % n]
    d.nset('STENTZ', [pts[0]])
    d.steps += explicit_step('Crimp', 0.005, body + ['*Boundary', 'STENTZ, 3, 3'] + general_contact(d, 0.1),
                             rot=True)
    register(d, 'explicit', 'B31 + S4R', 'stent crimping', 'Elastic-plastic zig-zag ring of radius 4 crimped to radius '
             '2 by a tube whose radius shrinks from 5 to 2 (beam-to-shell contact with friction).',
             [pts[0], pts[1], pts[6]], [], True, 'contact', category='feature', group='application', twin=False)

    for n in (20, 40, 80):
        d = Deck('par_perf_%d' % n, 'Clamped S4R plate %d x %d elements, pressure pulse (explicit, performance)' % (n, n))
        quads, ids = quad_plate(d, n, n, 40.0, 40.0)
        shell(d, 'S4R', quads, 0.5, 'STEEL')
        d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
        edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 40.0) or near(y, 40.0))
        d.title += '\n** validation: single'
        d.steps += explicit_step('Pulse', 1e-4, ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.5'], rot=True)
        register(d, 'explicit', 'S4R', 'performance (%d elements)' % (n * n), 'Explicit run time against the model '
                 'size, with Abaqus on 1 and 4 cpus and in single precision.', [ids[n // 2, n // 2]], edge, True,
                 'dynamic', category='feature', group='parallel', twin=False)


def coverage_features():
    """Cases added from the coverage report (keywords and element types read by asclepios without a case)."""
    d = Deck('dload_solid_pressure', 'C3D8R block compressed by a pressure on its top faces (*DLOAD P2)')
    hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
    el = solid(d, 'C3D8R', hexes, 'STEEL')
    top = faces_on(d, 'C3D8R', el, hexes, lambda x, y, z: near(z, 4.0))
    d.elsets['ETOP'] = sorted(set(e for e, f in top))
    fz = d.select(lambda x, y, z: near(z, 0))
    d.nset('FX', d.select(lambda x, y, z: near(x, 0))); d.nset('FY', d.select(lambda x, y, z: near(y, 0)))
    d.steps += static_step('Load', ['*Boundary', 'FIX, 3, 3', 'FX, 1, 1', 'FY, 2, 2', '*Dload', 'ETOP, P2, 100.'],
                           False, nframes=1)
    register(d, 'static', 'C3D8R', '*DLOAD on solid faces', 'Block under a 100 MPa pressure on its top faces '
             '(element-based distributed load), symmetry supports: homogeneous compression.',
             d.select(lambda x, y, z: near(z, 4.0))[:3], fz, False, 'exact', 'linear', reference='analytic',
             category='feature', group='load')

    d = Deck('connector_reference_length', 'CONN3D2 AXIAL spring with a reference length shorter than its length')
    a, b = d.node(0, 0, 0), d.node(10, 0, 0)
    d.elements('CONN3D2', 'CONN', [[a, b]])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL', '*Connector Behavior, name=SPRING',
                '*Connector Elasticity, component=1', '100.', '*Connector Constitutive Reference', '9.']
    d.steps += static_step('Release', ['*Boundary', '%d, 1, 3' % a, '%d, 2, 3' % b], True, inc=0.1)
    register(d, 'static', 'CONN3D2', 'connector reference length', 'Axial spring of length 10 with a free length of 9: '
             'the free end moves back by 1 mm.', [b], [a], False, 'exact', reference='analytic', category='feature',
             group='constraint', twin=False)

    for etype in ('S3', 'S3RS', 'S4RS', 'S4R5'):
        d = Deck('elem_%s' % etype.lower(), '%s clamped square plate under pressure' % etype)
        quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
        conns = tris_from_quads(quads) if etype.startswith('S3') else quads
        d.elements(etype, 'EALL', conns)
        if etype in ('S3RS', 'S4RS') and not MODE['twin']:  # Abaqus/Explicit elements: quasi-static explicit reference
            d.title += '\n** validation: explicit_static'
        d.materials.add('STEEL')
        d.model += ['*Shell Section, elset=EALL, material=STEEL', '0.2, 5', '*Surface, type=ELEMENT, name=STOP',
                    'EALL, SPOS']
        edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
        QS_DURATION[d.name] = 0.002
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.05'], False, rot=True,
                               nframes=1)
        register(d, 'static', etype, 'element %s' % etype, 'Clamped square plate under a small pressure (linear).',
                 [ids[4, 4], ids[2, 2]], edge, True, 'linear', 'linear', category='feature', group='element',
                 twin=etype not in ('S4R5',))

    # S3RS on thick elements (t > L / 5): full transverse shear stiffness
    d = Deck('elem_s3rs_thick', 'S3RS thick clamped square plate under pressure')
    if not MODE['twin']:
        d.title += '\n** validation: explicit_static'
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    d.elements('S3RS', 'EALL', tris_from_quads(quads))
    d.materials.add('STEEL')
    d.model += ['*Shell Section, elset=EALL, material=STEEL', '1.5, 5', '*Surface, type=ELEMENT, name=STOP',
                'EALL, SPOS']
    edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
    QS_DURATION[d.name] = 0.002
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 5.'], False, rot=True,
                           nframes=1)
    register(d, 'static', 'S3RS', 'element S3RS, thick', 'Thick clamped square plate (a / t = 13, elements thicker '
             'than a fifth of their shortest edge: full transverse shear stiffness) under pressure (linear).',
             [ids[4, 4], ids[2, 2]], edge, True, 'linear', 'linear', category='feature', group='element')

    # PIPE (thin-walled) and CIRC cantilevers, tip force and torque
    d = Deck('b31_pipe_sections', 'B31 PIPE and CIRC cantilevers, tip force and torque')
    tips, roots = [], []
    for k, (sec, dims) in enumerate((('PIPE', '3., 0.5'), ('PIPE', '3., 1.5'), ('PIPE', '3., 2.9'), ('CIRC', '3.'))):
        segs, ids = beam_line(d, 20, 100.0, x0=(0.0, 50.0 * k, 0.0))
        d.elements('B31', 'B%d' % (k + 1), segs)
        d.model += ['*Beam Section, elset=B%d, material=STEEL, section=%s' % (k + 1, sec), dims, '0., 0., -1.']
        tips.append(ids[-1]); roots.append(ids[0])
    d.materials.add('STEEL')
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 3, -10.' % n for n in tips] + ['%d, 4, 1000.' % n for n in tips], False, rot=True,
                           nframes=1)
    register(d, 'static', 'B31', 'PIPE and CIRC sections', 'Four cantilevers (PIPE with wall thickness 0.5, 1.5 and '
             '2.9 for an outer radius 3, CIRC of radius 3) under a tip force and a tip torque: thin-walled PIPE '
             'section (bending and torsion stiffness on the mean radius).', tips, roots, True, 'linear', 'linear',
             category='feature', group='element', twin=False)

    d = Deck('elem_b31h', 'B31H (hybrid) cantilever, tip loads')
    segs, ids = beam_line(d, 50, 100.0)
    d.elements('B31H', 'EALL', segs)
    d.materials.add('STEEL')
    d.model += ['*Beam Section, elset=EALL, material=STEEL, section=RECT', '2., 1.', '0., 0., -1.']
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload', '%d, 2, 0.5' % ids[-1],
                                    '%d, 3, -2.' % ids[-1]], False, rot=True, nframes=1)
    register(d, 'static', 'B31H', 'element B31H', 'Hybrid beam cantilever with tip loads in both planes (same solution '
             'as b31_bending).', [ids[-1], ids[25]], [ids[0]], True, 'theory', 'linear', reference='analytic',
             category='feature', group='element')


def generate():
    coverage_features()
    patch_features()
    benchmark_features()
    input_features()
    more_features()
    application_features()
