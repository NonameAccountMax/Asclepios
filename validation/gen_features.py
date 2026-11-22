"""Solver feature cases: contact between every pair of element types, constraints, connectors, element types,
loads and boundary conditions, materials, options specific to each procedure, parallel runs.

Small and fast models, one feature each. Cases that asclepios does not support yet are marked supported=False: they run
in Abaqus only and list what remains to be implemented. Static cases also get a quasi-static explicit twin and
explicit cases an implicit dynamic twin (see deck.py), except the procedure-specific ones.
"""
import math

from deck import (QS_DAMPING, QS_DURATION, Deck, _smooth_amplitude, beam_line, edge_shares, explicit_step, faces_on, fmt, hex_block, near, output_block,
                  quad_plate, register, static_step, tets_from_hexes)

NS = 'not supported by asclepios yet: '

# B31 sections of the rotary inertia cases: largest second moment of area (PIPE: thin wall at the mean radius)
BEAM_IMAX = {'CIRC': math.pi / 4, 'RECT': 8.0 / 12, 'PIPE': math.pi * 2.75 ** 3 * 0.5}
# twisted free elements: section, dimensions, end points, axial moment per node
BEAM_TWIST = [('RECT', '2., 1.', (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1e5),
              ('PIPE', '3., 0.5', (0.0, 50.0, 0.0), (2.0, 50.0, 0.0), 1e7),
              ('RECT', '2., 1.', (0.0, 100.0, 0.0), (1.0, 101.0, 1.0), 1e5)]


def shell(d, etype, conns, t, mat, elset='EALL'):
    d.elements(etype, elset, conns)
    d.materials.add(mat)
    d.model += ['*Shell Section, elset=%s, material=%s' % (elset, mat), '%s, 5' % fmt(t)]


def solid(d, etype, conns, mat, elset='EALL'):
    el = d.elements(etype, elset, conns)
    d.materials.add(mat)
    d.model += ['*Solid Section, elset=%s, material=%s' % (elset, mat), ',']
    return el


def beam(d, conns, section, dims, mat='STEEL', elset='EALL', etype='B31', n1='0., 0., -1.'):
    d.elements(etype, elset, conns)
    d.materials.add(mat)
    d.model += ['*Beam Section, elset=%s, material=%s, section=%s' % (elset, mat, section), dims, n1]


def surface(d, name, pairs):
    d.model += ['*Surface, type=ELEMENT, name=%s' % name] + ['%d, S%d' % ef for ef in pairs]


def general_contact(d, friction=None):
    d.model += ['*Surface Interaction, name=IP'] + (['*Friction', fmt(friction)] if friction else [])
    return ['*Contact', '*Contact Inclusions, ALL EXTERIOR', '*Contact Property Assignment', ' , , IP']


def clamped_plate(d, n=8, size=20.0, t=0.5, z=0.0):
    quads, ids = quad_plate(d, n, n, size, size, x0=(0.0, 0.0, z))
    shell(d, 'S4R', quads, t, 'STEEL', 'PLATE')
    edge = [ids[i, j] for (i, j) in ids if i in (0, n) or j in (0, n)]
    return ids, edge


# ------------------------------------------------------------------------------------------------ contact
def contact_features():
    d = Deck('contact_solid_on_solid', 'C3D8R block dropped on a fixed C3D8R base (explicit general contact)')
    base = hex_block(d, 4, 4, 1, 12.0, 12.0, 2.0)
    solid(d, 'C3D8R', base, 'STEEL', 'BASE')
    blk = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0, x0=(4.0, 4.0, 2.5))
    solid(d, 'C3D8R', blk, 'ALU', 'BLOCK')
    fix = d.select(lambda x, y, z: near(z, 0))
    moving = d.select(lambda x, y, z: z > 2.2)
    d.nset('MOVING', moving)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'MOVING, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE'] + general_contact(d))
    register(d, 'explicit', 'C3D8R', 'contact solid-solid', 'Aluminium block dropped at 2 m/s on a steel base: '
             'solid-to-solid impact (general contact).', moving[:3], fix, False, 'contact', category='feature',
             group='contact', supported=True)

    d = Deck('contact_beam_on_solid', 'B31 rod dropped on a fixed C3D8R base (explicit general contact)')
    base = hex_block(d, 4, 4, 1, 12.0, 12.0, 2.0)
    solid(d, 'C3D8R', base, 'STEEL', 'BASE')
    segs, rod = beam_line(d, 8, 10.0, x0=(1.0, 6.0, 3.0))
    beam(d, segs, 'CIRC', '0.5', elset='ROD')
    fix = d.select(lambda x, y, z: near(z, 0))
    d.nset('RODN', rod)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'RODN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE'] + general_contact(d), rot=True)
    register(d, 'explicit', 'B31 + C3D8R', 'contact beam-solid', 'Steel rod dropped at 2 m/s on a steel block '
             '(general contact).', [rod[0], rod[4]], fix, True, 'contact', category='feature', group='contact',
             supported=True)

    d = Deck('contact_beam_beam', 'B31 rod dropped across a fixed B31 rod (explicit general contact, edge to edge)')
    segs, low = beam_line(d, 9, 18.0, x0=(-9.0, 0.0, 0.0))
    beam(d, segs, 'CIRC', '0.5', elset='LOW')
    segs, top = beam_line(d, 9, 18.0, x0=(0.0, -9.0, 1.5), axis=1)
    beam(d, segs, 'CIRC', '0.5', elset='TOP', n1='0., 0., -1.')
    d.nset('TOPN', top)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'TOPN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', '%d, ENCASTRE' % low[0], '%d, ENCASTRE' % low[-1]] +
                             general_contact(d), rot=True)
    register(d, 'explicit', 'B31', 'contact beam-beam', 'Steel rod dropped at 2 m/s across a fixed perpendicular rod, '
             'crossing between nodes (edge-to-edge general contact).', [top[4], low[4]], [low[0], low[-1]], True,
             'contact', category='feature', group='contact', supported=False, reason=NS + 'beam-to-beam edge contact')

    d = Deck('contact_self_shell', 'S4R strip coiled by an end moment until it touches itself (self-contact)')
    quads, ids = quad_plate(d, 40, 2, 40.0, 4.0)
    shell(d, 'S4R', quads, 0.3, 'STEEL0')
    fix = [ids[0, j] for j in range(3)]
    tip = [ids[40, j] for j in range(3)]
    ei = 210000.0 * 4.0 * 0.3 ** 3 / 12
    m = 2.4 * math.pi * ei / 40.0  # 1.2 turns if free
    d.model += ['*Amplitude, name=COIL, definition=SMOOTH STEP', '0., 0., 0.008, 1.']
    # up to the end of the first contact (two presses, 5.74 to 5.79 ms): the later chatter of the coil on the
    # clamped edge goes straight into the reaction and Abaqus itself loses the symmetry of the strip
    d.steps += explicit_step('Coil', 0.0058, ['*Boundary', 'FIX, ENCASTRE', '*Cload, amplitude=COIL'] +
                             ['%d, 5, %s' % (n, fmt(-m * w)) for n, w in zip(tip, (0.25, 0.5, 0.25))] +
                             general_contact(d), rot=True, nframes=58, time_marks=True)
    register(d, 'explicit', 'S4R', 'self-contact', 'Strip coiled by an end moment (1.2 turns if free): the coil comes '
             'into contact with the clamped end of the strip (self-contact in general contact), up to the end of the '
             'first contact.', tip + [ids[20, 1]], fix, True, 'contact', category='feature', group='contact')

    d = Deck('contact_shell_edge', 'S4R strip dropped on its free edge onto a clamped S4R plate (general contact)')
    ids, edge = clamped_plate(d)
    squads, sids = quad_plate(d, 3, 2, 6.0, 4.0, x0=(7.0, 0.0, 0.0))
    for k, (x, y, z) in list(d.nodes.items()):
        if k in sids.values():  # strip in the plane y = 10, lower edge 1 above the plate
            d.nodes[k] = (x, 10.0, 1.0 + y)
    shell(d, 'S4R', squads, 0.6, 'STEEL', 'STRIP')
    d.nset('STRIPN', list(sids.values()))
    d.initial += ['*Initial Conditions, type=VELOCITY', 'STRIPN, 3, -2000.']
    d.steps += explicit_step('Impact', 4.3e-4, ['*Boundary', 'FIX, ENCASTRE'] + general_contact(d), rot=True)
    register(d, 'explicit', 'S4R', 'contact shell edge', 'Vertical strip (thickness 0.6) dropped at 2 m/s on its '
             'free edge onto a clamped plate, up to the end of the first impact: the free edge of a shell is rounded '
             'within the shell (moved back by the half thickness), which sets the instant of impact.', [sids[0, 0], sids[1, 0], sids[1, 2], ids[4, 4]],
             edge, True, 'contact', category='feature', group='contact')

    d = Deck('contact_penalty_soft_light', 'Soft and light S4R sheets dropped flat on clamped S4R plates (general contact penalty)')
    pq, pids = quad_plate(d, 2, 1, 2.0, 1.0)
    shell(d, 'S4R', pq, 0.2, 'STEEL', 'PLATE')
    sheets = []
    for k, (mat, x0, gap) in enumerate((('SOFT7', 0.3, 0.25), ('LIGHT7', 1.3, 0.27))):
        q, sid = quad_plate(d, 1, 1, 0.4, 0.4, x0=(x0, 0.3, gap))  # 0.05 / 0.07 above contact
        shell(d, 'S4R', q, 0.2, mat, 'SHEET%d' % (k + 1))
        sheets.append(sid)
    d.nset('SHEETN', [n for sid in sheets for n in sid.values()])
    d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -1000.']
    d.steps += explicit_step('Drop', 1.2e-4, ['*Boundary', 'FIX, ENCASTRE'] + general_contact(d), rot=True,
                             nframes=40)
    register(d, 'explicit', 'S4R', 'contact penalty', 'Two free single-element sheets dropped flat at 1 m/s on clamped '
             'plates: a soft one (E 7000) and a light one (density / 10). The penalty follows the stable increment of the '
             'elements of each node (stiffer spring for light nodes of soft elements), which sets the duration of the '
             'impact and the rebound velocity.', [sheets[0][0, 0], sheets[1][0, 0]], list(pids.values()), True,
             'contact', category='feature', group='contact')

    d = Deck('contact_pair_surface', 'S4R sheet dropped on a clamped S4R plate, surface-to-surface contact pair')
    ids, edge = clamped_plate(d)
    squads, sids = quad_plate(d, 3, 3, 6.0, 6.0, x0=(7.0, 7.0, 1.0))
    shell(d, 'S4R', squads, 0.5, 'ALU', 'SHEET')
    d.model += ['*Surface, type=ELEMENT, name=SPLATE', 'PLATE, SPOS', '*Surface, type=ELEMENT, name=SSHEET',
                'SHEET, SNEG', '*Surface Interaction, name=IP']
    d.nset('SHEETN', list(sids.values()))
    d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE', '*Contact Pair, interaction=IP',
                                              'SSHEET, SPLATE'], rot=True)
    register(d, 'explicit', 'S4R', 'contact pair', 'Sheet dropped at 2 m/s on a clamped plate with a '
             'surface-to-surface contact pair instead of general contact.', [ids[4, 4], sids[1, 1]], edge, True,
             'contact', category='feature', group='contact', supported=False, reason=NS + '*CONTACT PAIR')

    d = Deck('contact_rigid_surface', 'S4R sheet dropped on a rigid R3D4 surface (explicit general contact)')
    rquads, rids = quad_plate(d, 4, 4, 20.0, 20.0)
    d.elements('R3D4', 'RIGID', rquads)
    rp = d.node(10.0, 10.0, -1.0)
    d.model += ['*Rigid Body, ref node=%d, elset=RIGID' % rp]
    squads, sids = quad_plate(d, 3, 3, 6.0, 6.0, x0=(7.0, 7.0, 1.0))
    shell(d, 'S4R', squads, 0.5, 'ALU', 'SHEET')
    d.nset('SHEETN', list(sids.values()))
    d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', '%d, ENCASTRE' % rp] + general_contact(d), rot=True)
    register(d, 'explicit', 'R3D4 + S4R', 'contact with rigid surface', 'Sheet dropped at 2 m/s on a rigid R3D4 '
             'surface (rigid body with a reference node).', [sids[1, 1], sids[0, 0]], [rp], True, 'contact',
             category='feature', group='contact', supported=True)


# -------------------------------------------------------------------------------------------- constraints
def cantilever_bar(d, nx, ny, nz, lx, x0=0.0, elset='EALL', mat='STEEL'):
    hexes = hex_block(d, nx, ny, nz, lx, 4.0, 4.0, x0=(x0, 0.0, 0.0))
    el = solid(d, 'C3D8R', hexes, mat, elset)
    return el, hexes


def constraint_features():
    d = Deck('tie_solid', 'Two C3D8R bars with non-matching meshes tied at their interface, cantilever tip load')
    e1, h1 = cantilever_bar(d, 4, 2, 2, 10.0, elset='BAR1')
    e2, h2 = cantilever_bar(d, 5, 3, 3, 10.0, x0=10.0, elset='BAR2')
    surface(d, 'MASTER', faces_on(d, 'C3D8R', e1, h1, lambda x, y, z: near(x, 10.0)))
    surface(d, 'SLAVE', faces_on(d, 'C3D8R', e2, h2, lambda x, y, z: near(x, 10.0)))
    d.model += ['*Tie, name=TIE1, adjust=yes', 'SLAVE, MASTER']
    fix = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 20.0))
    QS_DURATION['tie_solid'] = 0.005
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 3, %s' % (n, fmt(-100.0 / len(tip))) for n in tip], False, nframes=1)
    register(d, 'static', 'C3D8R', '*TIE', 'Cantilever made of two bars with non-matching meshes tied at mid-length '
             '(surface-based tie), 100 N tip load.', tip[:2] + d.select(lambda x, y, z: near(x, 10.0) and near(z, 4.0)),
             fix, False, 'linear', 'linear', category='feature', group='constraint', supported=True)

    d = Deck('tie_shell', 'Two overlapping S4R strips tied face to face (lap joint), cantilever tip load')
    q1, i1 = quad_plate(d, 6, 2, 24.0, 6.0)
    shell(d, 'S4R', q1, 0.5, 'STEEL', 'STRIP1')
    q2, i2 = quad_plate(d, 5, 3, 20.0, 6.0, x0=(20.0, 0.0, 0.5))
    shell(d, 'S4R', q2, 0.5, 'STEEL', 'STRIP2')
    d.model += ['*Surface, type=ELEMENT, name=MASTER', 'STRIP1, SPOS', '*Surface, type=ELEMENT, name=SLAVE',
                'STRIP2, SNEG', '*Tie, name=TIE1, adjust=no, position tolerance=1.', 'SLAVE, MASTER']
    fix = [i1[0, j] for j in range(3)]
    tip = [i2[5, j] for j in range(4)]
    QS_DURATION['tie_shell'] = 0.05
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 3, %s' % (n, fmt(-0.5 / len(tip))) for n in tip], False, rot=True, nframes=1)
    register(d, 'static', 'S4R', '*TIE', 'Cantilever made of two overlapping strips tied face to face over 4 mm, '
             '0.5 N tip load.', tip[:2] + [i1[6, 1]], fix, True, 'linear', 'linear', category='feature',
             group='constraint', supported=True)

    for kind in ('kinematic', 'distributing'):
        d = Deck('coupling_%s' % kind, 'C3D8R bar, end face coupled to a reference node (%s coupling)' % kind)
        el, hexes = cantilever_bar(d, 8, 2, 2, 20.0)
        rp = d.node(20.0, 2.0, 2.0)
        d.nset('RP', [rp])
        surface(d, 'SEND', faces_on(d, 'C3D8R', el, hexes, lambda x, y, z: near(x, 20.0)))
        d.model += ['*Coupling, constraint name=CP1, ref node=RP, surface=SEND',
                    '*Kinematic' if kind == 'kinematic' else '*Distributing, weighting method=UNIFORM']
        fix = d.select(lambda x, y, z: near(x, 0))
        QS_DURATION[d.name] = 0.02
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload', 'RP, 3, -100.', 'RP, 4, 500.'], False,
                               rot=True, nframes=1)
        register(d, 'static', 'C3D8R', '*COUPLING %s' % kind, 'Bar loaded through a reference node coupled to its end '
                 'face: transverse force and torque.', [rp] + d.select(lambda x, y, z: near(x, 20.0) and near(z, 4.0)),
                 fix, True, 'linear', 'linear', category='feature', group='constraint', supported=True)

    d = Deck('mpc_beam', 'Two B31 cantilevers joined by a rigid MPC BEAM link')
    s1, b1 = beam_line(d, 40, 20.0)
    s2, b2 = beam_line(d, 40, 20.0, x0=(21.0, 0.0, 0.0))
    beam(d, s1 + s2, 'RECT', '2., 1.')
    d.model += ['*MPC', 'BEAM, %d, %d' % (b2[0], b1[-1])]
    QS_DURATION['mpc_beam'] = 0.05
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % b1[0], '*Cload', '%d, 3, -1.' % b2[-1]], False,
                           rot=True, nframes=1)
    register(d, 'static', 'B31', '*MPC BEAM', 'Two cantilever segments joined by a rigid beam MPC across a 1 mm gap, '
             'tip load.', [b2[-1], b1[-1], b2[0]], [b1[0]], True, 'theory', 'linear', category='feature',
             group='constraint', supported=True, reference='analytic')

    d = Deck('equation', 'Two parallel B31 cantilevers with equal tip deflections imposed by *EQUATION')
    s1, b1 = beam_line(d, 40, 20.0)
    s2, b2 = beam_line(d, 40, 20.0, x0=(0.0, 5.0, 0.0))
    beam(d, s1 + s2, 'RECT', '2., 1.')
    d.model += ['*Equation', '2', '%d, 3, 1., %d, 3, -1.' % (b1[-1], b2[-1])]
    QS_DURATION['equation'] = 0.05
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % b1[0], '%d, ENCASTRE' % b2[0], '*Cload',
                                    '%d, 3, -1.' % b1[-1]], False, rot=True, nframes=1)
    register(d, 'static', 'B31', '*EQUATION', 'Two parallel cantilevers, linear constraint u3(tip 1) = u3(tip 2), '
             'load on tip 1 only: each beam takes half.', [b1[-1], b2[-1]], [b1[0], b2[0]], True, 'theory', 'linear',
             category='feature', group='constraint', supported=True, reference='analytic')

    d = Deck('rigid_body', 'Rigid C3D8R block on four axial springs, eccentric load on its reference node')
    blk = hex_block(d, 2, 2, 1, 4.0, 4.0, 2.0, x0=(-2.0, -2.0, 10.0))
    solid(d, 'C3D8R', blk, 'STEEL', 'BLOCK')
    rp = d.node(1.0, 0.0, 11.0)
    d.model += ['*Rigid Body, ref node=%d, elset=BLOCK' % rp]
    corners = d.select(lambda x, y, z: near(z, 10.0) and abs(abs(x) - 2.0) < 1e-9 and abs(abs(y) - 2.0) < 1e-9)
    ground = [d.node(d.nodes[n][0], d.nodes[n][1], 0.0) for n in corners]
    d.elements('CONN3D2', 'CONN', [[g, n] for g, n in zip(ground, corners)])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL', '*Connector Behavior, name=SPRING',
                '*Connector Elasticity, component=1', '100.']
    d.nset('GROUND', ground)
    QS_DURATION['rigid_body'] = 0.02
    d.steps += static_step('Load', ['*Boundary', 'GROUND, 1, 3', '%d, 1, 2' % rp, '%d, 6, 6' % rp,
                                    '*Cload', '%d, 3, -4.' % rp], True, rot=True)
    register(d, 'static', 'C3D8R + CONN3D2', '*RIGID BODY', 'Rigid block (reference node, eccentric load) supported '
             'by four axial springs: translation and tilt.', [rp] + corners[:2], ground, True, 'exact',
             category='feature', group='constraint', supported=True, twin=False,
             reference='analytic')

    d = Deck('connector_join', 'Two B31 beams pinned by a CONN3D2 JOIN connector, far end on a roller')
    s1, b1 = beam_line(d, 40, 20.0)
    s2, b2 = beam_line(d, 40, 20.0, x0=(20.0, 0.0, 0.0))
    b2[0] = d.node(20.0, 0.0, 0.0)  # separate node at the joint
    s2[0][0] = b2[0]
    beam(d, s1 + s2, 'RECT', '2., 1.')
    d.elements('CONN3D2', 'JOINT', [[b1[-1], b2[0]]])
    d.model += ['*Connector Section, elset=JOINT', 'JOIN']
    QS_DURATION['connector_join'] = 0.05
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % b1[0], '%d, 3, 3' % b2[-1], '%d, 2, 2' % b2[-1],
                                    '%d, 4, 4' % b2[-1], '*Cload', '%d, 3, -1.' % b1[-1]], False, rot=True, nframes=1)
    register(d, 'static', 'CONN3D2 JOIN + B31', 'connector JOIN', 'Cantilever pinned to a second beam by a JOIN '
             'connector (translations tied, rotations free), load at the pin.', [b1[-1], b2[0], b2[20]],
             [b1[0], b2[-1]], True, 'theory', 'linear', category='feature', group='constraint', reference='analytic')

    d = Deck('connector_hinge', 'Two B31 beams joined by a CONN3D2 HINGE connector')
    s1, b1 = beam_line(d, 5, 20.0)
    s2, b2 = beam_line(d, 5, 20.0, x0=(20.0, 0.0, 0.0))
    b2[0] = d.node(20.0, 0.0, 0.0)
    s2[0][0] = b2[0]
    beam(d, s1 + s2, 'RECT', '2., 1.')
    d.elements('CONN3D2', 'JOINT', [[b1[-1], b2[0]]])
    d.model += ['*Orientation, name=HO', '1., 0., 0., 0., 1., 0.', '3, 0.',
                '*Connector Section, elset=JOINT', 'HINGE', 'HO']
    QS_DURATION['connector_hinge'] = 0.05
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % b1[0], '%d, 3, 4' % b2[-1], '*Cload',
                                    '%d, 2, 0.5' % b2[2], '%d, 3, -1.' % b1[-1]], False, rot=True, nframes=1)
    register(d, 'static', 'CONN3D2 HINGE + B31', 'connector HINGE', 'Two beams joined by a hinge about the global x '
             'axis.', [b1[-1], b2[0], b2[2]], [b1[0], b2[-1]], True, 'linear', 'linear', category='feature',
             group='constraint', supported=True)


# --------------------------------------------------------------------------------------- element types
def tet10(d, tets):
    mids, out = {}, []
    for t in tets:
        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in mids:
                pa, pb = d.nodes[a], d.nodes[b]
                mids[k] = d.node(*[(pa[i] + pb[i]) / 2 for i in range(3)])
            return mids[k]
        out.append(t + [mid(t[0], t[1]), mid(t[1], t[2]), mid(t[2], t[0]), mid(t[0], t[3]), mid(t[1], t[3]),
                        mid(t[2], t[3])])
    return out


def wedges(hexes):
    out = []
    for h in hexes:
        out += [[h[0], h[1], h[2], h[4], h[5], h[6]], [h[0], h[2], h[3], h[4], h[6], h[7]]]
    return out


def element_features():
    def cantilever(name, etype, conns_fn, title):
        d = Deck(name, title)
        hexes = hex_block(d, 10, 2, 2, 40.0, 4.0, 4.0)
        solid(d, etype, conns_fn(d, hexes), 'STEEL')
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, 40.0) and (near(y, 0) or near(y, 2.0) or near(y, 4.0))
                       and (near(z, 0) or near(z, 2.0) or near(z, 4.0)))
        QS_DURATION[name] = 0.03
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 3, %s' % (n, fmt(-200.0 / len(tip))) for n in tip], False, nframes=1)
        register(d, 'static', etype, 'element %s' % etype, 'Cantilever under a 200 N tip load.', tip[:3], fix, False,
                 'linear', 'linear', category='feature', group='element', supported=False,
                 reason=NS + '%s elements' % etype)
    cantilever('elem_c3d10m', 'C3D10M', lambda d, h: tet10(d, tets_from_hexes(d, h)), 'C3D10M cantilever, tip load')
    cantilever('elem_c3d6', 'C3D6', lambda d, h: wedges(h), 'C3D6 cantilever, tip load')

    # C3D8R with the enhanced hourglass control, one element through the depth (bending carried by the
    # hourglass modes), rectangular and distorted elements
    for name, skew, title in (('hourglass_enhanced_c3d8r', 0.0, 'rectangular'),
                              ('hourglass_enhanced_distorted_c3d8r', 1.0, 'distorted')):
        d = Deck(name, 'C3D8R cantilever, one element through the depth, enhanced hourglass control (%s)' % title)
        hexes = hex_block(d, 10, 1, 1, 40.0, 4.0, 4.0)
        if skew:  # shift the inner nodes along x alternately (trapezoidal elements)
            for n, (x, y, z) in list(d.nodes.items()):
                if 0 < x < 40.0:
                    d.nodes[n] = (x + skew * (1 if near(z, 4.0) else -1) * (1 if round(x / 4.0) % 2 else -1), y, z)
        d.elements('C3D8R', 'EALL', hexes)
        d.materials.add('STEEL')
        d.model += ['*Solid Section, elset=EALL, material=STEEL, controls=CENH', ',',
                    '*Section Controls, name=CENH, hourglass=ENHANCED']
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, 40.0))
        QS_DURATION[name] = 0.03
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 3, %s' % (n, fmt(-10.0 / len(tip))) for n in tip], False, nframes=1)
        register(d, 'static', 'C3D8R', '*SECTION CONTROLS, HOURGLASS=ENHANCED',
                 'Cantilever of ten %s C3D8R elements, one through the depth, tip load: bending carried by the '
                 'hourglass modes with the enhanced (physical) hourglass stiffness.' % title, tip, fix, False,
                 'linear', 'linear', category='feature', group='element', twin=False)

    d = Deck('elem_t3d2', 'T3D2 truss tripod loaded at its apex')
    base = [d.node(10.0 * math.cos(a), 10.0 * math.sin(a), 0.0) for a in (0.0, 2.0944, 4.18879)]
    apex = d.node(0.0, 0.0, 10.0)
    d.elements('T3D2', 'TRUSS', [[b, apex] for b in base])
    d.materials.add('STEEL')
    d.model += ['*Solid Section, elset=TRUSS, material=STEEL', '1.']
    d.nset('BASE', base)
    QS_DURATION['elem_t3d2'] = 0.002
    d.steps += static_step('Load', ['*Boundary', 'BASE, 1, 3', '*Cload', '%d, 3, -1000.' % apex, '%d, 1, 200.' % apex],
                           False, nframes=1)
    register(d, 'static', 'T3D2', 'element T3D2', 'Three-bar truss tripod loaded at its apex.', [apex], base, False,
             'linear', 'linear', category='feature', group='element', supported=True)

    d = Deck('elem_s4', 'S4 clamped square plate under pressure')
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    d.elements('S4', 'EALL', quads)
    d.materials.add('STEEL')
    d.model += ['*Shell Section, elset=EALL, material=STEEL', '0.2, 5', '*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
    edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
    QS_DURATION['elem_s4'] = 0.002
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.5'], True, rot=True)
    register(d, 'static', 'S4', 'element S4', 'Clamped square plate of fully integrated shells under follower '
             'pressure, large deflection.', [ids[4, 4], ids[2, 2]], edge, True, 'nonlinear', category='feature',
             group='element')

    d = Deck('elem_m3d4r', 'M3D4R membrane square under biaxial in-plane stretch')
    quads, ids = quad_plate(d, 4, 4, 20.0, 20.0)
    d.elements('M3D4R', 'EALL', quads)
    d.materials.add('STEEL')
    d.model += ['*Membrane Section, elset=EALL, material=STEEL', '0.2']
    x0 = d.select(lambda x, y, z: near(x, 0))
    d.nset('NALL', list(d.nodes)); d.nset('Y0', d.select(lambda x, y, z: near(y, 0)))
    d.nset('X1', d.select(lambda x, y, z: near(x, 20.0))); d.nset('Y1', d.select(lambda x, y, z: near(y, 20.0)))
    QS_DURATION['elem_m3d4r'] = 0.001
    d.steps += static_step('Stretch', ['*Boundary', 'NALL, 3, 3', 'FIX, 1, 1', 'Y0, 2, 2', '*Boundary', 'X1, 1, 1, 0.02',
                                       'Y1, 2, 2, 0.01'], False, nframes=1)
    register(d, 'static', 'M3D4R', 'element M3D4R', 'Membrane square under homogeneous biaxial stretch (0.1 % x '
             '0.05 %).', [ids[4, 4], ids[2, 2]], x0, False, 'linear', 'linear', category='feature', group='element',
             supported=False, reason=NS + 'M3D4R membrane elements')

    d = Deck('elem_b32', 'B32 cantilever, tip load (Abaqus/Standard only)')
    nodes = [d.node(5.0 * i, 0.0, 0.0) for i in range(11)]
    d.elements('B32', 'EALL', [[nodes[2 * i], nodes[2 * i + 1], nodes[2 * i + 2]] for i in range(5)])
    d.materials.add('STEEL')
    d.model += ['*Beam Section, elset=EALL, material=STEEL, section=RECT', '2., 1.', '0., 0., -1.']
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % nodes[0], '*Cload', '%d, 3, -2.' % nodes[-1]], False,
                           rot=True, nframes=1)
    register(d, 'static', 'B32', 'element B32', 'Quadratic beam cantilever under a tip load.', [nodes[-1], nodes[5]],
             [nodes[0]], True, 'linear', 'linear', category='feature', group='element', supported=False,
             reason=NS + 'B32 elements', twin=False)

    d = Deck('elem_spring_dashpot_mass', 'Point mass on a SPRINGA spring and a DASHPOTA damper (explicit)')
    g, m = d.node(0.0, 0.0, 0.0), d.node(10.0, 0.0, 0.0)
    d.elements('SPRINGA', 'ESP', [[g, m]])
    d.elements('DASHPOTA', 'EDP', [[g, m]])
    d.elements('MASS', 'EM', [[m], [g]])  # a mass at the fixed end too (Abaqus/Explicit checks both spring nodes)
    d.model += ['*Spring, elset=ESP', '', '100.', '*Dashpot, elset=EDP', '', '0.002', '*Mass, elset=EM', '1e-04']
    d.initial += ['*Initial Conditions, type=VELOCITY', '%d, 1, 100.' % m]
    d.steps += explicit_step('Oscillate', 0.02, ['*Boundary', '%d, 1, 3' % g, '%d, 2, 3' % m], fixed_dt=1e-05)
    register(d, 'explicit', 'SPRINGA + DASHPOTA + MASS', 'spring, dashpot, point mass', 'Damped mass-spring oscillator '
             'released with 100 mm/s.', [m], [g], False, 'dynamic', category='feature', group='element',
             supported=True)

    # single S4R elements, three corners clamped, the fourth loaded: each probe isolates one part of the
    # element stiffness (drilling, transverse shear with its stabilised mode, membrane hourglass, rotation
    # hourglass) on rectangles, a trapezoid and a distorted quad, thin and thick
    d = Deck('s4r_single_element_probes', 'Single S4R elements, one free corner loaded (stiffness probes)')
    probes = (  # corners, thickness, load on the free corner (dof, value)
        ([(0, 0), (2, 0), (2, 1), (0, 1)], 0.5, (6, 1.0)),       # drilling, thick rectangle
        ([(0, 0), (2, 0), (2, 1), (0, 1)], 0.05, (6, 0.001)),    # drilling, thin rectangle
        ([(0, 0), (2, 0), (1.5, 1), (0.5, 1)], 0.1, (3, 1.0)),   # transverse shear, trapezoid
        ([(0, 0), (2.2, -0.3), (2.5, 1.4), (-0.2, 1.1)], 1.0, (3, 100.0)),  # transverse shear, distorted
        ([(0, 0), (4, 0), (4, 2), (0, 2)], 0.2, (2, 10.0)),      # in-plane load across the element: hourglass
        ([(0, 0), (1, 0), (1, 4), (0, 4)], 0.1, (4, 0.1)),       # bending moment, slender rectangle
        ([(0, 0), (2, 0), (2.8, 1), (0.8, 1)], 0.3, (5, 1.0)))   # bending moment, parallelogram
    loads, mon, fix = [], [], []
    for k, (corners, t, (dof, val)) in enumerate(probes):
        ids = [d.node(x + 10.0 * k, y, 0.0) for (x, y) in corners]
        shell(d, 'S4R', [ids], t, 'STEEL', elset='P%d' % (k + 1))
        fix += [ids[0], ids[1], ids[3]]
        mon.append(ids[2])
        loads.append('%d, %d, %s' % (ids[2], dof, fmt(val)))
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] + loads, False, rot=True, nframes=1)
    register(d, 'static', 'S4R', 'element stiffness probes', 'Seven single S4R elements (rectangles, parallelogram, '
             'trapezoid, distorted quad; t/L from 0.025 to 1) with three corners clamped and a force or moment on '
             'the fourth: drilling, transverse shear, membrane and rotation hourglass stiffness.', mon, fix, True,
             'linear', 'linear', category='feature', group='element', twin=False)

    # in-plane dynamics of a single row of S4R elements: carried by the membrane hourglass and drilling
    # stiffnesses (explicit values of the membrane hourglass factor)
    d = Deck('s4r_inplane_cantilever_xpl', 'S4R strip, in-plane tip load (explicit): membrane hourglass')
    ids = [[d.node(4.0 * i, 2.0 * j, 0.0) for j in range(2)] for i in range(11)]
    shell(d, 'S4R', [[ids[i][0], ids[i + 1][0], ids[i + 1][1], ids[i][1]] for i in range(10)], 0.2, 'STEEL')
    d.nset('NALL', [n for row in ids for n in row])
    fix, tip = [ids[0][0], ids[0][1]], [ids[10][0], ids[10][1]]
    amp = _smooth_amplitude(d, 0.002)
    d.steps += explicit_step('Load', 0.004, ['*Boundary', 'FIX, ENCASTRE', 'NALL, 3, 5', '*Cload, amplitude=%s' % amp,
                                             '%d, 2, 0.05' % tip[0], '%d, 2, 0.05' % tip[1]], rot=True, nframes=40)
    register(d, 'explicit', 'S4R', 'in-plane bending, membrane hourglass', 'Strip of ten 4 x 2 S4R elements, clamped, '
             'in-plane tip load ramped over 2 ms then held: the in-plane bending of a single row of elements is carried '
             'by the membrane hourglass and drilling stiffnesses.', tip, fix, True, 'dynamic', category='feature',
             group='element', twin=False)

    # relaxed hourglass stiffness of C3D8R in explicit dynamics (Abaqus/Explicit default): a single layer of bricks
    # bends through its hourglass modes; the tip load is held while the stiffness relaxes to 0.1 of its initial value
    d = Deck('hourglass_relax_c3d8r_xpl', 'C3D8R single-layer cantilever, tip load held (explicit): relaxed hourglass')
    nx = 10
    ids = [[[d.node(2.0 * i, 2.0 * j, 2.0 * k) for k in range(2)] for j in range(2)] for i in range(nx + 1)]
    solid(d, 'C3D8R', [[ids[i][0][0], ids[i + 1][0][0], ids[i + 1][1][0], ids[i][1][0],
                        ids[i][0][1], ids[i + 1][0][1], ids[i + 1][1][1], ids[i][1][1]] for i in range(nx)], 'STEEL')
    fix = [ids[0][j][k] for j in range(2) for k in range(2)]
    tip = [ids[nx][j][k] for j in range(2) for k in range(2)]
    amp = _smooth_amplitude(d, 1e-4)
    d.steps += explicit_step('Load', 1.5e-3, ['*Boundary', 'FIX, ENCASTRE', '*Cload, amplitude=%s' % amp] +
                             ['%d, 3, -5.' % n for n in tip], nframes=60)
    register(d, 'explicit', 'C3D8R', 'relaxed hourglass stiffness', 'Cantilever of ten 2 mm C3D8R bricks, one element '
             'through the depth (bending carried by the hourglass modes), tip load ramped over 0.1 ms then held for '
             '1.4 ms (about six relaxation times): the deflection grows as the hourglass stiffness relaxes.', tip, fix,
             False, 'dynamic', category='feature', group='element', twin=False)

    # hourglass control of C3D8R on stretched bricks (explicit): every node follows a homogeneous stretch lambda
    # along z (lateral 1/sqrt(lambda)) except one corner, free along x and pushed sideways once the stretch is held:
    # its motion is resisted by the uniform-strain and hourglass stiffnesses of the current (stretched) brick. Every
    # other node prescribed, slow push: tolerance of the static Abaqus-referenced cases (the hourglass stiffness of
    # the undeformed brick gives 2.4 % on U)
    d = Deck('hourglass_stretched_c3d8r_xpl', 'Stretched C3D8R bricks, corner pushed sideways (explicit): hourglass')
    corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    variants = [(0.5, 'STEEL', 2500.), (0.75, 'STEEL', 4800.), (1.5, 'STEEL', 6900.), (0.5, 'COPPER', 600.)]
    body, mon, fix, conns = [], [], [], {'STEEL': [], 'COPPER': []}
    d.model += ['*Amplitude, name=PUSH, definition=SMOOTH STEP', '0., 0., 0.004, 0., 0.005, 1.']
    for k, (lam, mat, push) in enumerate(variants):
        ids = [d.node(3.0 * k + x, y, z) for x, y, z in corners]
        conns[mat].append(ids)
        lt = 1.0 / math.sqrt(lam)
        for n, (x, y, z) in zip(ids, corners):
            u = [(lt - 1) * (x - 0.5), (lt - 1) * (y - 0.5), (lam - 1) * z]
            free = (x, y, z) == (1, 1, 1)
            (mon if free else fix).append(n)
            for j in range(3):
                if free and j == 0:
                    continue
                name = 'S%dD%d' % (n, j + 1)
                d.model += ['*Amplitude, name=%s, definition=SMOOTH STEP' % name, '0., 0., 0.004, %s, 0.008, %s' % (fmt(u[j]), fmt(u[j]))]
                body += ['*Boundary, amplitude=%s' % name, '%d, %d, %d, 1.' % (n, j + 1, j + 1)]
        body += ['*Cload, amplitude=PUSH', '%d, 1, %s' % (ids[6], fmt(push))]
    solid(d, 'C3D8R', conns['STEEL'], 'STEEL', elset='EST')
    solid(d, 'C3D8R', conns['COPPER'], 'COPPER', elset='ECU')
    d.steps += explicit_step('Push', 0.008, body, nframes=40)
    register(d, 'explicit', 'C3D8R', 'hourglass of stretched bricks', 'Four 1 mm C3D8R bricks stretched homogeneously '
             'to lambda = 0.5, 0.75 and 1.5 (elastic steel) and 0.5 (elastic-plastic copper) over 4 ms, every node '
             'prescribed except one corner, free along x and pushed sideways (about 0.2 mm) once the stretch is held: the '
             'relaxed hourglass stiffness follows the current shape of the brick.', mon, fix, False, 'linear',
             category='feature', group='element', twin=False)

    # in-plane bending of a row of S3R triangles in explicit dynamics: drilling of the triangles (merged corner of
    # the degenerate quadrilateral counted twice, rotations against the element spin)
    d = Deck('s3r_inplane_cantilever_xpl', 'S3R strip, in-plane tip load (explicit): triangle drilling')
    ids = [[d.node(4.0 * i, 2.0 * j, 0.0) for j in range(2)] for i in range(11)]
    conns = []
    for i in range(10):
        conns += [[ids[i][0], ids[i + 1][0], ids[i + 1][1]], [ids[i][0], ids[i + 1][1], ids[i][1]]]
    shell(d, 'S3R', conns, 0.2, 'STEEL')
    d.nset('NALL', [n for row in ids for n in row])
    fix, tip = [ids[0][0], ids[0][1]], [ids[10][0], ids[10][1]]
    amp = _smooth_amplitude(d, 0.002)
    d.steps += explicit_step('Load', 0.004, ['*Boundary', 'FIX, ENCASTRE', 'NALL, 3, 5', '*Cload, amplitude=%s' % amp,
                                             '%d, 2, 0.5' % tip[0], '%d, 2, 0.5' % tip[1]], rot=True, nframes=40)
    register(d, 'explicit', 'S3R', 'in-plane bending, drilling', 'Strip of twenty S3R triangles, clamped, in-plane tip '
             'load ramped over 2 ms then held: the in-plane bending of the constant-strain triangles involves their '
             'drilling stiffness.', tip, fix, True, 'dynamic', category='feature', group='element', twin=False)

    # nodal masses and rotary inertias of shells: free single elements under constant nodal forces and moments, one
    # increment of 1e-9 (no internal force yet): u = a dt^2 / 2 with a = load / nodal mass or inertia
    d = Deck('shell_rotary_inertia_xpl', 'Free S3R and S4R elements, nodal forces and moments (explicit): lumped inertia')
    shapes = (('S3R', [(0, 0), (1, 0), (0, 1)]), ('S3R', [(0, 0), (1, 0), (0.5, 0.866)]), ('S3R', [(0, 0), (2, 0), (0, 1)]),
              ('S3R', [(0, 0), (4, 0), (1, 2)]), ('S3R', [(0, 0), (3, 0), (2.5, 0.5)]),
              ('S4R', [(0, 0), (1, 0), (1, 1), (0, 1)]), ('S4R', [(0, 0), (2, 0), (2, 1), (0, 1)]),
              ('S4R', [(0, 0), (3, 0), (2, 1), (0.5, 1)]), ('S4R', [(0, 0), (2, 0), (3, 1), (1, 1)]))
    loads, mon, k = [], [], 0
    for t in (0.01, 0.2, 2.0):
        for etype, corners in shapes:
            ids = [d.node(x + 10.0 * k, y, 0.0) for (x, y) in corners]
            shell(d, etype, [ids], t, 'STEEL', elset='P%d' % (k + 1))
            k += 1
            n = len(corners)
            area = 0.5 * abs(sum(corners[i][0] * corners[(i + 1) % n][1] - corners[(i + 1) % n][0] * corners[i][1]
                                 for i in range(n)))
            lmin2 = min((corners[(i + 1) % n][0] - corners[i][0]) ** 2 + (corners[(i + 1) % n][1] - corners[i][1]) ** 2
                        for i in range(n))
            mass = 7.85e-9 * t * area / n
            inertia = mass * (area / 12.0 + max(lmin2 / 16.0, t * t / 6.0))
            for nid in ids:  # loads giving displacements and rotations of about 0.01 after the increment
                loads += ['%d, 1, %s' % (nid, fmt(2e16 * mass)), '%d, 3, %s' % (nid, fmt(-2e16 * mass))]
                loads += ['%d, %d, %s' % (nid, dof, fmt(2e16 * inertia)) for dof in (4, 5, 6)]
            mon += ids
    d.steps += explicit_step('Load', 1e-9, ['*Cload'] + loads, rot=True, nframes=1, fixed_dt=1e-9)
    register(d, 'explicit', 'S4R + S3R', 'lumped mass and rotary inertia', 'Twenty-seven free single S3R and S4R '
             'elements (triangles and quadrilaterals of several shapes, thickness 0.01, 0.2 and 2) under constant '
             'nodal forces and moments, one increment of 1e-9 s: the motion is set by the lumped nodal masses and '
             'rotary inertias alone.', mon, [], True, 'linear', 'linear', category='feature', group='element',
             twin=False)

    # shell bulk viscosity: free single elements with all nodal velocities prescribed (in-plane dilatation or
    # uniform curvature rate), viscous dissipation ALLVD compared. Rates chosen so that every element dissipates
    # about the same energy (small strains: 1e-3 at most)
    def bulk_viscosity_case(name, title, items, element, description):
        d = Deck(name, title)
        bcs, mon = [], []
        for k, (etype, corners, mat, mode) in enumerate(items):
            n = len(corners)
            cx, cy = sum(c[0] for c in corners) / n, sum(c[1] for c in corners) / n
            area = 0.5 * abs(sum(corners[i][0] * corners[(i + 1) % n][1] - corners[(i + 1) % n][0] * corners[i][1]
                                 for i in range(n)))
            rho, modulus = {'STEEL': (7.85e-9, 2.3e5), 'STEEL49': (7.85e-9, 2.76e5), 'NH45': (1.1e-9, 729.0),
                            'NH499': (1.1e-9, 7896.0), 'FUNG': (1.06e-9, 7896.0)}[mat]
            rate = 10.0 * (0.0425 / math.sqrt(rho * modulus)) ** 0.5 / area ** 0.75
            if mode == 'bend':
                rate *= 35.0
            ids = [d.node(x - cx + 10.0 * k, y - cy, 0.0) for (x, y) in corners]
            elset = 'P%d' % (k + 1)
            shell(d, etype, [ids], 0.1, mat, elset=elset)
            if mat == 'FUNG':
                d.model[-2] += ', orientation=ORI'
                if '*Orientation, name=ORI, local directions=2' not in d.model:
                    d.model[:0] = ['*Orientation, name=ORI, local directions=2', '1., 0., 0., 0., 1., 0.']
            for nid, (x, y) in zip(ids, corners):
                x, y = x - cx, y - cy
                vals = (rate * x, rate * y, 0, 0, 0, 0) if mode == 'dil' else (0, 0, 0, rate * y, -rate * x, 0)
                bcs += ['%d, %d, %d, %s' % (nid, dof + 1, dof + 1, fmt(v)) for dof, v in enumerate(vals)]
            mon += ids
        d.model += ['*Amplitude, name=CONSTANT', '0., 1., 1., 1.']
        d.steps += explicit_step('Move', 5e-6, ['*Boundary, type=VELOCITY, amplitude=CONSTANT'] + bcs, rot=True,
                                 nframes=10)
        register(d, 'explicit', element, 'bulk viscosity', description, mon, [], True, 'linear', category='feature',
                 group='element', twin=False, energies=['ALLVD'])

    sq, rect2, para, rect4 = ([(0, 0), (1, 0), (1, 1), (0, 1)], [(0, 0), (2, 0), (2, 1), (0, 1)],
                              [(0, 0), (1, 0), (1.6, 1), (0.6, 1)], [(0, 0), (4, 0), (4, 1), (0, 1)])
    rt1, equi, obt, sliv = ([(0, 0), (1, 0), (0, 1)], [(0, 0), (1, 0), (0.5, 0.8660254)], [(0, 0), (3, 0), (1, 0.5)],
                            [(0, 0), (5, 0), (2.5, 0.3)])
    bulk_viscosity_case('shell_bulk_viscosity_s4r_xpl', 'Free S4R elements, prescribed velocities (explicit): bulk '
                        'viscosity', [('S4R', sq, 'STEEL', 'dil'), ('S4R', para, 'STEEL', 'dil'),
                                      ('S4R', rect2, 'STEEL', 'bend'), ('S4R', para, 'STEEL', 'bend')], 'S4R',
                        'Four free S4R elements (square, parallelogram, rectangle) with all nodal velocities '
                        'prescribed: constant in-plane dilatation rate or constant curvature rate. Viscous dissipation '
                        'of the linear bulk viscosity (membrane pressure and pressure moment) compared.')
    bulk_viscosity_case('shell_bulk_viscosity_s3r_xpl', 'Free S3R elements, prescribed velocities (explicit): bulk '
                        'viscosity', [('S3R', rt1, 'STEEL', 'dil'), ('S3R', sliv, 'STEEL', 'dil'),
                                      ('S3R', equi, 'STEEL', 'bend'), ('S3R', obt, 'STEEL', 'bend')], 'S3R',
                        'Four free S3R elements (right, slender, equilateral and obtuse triangles) with all nodal '
                        'velocities prescribed: constant in-plane dilatation rate or constant curvature rate. Viscous '
                        'dissipation of the linear bulk viscosity compared.')
    bulk_viscosity_case('shell_bulk_viscosity_hyper_xpl', 'Free S4R and S3R elements, compressible and nearly '
                        'incompressible materials (explicit): bulk viscosity',
                        [('S4R', sq, 'NH45', 'dil'), ('S4R', sq, 'NH499', 'bend'), ('S3R', rt1, 'FUNG', 'dil'),
                         ('S4R', rect4, 'STEEL49', 'bend')], 'S4R + S3R',
                        'Free S4R and S3R elements with all nodal velocities prescribed, neo-Hooke (nu = 0.45 and '
                        '0.499), Fung (nearly incompressible) and elastic (nu = 0.49) materials: wave speed of the '
                        'bulk viscosity, with the Poisson ratio limited to 0.49 at fixed bulk modulus.')

    # S3RS: same nodal masses as S3R, rotary inertia mass times max(A / 4, t^2 / 6)
    d = Deck('s3rs_rotary_inertia_xpl', 'Free S3RS elements, nodal forces and moments (explicit): lumped inertia')
    loads, mon, k = [], [], 0
    for t in (0.01, 0.2, 2.0):
        for etype, corners in shapes:
            if etype != 'S3R':
                continue
            ids = [d.node(x + 10.0 * k, y, 0.0) for (x, y) in corners]
            shell(d, 'S3RS', [ids], t, 'STEEL', elset='P%d' % (k + 1))
            k += 1
            area = 0.5 * abs((corners[1][0] - corners[0][0]) * (corners[2][1] - corners[0][1]) -
                             (corners[2][0] - corners[0][0]) * (corners[1][1] - corners[0][1]))
            mass = 7.85e-9 * t * area / 3
            inertia = mass * max(area / 4.0, t * t / 6.0)
            for nid in ids:
                loads += ['%d, 1, %s' % (nid, fmt(2e16 * mass)), '%d, 3, %s' % (nid, fmt(-2e16 * mass))]
                loads += ['%d, %d, %s' % (nid, dof, fmt(2e16 * inertia)) for dof in (4, 5, 6)]
            mon += ids
    d.steps += explicit_step('Load', 1e-9, ['*Cload'] + loads, rot=True, nframes=1, fixed_dt=1e-9)
    register(d, 'explicit', 'S3RS', 'lumped mass and rotary inertia', 'Fifteen free single S3RS elements (five '
             'triangles, thickness 0.01, 0.2 and 2) under constant nodal forces and moments, one increment of 1e-9 s: '
             'the motion is set by the lumped nodal masses and rotary inertias alone.', mon, [], True, 'linear',
             'linear', category='feature', group='element', twin=False)

    # nodal rotary inertia of B31 with ROTARY INERTIA=ISOTROPIC (one value about the three axes): free single
    # elements under constant nodal forces and moments, one increment of 1e-9
    d = Deck('beam_rotary_inertia_xpl', 'Free B31 elements, nodal forces and moments (explicit): lumped inertia')
    secs = (('CIRC', '1.'), ('RECT', '2., 1.'), ('PIPE', '3., 0.5'))
    loads, mon, k = [], [], 0
    for iso in (True,):
        for sec, dims in secs:
            for length in (1.0, 5.0, 20.0):
                a, b = d.node(0.0, 50.0 * k, 0.0), d.node(length, 50.0 * k, 0.0)
                d.elements('B31', 'B%d' % (k + 1), [[a, b]])
                d.model += ['*Beam Section, elset=B%d, material=STEEL, section=%s%s' % (
                    k + 1, sec, ', rotary inertia=ISOTROPIC' if iso else ''), dims, '0., 0., -1.']
                k += 1
                for nid in (a, b):
                    loads += ['%d, 2, 1e+07' % nid, '%d, 4, 1e+08' % nid, '%d, 5, -1e+08' % nid, '%d, 6, 1e+08' % nid]
                mon += [a, b]
    d.materials.add('STEEL')
    d.steps += explicit_step('Load', 1e-9, ['*Cload'] + loads, rot=True, nframes=1, fixed_dt=1e-9)
    register(d, 'explicit', 'B31', 'lumped rotary inertia', 'Nine free single B31 elements (CIRC, RECT, PIPE; '
             'length 1, 5 and 20; ROTARY INERTIA=ISOTROPIC) under constant nodal forces and moments, one increment of '
             '1e-9 s: the rotations are set by the nodal rotary inertias alone.', mon, [], True, 'linear',
             'linear', category='feature', group='element', twin=False)

    # default (exact) rotary inertia of B31, anisotropic: rho L / 2 (I11 + I22) about the beam axis. Free elements
    # twisted by equal moments about their axis at both nodes (rigid rotation, no strain): theta = M t^2 / (2 J)
    d = Deck('beam_rotary_inertia_torsion_xpl', 'Free B31 elements twisted by axial moments (explicit): default '
             'rotary inertia about the beam axis')
    loads, mon = [], []
    for k, (sec, dims, p0, p1, moment) in enumerate(BEAM_TWIST):
        a, b = d.node(*p0), d.node(*p1)
        d.elements('B31', 'B%d' % (k + 1), [[a, b]])
        d.model += ['*Beam Section, elset=B%d, material=STEEL, section=%s' % (k + 1, sec), dims, '0., 0., -1.']
        length = math.sqrt(sum((q - p) ** 2 for p, q in zip(p0, p1)))
        axis = [(q - p) / length for p, q in zip(p0, p1)]
        for nid in (a, b):
            loads += ['%d, %d, %s' % (nid, 4 + j, fmt(moment * axis[j])) for j in range(3) if abs(axis[j]) > 1e-12]
        mon += [a, b]
    d.materials.add('STEEL')
    d.steps += explicit_step('Load', 1e-8, ['*Cload'] + loads, rot=True, nframes=10, fixed_dt=1e-9)
    register(d, 'explicit', 'B31', 'anisotropic rotary inertia', 'Three free B31 elements (RECT along x, PIPE, RECT '
             'along a skew axis) twisted by equal moments about their axis at both nodes: rigid rotation set by the '
             'default rotary inertia about the beam axis, rho L / 2 (I11 + I22) per node.', mon, [], True, 'exact',
             'linear', category='feature', group='element', twin=False, reference='analytic')

    # the same nine elements as beam_rotary_inertia_xpl with the default inertia, rho L / 2 (I11 + I22, I11, I22)
    # about (axis, n1, n2): smooth moments (no load in the first increment, whose rotation Abaqus/Explicit does not
    # apply to these nodes), forty increments of 1e-9, rotations of a few hundredths of a radian
    d = Deck('beam_rotary_inertia_exact_xpl', 'Free B31 elements, smooth nodal moments (explicit): default anisotropic '
             'rotary inertia')
    amp = _smooth_amplitude(d, 2e-8)
    loads, mon, k = [], [], 0
    for sec, dims in secs:
        for length in (1.0, 5.0, 20.0):
            a, b = d.node(0.0, 50.0 * k, 0.0), d.node(length, 50.0 * k, 0.0)
            d.elements('B31', 'B%d' % (k + 1), [[a, b]])
            d.model += ['*Beam Section, elset=B%d, material=STEEL, section=%s' % (k + 1, sec), dims, '0., 0., -1.']
            k += 1
            # moments scaled with the largest flexural inertia: angular accelerations of about 5e13 rad/s^2
            mref = 5e13 * 7.85e-9 * length / 2 * BEAM_IMAX[sec]
            for nid in (a, b):
                loads += ['%d, 4, %s' % (nid, fmt(mref)), '%d, 5, %s' % (nid, fmt(-0.5 * mref)),
                          '%d, 6, %s' % (nid, fmt(0.7 * mref))]
            mon += [a, b]
    d.materials.add('STEEL')
    d.steps += explicit_step('Load', 4e-8, ['*Cload, amplitude=%s' % amp] + loads, rot=True, nframes=8,
                             fixed_dt=1e-9)
    register(d, 'explicit', 'B31', 'anisotropic rotary inertia', 'Nine free single B31 elements (CIRC, RECT, PIPE; '
             'length 1, 5 and 20; default rotary inertia) under smooth nodal moments about the three axes: '
             'anisotropic nodal inertia rho L / 2 (I11 + I22, I11, I22) carried by the beam frame, with the '
             'gyroscopic moment.', mon, [], True, 'linear', 'linear', category='feature', group='element', twin=False)


# ------------------------------------------------------------------------------------- loads and conditions
def cantilever_beam(name, title, n=10, section='RECT', dims='2., 1.'):  # noqa: E302
    d = Deck(name, title)
    segs, ids = beam_line(d, n, 100.0)
    beam(d, segs, section, dims)
    return d, ids


def load_features():
    d, ids = cantilever_beam('load_gravity', 'B31 cantilever under its own weight (*DLOAD GRAV)', n=50)
    QS_DURATION['load_gravity'] = 0.03
    d.steps += static_step('Weight', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Dload', 'EALL, GRAV, 98100., 0., 0., -1.'],
                           False, rot=True, nframes=1)
    register(d, 'static', 'B31', 'gravity', 'Cantilever under 10 g gravity.', [ids[-1], ids[25]], [ids[0]], True,
             'theory', 'linear', category='feature', group='load', supported=True, reference='analytic')

    d, ids = cantilever_beam('load_follower', 'B31 cantilever under a follower tip force (large deflection)')
    QS_DURATION['load_follower'] = 0.03
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload, follower', '%d, 3, -0.4' % ids[-1]],
                           True, rot=True)
    register(d, 'static', 'B31', 'follower force', 'Cantilever under a tip force that follows the rotation of the '
             'tip (large deflection).', [ids[-1], ids[5]], [ids[0]], True, 'nonlinear', category='feature',
             group='load', supported=True)

    d, ids = cantilever_beam('bc_transform', 'B31 cantilever, tip guided along an inclined direction (*TRANSFORM)',
                             n=50)
    c, s = math.cos(math.radians(30)), math.sin(math.radians(30))
    d.nset('TIPT', [ids[-1]])
    d.model += ['*Transform, nset=TIPT, type=R', '%s, 0., %s, 0., 1., 0.' % (fmt(c), fmt(s))]
    QS_DURATION['bc_transform'] = 0.03
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Boundary', 'TIPT, 1, 1', '*Cload',
                                    '%d, 3, -2.' % ids[-1]], False, rot=True, nframes=1)
    register(d, 'static', 'B31', 'local boundary conditions', 'Cantilever whose tip is held in a direction inclined at '
             '30 degrees in the x-z plane (nodal transformation), transverse load.', [ids[-1], ids[25]],
             [ids[0], ids[-1]], True, 'theory', 'linear', category='feature', group='load', supported=True,
             reference='analytic')

    d = Deck('amplitudes', 'S4R strip, tip displacement driven by TABULAR, SMOOTH STEP and total-time amplitudes')
    quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
    shell(d, 'S4R', quads, 0.5, 'STEEL', 'EALL')
    fix = [ids[0, j] for j in range(3)]
    tip = [ids[10, j] for j in range(3)]
    d.nset('TIP', tip); d.nset('NALL', list(d.nodes))
    d.model += ['*Amplitude, name=TAB', '0., 0., 0.3, 1., 0.6, 0.4, 1., 0.8',
                '*Amplitude, name=SMO, definition=SMOOTH STEP', '0., 0., 1., 1.',
                '*Amplitude, name=TOT, time=TOTAL TIME', '0., 0., 1., 0., 2., 1.']
    d.steps += static_step('Tabular', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 1', '%d, 2, 2' % ids[0, 0],
                                       '*Boundary, amplitude=TAB', 'TIP, 1, 1, 0.04'], False, rot=True, inc=0.05)
    # in-plane load spread over the tip edge (a point load in the plane of a membrane is singular)
    d.steps += static_step('Smooth', ['*Boundary, amplitude=SMO', 'TIP, 1, 1, 0.02', '*Cload, amplitude=TOT'] +
                           ['%d, 2, %s' % (n, fmt(100.0 * w)) for n, w in edge_shares(d, tip)], False, rot=True,
                           inc=0.05)
    register(d, 'static', 'S4R', 'amplitudes', 'Tip displacement following a tabular amplitude, then a smooth step; '
             'in-plane load following an amplitude in total time.', tip, fix, True, 'linear', 'linear',
             category='feature', group='load', twin=False)

    d = Deck('amp_periodic', 'S4R strip, tip displacement driven by a PERIODIC amplitude')
    quads, ids = quad_plate(d, 10, 2, 40.0, 10.0)
    shell(d, 'S4R', quads, 0.5, 'STEEL', 'EALL')
    fix = [ids[0, j] for j in range(3)]
    tip = [ids[10, j] for j in range(3)]
    d.nset('TIP', tip); d.nset('NALL', list(d.nodes))
    d.model += ['*Amplitude, name=PER, definition=PERIODIC', '1, 6.283185307, 0., 0.', '0., 1.']
    d.steps += static_step('Periodic', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 1', '%d, 2, 2' % ids[0, 0],
                                        '*Boundary, amplitude=PER', 'TIP, 1, 1, 0.04'], False, rot=True, inc=0.05)
    register(d, 'static', 'S4R', 'periodic amplitude', 'Tip displacement following a sine (periodic amplitude).', tip,
             fix, True, 'linear', 'linear', category='feature', group='load', supported=True, twin=False)

    d, ids = cantilever_beam('multistep_opnew', 'B31 cantilever: load, added support, then loads replaced (OP=NEW)')
    QS_DURATION['multistep_opnew'] = 0.03
    # the tip load removed by OP=NEW drops at once in explicit dynamics: damping (time constant 1 ms, first mode period
    # 0.9 ms) brings the twin to rest before the next frame, so its frames are static states
    QS_DAMPING['multistep_opnew'] = 2000.0
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload', '%d, 3, -1.' % ids[-1]], False,
                           rot=True, nframes=2)
    d.steps += static_step('Support', ['*Boundary', '%d, 3, 3' % ids[5]], False, rot=True, nframes=2)
    d.steps += static_step('Replace', ['*Cload, op=NEW', '%d, 2, 0.5' % ids[-1]], False, rot=True, nframes=2)
    register(d, 'static', 'B31', 'steps, OP=NEW', 'Three steps: tip load; a support added at mid-span (load '
             'propagated); loads replaced by a lateral load (OP=NEW).', [ids[-1], ids[5], ids[8]], [ids[0], ids[5]],
             True, 'linear', 'linear', category='feature', group='load')

    d = Deck('symmetry_bcs', 'Quarter of a clamped S4R plate under pressure (XSYMM / YSYMM)')
    quads, ids = quad_plate(d, 5, 5, 10.0, 10.0)
    shell(d, 'S4R', quads, 0.2, 'STEEL', 'EALL')
    d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
    clamp = [ids[i, j] for (i, j) in ids if i == 0 or j == 0]
    d.nset('XS', [ids[5, j] for j in range(6)]); d.nset('YS', [ids[i, 5] for i in range(6)])
    QS_DURATION['symmetry_bcs'] = 0.002
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', 'XS, XSYMM', 'YS, YSYMM', '*Dsload', 'STOP, P, 0.05'],
                           False, rot=True, nframes=1)
    register(d, 'static', 'S4R', 'symmetry conditions', 'Quarter model of a clamped square plate under pressure, '
             'XSYMM and YSYMM on the symmetry edges (linear).', [ids[5, 5], ids[3, 3], ids[5, 2]], clamp, True,
             'linear', 'linear', category='feature', group='load')

    d = Deck('initial_stress', 'C3D8R block released from an initial stress (*INITIAL CONDITIONS, TYPE=STRESS)')
    hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
    solid(d, 'C3D8R', hexes, 'STEEL', 'EALL')
    fx = d.select(lambda x, y, z: near(x, 0))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.initial += ['*Initial Conditions, type=STRESS', 'EALL, 100., 0., 0., 0., 0., 0.']
    QS_DURATION['initial_stress'] = 0.001
    d.steps += static_step('Release', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3'], False, nframes=1)
    register(d, 'static', 'C3D8R', 'initial stress', 'Block with an initial axial stress of 100 MPa, released on '
             'symmetry supports: it shortens by the corresponding strain.', d.select(lambda x, y, z: near(x, 4.0))[:3],
             fx, False, 'exact', 'linear', category='feature', group='load', supported=True, reference='analytic')

    d = Deck('temperature_dependence', 'C3D8R bar with temperature-dependent modulus, heated under load')
    hexes = hex_block(d, 4, 1, 1, 40.0, 4.0, 4.0)
    solid(d, 'C3D8R', hexes, 'STEEL_TDEP', 'EALL')
    fx = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 40.0))
    d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
    d.nset('NALL', list(d.nodes))
    d.initial += ['*Initial Conditions, type=TEMPERATURE', 'NALL, 0.']
    QS_DURATION['temperature_dependence'] = 0.003
    d.steps += static_step('Load', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Cload'] +
                           ['%d, 1, %s' % (n, fmt(1600.0 / len(tip))) for n in tip], False, nframes=1)
    d.steps += static_step('Heat', ['*Temperature', 'NALL, 200.'], False, nframes=4)
    register(d, 'static', 'C3D8R', 'temperature-dependent properties', 'Bar under a 100 MPa axial stress, then heated '
             'to 200 degrees: the modulus drops from 210 to 150 GPa while the bar expands.', tip, fx, False, 'exact',
             'linear', category='feature', group='load', supported=True, reference='analytic')

    d = Deck('connector_field', 'CONN3D2 AXIAL spring whose stiffness depends on a field variable')
    a, b = d.node(0, 0, 0), d.node(10, 0, 0)
    d.elements('CONN3D2', 'CONN', [[a, b]])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL', '*Connector Behavior, name=SPRING',
                '*Connector Elasticity, component=1, dependencies=1', '100., , , 0.', '400., , , 1.']
    d.nset('NALL', [a, b])
    d.initial += ['*Initial Conditions, type=FIELD, variable=1', 'NALL, 0.']
    d.steps += static_step('Load', ['*Boundary', '%d, 1, 3' % a, '%d, 2, 3' % b, '*Cload', '%d, 1, 50.' % b], False,
                           nframes=1)
    d.steps += static_step('Stiffen', ['*Field, variable=1', 'NALL, 1.'], False, nframes=4)
    register(d, 'static', 'CONN3D2', 'field-dependent connector', 'Axial spring under 50 N whose stiffness goes from '
             '100 to 400 N/mm with a field variable.', [b], [a], False, 'exact', 'linear', category='feature',
             group='load', twin=False, reference='analytic')


# ----------------------------------------------------------------------------------------------- materials
def material_features():
    def block(name, mat, title, stretch):
        d = Deck(name, title)
        hexes = hex_block(d, 2, 2, 2, 4.0, 4.0, 4.0)
        solid(d, 'C3D8R', hexes, mat, 'EALL')
        fx = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, 4.0))
        d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
        d.nset('TIP', tip)
        return d, fx, tip
    d, fx, tip = block('mat_mooney_rivlin', 'MOONEY', 'Mooney-Rivlin block, uniaxial stretch 1.5', 2.0)
    QS_DURATION['mat_mooney_rivlin'] = 0.01
    d.steps += static_step('Stretch', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary', 'TIP, 1, 1, 2.'],
                           True, inc=0.02, maxinc=0.02)
    register(d, 'static', 'C3D8R', 'Mooney-Rivlin', 'Homogeneous uniaxial stretch to 1.5 of a Mooney-Rivlin block.',
             tip[:3], fx, False, 'exact', category='feature', group='material', supported=True, reference='analytic')

    d, fx, tip = block('mat_kinematic_hardening', 'STEEL_KIN', 'Linear kinematic hardening, tension then compression', 0)
    QS_DURATION['mat_kinematic_hardening'] = 0.001
    d.steps += static_step('Tension', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary',
                                       'TIP, 1, 1, 0.08'], True, inc=0.02, maxinc=0.02)
    d.steps += static_step('Compression', ['*Boundary', 'TIP, 1, 1, -0.08'], True, inc=0.02, maxinc=0.02)
    register(d, 'static', 'C3D8R', 'kinematic hardening', 'Homogeneous tension to 2 % then compression to -2 % '
             '(Bauschinger effect).', tip[:3], fx, False, 'exact', category='feature', group='material',
             supported=True, reference='analytic')

    d, fx, tip = block('mat_viscoelastic', 'VISCO', 'Viscoelastic block (Prony series): step strain then relaxation', 0)
    d.model += ['*Amplitude, name=FAST, definition=SMOOTH STEP', '0., 0., 0.0005, 1.']
    d.steps += explicit_step('Relax', 0.01, ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3',
                                             '*Boundary, amplitude=FAST', 'TIP, 1, 1, 0.2'])
    register(d, 'explicit', 'C3D8R', 'viscoelasticity', 'Block stretched by 5 % in 0.5 ms, then held: stress '
             'relaxation with a 1 ms time constant (Prony series).', tip[:3], fx, False, 'dynamic', category='feature',
             group='material', supported=False, reason=NS + 'viscoelasticity')


# ------------------------------------------------------------------------------- procedure-specific options
def wave_bar(name, title):
    d = Deck(name, title)
    hexes = hex_block(d, 10, 2, 2, 40.0, 4.0, 4.0)
    solid(d, 'C3D8R', hexes, 'STEEL', 'EALL')
    fix = d.select(lambda x, y, z: near(x, 0))
    d.nset('MOVING', [n for n in d.nodes if n not in set(fix)])
    d.initial += ['*Initial Conditions, type=VELOCITY', 'MOVING, 1, -1000.']
    mon = d.select(lambda x, y, z: (near(x, 40.0) or near(x, 20.0)) and near(y, 0) and near(z, 0))
    return d, fix, mon


def explicit_features():
    for name, extra, title in (
            ('xpl_mass_scaling_fixed', ['*Fixed Mass Scaling, factor=4.'], 'fixed mass scaling (factor 4)'),
            ('xpl_mass_scaling_variable', ['*Variable Mass Scaling, dt=5e-07, type=below min, frequency=10'],
             'variable mass scaling (target increment 5e-7 s)'),
            ('xpl_bulk_viscosity', ['*Bulk Viscosity', '0.2, 2.4'], 'bulk viscosity b1 = 0.2, b2 = 2.4')):
        d, fix, mon = wave_bar(name, 'C3D8R bar impact, %s' % title)
        d.steps += explicit_step('Wave', 3.2e-05, ['*Boundary', 'FIX, ENCASTRE'], extra=extra)
        register(d, 'explicit', 'C3D8R', title, 'Bar moving at 1 m/s stopped by its clamped end, with %s.' % title,
                 mon, fix, False, 'dynamic', category='feature', group='explicit', twin=False)

    # stable increment of the elements from the uniform-gradient bound only: with variable mass scaling, the added
    # mass (and so the wave speed) depends on the element increment
    d, fix, mon = wave_bar('xpl_improved_dt_off', 'C3D8R bar impact, variable mass scaling, IMPROVED DT METHOD=NO')
    d.model = [l.replace('*Solid Section, elset=EALL, material=STEEL', '*Solid Section, elset=EALL, material=STEEL, '
                         'controls=CDT') for l in d.model] + ['*Section Controls, name=CDT, improved dt method=NO']
    d.steps += explicit_step('Wave', 3.2e-05, ['*Boundary', 'FIX, ENCASTRE'],
                             extra=['*Variable Mass Scaling, dt=5e-07, type=below min, frequency=10'])
    register(d, 'explicit', 'C3D8R', '*SECTION CONTROLS, IMPROVED DT METHOD=NO', 'Bar impact with variable mass '
             'scaling to a 5e-7 s increment, the element increments from the conventional estimate (uniform-gradient '
             'bound) instead of the default improved one.', mon, fix, False, 'dynamic', category='feature',
             group='explicit', twin=False)

    d, fix, mon = wave_bar('xpl_element_by_element', 'C3D8R bar impact, element-by-element stable increment')
    d.steps += ['*Step, name=Wave, nlgeom=YES', '*Dynamic, Explicit, element by element', ', 3.2e-05',
                '*Boundary', 'FIX, ENCASTRE'] + output_block(False, False, 20, 3.2e-05) + ['*End Step']
    register(d, 'explicit', 'C3D8R', 'element-by-element increment', 'Bar impact with the element-by-element estimate '
             'of the stable increment only (no global estimate).', mon, fix, False, 'dynamic', category='feature',
             group='explicit', twin=False)

    d, fix, mon = wave_bar('xpl_fixed_increment', 'C3D8R bar impact, fixed time increment (direct user control)')
    d.steps += ['*Step, name=Wave, nlgeom=YES', '*Dynamic, Explicit, direct user control', '2e-07, 3.2e-05',
                '*Boundary', 'FIX, ENCASTRE', '*Output, field, number interval=20', '*Node Output', 'U, RF, V',
                '*Output, history, variable=PRESELECT, time interval=1.6e-07', '*End Step']
    register(d, 'explicit', 'C3D8R', 'fixed time increment', 'Bar impact with a user-fixed time increment of 2e-7 s.',
             mon, fix, False, 'dynamic', category='feature', group='explicit', twin=False)

    d = Deck('xpl_damping_alpha', 'S4R clamped plate, pressure pulse, mass-proportional damping (alpha = 200 /s)')
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    shell(d, 'S4R', quads, 0.2, 'STEEL_DAMP', 'EALL')
    d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
    edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
    d.steps += explicit_step('Pulse', 1e-3, ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.1'], rot=True,
                             nframes=100)  # about 20 frames per period of the plate
    register(d, 'explicit', 'S4R', 'Rayleigh damping (alpha)', 'Clamped plate under a suddenly applied pressure, '
             'vibration damped by mass-proportional damping.', [ids[4, 4], ids[2, 2]], edge, True, 'dynamic',
             category='feature', group='explicit')


def standard_features():
    d = Deck('std_riks', 'Shallow B31 arch, snap-through under an apex load (arc-length method)')
    n = 20
    ids = []
    for i in range(n + 1):
        x = -20.0 + 40.0 * i / n
        ids.append(d.node(x, 0.0, 2.0 * (1 - (x / 20.0) ** 2)))
    beam(d, [[ids[i], ids[i + 1]] for i in range(n)], 'RECT', '1., 2.', n1='0., 1., 0.')
    d.steps += ['*Step, name=Riks, nlgeom=YES, inc=300', '*Static, riks', '0.02, 1., 1e-06, 0.05, 1.5',
                '*Boundary', '%d, 1, 6' % ids[0], '%d, 1, 6' % ids[-1], 'EALL_NODES, 2, 2',
                '*Cload', '%d, 3, -400.' % ids[n // 2], '*Output, field, frequency=5', '*Node Output', 'U, RF, UR',
                '*Output, history, variable=PRESELECT', '*End Step']
    d.nset('EALL_NODES', ids)
    register(d, 'riks', 'B31', 'arc length (Riks)', 'Clamped shallow arch under an apex load: snap-through traced '
             'with the arc-length method.', [ids[n // 2], ids[n // 4]], [ids[0], ids[-1]], True, 'nonlinear',
             category='feature', group='standard', supported=False, reason=NS + '*STATIC, RIKS', twin=False)

    d, ids = cantilever_beam('std_frequency', 'B31 cantilever, first five natural frequencies')
    d.steps += ['*Step, name=Modes, perturbation', '*Frequency, eigensolver=Lanczos', '5, , , , ,',
                '*Boundary', '%d, ENCASTRE' % ids[0], '*Output, field', '*Node Output', 'U', '*End Step']
    register(d, 'frequency', 'B31', 'natural frequencies', 'First five natural frequencies and mode shapes of a '
             'cantilever.', [ids[-1], ids[5]], [ids[0]], True, 'linear', 'linear', category='feature',
             group='standard', supported=False, reason=NS + '*FREQUENCY', twin=False)

    d, ids = cantilever_beam('std_buckle', 'B31 column, first three buckling loads')
    d.steps += ['*Step, name=Buckle, perturbation', '*Buckle', '3, , 18, 300', '*Boundary', '%d, ENCASTRE' % ids[0],
                '*Cload', '%d, 1, -1.' % ids[-1], '*Output, field', '*Node Output', 'U', '*End Step']
    register(d, 'buckle', 'B31', 'linear buckling', 'First three buckling loads of a cantilever column.', [ids[-1],
             ids[5]], [ids[0]], True, 'linear', 'linear', category='feature', group='standard', supported=False,
             reason=NS + '*BUCKLE', twin=False)


def parallel_features():
    d = Deck('par_plate_pulse', 'Clamped S4R plate 60 x 60 elements, pressure pulse (explicit, parallel run)')
    quads, ids = quad_plate(d, 60, 60, 60.0, 60.0)
    shell(d, 'S4R', quads, 0.5, 'STEEL', 'EALL')
    d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
    edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 60.0) or near(y, 60.0))
    d.steps += explicit_step('Pulse', 2e-4, ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.5'], rot=True)
    register(d, 'explicit', 'S4R', 'parallel explicit', '3600 shell elements: results must not depend on the number '
             'of threads; speed-up measured.', [ids[30, 30], ids[15, 15]], edge, True, 'dynamic', category='feature',
             group='parallel', twin=False)

    d = Deck('par_block_static', 'C3D8R cantilever 40 x 8 x 8 elements, static (parallel run)')
    hexes = hex_block(d, 40, 8, 8, 80.0, 8.0, 8.0)
    solid(d, 'C3D8R', hexes, 'STEEL', 'EALL')
    fix = d.select(lambda x, y, z: near(x, 0))
    tip = d.select(lambda x, y, z: near(x, 80.0))
    d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                           ['%d, 3, %s' % (n, fmt(-500.0 / len(tip))) for n in tip], False, nframes=1)
    register(d, 'static', 'C3D8R', 'parallel static', '2560 solid elements: results must not depend on the number of '
             'threads; speed-up measured.', tip[:3], fix, False, 'linear', 'linear', category='feature',
             group='parallel', twin=False)


def generate():
    contact_features()
    constraint_features()
    element_features()
    load_features()
    material_features()
    explicit_features()
    standard_features()
    parallel_features()
