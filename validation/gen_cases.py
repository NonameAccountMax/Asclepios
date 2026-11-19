#!/usr/bin/env python3
"""Validation decks: every element type supported by asclepios under the main deformation modes (this file), solver
features (gen_features.py), and the twin of every Abaqus-referenced case in the other procedure (static <-> explicit).

Each deck uses only keywords understood by both Abaqus and asclepios, unless the case is marked as not supported by
asclepios (it is then run in Abaqus only, to measure coverage). Writes cases/<name>.inp and cases/manifest.json.
Units: mm, N, tonne, s, MPa, mJ.
"""
import math
import os

from deck import (CASES, KNOWN, KNOWN_QUANTITIES, MODE, OUT, QS_DURATION, Deck, beam_line, boundary_stretch,
                  edge_shares, explicit_step, faces_on, fmt, hex_block, near, quad_plate, register, static_step,
                  tets_from_hexes, tris_from_quads, write_manifest)
import gen_features
import gen_more
import gen_literature

KNOWN.update({
    'contact_block_on_plate': 'general contact does not include solid element faces and nodes yet',
    'multistep_opnew_xpl': 'OP=NEW removes the tip load at once in explicit dynamics: free vibration of period 0.9 ms '
                           'sampled every 7.5 ms, phase not comparable (the kinetic energy agrees within 2 % over the '
                           'first periods)',
    's4r_twist': 'the closed form is the Kirchhoff plate: with transverse shear deformation (Reissner-Mindlin shell) '
                 'the corner point loads add a local deflection, 0.5 % of the twist',
})
ENERGY_NOTE = ('ALLIE and ALLWK are the work of the nodal forces integrated with the trapezoidal rule (thermal strain '
               'work included, increment-size error across the elastic-plastic transition); the reference is the '
               'strain energy plus the plastic dissipation, as Abaqus computes it')
KNOWN_QUANTITIES.update({name: (['E'], ENERGY_NOTE) for name in ('c3d8r_thermal', 'b31_thermal', 'c3d8r_plastic',
                                                                  'c3d8_plastic')})
KNOWN_QUANTITIES['symmetry_bcs_xpl'] = (['RM'], 'Reaction moments not compared: the corner nodes carry both the clamping '
                                        'and the symmetry conditions')
S3R_XPL_RM = ('S3R reaction moments about the clamped edge not compared in explicit dynamics (the reference values do '
              'not follow from equilibrium); asclepios gives the implicit value')
for _n in ('elem_s3_xpl', 's3r_pressure_xpl', 's3r_pressure_pulse'):
    KNOWN_QUANTITIES[_n] = (['RM'], S3R_XPL_RM)
KNOWN_QUANTITIES['s4r_single_element_probes'] = (['RM'], 'Reaction moments not compared: transverse shear of the strongly distorted '
                                                 'quad (probe 4) and drilling stiffness of distorted quads.')
# duration of the quasi-static explicit twins: a few periods of the lowest mode, loads following a smooth step
QS_DURATION.update({
    'c3d8r_bending': 0.02, 'c3d8_bending': 0.02, 'c3d4_bending': 0.01,
    'c3d8r_torsion': 0.01, 'c3d8_torsion': 0.01, 'c3d4_torsion': 0.005, 'c3d8r_pressure': 0.02,
    's4r_pressure': 0.002, 's3r_pressure': 0.002, 's3r_bending': 0.03, 's3r_rollup': 0.03,
    's4r_plastic_bending': 0.03, 'b31_plastic': 0.01,
})


# ------------------------------------------------------------------------------------------------ solids
L, B, H = 40.0, 4.0, 4.0


def solid_bar(name, etype, title, nx=10, ny=2, nz=2, lx=L, ly=B, lz=H, mat='STEEL'):
    d = Deck(name, title)
    hexes = hex_block(d, nx, ny, nz, lx, ly, lz)
    conns = tets_from_hexes(d, hexes) if etype == 'C3D4' else hexes
    elist = d.elements(etype, 'EALL', conns)
    d.materials.add(mat)
    d.model.append('*Solid Section, elset=EALL, material=%s' % mat)
    d.model.append(',')
    return d, elist, conns


def solid_cases(etype, modes):
    tag = etype.lower()
    if 'tension' in modes:
        d, el, cn = solid_bar('%s_tension' % tag, etype, '%s bar, uniaxial tension with free contraction (linear)' % etype)
        fx = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, L))
        d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
        d.nset('TIP', tip)
        d.steps += static_step('Load', ['*Boundary', 'FIX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3', '*Boundary', 'TIP, 1, 1, 0.04'],
                               False, nframes=1)
        register(d, 'static', etype, 'tension', 'Bar stretched by 0.1 % with free lateral contraction (homogeneous '
                 'linear elastic state).', tip + d.select(lambda x, y, z: near(x, L / 2) and near(z, H)), fx, False,
                 'exact', 'linear', reference='analytic')
    if 'bending' in modes:
        d, el, cn = solid_bar('%s_bending' % tag, etype, '%s cantilever, tip load (linear)' % etype)
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, L))
        d.nset('TIP', tip)
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 3, %s' % (n, fmt(-200.0 / len(tip))) for n in tip], False, nframes=1)
        register(d, 'static', etype, 'bending', 'Cantilever under a 200 N transverse tip load (linear elastic).',
                 tip + d.select(lambda x, y, z: near(x, L / 2) and near(z, H)), fix, False, 'linear', 'linear')
    if 'torsion' in modes:
        d, el, cn = solid_bar('%s_torsion' % tag, etype, '%s bar, end face twisted by 30 deg (NLGEOM)' % etype)
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, L))
        th = math.radians(30.0)
        body = ['*Boundary', 'FIX, ENCASTRE', '*Boundary']
        for n in tip:
            x, y, z = d.nodes[n]
            yc, zc = y - B / 2, z - H / 2
            body += ['%d, 2, 2, %s' % (n, fmt(yc * math.cos(th) - zc * math.sin(th) - yc)),
                     '%d, 3, 3, %s' % (n, fmt(yc * math.sin(th) + zc * math.cos(th) - zc))]
        d.steps += static_step('Twist', body, True)
        register(d, 'static', etype, 'torsion', 'Bar twisted by 30 degrees at its free end (large rotation).',
                 tip + d.select(lambda x, y, z: near(x, L / 2) and (near(y, 0) or near(y, B)) and near(z, H)),
                 fix, False, 'nonlinear')
    if 'shear' in modes:
        d, el, cn = solid_bar('%s_shear_hyper' % tag, etype, '%s rubber block, simple shear 50 %% (NLGEOM)' % etype,
                              nx=3, ny=3, nz=3, lx=6.0, ly=6.0, lz=6.0, mat='RUBBER')
        bnd = d.select(lambda x, y, z: min(x, y, z) < 1e-9 or max(x, y, z) > 6.0 - 1e-9)
        inner = [n for n in d.nodes if n not in set(bnd)]
        d.nset('BND', bnd)
        body = ['*Boundary', 'BND, 2, 3, 0.', '*Boundary'] + ['%d, 1, 1, %s' % (n, fmt(0.5 * d.nodes[n][2])) for n in bnd]
        d.steps += static_step('Shear', body, True, inc=0.02, maxinc=0.02)
        register(d, 'static', etype, 'shear', 'Neo-Hookean block, homogeneous simple shear of 50 % imposed on the whole '
                 'boundary (hyperelastic, large strain).', inner, [], False, 'exact', reference='analytic')
    if 'plastic' in modes:
        d, el, cn = solid_bar('%s_plastic' % tag, etype, '%s bar, tension into plasticity then reversal' % etype,
                              nx=4, ny=1, nz=1, mat='STEEL_PL')
        fx = d.select(lambda x, y, z: near(x, 0))
        fy = d.select(lambda x, y, z: near(y, 0))
        fz = d.select(lambda x, y, z: near(z, 0))
        tip = d.select(lambda x, y, z: near(x, L))
        d.nset('FX', fx); d.nset('FY', fy); d.nset('FZ', fz); d.nset('TIP', tip)
        bcs = ['*Boundary', 'FX, 1, 1', 'FY, 2, 2', 'FZ, 3, 3']
        d.steps += static_step('Load', bcs + ['*Boundary', 'TIP, 1, 1, 2.'], True, inc=0.02, maxinc=0.02)
        d.steps += static_step('Reverse', ['*Boundary', 'TIP, 1, 1, 1.'], True, inc=0.02, maxinc=0.02)
        register(d, 'static', etype, 'plasticity', 'Uniaxial tension to 5 % (J2, isotropic hardening), then reversed to '
                 '2.5 % (elastic unloading and reverse yielding).', tip, fx, False, 'exact', reference='analytic')
    if 'pressure' in modes:
        d, el, cn = solid_bar('%s_pressure' % tag, etype, '%s cantilever under follower pressure (NLGEOM)' % etype)
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, L))
        top = faces_on(d, etype, el, cn, lambda x, y, z: near(z, H))
        d.model.append('*Surface, type=ELEMENT, name=STOP')
        d.model += ['%d, S%d' % ef for ef in top]
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 20.'], True)
        register(d, 'static', etype, 'pressure', 'Cantilever under a 20 MPa follower pressure on its top face '
                 '(large deflection).', tip, fix, False, 'nonlinear')
    if 'impact' in modes:
        d, el, cn = solid_bar('%s_wave' % tag, etype, '%s bar with initial velocity against a clamped end' % etype)
        fix = d.select(lambda x, y, z: near(x, 0))
        free = [n for n in d.nodes if n not in set(fix)]
        d.nset('MOVING', free)
        d.initial += ['*Initial Conditions, type=VELOCITY', 'MOVING, 1, -1000.']
        d.steps += explicit_step('Wave', 3.2e-05, ['*Boundary', 'FIX, ENCASTRE'])
        register(d, 'explicit', etype, 'wave', 'Bar moving at 1 m/s stopped by its clamped end: compression wave and '
                 'reflections (explicit dynamics).', d.select(lambda x, y, z: (near(x, L) or near(x, L / 2)) and near(y, 0)
                                                              and near(z, 0)), fix, False, 'dynamic')
    if 'thermal' in modes:
        d, el, cn = solid_bar('%s_thermal' % tag, etype, '%s bar heated by 100 K between fixed ends' % etype)
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, L))
        d.nset('TIP', tip)
        d.nset('NALL', list(d.nodes))
        d.initial += ['*Initial Conditions, type=TEMPERATURE', 'NALL, 20.']
        d.nset('FY', d.select(lambda x, y, z: near(y, 0))); d.nset('FZ', d.select(lambda x, y, z: near(z, 0)))
        d.steps += static_step('Heat', ['*Boundary', 'FIX, 1, 1', 'TIP, 1, 1', 'FY, 2, 2', 'FZ, 3, 3',
                                        '*Temperature', 'NALL, 120.'], False, nframes=1)
        register(d, 'static', etype, 'thermal', 'Bar heated by 100 K with both ends held axially and free lateral '
                 'expansion (homogeneous thermal stress).', d.select(lambda x, y, z: near(x, L / 2)), fix, False, 'exact',
                 'linear', reference='analytic')


# ------------------------------------------------------------------------------------------------ shells
def shell_model(d, etype, conns, t, mat='STEEL', orient=None, npt=5):
    el = d.elements(etype, 'EALL', conns)
    d.materials.add(mat)
    if orient:
        d.model += ['*Orientation, name=ORI, local directions=2', orient]
        d.model += ['*Shell Section, elset=EALL, material=%s, orientation=ORI' % mat, '%s, %d' % (fmt(t), npt)]
    else:
        d.model += ['*Shell Section, elset=EALL, material=%s' % mat, '%s, %d' % (fmt(t), npt)]
    return el


def plate(name, etype, title, nx, ny, lx, ly, t, mat='STEEL', orient=None, x0=(0.0, 0.0, 0.0)):
    d = Deck(name, title)
    quads, ids = quad_plate(d, nx, ny, lx, ly, x0)
    shell_model(d, etype, tris_from_quads(quads) if etype == 'S3R' else quads, t, mat, orient)
    return d, ids


def shell_cases(etype, modes):
    tag = etype.lower()
    lx, ly, t = 40.0, 10.0, 0.5
    ei = 210000.0 * ly * t ** 3 / 12
    if 'membrane' in modes:
        d, ids = plate('%s_membrane' % tag, etype, '%s strip in uniaxial membrane tension' % etype, 10, 3, lx, ly, t)
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, lx))
        d.nset('TIP', tip); d.nset('NALL', list(d.nodes))
        d.steps += static_step('Load', ['*Boundary', 'NALL, 3, 5', 'FIX, 1, 1', '%d, 2, 2' % ids[0, 0],
                                        '*Boundary', 'TIP, 1, 1, 0.04'], False, rot=True, nframes=1)
        register(d, 'static', etype, 'membrane', 'Strip stretched in its plane by 0.1 % (homogeneous plane stress, '
                 'linear).', tip, fix, True, 'exact', 'linear', reference='analytic')
    if 'bending' in modes:
        d, ids = plate('%s_bending' % tag, etype, '%s cantilever strip (nu = 0), tip load (linear)' % etype,
                       40 if etype == 'S4R' else 10, 3, lx, ly, t, mat='STEEL0')
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, lx))
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 3, %s' % (n, fmt(-1.0 * w)) for n, w in edge_shares(d, tip)], False, rot=True,
                               nframes=1)
        register(d, 'static', etype, 'bending', 'Cantilever strip (Poisson ratio 0) under a 1 N tip load: Timoshenko '
                 'beam solution (linear).', tip + d.select(lambda x, y, z: near(x, lx / 2)), fix, True,
                 'theory' if etype == 'S4R' else 'linear', 'linear', reference='analytic' if etype == 'S4R' else 'abaqus')
    if 'rollup' in modes:
        d, ids = plate('%s_rollup' % tag, etype, '%s cantilever strip (nu = 0) rolled into a quarter circle by an end '
                       'moment' % etype, 16, 2, lx, ly, t, mat='STEEL0')
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, lx))
        m = math.pi * ei / (2 * lx)
        d.steps += static_step('Moment', ['*Boundary', 'FIX, ENCASTRE', '*Cload'] +
                               ['%d, 5, %s' % (n, fmt(-m * w)) for n, w in edge_shares(d, tip)], True, rot=True,
                               inc=0.02, maxinc=0.02)
        register(d, 'static', etype, 'large rotation', 'Cantilever strip (Poisson ratio 0) rolled into a quarter circle '
                 'by an end moment: nodes on the exact circular arc (90 degree rotation).',
                 tip + d.select(lambda x, y, z: near(x, lx / 2)), fix, True, 'theory' if etype == 'S4R' else 'nonlinear',
                 reference='analytic' if etype == 'S4R' else 'abaqus')
    if 'twist' in modes:
        d, ids = plate('%s_twist' % tag, etype, '%s square plate twisted by a corner load' % etype, 6, 6, 20.0, 20.0, t)
        c00, c10, c01, c11 = ids[0, 0], ids[6, 0], ids[0, 6], ids[6, 6]
        d.steps += static_step('Load', ['*Boundary', '%d, 1, 3' % c00, '%d, 2, 3' % c10, '%d, 3, 3' % c01,
                                        '*Cload', '%d, 3, -2.' % c11], False, rot=True, nframes=1)
        register(d, 'static', etype, 'twisting', 'Square plate on three corner supports, transverse load on the fourth '
                 '(uniform twisting moment, w = c x y, linear).', [c11, ids[3, 3], ids[6, 3], ids[3, 6]], [c00, c10, c01],
                 True, 'exact', 'linear', reference='analytic')
    if 'pressure' in modes:
        d, ids = plate('%s_pressure' % tag, etype, '%s clamped square plate under pressure (NLGEOM)' % etype,
                       8, 8, 20.0, 20.0, 0.2)
        edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
        d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.5'], True, rot=True)
        register(d, 'static', etype, 'pressure', 'Clamped square plate under 0.5 MPa follower pressure; deflection of '
                 'several thicknesses (membrane stiffening).', [ids[4, 4], ids[2, 2], ids[4, 2], ids[2, 4]], edge,
                 True, 'nonlinear')
    if 'fung' in modes:
        c, s = math.cos(math.radians(30)), math.sin(math.radians(30))
        d, ids = plate('%s_fung_biaxial' % tag, etype, '%s Fung anisotropic patch, biaxial stretch 1.15 x 1.10' % etype,
                       2, 2, 2.0, 2.0, 0.5, mat='FUNG', orient='%s, %s, 0., %s, %s, 0.' % (fmt(c), fmt(s), fmt(-s), fmt(c)))
        x0 = d.select(lambda x, y, z: near(x, 0))
        y0 = d.select(lambda x, y, z: near(y, 0))
        x1 = d.select(lambda x, y, z: near(x, 2.0))
        y1 = d.select(lambda x, y, z: near(y, 2.0))
        d.nset('X0', x0); d.nset('Y0', y0); d.nset('X1', x1); d.nset('Y1', y1); d.nset('NALL', list(d.nodes))
        bnd = x0 + y0 + x1 + y1
        d.steps += static_step('Stretch', ['*Boundary', 'NALL, 3, 5', '*Boundary'] + boundary_stretch(d, bnd, 0.15, 0.10),
                               True, rot=True, inc=0.01, maxinc=0.01)
        register(d, 'static', etype, 'hyperelastic (Fung)', 'Fung anisotropic membrane, material axes at 30 degrees, '
                 'homogeneous biaxial stretch 1.15 x 1.10 imposed on the boundary (plane stress, large strain).',
                 [ids[1, 1]], x0, True, 'exact', reference='analytic')
    if 'neohooke' in modes:
        d, ids = plate('%s_neohooke_biaxial' % tag, etype, '%s neo-Hookean membrane, biaxial stretch 1.5 x 1.3' % etype,
                       2, 2, 2.0, 2.0, 0.5, mat='RUBBER')
        x0 = d.select(lambda x, y, z: near(x, 0))
        y0 = d.select(lambda x, y, z: near(y, 0))
        d.nset('X0', x0); d.nset('Y0', y0); d.nset('NALL', list(d.nodes))
        d.nset('X1', d.select(lambda x, y, z: near(x, 2.0))); d.nset('Y1', d.select(lambda x, y, z: near(y, 2.0)))
        bnd = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 2.0) or near(y, 2.0))
        d.steps += static_step('Stretch', ['*Boundary', 'NALL, 3, 5', '*Boundary'] + boundary_stretch(d, bnd, 0.5, 0.3),
                               True, rot=True, inc=0.02, maxinc=0.02)
        register(d, 'static', etype, 'hyperelastic (neo-Hooke)', 'Nearly incompressible neo-Hookean membrane, homogeneous '
                 'biaxial stretch 1.5 x 1.3 imposed on the boundary (thickness reduction, large strain).', [ids[1, 1]], x0,
                 True, 'exact', reference='analytic')
    if 'plastic' in modes:
        d, ids = plate('%s_plastic_bending' % tag, etype, '%s elastic-plastic cantilever, tip displacement then '
                       'partial return' % etype, 10, 3, lx, ly, 1.0, mat='STEEL_PL')
        fix = d.select(lambda x, y, z: near(x, 0))
        tip = d.select(lambda x, y, z: near(x, lx))
        d.nset('TIP', tip)
        d.steps += static_step('Load', ['*Boundary', 'FIX, ENCASTRE', '*Boundary', 'TIP, 3, 3, -8.'], True, rot=True)
        d.steps += static_step('Return', ['*Boundary', 'TIP, 3, 3, -4.'], True, rot=True)
        register(d, 'static', etype, 'plasticity', 'Cantilever plate bent into plasticity by an 8 mm tip displacement, '
                 'then brought back to 4 mm (plastic hinge, springback).', tip + d.select(lambda x, y, z: near(x, lx / 2)),
                 fix, True, 'nonlinear')
    if 'pulse' in modes:
        d, ids = plate('%s_pressure_pulse' % tag, etype, '%s clamped plate, suddenly applied pressure (explicit)' % etype,
                       8, 8, 20.0, 20.0, 0.2)
        edge = d.select(lambda x, y, z: near(x, 0) or near(y, 0) or near(x, 20.0) or near(y, 20.0))
        d.model += ['*Surface, type=ELEMENT, name=STOP', 'EALL, SPOS']
        d.steps += explicit_step('Pulse', 5e-4, ['*Boundary', 'FIX, ENCASTRE', '*Dsload', 'STOP, P, 0.5'], rot=True,
                                 nframes=50)  # about 20 frames per period of the plate
        register(d, 'explicit', etype, 'dynamic pressure', 'Clamped square plate under a suddenly applied 0.5 MPa '
                 'pressure: large-amplitude vibration (explicit dynamics).', [ids[4, 4], ids[2, 2], ids[4, 2]], edge,
                 True, 'dynamic')
    if 'spin' in modes:
        d, ids = plate('%s_spin' % tag, etype, '%s free plate spinning in its plane (explicit)' % etype, 6, 6, 20.0,
                       20.0, 0.5, x0=(-10.0, -10.0, 0.0))
        w = 500.0
        d.nset('NALL', list(d.nodes))
        d.initial += ['*Initial Conditions, type=VELOCITY']
        for n, (x, y, z) in sorted(d.nodes.items()):
            d.initial += ['%d, 1, %s' % (n, fmt(-w * y)), '%d, 2, %s' % (n, fmt(w * x)), '%d, 6, %s' % (n, fmt(w))]
        d.steps += explicit_step('Spin', 2e-3, [], rot=True)
        register(d, 'explicit', etype, 'rigid rotation', 'Free plate spinning at 500 rad/s (about 9 degrees per frame, '
                 '57 degrees in total): rigid rotation, objectivity (explicit dynamics).',
                 [ids[6, 6], ids[6, 3], ids[3, 3]], [], True, 'exact', reference='analytic')


# ------------------------------------------------------------------------------------------------- beams
def beam_model(name, title, section, dims, n=10, length=100.0, mat='STEEL'):
    d = Deck(name, title)
    segs, ids = beam_line(d, n, length)
    d.elements('B31', 'EALL', segs)
    d.materials.add(mat)
    d.model += ['*Beam Section, elset=EALL, material=%s, section=%s' % (mat, section), dims, '0., 0., -1.']
    return d, ids


def beam_cases():
    d, ids = beam_model('b31_tension', 'B31 rectangular beam, axial displacement', 'RECT', '4., 2.')
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Boundary', '%d, 1, 1, 0.05' % ids[-1]],
                           False, rot=True, nframes=1)
    register(d, 'static', 'B31', 'tension', 'Clamped beam stretched by 0.05 % (linear).', [ids[-1], ids[5]], [ids[0]],
             True, 'exact', 'linear', reference='analytic')

    d, ids = beam_model('b31_bending', 'B31 rectangular cantilever, tip loads in both planes', 'RECT', '2., 1.', n=50)
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload', '%d, 2, 0.5' % ids[-1],
                                    '%d, 3, -2.' % ids[-1]], False, rot=True, nframes=1)
    register(d, 'static', 'B31', 'bending', 'Slender cantilever with tip loads in both principal planes (two bending '
             'stiffnesses, orientation of the section): Timoshenko beam solution (linear).', [ids[-1], ids[25]], [ids[0]],
             True, 'theory', 'linear', reference='analytic')

    d, ids = beam_model('b31_torsion', 'B31 circular cantilever, tip torque', 'CIRC', '2.')
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload', '%d, 4, 1000.' % ids[-1]],
                           False, rot=True, nframes=1)
    register(d, 'static', 'B31', 'torsion', 'Circular cantilever under a 1 N.m tip torque (linear).', [ids[-1], ids[5]],
             [ids[0]], True, 'exact', 'linear', reference='analytic')

    d, ids = beam_model('b31_rollup', 'B31 cantilever rolled into a quarter circle by an end moment', 'RECT', '1., 4.',
                        n=16)
    ei = 210000.0 * 1.0 ** 3 * 4.0 / 12
    d.steps += static_step('Moment', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Cload',
                                      '%d, 5, %s' % (ids[-1], fmt(-math.pi * ei / 200.0))], True, rot=True, inc=0.02,
                           maxinc=0.02)
    register(d, 'static', 'B31', 'large rotation', 'Cantilever rolled into a quarter circle by an end moment: nodes on '
             'the exact circular arc.', [ids[-1], ids[8]], [ids[0]], True, 'theory', reference='analytic')

    d, ids = beam_model('b31_plastic', 'B31 circular cantilever bent into plasticity, then partial return', 'CIRC', '2.',
                        mat='STEEL_PL')
    d.steps += static_step('Load', ['*Boundary', '%d, ENCASTRE' % ids[0], '*Boundary', '%d, 3, 3, -20.' % ids[-1]],
                           True, rot=True)
    d.steps += static_step('Return', ['*Boundary', '%d, 3, 3, -10.' % ids[-1]], True, rot=True)
    register(d, 'static', 'B31', 'plasticity', 'Elastic-plastic cantilever: 20 mm tip displacement, then back to 10 mm.',
             [ids[-1], ids[5]], [ids[0]], True, 'nonlinear')

    d, ids = beam_model('b31_thermal', 'B31 pipe heated by 100 K between fixed ends', 'PIPE', '3., 0.5')
    d.nset('NALL', ids)
    d.initial += ['*Initial Conditions, type=TEMPERATURE', 'NALL, 20.']
    d.steps += static_step('Heat', ['*Boundary', '%d, ENCASTRE' % ids[0], '%d, ENCASTRE' % ids[-1],
                                    '*Temperature', 'NALL, 120.'], False, rot=True, nframes=1)
    register(d, 'static', 'B31', 'thermal', 'Pipe heated by 100 K between clamped ends (thermal axial force, linear).',
             [ids[5], ids[3]], [ids[0]], True, 'exact', 'linear', reference='analytic')

    d, ids = beam_model('b31_vibration', 'B31 pipe cantilever, free vibration (explicit)', 'PIPE', '3., 0.5')
    d.initial += ['*Initial Conditions, type=VELOCITY']
    for i, n in enumerate(ids):
        s = i / (len(ids) - 1.0)
        d.initial.append('%d, 3, %s' % (n, fmt(100.0 * s * s)))
    # first period of the pipe about 1.8 ms: two periods, 20 frames per period
    d.steps += explicit_step('Vibrate', 3.6e-3, ['*Boundary', '%d, ENCASTRE' % ids[0]], rot=True, nframes=40)
    register(d, 'explicit', 'B31', 'vibration', 'Cantilever pipe released with a first-mode-like velocity field: two '
             'periods of free vibration (explicit dynamics).', [ids[-1], ids[5]], [ids[0]], True, 'dynamic')

    # same with the isotropic (scaled) rotary inertia instead of the default exact one
    d, ids = beam_model('b31_vibration_isotropic', 'B31 pipe cantilever, free vibration, ROTARY INERTIA=ISOTROPIC',
                        'PIPE', '3., 0.5')
    d.model = [l.replace('section=PIPE', 'section=PIPE, rotary inertia=ISOTROPIC') for l in d.model]
    d.initial += ['*Initial Conditions, type=VELOCITY']
    for i, n in enumerate(ids):
        s = i / (len(ids) - 1.0)
        d.initial.append('%d, 3, %s' % (n, fmt(100.0 * s * s)))
    d.steps += explicit_step('Vibrate', 3.6e-3, ['*Boundary', '%d, ENCASTRE' % ids[0]], rot=True, nframes=40)
    register(d, 'explicit', 'B31', '*BEAM SECTION, ROTARY INERTIA=ISOTROPIC', 'b31_vibration with the isotropic '
             'rotary inertia (scaled so that it does not limit the stable increment).', [ids[-1], ids[5]], [ids[0]],
             True, 'dynamic', twin=False)


# -------------------------------------------------------------------------------------------- connectors
def connector_cases():
    d = Deck('conn3d2_axial', 'CONN3D2 AXIAL spring: axial load, then large transverse motion (NLGEOM)')
    a, b = d.node(0, 0, 0), d.node(10, 0, 0)
    d.elements('CONN3D2', 'CONN', [[a, b]])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL',
                '*Connector Behavior, name=SPRING', '*Connector Elasticity, component=1', '100.']
    d.steps += static_step('Axial', ['*Boundary', '%d, 1, 3' % a, '%d, 2, 3' % b, '*Cload', '%d, 1, 50.' % b], True,
                           inc=0.02, maxinc=0.02)
    d.steps += static_step('Swing', ['*Boundary', '%d, 2, 2, 5.' % b], True, inc=0.02, maxinc=0.02)
    register(d, 'static', 'CONN3D2', 'axial', 'Axial connector: 50 N axial load, then 5 mm transverse displacement '
             '(rotation of the connector axis).', [b], [a], False, 'exact', reference='analytic')

    d = Deck('conn3d2_cartesian', 'CONN3D2 CARTESIAN spring loaded in three directions')
    a, b = d.node(0, 0, 0), d.node(10, 0, 0)
    d.elements('CONN3D2', 'CONN', [[a, b]])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'CARTESIAN',
                '*Connector Behavior, name=SPRING', '*Connector Elasticity, component=1', '100.',
                '*Connector Elasticity, component=2', '200.', '*Connector Elasticity, component=3', '300.']
    d.steps += static_step('Load', ['*Boundary', '%d, 1, 6' % a, '*Cload', '%d, 1, 10.' % b, '%d, 2, 20.' % b,
                                    '%d, 3, 30.' % b], False, nframes=1)
    register(d, 'static', 'CONN3D2', 'cartesian', 'Cartesian connector with three different stiffnesses (linear).',
             [b], [a], False, 'exact', 'linear', reference='analytic')

    d = Deck('conn3d2_oscillator', 'Mass on four axial connectors (explicit)')
    hexes = hex_block(d, 1, 1, 1, 2.0, 2.0, 2.0, x0=(10.0, -1.0, -1.0))
    d.elements('C3D8R', 'MASS', hexes)
    d.materials.add('HEAVY')
    d.model += ['*Solid Section, elset=MASS, material=HEAVY', ',']
    face = d.select(lambda x, y, z: near(x, 10.0))
    ground = [d.node(0.0, d.nodes[n][1], d.nodes[n][2]) for n in face]
    d.elements('CONN3D2', 'CONN', [[g, n] for g, n in zip(ground, face)])
    d.model += ['*Connector Section, elset=CONN, behavior=SPRING', 'AXIAL',
                '*Connector Behavior, name=SPRING', '*Connector Elasticity, component=1', '25.']
    cube = [n for n in d.nodes if n not in set(ground)]
    d.nset('CUBE', cube); d.nset('GROUND', ground)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'CUBE, 1, 100.']
    d.steps += explicit_step('Oscillate', 0.01, ['*Boundary', 'GROUND, 1, 3', 'CUBE, 2, 3'], time_marks=True)
    register(d, 'explicit', 'CONN3D2', 'oscillation', 'Steel-density block on four axial springs, released with '
             '100 mm/s: two periods of oscillation (explicit dynamics).', cube[:2], ground, False, 'exact',
             reference='analytic')


# ----------------------------------------------------------------------------------------------- contact
def contact_block(d, mat):
    d.model += ['*Surface Interaction, name=IP']
    return ['*Contact', '*Contact Inclusions, ALL EXTERIOR', '*Contact Property Assignment', ' , , IP']


def contact_cases():
    d = Deck('contact_block_on_plate', 'C3D8R block dropped on a clamped S4R plate (explicit general contact)')
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    shell_model(d, 'S4R', quads, 0.5)
    hexes = hex_block(d, 3, 3, 2, 6.0, 6.0, 4.0, x0=(7.0, 7.0, 0.75))
    d.elements('C3D8R', 'BLOCK', hexes)
    d.materials.add('ALU')
    d.model += ['*Solid Section, elset=BLOCK, material=ALU', ',']
    edge = [ids[i, j] for (i, j) in ids if i in (0, 8) or j in (0, 8)]
    block = d.select(lambda x, y, z: z > 0.5)
    d.nset('BLOCKN', block)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'BLOCKN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE'] + contact_block(d, 'ALU'), rot=True,
                             nframes=200)  # plate period about 0.09 ms: 20 frames gave 2 per period
    register(d, 'explicit', 'C3D8R + S4R', 'contact (impact)', 'Aluminium block dropped at 2 m/s on a clamped steel '
             'plate: impact, rebound and plate vibration (general contact, frictionless).',
             [ids[4, 4], ids[2, 4]] + d.select(lambda x, y, z: near(x, 10.0) and near(y, 10.0)), edge, True, 'contact')

    d = Deck('contact_shell_on_plate', 'S4R sheet dropped flat on a clamped S4R plate (explicit general contact)')
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    shell_model(d, 'S4R', quads, 0.5)
    squads, sids = quad_plate(d, 3, 3, 6.0, 6.0, x0=(7.0, 7.0, 1.0))
    d.elements('S4R', 'SHEET', squads)
    d.materials.add('ALU')
    d.model += ['*Shell Section, elset=SHEET, material=ALU', '0.5, 5']
    edge = [ids[i, j] for (i, j) in ids if i in (0, 8) or j in (0, 8)]
    d.nset('SHEETN', list(sids.values()))
    d.initial += ['*Initial Conditions, type=VELOCITY', 'SHEETN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE'] + contact_block(d, 'ALU'), rot=True,
                             nframes=200)  # plate period about 0.09 ms: 20 frames gave 2 per period
    register(d, 'explicit', 'S4R', 'contact (shell on shell)', 'Aluminium sheet dropped flat at 2 m/s on a clamped steel '
             'plate: shell-to-shell impact and rebound (general contact with thickness offsets, frictionless).',
             [ids[4, 4], ids[2, 4], sids[1, 1], sids[0, 0]], edge, True, 'contact')

    d = Deck('contact_friction_sliding', 'S4R sheet pressed and dragged on a fixed S4R plate (Coulomb friction 0.3)')
    quads, ids = quad_plate(d, 12, 4, 30.0, 10.0)
    shell_model(d, 'S4R', quads, 0.5)
    # sled 0.001 mm above the contact position (half thicknesses 0.25 for the plate and 0.3 x 0.707 for the sled
    # facets), so that it is pressed without dropping onto the plate
    squads, sids = quad_plate(d, 4, 4, 4.0, 4.0, x0=(5.0, 3.0, 0.463))
    sel = d.elements('S4R', 'SLED', squads)
    d.model += ['*Shell Section, elset=SLED, material=STEEL', '0.5, 5']
    plate_nodes = list(ids.values())
    sled = list(sids.values())
    d.nset('SLEDN', sled)
    d.model += ['*Surface Interaction, name=IP', '*Friction', '0.3']
    d.model += ['*Amplitude, name=RAMP, definition=SMOOTH STEP', '0., 0., 0.001, 1.']
    d.model += ['*Amplitude, name=DRAG', '0., 0., 0.001, 0., 0.0012, 1., 0.003, 1.']
    d.model += ['*Surface, type=ELEMENT, name=STOP', 'SLED, SPOS']
    d.steps += explicit_step('Slide', 3e-3, ['*Boundary', 'FIX, ENCASTRE',
                                            '*Boundary, type=VELOCITY, amplitude=DRAG', 'SLEDN, 1, 1, 100.',
                                            '*Dsload, amplitude=RAMP', 'STOP, P, 1.',
                                            '*Contact', '*Contact Inclusions, ALL EXTERIOR',
                                            '*Contact Property Assignment', ' , , IP'], rot=True)
    register(d, 'explicit', 'S4R', 'contact (friction)', 'Steel sheet pressed with 16 N, then dragged at 100 mm/s on a '
             'fixed plate: sticking then sliding, friction force 0.3 x 16 N (general contact).',
             [sids[2, 2], sids[0, 0]], plate_nodes, True, 'contact')

    d = Deck('contact_beam_on_plate', 'B31 rod dropped on a clamped S4R plate (explicit general contact)')
    quads, ids = quad_plate(d, 8, 8, 20.0, 20.0)
    shell_model(d, 'S4R', quads, 0.5)
    segs, bids = beam_line(d, 10, 16.0, x0=(2.0, 10.0, 1.5))
    d.elements('B31', 'ROD', segs)
    d.model += ['*Beam Section, elset=ROD, material=STEEL, section=CIRC', '1.', '0., 0., -1.']
    edge = [ids[i, j] for (i, j) in ids if i in (0, 8) or j in (0, 8)]
    d.nset('RODN', bids)
    d.initial += ['*Initial Conditions, type=VELOCITY', 'RODN, 3, -2000.']
    d.steps += explicit_step('Impact', 1e-3, ['*Boundary', 'FIX, ENCASTRE'] + contact_block(d, 'STEEL'), rot=True,
                             nframes=200)  # plate period about 0.09 ms: 20 frames gave 2 per period
    register(d, 'explicit', 'B31 + S4R', 'contact (beam on shell)', 'Steel rod dropped at 2 m/s across a clamped '
             'plate: beam-to-shell impact (general contact, frictionless).', [ids[4, 4], bids[5], bids[0]], edge, True,
             'contact')


def generate():
    solid_cases('C3D8R', ['tension', 'bending', 'torsion', 'shear', 'plastic', 'pressure', 'impact', 'thermal'])
    solid_cases('C3D8', ['tension', 'bending', 'torsion', 'shear', 'plastic', 'impact'])
    solid_cases('C3D4', ['tension', 'bending', 'torsion', 'shear', 'impact'])
    shell_cases('S4R', ['membrane', 'bending', 'rollup', 'twist', 'pressure', 'fung', 'neohooke', 'plastic', 'pulse',
                        'spin'])
    shell_cases('S3R', ['membrane', 'bending', 'rollup', 'pressure', 'fung', 'pulse'])
    beam_cases()
    connector_cases()
    contact_cases()
    gen_features.generate()
    gen_more.generate()
    gen_literature.generate()


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(OUT):
        if f.endswith(('.inp', '.inc')):
            os.remove(os.path.join(OUT, f))
    MODE['twin'] = False
    generate()
    MODE['twin'] = True
    generate()
    write_manifest()
    n_twin = sum(1 for c in CASES if c['twin_of'])
    n_sup = sum(1 for c in CASES if c['supported'])
    print('%d cases (%d twins), %d run by asclepios, written to %s' % (len(CASES), n_twin, n_sup, OUT))


if __name__ == '__main__':
    main()
