#!/usr/bin/env python3
"""Verification suite: generates small decks, runs asclepios, checks results against closed-form values."""
import base64, glob, math, os, re, struct, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gen_bar import bar
from gen_plate import plate
from gen_beam import beam
from gen_fung import fung, analytic as fung_ref
from gen_contact import contact
BIN = os.path.join(HERE, '..', 'bin', 'asclepios')
OUT = os.path.join(HERE, 'out')
os.makedirs(OUT, exist_ok=True)

def run(name):
    r = subprocess.run([BIN, os.path.join(OUT, name + '.inp'), '-q', '-j', '4'], cwd=OUT, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.strip())

def last(name, field):
    f = sorted(glob.glob(os.path.join(OUT, name + '_vtk', '*.vtu')))[-1]
    s = open(f).read()
    m = re.search(r'Name="%s"[^>]*>\n(.*?)\n</DataArray>' % field, s, re.S)
    b = base64.b64decode(m.group(1)); n = struct.unpack('<I', b[:4])[0]
    return struct.unpack('<%df' % (n // 4), b[4:4 + n])

def history(name):
    rows = open(os.path.join(OUT, name + '_history.csv')).read().strip().splitlines()
    keys = rows[0].split(','); vals = list(map(float, rows[-1].split(',')))
    return dict(zip(keys, vals))

cases = []
def case(fn): cases.append(fn); return fn

def check(label, got, ref, tol):
    err = abs(got - ref) / max(abs(ref), 1e-30)
    ok = err <= tol
    print('  %-46s got %12.6g  ref %12.6g  err %6.2f%%  %s' % (label, got, ref, 100 * err, 'ok' if ok else 'FAIL'))
    return ok

@case
def bar_explicit():
    ok = True
    for et in ('C3D8R', 'C3D8'):
        n = 'bar_' + et.lower(); bar(os.path.join(OUT, n + '.inp'), et); run(n)
        h = history(n)
        ok &= check(et + ' explicit: strain energy', h['ALLIE'], 1.0, 0.01)
        ok &= check(et + ' explicit: energy balance |ETOTAL|/ALLWK', 1 + abs(h['ETOTAL']) / h['ALLWK'], 1.0, 1e-3)
    return ok

@case
def bar_static():
    ok = True
    bar(os.path.join(OUT, 'bars.inp'), 'C3D8', static=True, nlgeom='NO'); run('bars')
    RF = last('bars', 'RF'); fx = sum(RF[3 * i] for i in range(len(RF) // 3) if i % 21 == 20)
    ok &= check('C3D8 static: reaction', fx, 200.0, 1e-4)
    bar(os.path.join(OUT, 'barsnl.inp'), 'C3D8R', static=True, nlgeom='YES', disp=1.0); run('barsnl')
    S = last('barsnl', 'S11')
    ok &= check('C3D8R static nlgeom: Cauchy S11 = E ln(1.1)', sum(S) / len(S), 200000 * math.log(1.1), 2e-3)
    return ok

@case
def plates():
    ok = True
    # beam theory for a plate strip with Poisson ratio 0.3 (anticlastic curvature restrained at the clamp), equal
    # nodal shares of the tip load, explicit with damping: about 1 % from the reference; the exact version of this
    # test is validation/s4r_bending (0.02 %)
    for et, tol in (('S4R', 0.02), ('S3R', 0.06)):
        n = 'plate_' + et.lower(); plate(os.path.join(OUT, n + '.inp'), et); run(n)
        U = last(n, 'U'); tip = U[3 * 20 + 2]
        ok &= check(et + ' cantilever tip deflection', tip, 0.01 * 10 ** 3 / (3 * 200000 * 0.1 ** 3 / 12), tol)
    return ok

@case
def beams():
    ok = True
    beam(os.path.join(OUT, 'bc.inp')); run('bc')
    ok &= check('B31 CIRC cantilever', last('bc', 'U')[3 * 20 + 2], 0.001 * 1000 / (3 * 200000 * math.pi * 0.1 ** 4 / 4), 0.01)
    beam(os.path.join(OUT, 'br.inp'), sec='RECT', dims='0.2, 0.1'); run('br')
    ok &= check('B31 RECT cantilever', last('br', 'U')[3 * 20 + 2], 0.001 * 1000 / (3 * 200000 * 0.2 * 0.1 ** 3 / 12), 0.01)
    beam(os.path.join(OUT, 'bt.inp'), P=0, dT=-10.); run('bt')
    ok &= check('B31 thermal contraction', last('bt', 'U')[3 * 20], -0.1, 1e-3)
    return ok

@case
def fung_membrane():
    ok = True
    l1, l2, t = 1.15, 1.10, 0.5
    for et, ang in (('S4R', 0.0), ('S4R', 30.0), ('S3R', 30.0)):
        n = 'fung_%s_%d' % (et.lower(), ang); fung(os.path.join(OUT, n + '.inp'), et, l1, l2, ang, t); run(n)
        RF = last(n, 'RF'); sxx, syy, l3 = fung_ref(l1, l2, ang)
        fx = RF[3 * 1] + RF[3 * 2]; fy = RF[3 * 2 + 1] + RF[3 * 3 + 1]
        ok &= check('%s Fung biaxial, material at %g deg: sigma_xx' % (et, ang), fx / (t * l3 * l2), sxx, 2e-3)
        ok &= check('%s Fung biaxial, material at %g deg: sigma_yy' % (et, ang), fy / (t * l3 * l1), syy, 2e-3)
    # orientation without its second line: rotation axis = local 1-axis (here in the shell plane), so the projected
    # local 2-axis is the material 1-direction: material axes at 30 + 90 degrees
    n = 'fung_s4r_noaxis'; fung(os.path.join(OUT, n + '.inp'), 'S4R', l1, l2, 30.0, t, rot_line=False); run(n)
    RF = last(n, 'RF'); sxx, syy, l3 = fung_ref(l1, l2, 120.0)
    fx = RF[3 * 1] + RF[3 * 2]; fy = RF[3 * 2 + 1] + RF[3 * 3 + 1]
    ok &= check('S4R Fung biaxial, orientation without rotation line (axes at 120 deg): sigma_xx', fx / (t * l3 * l2), sxx, 2e-3)
    ok &= check('S4R Fung biaxial, orientation without rotation line (axes at 120 deg): sigma_yy', fy / (t * l3 * l1), syy, 2e-3)
    return ok

@case
def contact_impact():
    # plate thrown at a fixed plate (general contact, penalty 0.3 m/dt^2, normal damping 3 % of critical): rebound
    # without energy creation; the damped penalty contact (resolved in 2-3 increments) keeps 70-90 % of the kinetic
    # energy (the shell-to-shell and beam-to-shell impacts of the validation suite agree with Abaqus with this damping)
    contact(os.path.join(OUT, 'impact.inp')); run('impact')
    h = history('impact'); ke0 = 0.5 * 7.8e-9 * 0.1 * 0.8 * 0.8 * 1000.0 ** 2
    ok = check('impact: energy balance |ETOTAL - KE0| / KE0', 1 + abs(h['ETOTAL'] - ke0) / ke0, 1.0, 5e-3)
    r = h['ALLKE'] / ke0
    ok &= check('impact: kinetic energy kept (0.70-1.0)', min(max(r, 0.70), 1.0), r, 1e-12)
    # the energy taken by the contact damper is reported in ALLVD (as the contact damping energy of Abaqus, 10-20 %
    # of the kinetic energy per rebound in single-element Abaqus/Explicit drops); no energy is created
    vd = h['ALLVD'] / ke0
    ok &= check('impact: contact damping in ALLVD (0.05-0.25 KE0)', min(max(vd, 0.05), 0.25), vd, 1e-12)
    m = (h['ALLKE'] + h['ALLIE'] + h['ALLVD']) / ke0
    ok &= check('impact: KE + IE + VD <= KE0', min(m, 1.0), m, 1e-3)
    return ok

@case
def shell_stable_increment():
    # S4R strip of 1 x 2 elements (E = 210000, nu = 0, bulk viscosity 0.06): thin shells, the membrane sets the
    # stable increment whatever the thickness (bending and transverse shear never control: the uniform transverse
    # shear and the curvature / hourglass modes are orthogonal, the bound takes the larger of the two groups), with
    # the element estimate of Abaqus/Explicit L_m / c (sqrt(1 + 0.06^2) - 0.06), L_m^2 = 0.92 L_e^2; Abaqus/Explicit
    # increment for the strip of the self-contact case (same elements, t = 0.3): 1.7462e-7
    def dt_of(t):
        name = 'shelldt_%g' % t
        nodes = ['%d, %g, %g, 0.' % (1 + i + 9 * j, float(i), 2.0 * j) for j in range(3) for i in range(9)]
        els = ['%d, %d, %d, %d, %d' % (1 + i + 8 * j, 1 + i + 9 * j, 2 + i + 9 * j, 11 + i + 9 * j, 10 + i + 9 * j)
               for j in range(2) for i in range(8)]
        deck = ['*Heading', '*Node'] + nodes + ['*Element, type=S4R, elset=ALL'] + els + [
            '*Nset, nset=NALL, elset=ALL', '*Shell Section, elset=ALL, material=STEEL', '%g, 5' % t,
            '*Material, name=STEEL', '*Elastic', '210000., 0.', '*Density', '7.85e-9',
            '*Initial Conditions, type=VELOCITY', 'NALL, 3, 100.',
            '*Step, nlgeom=YES', '*Dynamic, Explicit', ', 4e-6', '*Output, field, number interval=1', '*Node Output', 'U',
            '*Output, history, time interval=1e-6', '*Energy Output', 'ALLKE', '*End Step']
        open(os.path.join(OUT, name + '.inp'), 'w').write('\n'.join(deck) + '\n')
        run(name)
        rows = open(os.path.join(OUT, name + '_history.csv')).read().strip().splitlines()
        keys = rows[0].split(',')
        return float(rows[1].split(',')[keys.index('dt')])  # first increment (later ones cut at the output times)
    d3, d03 = dt_of(0.3), dt_of(0.03)
    ok = check('S4R t = 0.3: increment (Abaqus 1.7462e-7)', d3, 1.7462e-7, 3e-3)
    ok &= check('S4R t = 0.03 / t = 0.3: same increment', d03, d3, 5e-3)
    return ok

@case
def contact_facet_weights():
    # general contact on one clamped S4R plate facet 1 x 1 (t = 0.1), small S4R sheet 0.3 x 0.3 (t = 0.1) dropped
    # flat at 1000 off-centre: the reactions of the four plate nodes follow the bilinear weights of the contact
    # points (Abaqus/Explicit: 0.188 0.063 0.188 0.562 of the impulse, peak reaction 0.5111); light plate (density
    # / 10, the master side limits the penalty, stiffness and sharing from the triangles of the facet): Abaqus
    # peak 0.7013
    def drop(name, rp):
        deck = ['*Heading', '*Node', '1, 0., 0., 0.', '2, 1., 0., 0.', '3, 1., 1., 0.', '4, 0., 1., 0.',
                '11, 0.1, 0.6, 0.12', '12, 0.4, 0.6, 0.12', '13, 0.4, 0.9, 0.12', '14, 0.1, 0.9, 0.12',
                '*Element, type=S4R, elset=PLATE', '1, 1, 2, 3, 4', '*Element, type=S4R, elset=SHEET', '2, 11, 12, 13, 14',
                '*Shell Section, elset=PLATE, material=PM', '0.1, 5', '*Shell Section, elset=SHEET, material=SM', '0.1, 5',
                '*Nset, nset=PN', '1, 2, 3, 4', '*Nset, nset=SN', '11, 12, 13, 14',
                '*Material, name=PM', '*Density', '%g' % rp, '*Elastic', '210000., 0.3',
                '*Material, name=SM', '*Density', '7.85e-09', '*Elastic', '70000., 0.33',
                '*Initial Conditions, type=VELOCITY', 'SN, 3, -1000.', '*Surface Interaction, name=NOFRIC',
                '*Step, name=Drop, nlgeom=YES', '*Dynamic, Explicit', ', 4e-05', '*Boundary', 'PN, ENCASTRE',
                '*Contact', '*Contact Inclusions, ALL EXTERIOR', '*Contact Property Assignment', ' , , NOFRIC',
                '*Output, field, time interval=1e-8', '*Node Output', 'U, RF', '*End Step']
        open(os.path.join(OUT, name + '.inp'), 'w').write('\n'.join(deck) + '\n')
        run(name)
        imp, peak = [0.0] * 4, 0.0
        for f in sorted(glob.glob(os.path.join(OUT, name + '_vtk', '*.vtu'))):
            s = open(f).read()
            def arr(field, fmt):
                m = re.search(r'Name="%s"[^>]*>\n(.*?)\n</DataArray>' % field, s, re.S)
                b = base64.b64decode(m.group(1)); n = struct.unpack('<I', b[:4])[0]
                return struct.unpack('<%d%s' % (n // 4, fmt), b[4:4 + n])
            lab, rf = arr('NodeLabel', 'i'), arr('RF', 'f')
            rz = [rf[3 * lab.index(k) + 2] for k in (1, 2, 3, 4)]
            imp = [a + b for a, b in zip(imp, rz)]
            peak = max(peak, abs(sum(rz)))
        return [a / sum(imp) for a in imp], peak
    w, pk = drop('cfacet', 7.85e-9)
    ok = True
    for k, ref in enumerate((0.1875, 0.0625, 0.1875, 0.5625)):
        ok &= check('plate node %d: share of the impulse' % (k + 1), w[k], ref, 0.03)
    ok &= check('peak reaction (Abaqus 0.5111)', pk, 0.5111, 0.01)
    _, pkl = drop('cfacet_light', 7.85e-10)
    ok &= check('light plate: peak reaction (Abaqus 0.7013)', pkl, 0.7013, 0.03)
    return ok

@case
def contact_shell_beam():
    # general contact between B31 rods and S4R shells (t = 0.5), against Abaqus/Explicit:
    #  - rod (CIRC r = 0.4, L = 1.6) dropped flat at 1000 on a clamped facet 2.5 x 2.5: peak reaction 18.504;
    #  - free 2 x 2 sheet dropped on a clamped rod (CIRC r = 0.4) lying under its middle node line, or under the
    #    middle of two perimeter edges (no node above it): the shell nodes do not contact the rod, the perimeter
    #    edges do (edge to segment, on the edges moved back by the rounding of free edges); force of the first
    #    contact increment 6.191 and 3.810
    def rf_hist(name, deck):
        open(os.path.join(OUT, name + '.inp'), 'w').write('\n'.join(deck) + '\n')
        run(name)
        rows = open(os.path.join(OUT, name + '_rf.csv')).read().strip().splitlines()[1:]
        return [float(r.split(',')[3]) for r in rows]
    def common(fixed, moving):
        return ['*Element, type=B31, elset=ROD', '10, 11, 12', '*Shell Section, elset=PLATE, material=M', '0.5, 5',
                '*Beam Section, elset=ROD, material=M, section=CIRC', '0.4', '0., 0., -1.',
                '*Nset, nset=FIXED', fixed, '*Nset, nset=MOVING', moving,
                '*Material, name=M', '*Density', '7.85e-09', '*Elastic', '210000., 0.3',
                '*Initial Conditions, type=VELOCITY', 'MOVING, 3, -1000.', '*Surface Interaction, name=NOFRIC',
                '*Step, name=Drop, nlgeom=YES', '*Dynamic, Explicit', ', 4e-05', '*Boundary', 'FIXED, ENCASTRE',
                '*Contact', '*Contact Inclusions, ALL EXTERIOR', '*Contact Property Assignment', ' , , NOFRIC',
                '*Output, field, number interval=1', '*Node Output', 'U', '*Output, history, frequency=1',
                '*Node Output, nset=FIXED', 'RF1, RF2, RF3', '*End Step']
    rod = ['*Heading', '*Node', '1, 0., 0., 0.', '2, 2.5, 0., 0.', '3, 2.5, 2.5, 0.', '4, 0., 2.5, 0.',
           '11, 0.45, 1.25, 0.67', '12, 2.05, 1.25, 0.67', '*Element, type=S4R, elset=PLATE', '1, 1, 2, 3, 4']
    rf = rf_hist('cshb_rod', rod + common('1, 2, 3, 4', '11, 12'))
    ok = check('rod on a facet: peak reaction (Abaqus 18.504)', max(abs(x) for x in rf), 18.504, 0.01)
    for name, y, ref in (('cshb_node', 1.25, 6.191), ('cshb_edge', 0.625, 3.810)):
        deck = ['*Heading', '*Node', '11, -0.8, %g, 0.' % y, '12, 3.3, %g, 0.' % y]
        deck += ['%d, %g, %g, 0.67' % (1 + i + 3 * j, 1.25 * i, 1.25 * j) for j in range(3) for i in range(3)]
        deck += ['*Element, type=S4R, elset=PLATE'] + ['%d, %d, %d, %d, %d' % (1 + i + 2 * j, 1 + i + 3 * j, 2 + i + 3 * j,
                                                                            5 + i + 3 * j, 4 + i + 3 * j)
                                                     for j in range(2) for i in range(2)]
        rf = rf_hist(name, deck + common('11, 12', '1, 2, 3, 4, 5, 6, 7, 8, 9'))
        first = next((abs(x) for x in rf if abs(x) > 1e-9), 0.0)
        ok &= check('sheet on a rod, rod at y = %g: first force' % y, first, ref, 0.02)
    return ok

@case
def contact_solid_faces():
    # general contact with the exterior faces of C3D8R solids, against Abaqus/Explicit: free alu cube 2 x 2 x 2
    # dropped at 1000 on the top face (4 x 4) of a steel C3D8R base,
    #  - base fully clamped: contact force of the first increments 11.43, 22.19, 29.65, peak 32.97;
    #  - base clamped at its bottom only, cube nodes on the diagonal of the base face (two triangles): no contact
    #    lost at the third increment (28.84), peak 32.00;
    #  - global estimate of the increment with the boundary conditions ignored (fully clamped base, heavy cube):
    #    2.728e-7 before contact
    def deck(name, fixb, rho):
        nodes = ['1, 0., 0., 0.', '2, 4., 0., 0.', '3, 4., 4., 0.', '4, 0., 4., 0.', '5, 0., 0., 2.', '6, 4., 0., 2.',
                 '7, 4., 4., 2.', '8, 0., 4., 2.']
        nodes += ['%d, %g, %g, %g' % (11 + i + 4 * k, x, y, 2.01 + 2 * k) for k in range(2)
                  for i, (x, y) in enumerate(((1, 1), (3, 1), (3, 3), (1, 3)))]
        d = ['*Heading', '*Node'] + nodes + [
            '*Element, type=C3D8R, elset=BASE', '1, 1, 2, 3, 4, 5, 6, 7, 8',
            '*Element, type=C3D8R, elset=CUBE', '2, 11, 12, 13, 14, 15, 16, 17, 18',
            '*Solid Section, elset=BASE, material=BM', ',', '*Solid Section, elset=CUBE, material=CM', ',',
            '*Nset, nset=BN', '1, 2, 3, 4, 5, 6, 7, 8' if fixb else '1, 2, 3, 4', '*Nset, nset=CN',
            '11, 12, 13, 14, 15, 16, 17, 18', '*Material, name=BM', '*Density', '7.85e-09', '*Elastic', '210000., 0.3',
            '*Material, name=CM', '*Density', '%g' % rho, '*Elastic', '70000., 0.33',
            '*Initial Conditions, type=VELOCITY', 'CN, 3, -1000.', '*Surface Interaction, name=NOFRIC',
            '*Step, name=Drop, nlgeom=YES', '*Dynamic, Explicit', ', 1.4e-05', '*Boundary', 'BN, ENCASTRE',
            '*Contact', '*Contact Inclusions, ALL EXTERIOR', '*Contact Property Assignment', ' , , NOFRIC',
            '*Output, field, number interval=1', '*Node Output', 'U', '*Output, history, frequency=1',
            '*Node Output, nset=CN', 'U1, U2, U3', '*End Step']
        open(os.path.join(OUT, name + '.inp'), 'w').write('\n'.join(d) + '\n')
        run(name)
        by = {}
        for r in open(os.path.join(OUT, name + '_mon.csv')).read().strip().splitlines()[1:]:
            f = r.split(','); by.setdefault(float(f[0]), {})[int(f[1])] = float(f[4])
        t = sorted(by)
        m = rho * 8 / 8
        F = []  # contact force on the cube: momentum change, central differences of the displacements
        for i in range(1, len(t) - 1):
            a = sum((by[t[i + 1]][n] - by[t[i]][n]) / (t[i + 1] - t[i]) - (by[t[i]][n] - by[t[i - 1]][n]) / (t[i] - t[i - 1])
                    for n in range(11, 19)) / (0.5 * (t[i + 1] - t[i - 1]))
            F.append(m * a)
        rows = open(os.path.join(OUT, name + '_history.csv')).read().strip().splitlines()
        dt = float(rows[5].split(',')[rows[0].split(',').index('dt')])
        i0 = next(i for i, x in enumerate(F) if x > 1e-3)
        return F[i0:i0 + 10], dt
    F, _ = deck('csol_fixed', 1, 2.7e-9)
    ok = check('cube on a clamped base: 1st increment force', F[0], 11.43, 0.01)
    ok &= check('cube on a clamped base: 3rd increment force', F[2], 29.65, 0.01)
    ok &= check('cube on a clamped base: peak force', max(F), 32.97, 0.01)
    F, _ = deck('csol_diag', 0, 2.7e-9)
    ok &= check('cube on the base diagonal: 3rd increment force', F[2], 28.84, 0.01)
    ok &= check('cube on the base diagonal: peak force', max(F), 32.00, 0.01)
    _, dt = deck('csol_heavy', 1, 2.7e-8)
    ok &= check('global increment, clamped base ignored (Abaqus)', dt, 2.728e-7, 0.01)
    return ok

@case
def global_increment_mixed():
    # global estimate of the increment in a model with solids and shells: made on the solids alone, applied to their
    # element bound (alu C3D8R block 6 x 6 x 4 of 2 mm bricks, element bound 2.158e-7, next to a steel S4R plate:
    # Abaqus/Explicit 2.585e-7)
    nodes, ids = [], {}
    for k in range(3):
        for j in range(4):
            for i in range(4):
                ids[i, j, k] = 1 + len(nodes); nodes.append('%d, %g, %g, %g' % (ids[i, j, k], 2. * i, 2. * j, 1. + 2. * k))
    els = []
    for k in range(2):
        for j in range(3):
            for i in range(3):
                c = [ids[i, j, k], ids[i + 1, j, k], ids[i + 1, j + 1, k], ids[i, j + 1, k], ids[i, j, k + 1],
                     ids[i + 1, j, k + 1], ids[i + 1, j + 1, k + 1], ids[i, j + 1, k + 1]]
                els.append('%d, %s' % (1 + len(els), ', '.join(map(str, c))))
    d = ['*Heading', '*Node'] + nodes + ['101, 20., 0., 0.', '102, 40., 0., 0.', '103, 40., 20., 0.', '104, 20., 20., 0.',
         '*Element, type=C3D8R, elset=BLOCK'] + els + ['*Element, type=S4R, elset=PLATE', '100, 101, 102, 103, 104',
         '*Solid Section, elset=BLOCK, material=ALU', ',', '*Shell Section, elset=PLATE, material=ST', '0.5, 5',
         '*Nset, nset=BN, elset=BLOCK', '*Material, name=ALU', '*Density', '2.7e-09', '*Elastic', '70000., 0.33',
         '*Material, name=ST', '*Density', '7.85e-09', '*Elastic', '210000., 0.3',
         '*Initial Conditions, type=VELOCITY', 'BN, 3, -1000.', '*Step, nlgeom=YES', '*Dynamic, Explicit', ', 5e-06',
         '*Output, field, number interval=1', '*Node Output', 'U', '*Output, history, frequency=1', '*Energy Output',
         'ALLKE', '*End Step']
    open(os.path.join(OUT, 'gdt_mixed.inp'), 'w').write('\n'.join(d) + '\n')
    run('gdt_mixed')
    rows = open(os.path.join(OUT, 'gdt_mixed_history.csv')).read().strip().splitlines()
    dt = float(rows[10].split(',')[rows[0].split(',').index('dt')])
    return check('C3D8R block next to an S4R plate: increment', dt, 2.585e-7, 0.01)

@case
def beam_stable_increment():
    # free single B31 elements (E = 210000, nu = 0.3, default bulk viscosity): element estimate of Abaqus/Explicit,
    # larger of the opposite-rotation and coupled-shear modes with the cross pairing of the bending axes (RECT
    # sections), axial wave with ISOTROPIC rotary inertia
    cases = [('CIRC', '1.', 0.5, '', 8.491654e-08), ('CIRC', '1.', 10.0, '', 8.939024e-07),
             ('PIPE', '3., 0.5', 10.0, '', 1.516126e-06), ('RECT', '4., 0.5', 0.5, '', 1.073992e-08),
             ('RECT', '4., 0.5', 50.0, '', 5.497078e-07), ('RECT', '2., 1.', 1.0, ', rotary inertia=ISOTROPIC', 1.802679e-07)]
    ok = True
    for k, (sec, dims, L, opt, ref) in enumerate(cases):
        name = 'beamdt_%d' % k
        deck = ['*Heading', '*Node', '1, 0., 0., 0.', '2, %g, 0., 0.' % L, '*Element, type=B31, elset=B', '1, 1, 2',
                '*Nset, nset=NALL, elset=B', '*Beam Section, elset=B, material=M, section=%s%s' % (sec, opt), dims,
                '0., 0., -1.', '*Material, name=M', '*Density', '7.85e-09', '*Elastic', '210000., 0.3',
                '*Initial Conditions, type=VELOCITY', 'NALL, 2, 1.',
                '*Step, nlgeom=YES', '*Dynamic, Explicit', ', 1e-4', '*Output, field, number interval=1', '*Node Output',
                'U', '*Output, history, time interval=5e-5', '*Energy Output', 'ALLKE', '*End Step']
        open(os.path.join(OUT, name + '.inp'), 'w').write('\n'.join(deck) + '\n')
        run(name)
        rows = open(os.path.join(OUT, name + '_history.csv')).read().strip().splitlines()
        dt = float(rows[1].split(',')[rows[0].split(',').index('dt')])
        ok &= check('B31 %s L = %g%s: increment' % (sec + ' ' + dims.replace(' ', ''), L, ' iso' if opt else ''), dt, ref, 3e-3)
    return ok

@case
def output_frequency_and_paths():
    # *OUTPUT, FIELD, FREQUENCY=n in a static step (a frame every n increments and at the end of the step; ignored
    # in explicit steps, as Abaqus/Explicit does), job written through an absolute -o path in another directory:
    # every frame listed in the .pvd exists next to it
    src = os.path.join(OUT, 'freq.inp')
    bar(src, 'C3D8R', static=True, nlgeom='YES', disp=1.0)
    deck = open(src).read()
    deck = re.sub(r'\*Output, field[^\n]*', '*Output, field, frequency=3', deck)
    deck = re.sub(r'(\*Static[^\n]*\n)[^\n]*', r'\g<1>0.05, 1.', deck)
    open(src, 'w').write(deck)
    job = os.path.join(OUT, 'sub', 'freqjob')
    os.makedirs(os.path.dirname(job), exist_ok=True)
    r = subprocess.run([BIN, src, '-q', '-o', job], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.strip())
    files = re.findall(r'file="([^"]+)"', open(job + '.pvd').read())
    ninc = int(re.search(r'\((\d+) increments', open(job + '.dat').read()).group(1))
    ok = check('field frequency 3: number of frames', len(files), 1 + ninc // 3 + (1 if ninc % 3 else 0), 1e-12)
    present = sum(os.path.exists(os.path.join(os.path.dirname(job), f)) for f in files)
    ok &= check('frames found next to the .pvd (absolute -o)', present, len(files), 1e-12)
    return ok

@case
def reaction_history():
    # *OUTPUT, HISTORY, FREQUENCY=1 with *NODE OUTPUT, NSET=, RF1..RF3 in an explicit step: total reaction of the set
    # at every increment in <job>_rf.csv, equal to the field output reaction at the end of the step
    src = os.path.join(OUT, 'rfhist.inp')
    bar(src, 'C3D8R')
    deck = open(src).read().replace('*End Step', '*Output, history, frequency=1\n*Node Output, nset=XL\nRF1, RF2, RF3\n*End Step')
    open(src, 'w').write(deck)
    run('rfhist')
    rows = open(os.path.join(OUT, 'rfhist_rf.csv')).read().strip().splitlines()
    ninc = int(re.search(r'\((\d+) increments', open(os.path.join(OUT, 'rfhist.dat')).read()).group(1))
    ok = check('reaction history: one line per increment (+ initial)', len(rows) - 1, ninc + 1, 1e-12)
    RF = last('rfhist', 'RF')
    xl = [int(x) for x in re.search(r'\*Nset, nset=XL\n([^\n]+)', deck).group(1).split(',')]
    f = sorted(glob.glob(os.path.join(OUT, 'rfhist_vtk', '*.vtu')))[-1]
    m = re.search(r'Name="NodeLabel"[^>]*>\n(.*?)\n</DataArray>', open(f).read(), re.S)
    bl = base64.b64decode(m.group(1)); nb = struct.unpack('<I', bl[:4])[0]
    lab = list(struct.unpack('<%di' % (nb // 4), bl[4:4 + nb]))
    fx = sum(RF[3 * lab.index(n)] for n in xl)
    ok &= check('reaction history: last value = field reaction', float(rows[-1].split(',')[1]), fx, 1e-5)
    return ok

def tied_bars(path, static):
    # two C3D8 bars with non-matching meshes (2 x 2 x 2 and 3 x 3 x 3 elements, section 4 x 4) tied at x = 10;
    # static: uniform tension 10 MPa on the x = 20 face with symmetry conditions; explicit: free rigid motion
    nodes, els, ids = [], [], {}
    def block(tag, x0, nx, ny):
        for k in range(ny + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    ids[tag, i, j, k] = 1 + len(nodes)
                    nodes.append('%d, %.10g, %.10g, %.10g' % (ids[tag, i, j, k], x0 + 10. * i / nx, 4. * j / ny, 4. * k / ny))
        for k in range(ny):
            for j in range(ny):
                for i in range(nx):
                    c = [ids[tag, i, j, k], ids[tag, i + 1, j, k], ids[tag, i + 1, j + 1, k], ids[tag, i, j + 1, k],
                         ids[tag, i, j, k + 1], ids[tag, i + 1, j, k + 1], ids[tag, i + 1, j + 1, k + 1], ids[tag, i, j + 1, k + 1]]
                    els.append((tag, i, j, k, '%d, %s' % (1 + len(els), ', '.join(map(str, c)))))
    block('A', 0., 2, 2)
    block('B', 10., 3, 3)
    ea = [e[4] for e in els if e[0] == 'A']; eb = [e[4] for e in els if e[0] == 'B']
    master = ['%s, S4' % e[4].split(',')[0] for e in els if e[0] == 'A' and e[1] == 1]
    slave = ['%s, S6' % e[4].split(',')[0] for e in els if e[0] == 'B' and e[1] == 0]
    d = ['*Heading', '*Node'] + nodes + ['*Element, type=C3D8, elset=BA'] + ea + ['*Element, type=C3D8, elset=BB'] + eb
    d += ['*Solid Section, elset=BA, material=ST', ',', '*Solid Section, elset=BB, material=ST', ',',
          '*Surface, type=ELEMENT, name=M'] + master + ['*Surface, type=ELEMENT, name=S'] + slave
    d += ['*Tie, name=T, adjust=yes', 'S, M', '*Material, name=ST', '*Density', '7.85e-09', '*Elastic', '210000., 0.3']
    if static:
        x0 = [str(ids['A', 0, j, k]) for j in range(3) for k in range(3)]
        y0 = [str(n) for key, n in ids.items() if key[2] == 0 and key[0] == 'A'] + [str(n) for key, n in ids.items() if key[2] == 0 and key[0] == 'B' and key[1] > 0]
        z0 = [str(n) for key, n in ids.items() if key[3] == 0 and key[0] == 'A'] + [str(n) for key, n in ids.items() if key[3] == 0 and key[0] == 'B' and key[1] > 0]
        load = []
        for j in range(4):
            for k in range(4):
                share = (1 if j in (0, 3) else 2) * (1 if k in (0, 3) else 2) / 36.
                load.append('%d, 1, %.10g' % (ids['B', 3, j, k], 160. * share))
        d += ['*Nset, nset=X0'] + [', '.join(x0)] + ['*Nset, nset=Y0'] + [', '.join(y0[i:i + 12]) for i in range(0, len(y0), 12)]
        d += ['*Nset, nset=Z0'] + [', '.join(z0[i:i + 12]) for i in range(0, len(z0), 12)]
        d += ['*Step, nlgeom=NO', '*Static', '1., 1.', '*Boundary', 'X0, 1', 'Y0, 2', 'Z0, 3', '*Cload'] + load
        d += ['*Output, field, number interval=1', '*Node Output', 'U, RF', '*Element Output', 'S', '*End Step']
    else:
        d += ['*Nset, nset=ALL, elset=BA', '*Nset, nset=ALLB, elset=BB', '*Initial Conditions, type=VELOCITY',
              'ALL, 1, 100.', 'ALL, 2, -50.', 'ALLB, 1, 100.', 'ALLB, 2, -50.', '*Step, nlgeom=YES', '*Dynamic, Explicit',
              ', 2e-05', '*Output, field, number interval=1', '*Node Output', 'U, V', '*End Step']
    open(path, 'w').write('\n'.join(d) + '\n')
    return ids

@case
def surface_tie():
    # *TIE of non-matching meshes: patch test (uniform tension passes the tied interface exactly, implicit surface
    # formulation) and free rigid motion of the tied bars (explicit point formulation: no strain energy, momentum kept)
    ok = True
    ids = tied_bars(os.path.join(OUT, 'tie_patch.inp'), True)
    run('tie_patch')
    S = last('tie_patch', 'S11')
    ok &= check('tied non-matching bars: min S11', min(S), 10.0, 1e-6)
    ok &= check('tied non-matching bars: max S11', max(S), 10.0, 1e-6)
    U = last('tie_patch', 'U')
    tip = ids['B', 3, 3, 3] - 1
    ok &= check('tied non-matching bars: tip u1 = sigma L / E', U[3 * tip], 10. * 20. / 210000., 1e-6)
    tied_bars(os.path.join(OUT, 'tie_rigid.inp'), False)
    run('tie_rigid')
    h = history('tie_rigid')
    m = 7.85e-9 * 16 * 20
    ok &= check('tied bars, rigid motion: kinetic energy kept', h['ALLKE'], 0.5 * m * (100. ** 2 + 50. ** 2), 1e-9)
    ok &= check('tied bars, rigid motion: strain energy / kinetic', 1 + h['ALLIE'] / h['ALLKE'], 1.0, 1e-9)
    return ok

def coupled_bar(path, kind, static):
    # C3D8 bar 10 x 4 x 4 (5 x 2 x 2 elements), end face x = 10 coupled to a reference node at its centre
    nodes, ids = [], {}
    for k in range(3):
        for j in range(3):
            for i in range(6):
                ids[i, j, k] = 1 + len(nodes); nodes.append('%d, %g, %g, %g' % (ids[i, j, k], 2. * i, 2. * j, 2. * k))
    els = []
    for k in range(2):
        for j in range(2):
            for i in range(5):
                c = [ids[i, j, k], ids[i + 1, j, k], ids[i + 1, j + 1, k], ids[i, j + 1, k], ids[i, j, k + 1],
                     ids[i + 1, j, k + 1], ids[i + 1, j + 1, k + 1], ids[i, j + 1, k + 1]]
                els.append('%d, %s' % (1 + len(els), ', '.join(map(str, c))))
    endf = ['%d, S4' % (1 + i + 5 * (j + 2 * k)) for i in (4,) for j in range(2) for k in range(2)]
    d = ['*Heading', '*Node'] + nodes + ['100, 10., 2., 2.', '*Element, type=C3D8, elset=B'] + els
    d += ['*Solid Section, elset=B, material=ST', ',', '*Surface, type=ELEMENT, name=E'] + endf
    d += ['*Coupling, constraint name=C, ref node=100, surface=E', '*Kinematic' if kind == 'kinematic' else '*Distributing']
    d += ['*Material, name=ST', '*Density', '7.85e-09', '*Elastic', '210000., 0.3']
    if static:
        x0 = [str(ids[0, j, k]) for j in range(3) for k in range(3)]
        y0 = [str(ids[i, 0, k]) for i in range(5) for k in range(3)]
        z0 = [str(ids[i, j, 0]) for i in range(5) for j in range(3)]
        d += ['*Nset, nset=X0', ', '.join(x0), '*Nset, nset=Y0', ', '.join(y0), '*Nset, nset=Z0', ', '.join(z0)]
        d += ['*Step, nlgeom=NO', '*Static', '1., 1.', '*Boundary', 'X0, 1', 'Y0, 2', 'Z0, 3', '*Cload', '100, 1, 160.',
              '*Output, field, number interval=1', '*Node Output', 'U, RF', '*Element Output', 'S', '*End Step']
    else:
        d += ['*Nset, nset=ALL, elset=B', '*Initial Conditions, type=VELOCITY', 'ALL, 1, 100.', 'ALL, 3, 40.',
              '100, 1, 100.', '100, 3, 40.',
              '*Step, nlgeom=YES', '*Dynamic, Explicit', ', 2e-05', '*Output, field, number interval=1', '*Node Output',
              'U, V', '*End Step']
    open(path, 'w').write('\n'.join(d) + '\n')

@case
def node_couplings():
    # *COUPLING: a force on the reference node of a distributing coupling is spread as a uniform traction (end face
    # weighted by the tributary areas: uniform stress, exact); a bar moving freely with a massless kinematic
    # reference node keeps its kinetic energy without straining
    ok = True
    coupled_bar(os.path.join(OUT, 'cpl_dist.inp'), 'distributing', True)
    run('cpl_dist')
    S = last('cpl_dist', 'S11')
    ok &= check('distributing coupling, axial force: min S11', min(S), 10.0, 1e-6)
    ok &= check('distributing coupling, axial force: max S11', max(S), 10.0, 1e-6)
    U = last('cpl_dist', 'U')
    ok &= check('distributing coupling: reference node u1 = sigma L / E', U[3 * 54], 10. * 10. / 210000., 1e-6)
    coupled_bar(os.path.join(OUT, 'cpl_kin.inp'), 'kinematic', False)
    run('cpl_kin')
    h = history('cpl_kin')
    m = 7.85e-9 * 160
    ok &= check('kinematic coupling, rigid motion: kinetic energy kept', h['ALLKE'], 0.5 * m * (100. ** 2 + 40. ** 2), 1e-9)
    ok &= check('kinematic coupling, rigid motion: strain energy / kinetic', 1 + h['ALLIE'] / h['ALLKE'], 1.0, 1e-9)
    return ok

def rigid_deck(path, kind, ccw=True):
    # S4R sheet 4 x 4 (2 x 2 elements) falling at 1 m/s from 0.5 mm above a rigid surface: R3D4 plane at z = 0, or the
    # upper half of an analytical cylinder of radius 5 (axis y, top at z = 0) described counterclockwise (normal into
    # the cylinder: the sheet comes from behind) or clockwise (normal outwards)
    nodes = []
    for j in range(3):
        for i in range(3):
            nodes.append('%d, %g, %g, 0.5' % (1 + i + 3 * j, -2. + 2. * i, -2. + 2. * j))
    els = ['%d, %d, %d, %d, %d' % (1 + i + 2 * j, 1 + i + 3 * j, 2 + i + 3 * j, 5 + i + 3 * j, 4 + i + 3 * j)
           for j in range(2) for i in range(2)]
    d = ['*Heading', '*Node'] + nodes + ['100, 0., 0., -5.']
    if kind == 'r3d4':
        d += ['101, -10., -10., 0.', '102, 10., -10., 0.', '103, 10., 10., 0.', '104, -10., 10., 0.',
              '*Element, type=R3D4, elset=RIG', '50, 101, 102, 103, 104', '*Rigid Body, ref node=100, elset=RIG']
    d += ['*Element, type=S4R, elset=SH'] + els + ['*Shell Section, elset=SH, material=AL', '0.5, 5']
    if kind == 'cyl':
        prof = ['START, 5., 0.', 'CIRCL, 0., 5., 0., 0.', 'CIRCL, -5., 0., 0., 0.'] if ccw else \
               ['START, -5., 0.', 'CIRCL, 0., 5., 0., 0.', 'CIRCL, 5., 0., 0., 0.']
        d += ['*Surface, type=CYLINDER, name=CY', '0., 0., -5., 1., 0., -5.', '0., 1., -5.'] + prof
        d += ['*Rigid Body, ref node=100, analytical surface=CY']
    d += ['*Nset, nset=ALL, elset=SH', '*Material, name=AL', '*Density', '2.7e-09', '*Elastic', '70000., 0.33',
          '*Initial Conditions, type=VELOCITY', 'ALL, 3, -1000.', '*Step, nlgeom=YES', '*Dynamic, Explicit', ', 0.002',
          '*Boundary', '100, ENCASTRE', '*Contact', '*Contact Inclusions, ALL EXTERIOR',
          '*Output, field, number interval=4', '*Node Output', 'U, V', '*End Step']
    open(path, 'w').write('\n'.join(d) + '\n')

@case
def rigid_surfaces():
    # rigid surfaces in general contact: a sheet bounces on an R3D4 plane held by its reference node (work of the
    # contact returned, energy balance); one-sided analytical cylinder: passed through from behind without any
    # contact, bounced on from the front
    ok = True
    m = 2.7e-9 * 16 * 0.5
    ke0 = 0.5 * m * 1000. ** 2
    rigid_deck(os.path.join(OUT, 'rig_plane.inp'), 'r3d4'); run('rig_plane')
    h = history('rig_plane')
    ok &= check('sheet on a rigid R3D4 plane: energy kept by the bounce > 75 %', 1.0 if h['ALLKE'] + h['ALLIE'] > 0.75 * ke0 else 0.0, 1.0, 1e-12)
    ok &= check('sheet on a rigid R3D4 plane: ETOTAL = KE0', h['ETOTAL'], ke0, 1e-3)
    U = last('rig_plane', 'U')
    ok &= check('sheet on a rigid R3D4 plane: bounced (u3 of the centre > 0)', 1.0 if U[3 * 4 + 2] > 0 else 0.0, 1.0, 1e-12)
    rigid_deck(os.path.join(OUT, 'rig_cyl_behind.inp'), 'cyl', True); run('rig_cyl_behind')
    U = last('rig_cyl_behind', 'U')
    ok &= check('analytical cylinder seen from behind: free fall', U[3 * 4 + 2], -2.0, 1e-9)
    rigid_deck(os.path.join(OUT, 'rig_cyl_front.inp'), 'cyl', False); run('rig_cyl_front')
    U = last('rig_cyl_front', 'U')
    ok &= check('analytical cylinder seen from the front: bounced', 1.0 if U[3 * 4 + 2] > 0 else 0.0, 1.0, 1e-12)
    return ok

def hgo_energy(F, C10, D, k1, k2, ka, fibres):
    # Holzapfel-Gasser-Ogden strain energy (pure Python, for the numerical derivative of the test below)
    C = [[sum(F[k][i] * F[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    J = (F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1]) - F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0]) +
         F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]))
    j23 = J ** (-2. / 3.)
    I1 = j23 * (C[0][0] + C[1][1] + C[2][2])
    W = C10 * (I1 - 3) + ((J * J - 1) / 2 - math.log(J)) / D
    for a in fibres:
        I4 = j23 * sum(a[i] * C[i][j] * a[j] for i in range(3) for j in range(3))
        E = ka * (I1 - 3) + (1 - 3 * ka) * (I4 - 1)
        if E > 0:
            W += k1 / (2 * k2) * (math.exp(k2 * E * E) - 1)
    return W, J

@case
def hgo_homogeneous():
    # Holzapfel-Gasser-Ogden C3D8R cube under a prescribed homogeneous deformation (stretch, shear), fibres at +-30
    # degrees in the (x, y) plane: Cauchy stress = (1/J) dW/dF F^T from a central-difference derivative of the energy
    Fm = [[1.2, 0.1, 0.0], [0.05, 0.95, 0.0], [0.0, 0.0, 0.9]]
    pars = (0.05, 0.02, 2.0, 10.0, 0.1)
    c, s_ = math.cos(math.radians(30)), math.sin(math.radians(30))
    fib = [(c, s_, 0.0), (c, -s_, 0.0)]
    X = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    d = ['*Heading', '*Node'] + ['%d, %g, %g, %g' % (i + 1, *x) for i, x in enumerate(X)]
    d += ['*Element, type=C3D8R, elset=E', '1, 1, 2, 3, 4, 5, 6, 7, 8', '*Orientation, name=O, local directions=2',
          '1., 0., 0., 0., 1., 0.', '3, 0.', '%.12g, %.12g, 0.' % fib[0][:2], '%.12g, %.12g, 0.' % fib[1][:2],
          '*Solid Section, elset=E, material=H, orientation=O', ',', '*Material, name=H', '*Density', '1e-09',
          '*Anisotropic Hyperelastic, holzapfel, local directions=2', '%g, %g, %g, %g, %g' % pars,
          '*Step, nlgeom=YES, inc=1000', '*Static', '0.1, 1., 1e-05, 0.1', '*Boundary']
    for i, x in enumerate(X):
        for k in range(3):
            d.append('%d, %d, %d, %.12g' % (i + 1, k + 1, k + 1, sum(Fm[k][j] * x[j] for j in range(3)) - x[k]))
    d += ['*Output, field, number interval=1', '*Node Output', 'U', '*Element Output', 'S', '*End Step']
    open(os.path.join(OUT, 'hgo_hom.inp'), 'w').write('\n'.join(d) + '\n')
    run('hgo_hom')
    h = 1e-6
    P = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            Fp = [r[:] for r in Fm]; Fp[i][j] += h
            Fn = [r[:] for r in Fm]; Fn[i][j] -= h
            P[i][j] = (hgo_energy(Fp, *pars, fib)[0] - hgo_energy(Fn, *pars, fib)[0]) / (2 * h)
    J = hgo_energy(Fm, *pars, fib)[1]
    sig = [[sum(P[i][k] * Fm[j][k] for k in range(3)) / J for j in range(3)] for i in range(3)]
    ok = True
    for lab, (i, j) in (('S11', (0, 0)), ('S22', (1, 1)), ('S33', (2, 2)), ('S12', (0, 1))):
        ok &= check('HGO homogeneous deformation: ' + lab, last('hgo_hom', lab)[0], sig[i][j], 2e-4)
    return ok

def sma_uniaxial_ref(eps_path, EA=60000., EM=30000., eL=0.05, Rs=520., Rf=600., RsR=300., RfR=200.):
    # uniaxial response of the superelastic model (smooth cubic kinetics between the start and finish stresses),
    # integrated in closed form (log strains)
    xi, sig, out = 0.0, 0.0, []
    E = lambda x: EA + x * (EM - EA)
    def z(F, s, f):
        t = max(-1.0, min(1.0, 2 * (F - 0.5 * (s + f)) / (f - s)))
        return 0.5 + 0.75 * t - 0.25 * t ** 3
    def bis(g, lo, hi):
        glo = g(lo)
        for _ in range(200):
            m = 0.5 * (lo + hi); gm = g(m)
            if (gm > 0) == (glo > 0): lo, glo = m, gm
            else: hi = m
        return 0.5 * (lo + hi)
    for e in eps_path:
        S = lambda x: E(x) * (e - eL * x)
        st = S(xi)
        if st > sig and st > Rs and xi < 1:
            z0, xn = z(sig, Rs, Rf), xi
            g = lambda x: (1 - x) * (1 - z0) - (1 - xn) * (1 - z(S(x), Rs, Rf))
            xi = bis(g, xn, 1.0) if g(1.0) < 0 else 1.0
        elif st < sig and st < RsR and xi > 0:
            z0, xn = z(sig, RsR, RfR), xi
            g = lambda x: x * (1 - z0) - xn * (1 - z(S(x), RsR, RfR))
            xi = bis(g, 0.0, xn) if g(0.0) < 0 else 0.0
        sig = S(xi)
        out.append((sig, xi))
    return out

@case
def superelastic_uniaxial():
    # superelastic NiTi cube (ABQ_SUPER_ELASTIC constants) stretched to 7 % and unloaded with homogeneous uniaxial
    # stress: transformation plateaux of the closed-form uniaxial law, lateral contraction of the Drucker-Prager flow
    X = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    d = ['*Heading', '*Node'] + ['%d, %g, %g, %g' % (i + 1, *x) for i, x in enumerate(X)]
    d += ['*Element, type=C3D8R, elset=E', '1, 1, 2, 3, 4, 5, 6, 7, 8', '*Nset, nset=ALL, elset=E',
          '*Solid Section, elset=E, material=ABQ_SUPER_ELASTIC_N3D', ',', '*Material, name=ABQ_SUPER_ELASTIC_N3D',
          '*Density', '6.5e-09', '*Depvar', '31', '*User Material, constants=15',
          '60000., 0.3, 30000., 0.3, 0.05, 6.7, 520., 600.', '37., 6.7, 300., 200., 600., 0.05, 0.',
          '*Initial Conditions, type=TEMPERATURE', 'ALL, 37.']
    fixb = ['*Boundary', '1, 1, 3', '4, 1, 1', '4, 3, 3', '5, 1, 2', '8, 1, 1', '2, 2, 3', '3, 3, 3', '6, 2, 2']
    out = ['*Output, field, number interval=10', '*Node Output', 'U', '*Element Output', 'S', '*End Step']
    for name, target in (('Load', 0.07), ('Unload', 0.0)):
        d += ['*Step, name=%s, nlgeom=YES, inc=1000' % name, '*Static', '0.01, 1., 1e-05, 0.01']
        d += (fixb if name == 'Load' else ['*Boundary'])
        d += ['%d, 1, 1, %g' % (n, target) for n in (2, 3, 6, 7)] + out
    open(os.path.join(OUT, 'sma_uni.inp'), 'w').write('\n'.join(d) + '\n')
    run('sma_uni')
    frames = sorted(glob.glob(os.path.join(OUT, 'sma_uni_vtk', '*.vtu')))
    def frame(f, field):
        s_ = open(f).read()
        m = re.search(r'Name="%s"[^>]*>\n(.*?)\n</DataArray>' % field, s_, re.S)
        b = base64.b64decode(m.group(1)); n = struct.unpack('<I', b[:4])[0]
        return struct.unpack('<%df' % (n // 4), b[4:4 + n])
    path = [math.log(1 + 0.07 * k / 10) for k in range(11)] + [math.log(1 + 0.07 * (10 - k) / 10) for k in range(1, 11)]
    ref = sma_uniaxial_ref(path)
    ok = True
    for k in (3, 10, 13, 17, 19):  # plateau, end of loading, elastic unloading, reverse plateau
        ok &= check('superelastic uniaxial: S11 at frame %d' % k, frame(frames[k], 'S11')[0], ref[k][0], 2e-3)
    tb = 3 * (600 - 520) / (600 + 520)
    sig, xi = ref[10]
    lat = -0.3 * sig / (60000 + xi * (30000 - 60000)) + (-0.5 + tb / 3) / (1 + tb / 3) * 0.05 * xi
    ok &= check('superelastic uniaxial: lateral strain at 7 %', math.log(1 + frame(frames[10], 'U')[3 * 6 + 1]), lat, 2e-3)
    return ok

if __name__ == '__main__':
    allok = True
    for c in cases:
        print(c.__name__)
        try: allok &= c()
        except Exception as e: print('  ERROR', e); allok = False
    print('ALL PASSED' if allok else 'SOME TESTS FAILED')
    sys.exit(0 if allok else 1)
