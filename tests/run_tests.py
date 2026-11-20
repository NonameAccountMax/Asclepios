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
    return ok

@case
def contact_impact():
    # plate thrown at a fixed plate (general contact, penalty 0.3 m/dt^2, normal damping 3 % of critical): rebound
    # without energy creation; the damped penalty contact (resolved in 2-3 increments) keeps 70-90 % of the kinetic
    # energy (the shell-to-shell and beam-to-shell impacts of the validation suite agree with Abaqus with this damping)
    contact(os.path.join(OUT, 'impact.inp')); run('impact')
    h = history('impact'); ke0 = 0.5 * 7.8e-9 * 0.1 * 0.8 * 0.8 * 1000.0 ** 2
    ok = check('impact: energy balance |ETOTAL - KE0| / KE0', 1 + abs(h['ETOTAL'] - ke0) / ke0, 1.0, 5e-3)
    r = (h['ALLKE'] + h['ALLIE'] + h['ALLVD']) / ke0
    ok &= check('impact: mechanical energy kept (0.70-1.0)', min(max(r, 0.70), 1.0), r, 1e-12)
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

if __name__ == '__main__':
    allok = True
    for c in cases:
        print(c.__name__)
        try: allok &= c()
        except Exception as e: print('  ERROR', e); allok = False
    print('ALL PASSED' if allok else 'SOME TESTS FAILED')
    sys.exit(0 if allok else 1)
