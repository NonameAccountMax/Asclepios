"""Cantilever beam (B31), clamped at x=0, tip load in z; section CIRC or RECT."""
import sys
def beam(path, n=20, L=10.0, sec='CIRC', dims='0.1', P=0.001, T=2e-2, alpha=1000., nlgeom='NO', plastic=False, dT=None):
    with open(path, 'w') as f:
        f.write('*Node\n')
        for i in range(n + 1): f.write(f'{i+1}, {L*i/n}, 0., 0.\n')
        f.write('*Element, type=B31, elset=B\n')
        for i in range(n): f.write(f'{i+1}, {i+1}, {i+2}\n')
        f.write(f'*Nset, nset=ROOT\n1\n*Nset, nset=TIP\n{n+1}\n*Nset, nset=ALL, generate\n1, {n+1}, 1\n')
        f.write(f'*Beam Section, elset=B, material=M, section={sec}\n{dims}\n0., 1., 0.\n')
        f.write('*Material, name=M\n*Elastic\n200000., 0.3\n*Density\n7.8e-9\n' + f'*Damping, alpha={alpha}\n')
        if plastic: f.write('*Plastic\n250., 0.\n300., 0.2\n')
        if dT is not None: f.write('*Expansion\n1e-3\n')
        f.write('*Boundary\nROOT, ENCASTRE\n')
        f.write(f'*Amplitude, name=SM, definition=SMOOTH STEP\n0., 0., {T/3}, 1.\n')
        f.write(f'*Step, nlgeom={nlgeom}\n*Dynamic, Explicit\n, {T}\n')
        if P: f.write(f'*Cload, amplitude=SM\nTIP, 3, {P}\n')
        if dT is not None: f.write(f'*Temperature, amplitude=SM\nALL, {dT}\n')
        f.write('*Output, field, number interval=5\n*Node Output\nU\n*End Step\n')
if __name__ == '__main__':
    beam(sys.argv[1])
