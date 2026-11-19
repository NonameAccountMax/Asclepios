"""Cantilever plate (S4R or S3R) in the xy plane, clamped at x=0, tip load in z (explicit, damped)."""
import sys
def plate(path, etype='S4R', nx=20, ny=2, L=10.0, W=1.0, t=0.1, P=0.01, T=6e-3, alpha=6000., nlgeom='NO', load_dir=3, mat=None):
    def nid(i, j): return 1 + i + (nx + 1) * j
    with open(path, 'w') as f:
        f.write('*Heading\ncantilever plate\n*Node\n')
        for j in range(ny + 1):
            for i in range(nx + 1):
                f.write(f'{nid(i,j)}, {L*i/nx:.6f}, {W*j/ny:.6f}, 0.\n')
        f.write(f'*Element, type={etype}, elset=PLATE\n')
        e = 1
        for j in range(ny):
            for i in range(nx):
                a, b, c, d = nid(i,j), nid(i+1,j), nid(i+1,j+1), nid(i,j+1)
                if etype == 'S4R':
                    f.write(f'{e}, {a}, {b}, {c}, {d}\n'); e += 1
                else:
                    f.write(f'{e}, {a}, {b}, {c}\n'); e += 1
                    f.write(f'{e}, {a}, {c}, {d}\n'); e += 1
        f.write('*Nset, nset=ROOT\n' + ', '.join(str(nid(0,j)) for j in range(ny+1)) + '\n')
        f.write('*Nset, nset=TIP\n' + ', '.join(str(nid(nx,j)) for j in range(ny+1)) + '\n')
        f.write(f'*Shell Section, elset=PLATE, material=MAT\n{t}, 5\n')
        f.write('*Material, name=MAT\n' + (mat or '*Elastic\n200000., 0.3\n') + f'*Density\n7.8e-9\n*Damping, alpha={alpha}\n')
        f.write('*Boundary\nROOT, ENCASTRE\n')
        f.write(f'*Amplitude, name=SM, definition=SMOOTH STEP\n0., 0., {T/3}, 1.\n')
        f.write(f'*Step, nlgeom={nlgeom}\n*Dynamic, Explicit\n, {T}\n*Cload, amplitude=SM\nTIP, {load_dir}, {P/(ny+1)}\n')
        f.write('*Output, field, number interval=10\n*Node Output\nU, RF\n*End Step\n')
if __name__ == '__main__':
    plate(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'S4R')
