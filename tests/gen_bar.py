"""Hex bar along x, fixed at x=0, displacement ramp at x=L. Writes an Abaqus-style deck."""
import sys
def bar(path, etype='C3D8R', nx=20, ny=2, nz=2, L=10.0, W=1.0, H=1.0, disp=0.01, T=1e-4, plastic=False, nlgeom='YES', static=False):
    nodes = {}
    def nid(i, j, k): return 1 + i + (nx + 1) * (j + (ny + 1) * k)
    with open(path, 'w') as f:
        f.write('*Heading\nbar test\n*Node\n')
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    f.write(f'{nid(i,j,k)}, {L*i/nx:.6f}, {W*j/ny:.6f}, {H*k/nz:.6f}\n')
        f.write(f'*Element, type={etype}, elset=BAR\n')
        e = 1
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    n = [nid(i,j,k), nid(i+1,j,k), nid(i+1,j+1,k), nid(i,j+1,k), nid(i,j,k+1), nid(i+1,j,k+1), nid(i+1,j+1,k+1), nid(i,j+1,k+1)]
                    f.write(f'{e}, ' + ', '.join(map(str, n)) + '\n'); e += 1
        f.write('*Nset, nset=X0\n' + ', '.join(str(nid(0,j,k)) for k in range(nz+1) for j in range(ny+1)) + '\n')
        f.write('*Nset, nset=XL\n' + ', '.join(str(nid(nx,j,k)) for k in range(nz+1) for j in range(ny+1)) + '\n')
        f.write(f'*Nset, nset=ORIG\n{nid(0,0,0)}\n')
        f.write('*Solid Section, elset=BAR, material=STEEL\n,\n*Material, name=STEEL\n*Elastic\n200000., 0.3\n*Density\n7.8e-9\n')
        if plastic: f.write('*Plastic\n250., 0.\n350., 0.1\n')
        f.write('*Amplitude, name=RAMP, definition=SMOOTH STEP\n0., 0., %g, 1.\n' % T)
        f.write('*Boundary\nX0, 1, 1\nORIG, 2, 3\n')
        f.write(f'*Nset, nset=Y0\n' + ', '.join(str(nid(0,0,k)) for k in range(nz+1)) + '\n')
        f.write('*Boundary\nY0, 2, 2\n')
        f.write(f'*Nset, nset=Z0\n' + ', '.join(str(nid(0,j,0)) for j in range(ny+1)) + '\n')
        f.write('*Boundary\nZ0, 3, 3\n')
        if static:
            f.write(f'*Step, nlgeom={nlgeom}\n*Static\n0.1, 1.0, 1e-5, 0.25\n*Boundary\nXL, 1, 1, {disp}\n')
        else:
            f.write(f'*Step, nlgeom={nlgeom}\n*Dynamic, Explicit\n, {T}\n*Boundary, amplitude=RAMP\nXL, 1, 1, {disp}\n')
        f.write('*Output, field, number interval=10\n*Node Output\nU, RF\n*Element Output\nS\n*End Step\n')
if __name__ == '__main__':
    bar(sys.argv[1], *sys.argv[2:3])
