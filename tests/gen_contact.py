"""Plate (S4R, 2x2) thrown at a fixed plate (S4R, 4x4) with general contact: elastic rebound, energy balance."""
import sys
def contact(path, v=-1000.0, gap=0.3, T=6e-4, mu=0.0):
    with open(path, 'w') as f:
        f.write('*Heading\nplate impact\n*Node\n')
        nid = 1; low = {}; up = {}
        for j in range(5):
            for i in range(5):
                low[i, j] = nid; f.write(f'{nid}, {-0.5 + 0.5 * i}, {-0.5 + 0.5 * j}, 0.\n'); nid += 1
        for j in range(3):
            for i in range(3):
                up[i, j] = nid; f.write(f'{nid}, {0.1 + 0.4 * i}, {0.1 + 0.4 * j}, {gap}\n'); nid += 1
        f.write('*Element, type=S4R, elset=LOW\n')
        e = 1
        for j in range(4):
            for i in range(4):
                f.write(f'{e}, {low[i,j]}, {low[i+1,j]}, {low[i+1,j+1]}, {low[i,j+1]}\n'); e += 1
        f.write('*Element, type=S4R, elset=UP\n')
        for j in range(2):
            for i in range(2):
                f.write(f'{e}, {up[i,j]}, {up[i+1,j]}, {up[i+1,j+1]}, {up[i,j+1]}\n'); e += 1
        f.write('*Nset, nset=NLOW\n' + ', '.join(str(n) for n in low.values()) + '\n')
        f.write('*Nset, nset=NUP\n' + ', '.join(str(n) for n in up.values()) + '\n')
        f.write('*Shell Section, elset=LOW, material=M\n0.1, 5\n*Shell Section, elset=UP, material=M\n0.1, 5\n')
        f.write('*Material, name=M\n*Elastic\n200000., 0.3\n*Density\n7.8e-9\n')
        f.write('*Surface Interaction, name=IP\n*Friction\n%g\n' % mu)
        f.write(f'*Initial Conditions, type=velocity\nNUP, 3, {v}\n')
        f.write('*Boundary\nNLOW, ENCASTRE\n')
        f.write(f'*Step, nlgeom=YES\n*Dynamic, Explicit\n, {T}\n*Contact\n*Contact Inclusions, ALL EXTERIOR\n*Contact property assignment\n, , IP\n')
        f.write('*Output, field, number interval=10\n*Node Output\nU, V\n*Output, history\n*End Step\n')
if __name__ == '__main__':
    contact(sys.argv[1])
