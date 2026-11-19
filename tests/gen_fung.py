"""Single S4R (or 2 x S3R) with Fung anisotropic material under homogeneous biaxial stretch (static, NLGEOM).
The orientation is rotated by `ang` degrees about the normal; the stretch is applied in the global axes."""
import math, sys
FUNG_MAT = ('*Anisotropic Hyperelastic, fung-anisotropic, local directions=2\n'
         '8.00000, 0.00000, 100.000, 0.00000, 0.00000, 8.00000, 0.00000, 0.00000\n'
         '0.00000, 12.0000, 0.00000, 0.00000, 0.00000, 0.00000, 8.00000, 0.00000\n'
         '0.00000, 0.00000, 0.00000, 0.00000, 12.0000, 2.00E-05, 2.00E-05\n')
def fung(path, etype='S4R', l1=1.15, l2=1.10, ang=0.0, t=0.5, static=True):
    c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    with open(path, 'w') as f:
        f.write('*Heading\nFung patch\n*Node\n1, 0., 0., 0.\n2, 1., 0., 0.\n3, 1., 1., 0.\n4, 0., 1., 0.\n')
        f.write(f'*Element, type={etype}, elset=P\n')
        f.write('1, 1, 2, 3, 4\n' if etype == 'S4R' else '1, 1, 2, 3\n2, 1, 3, 4\n')
        f.write(f'*Orientation, name=OR, local directions=2\n{c:.12f}, {s:.12f}, 0., {-s:.12f}, {c:.12f}, 0.\n')
        f.write(f'*Shell Section, elset=P, material=M, orientation=OR\n{t}, 5\n')
        f.write('*Material, name=M\n' + FUNG_MAT + '*Density\n1.06e-9\n')
        f.write('*Amplitude, name=R\n0., 0., 1., 1.\n')
        f.write('*Boundary\n1, 1, 6\n2, 2, 6\n3, 3, 6\n4, 1, 1\n4, 3, 6\n')
        f.write('*Step, nlgeom=YES\n*Static\n0.05, 1.\n*Boundary, amplitude=R\n')
        f.write(f'2, 1, 1, {l1-1}\n3, 1, 1, {l1-1}\n3, 2, 2, {l2-1}\n4, 2, 2, {l2-1}\n')
        f.write('*Output, field, number interval=1\n*Node Output\nU, RF\n*End Step\n')

def analytic(l1, l2, ang=0.0):
    """Plane-stress Fung (FUNG_MAT constants: generic, strongly anisotropic) under F = diag(l1, l2, l3) in global axes, material axes rotated by
    `ang` about z: returns (sigma_xx, sigma_yy, l3). Independent of the solver: energy written in closed form
    and differentiated numerically in 40-digit arithmetic (mpmath), S33 = 0 solved for l3."""
    from mpmath import mp, mpf, matrix, exp, log, sqrt, diff, findroot, cos, sin, pi
    mp.dps = 40
    bd = {(0, 0): '8', (1, 1): '100', (2, 2): '8', (3, 3): '12', (4, 4): '8', (5, 5): '12'}
    cc, D = mpf('2.0e-5'), mpf('2.0e-5')
    a = mpf(ang) * pi / 180
    R = matrix([[cos(a), -sin(a), 0], [sin(a), cos(a), 0], [0, 0, 1]])
    idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    w = [1, 1, 1, 2, 2, 2]
    def U(C):
        J = sqrt(mp.det(C)); Cl = R.T * C * R; f = J ** (mpf(-2) / 3)
        E = [(f * Cl[i, j] - (1 if i == j else 0)) / 2 for i, j in idx]
        Q = sum(mpf(v) * w[I] * w[K] * E[I] * E[K] for (I, K), v in bd.items())
        return cc / 2 * (exp(Q) - 1) + ((J * J - 1) / 2 - log(J)) / D
    def S(C, i, j):
        def f(x):
            Cp = C.copy(); Cp[i, j] = x
            if i != j: Cp[j, i] = x
            return U(Cp)
        d = diff(f, C[i, j])
        return 2 * d if i == j else d
    def Cof(l3): return matrix([[mpf(l1) ** 2, 0, 0], [0, mpf(l2) ** 2, 0], [0, 0, l3 ** 2]])
    l3 = findroot(lambda x: S(Cof(x), 2, 2), mpf(1) / (mpf(l1) * mpf(l2)))
    C = Cof(l3); J = mpf(l1) * mpf(l2) * l3
    return float(mpf(l1) ** 2 * S(C, 0, 0) / J), float(mpf(l2) ** 2 * S(C, 1, 1) / J), float(l3)
if __name__ == '__main__':
    fung(sys.argv[1])
    print(analytic(1.15, 1.10))
