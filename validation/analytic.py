#!/usr/bin/env python3
"""Closed-form reference solutions of the validation cases whose nodal solution is known by hand.

Used for the cases flagged reference = 'analytic' in cases/manifest.json; the other cases are compared with Abaqus.
Each function returns the same layout as the Abaqus extraction (frames with t, U, UR, RF; energy histories), built
from the node coordinates of the deck. Two kinds of cases:
  - exact: homogeneous states, springs, rigid motion, constant moment: the closed form is the exact solution of
    the discrete model;
  - theory: slender cantilevers under tip load (Timoshenko beam), with a discretisation error below 1 %.

usage: analytic.py [case ...]     prints the reference values at the last frame
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
E, NU, ALPHA = 210000.0, 0.3, 1.2e-05
G = E / (2 * (1 + NU))


def deck_nodes(name, path=None):
    nodes, on = {}, False
    for line in open(path or os.path.join(HERE, 'cases', name + '.inp')):
        if line.startswith('*'):
            low = line.lower()
            if low.startswith('*include'):
                inc = line.split('=', 1)[1].strip()
                nodes.update(deck_nodes(None, os.path.join(HERE, 'cases', inc)))
            on = low.startswith('*node') and not low.startswith('*node output')
            continue
        if on and line.strip():
            p = [float(v) for v in line.split(',')]
            nodes[int(p[0])] = tuple(p[1:4])
    return nodes


def ramp(n, t_end=1.0, t0=0.0):
    return [t0 + (t_end - t0) * k / n for k in range(1, n + 1)]


def frames_of(case, nodes, times, disp, rot=None, rf=None, elem=None):
    """Reference frames; elem(t) gives the element quantities of every monitored element (homogeneous states)."""
    fr = []
    for t in times:
        f = {'t': t, 'U': {str(n): list(disp(t, *nodes[n])) for n in case['mon']}}
        if rot:
            f['UR'] = {str(n): list(rot(t, *nodes[n])) for n in case['mon']}
        if rf:
            f['RF'] = list(rf(t))
        if elem and case.get('mone'):
            e = elem(t)
            f['E'] = {str(lab): dict(e) for lab in case['mone']}
        fr.append(f)
    return fr


def mises_plane(sx, sy, txy):
    return math.sqrt(sx * sx - sx * sy + sy * sy + 3 * txy * txy)


def mises3(s):
    a, b, c, d, e, f = s
    return math.sqrt(0.5 * ((a - b) ** 2 + (b - c) ** 2 + (c - a) ** 2) + 3 * (d * d + e * e + f * f))


def zero_rot(t, x, y, z):
    return (0.0, 0.0, 0.0)


# ------------------------------------------------------------------------------------------------ solids
def solid_tension(case, nodes):
    eps, A, V = 0.04 / 40.0, 16.0, 640.0
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (eps * x, -NU * eps * y, -NU * eps * z),
                   rf=lambda t: (-E * A * eps, 0.0, 0.0), elem=lambda t: {'mises': E * eps, 'smem': E * eps})
    w = 0.5 * E * eps ** 2 * V
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def solid_thermal(case, nodes):
    dT, A, V = 100.0, 16.0, 640.0
    sig = -E * ALPHA * dT
    elat = ALPHA * dT * (1 + NU)
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, elat * y, elat * z), rf=lambda t: (-sig * A, 0.0, 0.0),
                   elem=lambda t: {'mises': abs(sig), 'smem': abs(sig)})
    return fr, {'ALLIE': [[1.0, 0.5 * sig ** 2 / E * V]]}


def solid_shear(case, nodes):
    c10, V = 0.5, 216.0  # simple shear is isochoric: W = C10 (I1 - 3) = C10 gamma^2
    times = ramp(10)
    def elem(t):
        g = 0.5 * t
        m = 2 * c10 * math.sqrt(g ** 4 + 3 * g * g)
        return {'mises': m, 'smem': m}
    fr = frames_of(case, nodes, times, lambda t, x, y, z: (0.5 * t * z, 0.0, 0.0), elem=elem)
    hist = [[t, c10 * (0.5 * t) ** 2 * V] for t in [0.0] + ramp(200)]
    return fr, {'ALLIE': hist, 'ALLWK': hist}


HARD = [(250.0, 0.0), (300.0, 0.05), (340.0, 0.2), (360.0, 1.0)]


def flow_stress(ep):
    if ep >= HARD[-1][1]:
        return HARD[-1][0]
    for (s0, e0), (s1, e1) in zip(HARD, HARD[1:]):
        if ep <= e1:
            return s0 + (s1 - s0) * (ep - e0) / (e1 - e0)


def uniaxial_plastic(path, times, substeps=500):
    """1D J2 with isotropic hardening in logarithmic strain (uniaxial stress state of the homogeneous bar).
    path(t) = tip displacement; returns {t: (stress, plastic strain, log strain, work)} at the requested times."""
    L, A0 = 40.0, 16.0
    V0 = L * A0
    sig, ep, ebar, work = 0.0, 0.0, 0.0, 0.0
    t_prev, eps_prev = 0.0, 0.0
    out = {}
    for t_target in times:
        for k in range(1, substeps + 1):
            t = t_prev + (t_target - t_prev) * k / substeps
            eps = math.log(1.0 + path(t) / L)
            de = eps - eps_prev
            s_tr = sig + E * de
            f = abs(s_tr) - flow_stress(ebar)
            s_old = sig
            if f > 0:
                lo, hi = 0.0, abs(s_tr) / E
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    if abs(s_tr) - E * mid - flow_stress(ebar + mid) > 0:
                        lo = mid
                    else:
                        hi = mid
                dg = 0.5 * (lo + hi)
                sgn = 1.0 if s_tr > 0 else -1.0
                sig = s_tr - sgn * E * dg
                ep += sgn * dg
                ebar += dg
            else:
                sig = s_tr
            vol = V0 * math.exp((1 - 2 * NU) * 0.5 * (sig + s_old) / E)
            work += 0.5 * (sig + s_old) * de * vol
            eps_prev = eps
        t_prev = t_target
        out[t_target] = (sig, ep, eps_prev, work)
    return out


def solid_plastic(case, nodes):
    def path(t):
        return 2.0 * t if t <= 1.0 else 2.0 - (t - 1.0)
    times = ramp(10) + ramp(10, 2.0, 1.0)
    st = uniaxial_plastic(path, ramp(40) + ramp(40, 2.0, 1.0))

    def elat(t):
        s, ep, eps, w = st[min(st, key=lambda k: abs(k - t))]
        return -NU * s / E - 0.5 * ep

    def disp(t, x, y, z):
        f = math.exp(elat(t)) - 1.0
        return (path(t) * x / 40.0, f * y, f * z)

    def rf(t):
        s = st[min(st, key=lambda k: abs(k - t))][0]
        return (-s * 16.0 * math.exp(2 * elat(t)), 0.0, 0.0)
    def elem(t):
        s_, ep_ = st[min(st, key=lambda k: abs(k - t))][:2]
        return {'mises': abs(s_), 'smem': abs(s_), 'peeq': plastic_bar_ebar(t)}
    ebar_hist = {}

    def plastic_bar_ebar(t):  # equivalent plastic strain along the path (1D model of uniaxial_plastic)
        if not ebar_hist:
            sig, ep, ebar, eps_prev = 0.0, 0.0, 0.0, 0.0
            for k in range(1, 8001):
                tt = 2.0 * k / 8000
                eps = math.log(1.0 + path(tt) / 40.0)
                s_tr = sig + E * (eps - eps_prev)
                if abs(s_tr) > flow_stress(ebar):
                    lo, hi = 0.0, abs(s_tr) / E
                    for _ in range(60):
                        mid = 0.5 * (lo + hi)
                        lo, hi = (mid, hi) if abs(s_tr) - E * mid - flow_stress(ebar + mid) > 0 else (lo, mid)
                    dg = 0.5 * (lo + hi)
                    sig = s_tr - (1 if s_tr > 0 else -1) * E * dg
                    ebar += dg
                else:
                    sig = s_tr
                eps_prev = eps
                ebar_hist[round(tt, 6)] = ebar
        return ebar_hist[min(ebar_hist, key=lambda k: abs(k - t))]
    fr = frames_of(case, nodes, times, disp, rf=rf, elem=elem)
    hist = [[t, v[3]] for t, v in sorted(st.items())]
    return fr, {'ALLIE': hist, 'ALLWK': hist}


# ------------------------------------------------------------------------------------------------ shells
def shell_membrane(case, nodes):
    eps, t0, W, V = 0.001, 0.5, 10.0, 200.0
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (eps * x, -NU * eps * y, 0.0), zero_rot,
                   lambda t: (-E * t0 * W * eps, 0.0, 0.0), elem=lambda t: {'mises': E * eps, 'smem': E * eps})
    w = 0.5 * E * eps ** 2 * V
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def cantilever(P, L, EI, kGA):
    """Timoshenko cantilever under a tip load P: deflection and bending slope."""
    def w(x):
        return P / EI * (L * x * x / 2 - x ** 3 / 6) + P * x / kGA

    def slope(x):
        return P / EI * (L * x - x * x / 2)
    return w, slope


def shell_bending(case, nodes):
    P, L, W, t0 = -1.0, 40.0, 10.0, 0.5
    w, slope = cantilever(P, L, E * W * t0 ** 3 / 12, 5.0 / 6.0 * (E / 2) * W * t0)  # Poisson ratio 0
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, w(x)), lambda t, x, y, z: (0.0, -slope(x), 0.0),
                   lambda t: (0.0, 0.0, -P))
    e = 0.5 * P * w(L)
    return fr, {'ALLIE': [[1.0, e]], 'ALLWK': [[1.0, e]]}


def rollup(case, nodes, L, EI):
    """End moment M(t) = t pi EI / (2 L): constant curvature, nodes on a circular arc curling towards +z."""
    k1 = math.pi / (2 * L)

    def disp(t, x, y, z):
        k = t * k1
        if k == 0:
            return (0.0, 0.0, 0.0)
        th = k * x
        return (math.sin(th) / k - x, 0.0, (1 - math.cos(th)) / k)

    def rot(t, x, y, z):
        return (0.0, -t * k1 * x, 0.0)
    fr = frames_of(case, nodes, ramp(10), disp, rot, lambda t: (0.0, 0.0, 0.0))
    m1 = k1 * EI
    hist = [[t, 0.5 * (t * m1) * (t * k1 * L)] for t in [0.0] + ramp(200)]
    return fr, {'ALLIE': hist, 'ALLWK': hist}


def shell_twist(case, nodes):
    P, a, t0 = -2.0, 20.0, 0.5
    D = E * t0 ** 3 / (12 * (1 - NU ** 2))
    c = P / (2 * D * (1 - NU))  # w = c x y, uniform twisting moment
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, c * x * y),
                   lambda t, x, y, z: (c * x, -c * y, 0.0), lambda t: (0.0, 0.0, -P))
    e = 0.5 * P * c * a * a
    return fr, {'ALLIE': [[1.0, e]], 'ALLWK': [[1.0, e]]}


FUNG_B = {(0, 0): '8', (1, 1): '100', (2, 2): '8', (3, 3): '12', (4, 4): '8', (5, 5): '12'}  # material FUNG
FUNG_C, FUNG_D = '2e-05', '2e-05'


def fung_state(l1, l2, ang):
    """Plane-stress Fung membrane under F = diag(l1, l2, l3), material axes rotated by ang (degrees) about z:
    Cauchy stresses (xx, yy, xy), l3 and strain energy density. Energy in closed form, derivatives in 40-digit
    arithmetic (mpmath), S33 = 0 solved for l3."""
    from mpmath import mp, mpf, matrix, exp, log, sqrt, diff, findroot, cos, sin, pi
    mp.dps = 40
    a = mpf(ang) * pi / 180
    R = matrix([[cos(a), -sin(a), 0], [sin(a), cos(a), 0], [0, 0, 1]])
    idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    wgt = [1, 1, 1, 2, 2, 2]
    cc, D = mpf(FUNG_C), mpf(FUNG_D)

    def U(C):
        J = sqrt(mp.det(C))
        Cl = R.T * C * R
        f = J ** (mpf(-2) / 3)
        Eg = [(f * Cl[i, j] - (1 if i == j else 0)) / 2 for i, j in idx]
        Q = sum(mpf(v) * wgt[I] * wgt[K] * Eg[I] * Eg[K] for (I, K), v in FUNG_B.items())
        return cc / 2 * (exp(Q) - 1) + ((J * J - 1) / 2 - log(J)) / D

    def S(C, i, j):
        def f(x):
            Cp = C.copy()
            Cp[i, j] = x
            if i != j:
                Cp[j, i] = x
            return U(Cp)
        d = diff(f, C[i, j])
        return 2 * d if i == j else d

    def Cof(l3):
        return matrix([[mpf(l1) ** 2, 0, 0], [0, mpf(l2) ** 2, 0], [0, 0, l3 ** 2]])
    l3 = findroot(lambda x: S(Cof(x), 2, 2), mpf(1) / (mpf(l1) * mpf(l2)))
    C = Cof(l3)
    J = mpf(l1) * mpf(l2) * l3
    sxx = mpf(l1) ** 2 * S(C, 0, 0) / J
    syy = mpf(l2) ** 2 * S(C, 1, 1) / J
    sxy = mpf(l1) * mpf(l2) * S(C, 0, 1) / J
    return float(sxx), float(syy), float(sxy), float(l3), float(U(C))


def neohooke_state(l1, l2, c10=0.5, d1=0.04):
    """Plane-stress compressible neo-Hooke (W = C10 (I1bar - 3) + (J - 1)^2 / D1): sxx, syy, l3, W."""
    def sig(l3, li):
        J = l1 * l2 * l3
        f = J ** (-2.0 / 3.0)
        i1 = f * (l1 * l1 + l2 * l2 + l3 * l3)
        return 2 * c10 / J * (f * li * li - i1 / 3) + 2 / d1 * (J - 1)
    lo, hi = 0.05, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sig(mid, mid) > 0:
            hi = mid
        else:
            lo = mid
    l3 = 0.5 * (lo + hi)
    J = l1 * l2 * l3
    i1 = J ** (-2.0 / 3.0) * (l1 * l1 + l2 * l2 + l3 * l3)
    return sig(l3, l1), sig(l3, l2), l3, c10 * (i1 - 3) + (J - 1) ** 2 / d1


def membrane_patch(case, nodes, e1, e2, state):
    t0, W, V0 = 0.5, 2.0, 2.0
    times = ramp(10)
    st = {t: state(1 + e1 * t, 1 + e2 * t) for t in sorted(set(times + ramp(100)))}

    def rf(t):
        sxx, sxy, l3 = st[t][0], st[t][2], st[t][3]
        area = W * (1 + e2 * t) * t0 * l3
        return (-sxx * area, -sxy * area, 0.0)
    def elem(t):
        m = mises_plane(st[t][0], st[t][1], st[t][2])
        return {'mises': m, 'smem': m}
    fr = frames_of(case, nodes, times, lambda t, x, y, z: (e1 * t * x, e2 * t * y, 0.0), zero_rot, rf, elem)
    hist = [[0.0, 0.0]] + [[t, st[t][4] * V0] for t in ramp(100)]
    return fr, {'ALLIE': hist, 'ALLWK': hist}


def shell_fung(case, nodes):
    def state(l1, l2):
        sxx, syy, sxy, l3, w = fung_state(l1, l2, 30.0)
        return (sxx, syy, sxy, l3, w)
    return membrane_patch(case, nodes, 0.15, 0.10, state)


def shell_neohooke(case, nodes):
    def state(l1, l2):
        sxx, syy, l3, w = neohooke_state(l1, l2)
        return (sxx, syy, 0.0, l3, w)
    return membrane_patch(case, nodes, 0.5, 0.3, state)


def shell_spin(case, nodes):
    w = 500.0

    def disp(t, x, y, z):
        c, s = math.cos(w * t), math.sin(w * t)
        return (c * x - s * y - x, s * x + c * y - y, 0.0)
    fr = frames_of(case, nodes, ramp(20, 2e-3), disp, lambda t, x, y, z: (0.0, 0.0, w * t))
    return fr, {}


# ------------------------------------------------------------------------------------------------- beams
def b31_tension(case, nodes):
    eps, A, L = 0.05 / 100.0, 8.0, 100.0
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (eps * x, 0.0, 0.0), zero_rot,
                   lambda t: (-E * A * eps, 0.0, 0.0))
    w = 0.5 * E * A * eps ** 2 * L
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def b31_bending(case, nodes):
    # RECT 2 x 1 with n1 = -z: 2 mm along z, 1 mm along y
    L, Fy, Fz, A = 100.0, 0.5, -2.0, 2.0
    kga = 5.0 / 6.0 * G * A
    wy, sy = cantilever(Fy, L, E * 2.0 * 1.0 ** 3 / 12, kga)
    wz, sz = cantilever(Fz, L, E * 1.0 * 2.0 ** 3 / 12, kga)
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, wy(x), wz(x)),
                   lambda t, x, y, z: (0.0, -sz(x), sy(x)), lambda t: (0.0, -Fy, -Fz))
    e = 0.5 * (Fy * wy(L) + Fz * wz(L))
    return fr, {'ALLIE': [[1.0, e]], 'ALLWK': [[1.0, e]]}


def b31_torsion(case, nodes):
    M, r, L = 1000.0, 2.0, 100.0
    GJ = G * math.pi * r ** 4 / 2
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, 0.0), lambda t, x, y, z: (M * x / GJ, 0.0, 0.0),
                   lambda t: (0.0, 0.0, 0.0))
    e = 0.5 * M * M * L / GJ
    return fr, {'ALLIE': [[1.0, e]], 'ALLWK': [[1.0, e]]}


def b31_thermal(case, nodes):
    A, L, dT = math.pi * (3.0 ** 2 - 2.5 ** 2), 100.0, 100.0
    N = E * A * ALPHA * dT
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, 0.0), zero_rot, lambda t: (N, 0.0, 0.0))
    return fr, {'ALLIE': [[1.0, 0.5 * N * N * L / (E * A)]]}


# -------------------------------------------------------------------------------------------- connectors
def conn_axial(case, nodes):
    k, L0, F = 100.0, 10.0, 50.0

    def state(t):
        if t <= 1.0:
            return 0.5 * t, 0.0
        u2 = 5.0 * (t - 1.0)
        lo, hi = -L0 / 2, 1.0  # axial equilibrium k (l - L0) (L0 + u1) / l = F
        for _ in range(200):
            u1 = 0.5 * (lo + hi)
            l = math.hypot(L0 + u1, u2)
            if k * (l - L0) * (L0 + u1) / l > F:
                hi = u1
            else:
                lo = u1
        return 0.5 * (lo + hi), u2

    def rf(t):
        u1, u2 = state(t)
        l = math.hypot(L0 + u1, u2)
        f = k * (l - L0)
        return (-f * (L0 + u1) / l, -f * u2 / l, 0.0)

    def energy(t):
        u1, u2 = state(t)
        return 0.5 * k * (math.hypot(L0 + u1, u2) - L0) ** 2
    times = ramp(10) + ramp(10, 2.0, 1.0)
    fr = frames_of(case, nodes, times, lambda t, x, y, z: (state(t)[0], state(t)[1], 0.0), rf=rf)
    return fr, {'ALLIE': [[t, energy(t)] for t in ramp(40) + ramp(40, 2.0, 1.0)]}


def conn_cartesian(case, nodes):
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.1, 0.1, 0.1), rf=lambda t: (-10.0, -20.0, -30.0))
    return fr, {'ALLIE': [[1.0, 3.0]], 'ALLWK': [[1.0, 3.0]]}


def conn_oscillator(case, nodes):
    m, k, v0 = 7.85e-06 * 8.0, 100.0, 100.0
    w = math.sqrt(k / m)

    def u(t):
        return v0 / w * math.sin(w * t)
    fr = frames_of(case, nodes, ramp(20, 0.01), lambda t, x, y, z: (u(t), 0.0, 0.0), rf=lambda t: (-k * u(t), 0.0, 0.0))
    ts = [0.0] + ramp(4000, 0.01)
    return fr, {'ALLKE': [[t, 0.5 * m * (v0 * math.cos(w * t)) ** 2] for t in ts],
                'ALLIE': [[t, 0.5 * k * u(t) ** 2] for t in ts]}


# ------------------------------------------------------------------------------- beams with constraints
def beam_props(height_z=2.0, width_y=1.0):
    A = height_z * width_y
    return E * width_y * height_z ** 3 / 12, 5.0 / 6.0 * G * A, E * A


def integrate_beam(L, segments, P, n=200000):
    """Cantilever with a tip load P at x = L; segments = list of (x0, x1, flexible). Returns w(x), theta_y(x)
    by integrating curvature P (L - x) / EI and shear strain P / kGA on the flexible parts."""
    EI, kGA, _ = beam_props()
    xs = [L * k / n for k in range(n + 1)]
    phi, w, ws, ph = 0.0, 0.0, [0.0], [0.0]
    for k in range(n):
        x0, x1 = xs[k], xs[k + 1]
        xm = 0.5 * (x0 + x1)
        flex = any(a <= xm <= b and f for a, b, f in segments)
        dx = x1 - x0
        kap = P * (L - xm) / EI if flex else 0.0
        phi_new = phi + kap * dx
        w += 0.5 * (phi + phi_new) * dx + (P / kGA * dx if flex else 0.0)
        phi = phi_new
        ws.append(w)
        ph.append(phi)

    def at(arr, x):
        k = min(int(x / L * n), n - 1)
        a = (x - xs[k]) / (xs[k + 1] - xs[k])
        return arr[k] + a * (arr[k + 1] - arr[k])
    return (lambda x: at(ws, x)), (lambda x: -at(ph, x))


def mpc_beam(case, nodes):
    w, th = integrate_beam(41.0, [(0.0, 20.0, True), (20.0, 21.0, False), (21.0, 41.0, True)], -1.0)
    return frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, w(x)), lambda t, x, y, z: (0.0, th(x), 0.0),
                     lambda t: (0.0, 0.0, 1.0)), {}


def equation(case, nodes):
    w, th = integrate_beam(20.0, [(0.0, 20.0, True)], -0.5)
    return frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, w(x)), lambda t, x, y, z: (0.0, th(x), 0.0),
                     lambda t: (0.0, 0.0, 1.0)), {}


def connector_join(case, nodes):
    w1, th1 = integrate_beam(20.0, [(0.0, 20.0, True)], -1.0)
    wp = w1(20.0)
    first = case['mon'][0]

    def disp(t, x, y, z):
        return (0.0, 0.0, w1(x) if x <= 20.0 + 1e-9 else wp * (40.0 - x) / 20.0)

    def rot(t, x, y, z):
        return (0.0, th1(x), 0.0) if x < 20.0 - 1e-9 else (0.0, wp / 20.0, 0.0)
    fr = frames_of(case, nodes, [1.0], disp, rot, lambda t: (0.0, 0.0, 1.0))
    fr[0]['UR'][str(first)] = [0.0, th1(20.0), 0.0]  # pin node of beam 1: rotation of the cantilever
    return fr, {}


def load_gravity(case, nodes):
    EI, kGA, _ = beam_props()
    L, q = 100.0, -7.85e-09 * 2.0 * 98100.0

    def w(x):
        return q * x * x * (6 * L * L - 4 * L * x + x * x) / (24 * EI) + q * (L * x - x * x / 2) / kGA
    return frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, w(x)),
                     lambda t, x, y, z: (0.0, -q * (3 * L * L * x - 3 * L * x * x + x ** 3) / (6 * EI), 0.0),
                     lambda t: (0.0, 0.0, -q * L)), {}


def bc_transform(case, nodes):
    # tip node in the axes e1 = (c, 0, s), e2 = y, e3 = (-s, 0, c): its load (-2 along local 3) and its condition
    # (no motion along e1) are both given in these axes; reaction R along e1
    EI, kGA, EA = beam_props()
    L, c, s_ = 100.0, math.cos(math.radians(30)), math.sin(math.radians(30))
    fw, ka = L ** 3 / (3 * EI) + L / kGA, EA / L
    Px, Pz = 2.0 * s_, -2.0 * c
    R = -(c * Px / ka + s_ * Pz * fw) / (c * c / ka + s_ * s_ * fw)
    V, N = Pz + R * s_, Px + R * c
    return frames_of(case, nodes, [1.0],
                     lambda t, x, y, z: (N * x / EA, 0.0, V / EI * (L * x * x / 2 - x ** 3 / 6) + V * x / kGA),
                     lambda t, x, y, z: (0.0, -V / EI * (L * x - x * x / 2), 0.0), lambda t: (-Px, 0.0, -Pz)), {}


def rigid_body(case, nodes):
    k, F, xl = 100.0, -4.0, 1.0
    a = F / (4 * k)
    b = F * xl / (16 * k)
    th = -b
    rp = max(case['mon'])  # the reference node is created after the block nodes
    X0 = nodes[rp]

    def disp(t, x, y, z):
        return (th * (z - X0[2]), 0.0, a + b * x)
    fr = frames_of(case, nodes, [1.0], disp, rf=lambda t: (0.0, 0.0, -F))
    fr[0]['UR'] = {str(rp): [0.0, th, 0.0]}
    return fr, {}


# ---------------------------------------------------------------------------------- homogeneous solids
def uniaxial_hyper(W, stretch_end=1.5, L=4.0, A0=16.0):
    """Homogeneous uniaxial stretch with free lateral faces for a strain energy W(l1, l2, l3): Cauchy stress
    s1 = l1 / J dW/dl1, lateral stretch from s2 = 0."""
    def sig(l1, lt, i):
        h = 1e-6
        l = [l1, lt, lt]
        J = l1 * lt * lt
        lp, lm = list(l), list(l)
        lp[i] += h
        lm[i] -= h
        return l[i] / J * (W(*lp) - W(*lm)) / (2 * h)

    def state(t):
        l1 = 1.0 + (stretch_end - 1.0) * t
        lo, hi = 0.3, 1.5
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if sig(l1, mid, 1) > 0:
                hi = mid
            else:
                lo = mid
        lt = 0.5 * (lo + hi)
        return l1, lt, sig(l1, lt, 0), W(l1, lt, lt)
    return state


def iso(l1, l2, l3):
    J = l1 * l2 * l3
    f = J ** (-1.0 / 3.0)
    return J, f * l1, f * l2, f * l3


def W_mooney(l1, l2, l3, c10=0.3, c01=0.2, d1=0.04):
    J, a, b, c = iso(l1, l2, l3)
    return c10 * (a * a + b * b + c * c - 3) + c01 * (a * a * b * b + b * b * c * c + c * c * a * a - 3) + (J - 1) ** 2 / d1


def W_ogden(l1, l2, l3):
    J, a, b, c = iso(l1, l2, l3)
    w = 0.0
    for mu, al in ((0.6, 2.5), (-0.1, -2.0)):
        w += 2 * mu / al ** 2 * (a ** al + b ** al + c ** al - 3)
    return w + (J - 1) ** 2 / 0.05


def W_yeoh(l1, l2, l3):
    J, a, b, c = iso(l1, l2, l3)
    i1 = a * a + b * b + c * c - 3
    return 0.5 * i1 - 0.02 * i1 ** 2 + 0.003 * i1 ** 3 + (J - 1) ** 2 / 0.04


def hyper_block(W):
    def solve(case, nodes):
        st = uniaxial_hyper(W)
        cache = {t: st(t) for t in sorted(set(ramp(10) + ramp(100)))}

        def disp(t, x, y, z):
            l1, lt = cache[t][:2]
            return ((l1 - 1) * x, (lt - 1) * y, (lt - 1) * z)
        fr = frames_of(case, nodes, ramp(10), disp, rf=lambda t: (-cache[t][2] * 16.0 * cache[t][1] ** 2, 0.0, 0.0),
                       elem=lambda t: {'mises': abs(cache[t][2]), 'smem': abs(cache[t][2])})
        hist = [[0.0, 0.0]] + [[t, cache[t][3] * 64.0] for t in ramp(100)]
        return fr, {'ALLIE': hist, 'ALLWK': hist}
    return solve


def kinematic_hardening(case, nodes):
    s0, C, L, A0 = 250.0, 1000.0, 4.0, 16.0

    def path(t):
        return 0.08 * t if t <= 1.0 else 0.08 - 0.16 * (t - 1.0)
    times = ramp(10) + ramp(10, 2.0, 1.0)
    sig, a, ep, eps_prev, work = 0.0, 0.0, 0.0, 0.0, 0.0
    out, wk = {}, []
    n = 4000
    for k in range(1, 2 * n + 1):
        t = 2.0 * k / (2 * n)
        eps = math.log(1.0 + path(t) / L)
        de = eps - eps_prev
        s_tr = sig + E * de
        f = abs(s_tr - a) - s0
        s_old = sig
        if f > 0:
            sg = 1.0 if s_tr - a > 0 else -1.0
            dg = f / (E + C)
            sig = s_tr - sg * E * dg
            a += sg * C * dg
            ep += sg * dg
        else:
            sig = s_tr
        vol = L * A0 * math.exp((1 - 2 * NU) * 0.5 * (sig + s_old) / E)
        work += 0.5 * (sig + s_old) * de * vol
        eps_prev = eps
        out[round(t, 9)] = (sig, ep)
        if k % 20 == 0:
            wk.append([t, work])

    def get(t):
        return out[min(out, key=lambda q: abs(q - t))]

    def elat(t):
        s_, e_ = get(t)
        return -NU * s_ / E - 0.5 * e_

    def disp(t, x, y, z):
        f = math.exp(elat(t)) - 1.0
        return (path(t) * x / L, f * y, f * z)
    fr = frames_of(case, nodes, times, disp, rf=lambda t: (-get(t)[0] * A0 * math.exp(2 * elat(t)), 0.0, 0.0),
                   elem=lambda t: {'mises': abs(get(t)[0]), 'smem': abs(get(t)[0])})
    return fr, {'ALLIE': wk, 'ALLWK': wk}


def orthotropic(case, nodes):
    E1, n12, n13, eps = 150000.0, 0.3, 0.3, 0.004 / 4.0
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (eps * x, -n12 * eps * y, -n13 * eps * z),
                   rf=lambda t: (-E1 * eps * 16.0, 0.0, 0.0), elem=lambda t: {'mises': E1 * eps, 'smem': E1 * eps})
    w = 0.5 * E1 * eps ** 2 * 64.0
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def initial_stress(case, nodes):
    s0 = 100.0
    return frames_of(case, nodes, [1.0], lambda t, x, y, z: (-s0 / E * x, NU * s0 / E * y, NU * s0 / E * z)), {}


def temperature_dependence(case, nodes):
    def state(t):
        if t <= 1.0:
            return 100.0 * t, 0.0
        return 100.0, 200.0 * (t - 1.0)

    def disp(t, x, y, z):
        s_, T = state(t)
        Et = 210000.0 - 300.0 * T
        ex, el = s_ / Et + ALPHA * T, -NU * s_ / Et + ALPHA * T
        return (ex * x, el * y, el * z)
    times = [1.0] + ramp(4, 2.0, 1.0)
    return frames_of(case, nodes, times, disp, rf=lambda t: (-16.0 * state(t)[0], 0.0, 0.0),
                     elem=lambda t: {'mises': state(t)[0], 'smem': state(t)[0]}), {}


def visco_long_term(case, nodes):
    # instantaneous moduli E0 = 20, nu0 = 0.3; Prony sums g = 0.5, k = 0.2 -> long-term G = 0.5 G0, K = 0.8 K0
    G0, K0, eps = 20.0 / 2.6, 20.0 / 1.2, 0.01

    def moduli(t):
        G, K = (G0, K0) if t <= 1.0 else (0.5 * G0, 0.8 * K0)
        return 9 * K * G / (3 * K + G), (3 * K - 2 * G) / (2 * (3 * K + G))
    return frames_of(case, nodes, [1.0, 2.0],
                     lambda t, x, y, z: (eps * x, -moduli(t)[1] * eps * y, -moduli(t)[1] * eps * z),
                     rf=lambda t: (-16.0 * moduli(t)[0] * eps, 0.0, 0.0),
                     elem=lambda t: {'mises': moduli(t)[0] * eps, 'smem': moduli(t)[0] * eps}), {}


def connector_field(case, nodes):
    def u(t):
        return 0.5 * t if t <= 1.0 else 50.0 / (100.0 + 300.0 * (t - 1.0))
    times = [1.0] + ramp(4, 2.0, 1.0)
    return frames_of(case, nodes, times, lambda t, x, y, z: (u(t), 0.0, 0.0),
                     rf=lambda t: (-50.0 * min(t, 1.0), 0.0, 0.0)), {}


def connector_nonlinear(case, nodes):
    def u(F):
        return F / 200.0 if F <= 100.0 else 0.5 + (F - 100.0) / 400.0

    def W(x):
        return 100.0 * x * x if x <= 0.5 else 25.0 + 100.0 * (x - 0.5) + 200.0 * (x - 0.5) ** 2
    fr = frames_of(case, nodes, ramp(4), lambda t, x, y, z: (u(200.0 * t), 0.0, 0.0),
                   rf=lambda t: (-200.0 * t, 0.0, 0.0))
    return fr, {'ALLIE': [[t, W(u(200.0 * t))] for t in ramp(100)]}


# ------------------------------------------------------------------------------------------ patch tests
def membrane_patch_test(case, nodes):
    ex = ey = g = 1e-3
    c = E / (1 - NU * NU)
    sx, sy, txy = c * (ex + NU * ey), c * (ey + NU * ex), G * g
    m = mises_plane(sx, sy, txy)
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (1e-3 * (x + y / 2), 1e-3 * (y + x / 2), 0.0), zero_rot,
                   elem=lambda t: {'mises': m, 'smem': m})
    w = 0.5 * (sx * ex + sy * ey + txy * g) * 24.0 * 12.0 * 0.1
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def plate_patch_test(case, nodes):
    t0 = 0.5
    kx = ky = 1e-3
    chi = 1e-3
    c = E / (1 - NU * NU)
    z = t0 / 2
    sx, sy, txy = c * z * (kx + NU * ky), c * z * (ky + NU * kx), G * z * chi
    D = E * t0 ** 3 / (12 * (1 - NU * NU))
    u = D / 2 * (kx * kx + ky * ky + 2 * NU * kx * ky + (1 - NU) / 2 * chi * chi) * 24.0 * 12.0
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z_: (0.0, 0.0, 1e-3 * (x * x + x * y + y * y) / 2),
                   lambda t, x, y, z_: (1e-3 * (x / 2 + y), -1e-3 * (x + y / 2), 0.0),
                   elem=lambda t: {'mises': mises_plane(sx, sy, txy)})
    return fr, {'ALLIE': [[1.0, u]], 'ALLWK': [[1.0, u]]}


def solid_patch_test(case, nodes):
    e = [1e-3, 1e-3, 1e-3, 0.5e-3, 0.5e-3, 0.5e-3]  # tensor components 11, 22, 33, 12, 13, 23
    lam, mu = E * NU / ((1 + NU) * (1 - 2 * NU)), G
    tr = e[0] + e[1] + e[2]
    sg = [lam * tr + 2 * mu * e[i] for i in range(3)] + [2 * mu * e[i] for i in range(3, 6)]
    m = mises3(sg)
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (1e-3 * (2 * x + y + z) / 2, 1e-3 * (x + 2 * y + z) / 2,
                                                            1e-3 * (x + y + 2 * z) / 2),
                   elem=lambda t: {'mises': m, 'smem': m})
    w = 0.5 * (sum(sg[i] * e[i] for i in range(3)) + 2 * sum(sg[i] * e[i] for i in range(3, 6))) * 1000.0
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def thin_shell_moment(case, nodes):
    t0, L = 0.01, 40.0
    k = 0.01 / L
    ei = E * 4.0 * t0 ** 3 / 12
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (0.0, 0.0, k * x * x / 2),
                   lambda t, x, y, z: (0.0, -k * x, 0.0), elem=lambda t: {'mises': E * k * t0 / 2})
    w = 0.5 * k * ei * 0.01
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def spin_two_turns(case, nodes):
    w = 500.0

    def disp(t, x, y, z):
        c, s_ = math.cos(w * t), math.sin(w * t)
        return (c * x - s_ * y - x, s_ * x + c * y - y, 0.0)
    return frames_of(case, nodes, ramp(24, 4 * math.pi / w), disp), {}


def parts_instances(case, nodes):
    eps = 1e-3
    fr = [{'t': 1.0, 'U': {}}]
    for n in case['mon']:
        x, y, z = nodes[n]  # part coordinates
        if n < 1000:
            fr[0]['U'][str(n)] = [eps * x, -NU * eps * y, -NU * eps * z]
        else:
            fr[0]['U'][str(n)] = [NU * eps * y, eps * x, -NU * eps * z]
    fr[0]['RF'] = [-E * 16.0 * eps, -E * 16.0 * eps, 0.0]
    return fr, {}


def include_parameter(case, nodes):
    return solid_tension(case, nodes)


def amplitude_strip(amp):
    def solve(case, nodes):
        eps = 1e-3
        fr = frames_of(case, nodes, ramp(8), lambda t, x, y, z: (eps * amp(t) * x, -NU * eps * amp(t) * y, 0.0),
                       zero_rot, lambda t: (-E * 0.5 * 10.0 * eps * amp(t), 0.0, 0.0),
                       elem=lambda t: {'mises': E * eps * abs(amp(t)), 'smem': E * eps * abs(amp(t))})
        return fr, {}
    return solve


def equally_spaced(t):
    v = [0.0, 0.5, 1.0, 0.8, 0.6]
    k = min(int(t / 0.25), 3)
    return v[k] + (v[k + 1] - v[k]) * (t - 0.25 * k) / 0.25


def bc_acceleration(case, nodes):
    a, m = 1e5, 7.85e-09 * 64.0
    fr = frames_of(case, nodes, ramp(20, 1e-3), lambda t, x, y, z: (0.5 * a * t * t, 0.0, 0.0),
                   rf=lambda t: (m * a, 0.0, 0.0))
    return fr, {'ALLKE': [[t, 0.5 * m * (a * t) ** 2] for t in ramp(200, 1e-3)]}


def heat_transfer(case, nodes):
    fr = [{'t': 1.0, 'U': {}, 'NT': {str(n): [20.0 + 2.0 * nodes[n][0]] for n in case['mon']}}]
    return fr, {}


def beam_twist(case, nodes):
    """Free B31 elements twisted about their axis by a moment M at both nodes: rigid rotation with the default
    rotary inertia J = rho L / 2 (I11 + I22) per node (RECT: (b h^3 + h b^3) / 12, PIPE: thin wall at the mean
    radius, 2 pi r^3 t), theta = M t^2 / (2 J)."""
    from gen_features import BEAM_TWIST
    rho = 7.85e-9
    polar = {'RECT': (2.0 * 1.0 ** 3 + 1.0 * 2.0 ** 3) / 12, 'PIPE': 2 * math.pi * (3.0 - 0.25) ** 3 * 0.5}
    elems = []
    for sec, dims, p0, p1, moment in BEAM_TWIST:
        length = math.dist(p0, p1)
        axis = [(q - p) / length for p, q in zip(p0, p1)]
        elems.append((p0, p1, axis, moment, rho * length / 2 * polar[sec]))

    def owner(x, y, z):
        for p0, p1, axis, moment, J in elems:
            if math.dist((x, y, z), p0) < 1e-9 or math.dist((x, y, z), p1) < 1e-9:
                return axis, moment, J
        raise KeyError((x, y, z))

    def rot(t, x, y, z):
        axis, moment, J = owner(x, y, z)
        th = moment * t * t / (2 * J)
        return tuple(th * a for a in axis)
    times = ramp(10, 1e-8)
    fr = frames_of(case, nodes, times, lambda t, x, y, z: (0.0, 0.0, 0.0), rot)
    ts = [0.0] + times
    ke = [[t, sum(2 * 0.5 * (moment * t) ** 2 / J for _, _, _, moment, J in elems)] for t in ts]
    return fr, {'ALLKE': ke, 'ALLWK': ke}


SOLUTIONS = {
    'tension': solid_tension, 'thermal': solid_thermal, 'shear_hyper': solid_shear, 'plastic': solid_plastic,
    's4r_membrane': shell_membrane, 's3r_membrane': shell_membrane, 's4r_bending': shell_bending,
    's3r_bending': shell_bending, 's4r_rollup': lambda c, n: rollup(c, n, 40.0, E * 10.0 * 0.5 ** 3 / 12),
    's3r_rollup': lambda c, n: rollup(c, n, 40.0, E * 10.0 * 0.5 ** 3 / 12), 's4r_twist': shell_twist,
    's4r_fung_biaxial': shell_fung, 's3r_fung_biaxial': shell_fung, 's4r_neohooke_biaxial': shell_neohooke,
    's4r_spin': shell_spin, 'b31_tension': b31_tension, 'b31_bending': b31_bending, 'b31_torsion': b31_torsion,
    'b31_rollup': lambda c, n: rollup(c, n, 100.0, E * 4.0 * 1.0 ** 3 / 12), 'b31_thermal': b31_thermal,
    'conn3d2_axial': conn_axial, 'conn3d2_cartesian': conn_cartesian, 'conn3d2_oscillator': conn_oscillator,
    'mpc_beam': mpc_beam, 'equation': equation, 'rigid_body': rigid_body, 'connector_join': connector_join,
    'load_gravity': load_gravity, 'bc_transform': bc_transform, 'initial_stress': initial_stress,
    'temperature_dependence': temperature_dependence, 'connector_field': connector_field,
    'mat_mooney_rivlin': hyper_block(W_mooney), 'mat_ogden': hyper_block(W_ogden), 'mat_yeoh': hyper_block(W_yeoh),
    'mat_kinematic_hardening': kinematic_hardening, 'mat_orthotropic': orthotropic,
    'visco_long_term': visco_long_term,
    'patch_membrane_s4r': membrane_patch_test, 'patch_membrane_s3r': membrane_patch_test,
    'patch_membrane_s4': membrane_patch_test, 'patch_plate_s4': plate_patch_test,
    'patch_plate_s4r': plate_patch_test, 'patch_solid_c3d8r': solid_patch_test, 'patch_solid_c3d8': solid_patch_test, 'patch_solid_c3d8i': solid_patch_test,
    'patch_solid_c3d4': solid_patch_test, 'thin_shell_moment_s4r': thin_shell_moment,
    'thin_shell_moment_s3r': thin_shell_moment, 'spin_two_turns_s4r': spin_two_turns,
    'input_parts_instances': parts_instances, 'input_include_parameter': include_parameter,
    'connector_nonlinear': connector_nonlinear, 'amp_decay': amplitude_strip(lambda t: 1 - math.exp(-t / 0.3)),
    'amp_equally_spaced': amplitude_strip(equally_spaced), 'bc_acceleration': bc_acceleration,
    'heat_transfer_steady': heat_transfer, 'dload_solid_pressure': lambda c, n: dload_solid(c, n),
    # the initial connector force is not a load: equilibrium (free length reached) from the first increment
    'connector_reference_length': lambda c, n: (frames_of(c, n, ramp(10), lambda t, x, y, z: (-1.0, 0.0, 0.0),
                                                          rf=lambda t: (0.0, 0.0, 0.0)), {}),
    'elem_b31h': b31_bending, 'beam_rotary_inertia_torsion_xpl': beam_twist,
}


def dload_solid(case, nodes):
    p = 100.0
    ez, el = -p / E, NU * p / E
    fr = frames_of(case, nodes, [1.0], lambda t, x, y, z: (el * x, el * y, ez * z), rf=lambda t: (0.0, 0.0, p * 16.0),
                   elem=lambda t: {'mises': p, 'smem': p})
    w = 0.5 * p * p / E * 64.0
    return fr, {'ALLIE': [[1.0, w]], 'ALLWK': [[1.0, w]]}


def reference(case):
    name = case['name']
    fn = SOLUTIONS.get(name) or SOLUTIONS.get(name.split('_', 1)[1] if '_' in name else name)
    if fn is None:
        raise KeyError('no closed-form solution for ' + name)
    frames, history = fn(case, deck_nodes(name))
    return {'case': name, 'source': 'analytic', 'frames': frames, 'history': history}


def main():
    cases = [c for c in json.load(open(os.path.join(HERE, 'cases', 'manifest.json'))) if c['reference'] == 'analytic']
    names = sys.argv[1:]
    for c in cases:
        if names and c['name'] not in names:
            continue
        r = reference(c)
        f = r['frames'][-1]
        print(c['name'], 't=%g' % f['t'], 'U', {k: [round(x, 5) for x in v] for k, v in list(f['U'].items())[:1]},
              'RF', [round(x, 4) for x in f.get('RF', [])])


if __name__ == '__main__':
    main()
