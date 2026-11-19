#include "material/material.hpp"

#include <cmath>

namespace abq {

namespace {
inline double ddot(const Sym6& a, const Sym6& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]);
}
}  // namespace

double dilatational_modulus(const Material& m) {
  if (m.elastic) {
    double E = m.E, nu = m.nu;
    return E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
  }
  double G = m.shear_modulus(), K = m.bulk_modulus();
  return K + 4.0 / 3.0 * G;
}

double plane_stress_modulus(const Material& m) {
  if (m.elastic) return m.E / (1 - m.nu * m.nu);
  double G = m.shear_modulus(), K = m.bulk_modulus();
  // E/(1-nu^2) from (K, G)
  double E = 9 * K * G / (3 * K + G), nu = (3 * K - 2 * G) / (2 * (3 * K + G));
  return E / (1 - nu * nu);
}

double j2_update_3d(const Material& m, const Sym6& d, Sym6& s, double& ep) {
  const double E = m.E, nu = m.nu;
  const double G = E / (2 * (1 + nu)), lam = E * nu / ((1 + nu) * (1 - 2 * nu));
  double tr = d[0] + d[1] + d[2];
  Sym6 t = s;
  for (int i = 0; i < 3; ++i) t[i] += lam * tr + 2 * G * d[i];
  for (int i = 3; i < 6; ++i) t[i] += 2 * G * d[i];
  if (!m.plastic()) { s = t; return 0; }
  double p = (t[0] + t[1] + t[2]) / 3.0;
  Sym6 dev = t;
  for (int i = 0; i < 3; ++i) dev[i] -= p;
  double q = std::sqrt(1.5 * ddot(dev, dev));
  double sy0 = m.yield(ep);
  if (q <= sy0) { s = t; return 0; }
  double dp = 0;
  for (int it = 0; it < 50; ++it) {
    double H;
    double sy = m.yield(ep + dp, &H);
    double r = q - 3 * G * dp - sy;
    double dr = -3 * G - H;
    double ddp = -r / dr;
    dp += ddp;
    if (dp < 0) dp = 0;
    if (std::fabs(ddp) <= 1e-12 * (ep + dp + 1e-12)) break;
  }
  double f = 1 - 3 * G * dp / q;
  for (int i = 0; i < 6; ++i) s[i] = dev[i] * f;
  for (int i = 0; i < 3; ++i) s[i] += p;
  double sy1 = m.yield(ep + dp);
  ep += dp;
  return 0.5 * (sy0 + sy1) * dp;
}

double j2_update_ps(const Material& m, double d11, double d22, double d12, double sig[3], double& ep, double& d33) {
  const double E = m.E, nu = m.nu;
  if (!m.plastic()) {
    double c = E / (1 - nu * nu);
    sig[0] += c * (d11 + nu * d22);
    sig[1] += c * (d22 + nu * d11);
    sig[2] += E / (1 + nu) * d12;
    d33 = -nu / (1 - nu) * (d11 + d22);
    return 0;
  }
  // generic: 3D return mapping, thickness strain found by secant iterations on sigma33 = 0
  const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), G = E / (2 * (1 + nu));
  Sym6 s0;
  s0[0] = sig[0]; s0[1] = sig[1]; s0[3] = sig[2];
  auto eval = [&](double x, Sym6& so, double& epo, double& w) {
    Sym6 d;
    d[0] = d11; d[1] = d22; d[2] = x; d[3] = d12;
    so = s0; epo = ep;
    w = j2_update_3d(m, d, so, epo);
    return so[2];
  };
  double x0 = -nu / (1 - nu) * (d11 + d22);
  Sym6 so; double epo, w;
  double f0 = eval(x0, so, epo, w);
  double scale = std::fabs(so[0]) + std::fabs(so[1]) + std::fabs(so[3]) + m.yield(ep) * 1e-3 + 1e-30;
  double x1 = x0 - f0 / (lam + 2 * G);
  for (int it = 0; it < 30 && std::fabs(f0) > 1e-10 * scale; ++it) {
    double f1 = eval(x1, so, epo, w);
    if (std::fabs(f1) <= 1e-10 * scale) { x0 = x1; f0 = f1; break; }
    double den = f1 - f0;
    double x2 = std::fabs(den) > 1e-300 ? x1 - f1 * (x1 - x0) / den : x1 - f1 / (lam + 2 * G);
    x0 = x1; f0 = f1; x1 = x2;
  }
  eval(x0, so, epo, w);
  sig[0] = so[0]; sig[1] = so[1]; sig[2] = so[3];
  ep = epo;
  d33 = x0;
  return w;
}

double j2_update_beam(const Material& m, double G, double d11, double dg12, double dg13, double& s11, double& t12,
                      double& t13, double& ep) {
  const double E = m.E;
  double s = s11 + E * d11, a = t12 + G * dg12, b = t13 + G * dg13;
  if (!m.plastic()) { s11 = s; t12 = a; t13 = b; return 0; }
  double qtr = std::sqrt(s * s + 3 * (a * a + b * b));
  double sy0 = m.yield(ep);
  if (qtr <= sy0) { s11 = s; t12 = a; t13 = b; return 0; }
  auto resid = [&](double mu, double& q) {
    double ss = s / (1 + E * mu), ta = a / (1 + 3 * G * mu), tb = b / (1 + 3 * G * mu);
    q = std::sqrt(ss * ss + 3 * (ta * ta + tb * tb));
    return q - m.yield(ep + mu * q);
  };
  double lo = 0, hi = 1.0 / std::min(E, 3 * G), q;
  while (resid(hi, q) > 0) hi *= 4;
  double flo = qtr - sy0, fhi = resid(hi, q);
  double mu = hi;
  for (int it = 0; it < 100; ++it) {  // Illinois regula falsi
    mu = (lo * fhi - hi * flo) / (fhi - flo);
    double f = resid(mu, q);
    if (std::fabs(f) < 1e-12 * sy0) break;
    if (f > 0) { lo = mu; flo = f; fhi *= 0.5; }
    else { hi = mu; fhi = f; flo *= 0.5; }
    if (hi - lo < 1e-15 * hi) break;
  }
  resid(mu, q);
  s11 = s / (1 + E * mu);
  t12 = a / (1 + 3 * G * mu);
  t13 = b / (1 + 3 * G * mu);
  double dp = mu * q;
  double sy1 = m.yield(ep + dp);
  ep += dp;
  return 0.5 * (sy0 + sy1) * dp;
}

// ------------------------------------------------------------------ hyperelasticity
namespace {

// 2nd Piola-Kirchhoff stress (tensor components, material basis) from C (Sym6).
void hyper_pk2(const Material& m, const Sym6& C, Sym6& S, double* energy) {
  Mat3 Cm = C.full();
  double J2 = Cm.det();
  double J = std::sqrt(std::max(J2, 1e-300));
  Mat3 Ci = Cm.inverse();
  double Jm23 = std::pow(J, -2.0 / 3.0);
  Sym6 Sb;  // dU/dEbar (= 2 dU/dCbar)
  double Udev = 0, dUdJ = 0, Uvol = 0;
  if (m.fung) {
    Sym6 Eb;
    for (int i = 0; i < 3; ++i) Eb[i] = 0.5 * (Jm23 * C[i] - 1.0);
    for (int i = 3; i < 6; ++i) Eb[i] = 0.5 * Jm23 * C[i];
    const double w[6] = {1, 1, 1, 2, 2, 2};
    double bE[6];
    double Q = 0;
    for (int I = 0; I < 6; ++I) {
      double sI = 0;
      for (int Jx = 0; Jx < 6; ++Jx) sI += m.fung_b[I][Jx] * w[Jx] * Eb[Jx];
      bE[I] = sI;
      Q += w[I] * Eb[I] * sI;
    }
    Q = std::min(Q, 600.0);
    double eQ = std::exp(Q);
    for (int I = 0; I < 6; ++I) Sb[I] = m.fung_c * eQ * bE[I];
    Udev = 0.5 * m.fung_c * (eQ - 1);
    double D = m.fung_D > 0 ? m.fung_D : 1.0 / (10.0 * std::max(m.shear_modulus(), 1e-30));
    dUdJ = (J - 1.0 / J) / D;
    Uvol = ((J * J - 1) * 0.5 - std::log(J)) / D;
  } else if (m.neohooke) {
    for (int i = 0; i < 3; ++i) Sb[i] = 2 * m.C10;
    double I1b = Jm23 * (C[0] + C[1] + C[2]);
    Udev = m.C10 * (I1b - 3);
    double D1 = m.D1 > 0 ? m.D1 : 1.0 / (10.0 * std::max(m.C10, 1e-30));
    dUdJ = 2.0 / D1 * (J - 1);
    Uvol = (J - 1) * (J - 1) / D1;
  }
  // DEV(Sb) = Sb - 1/3 (Sb:C) C^-1
  double SC = Sb[0] * C[0] + Sb[1] * C[1] + Sb[2] * C[2] + 2 * (Sb[3] * C[3] + Sb[4] * C[4] + Sb[5] * C[5]);
  Sym6 Cis = Sym6::from(Ci);
  for (int i = 0; i < 6; ++i) S[i] = Jm23 * (Sb[i] - SC / 3.0 * Cis[i]) + J * dUdJ * Cis[i];
  if (energy) *energy = Udev + Uvol;
}

// Fung, plane-stress layout C = [[C11 C12 0][C12 C22 0][0 0 x]]: S11, S22, S12, S33 and dS33/dx
// (block algebra, exact derivative; used by the plane-stress Newton on the thickness stretch).
struct FungPS {
  double S[4];   // 11, 22, 12, 33
  double dS33;   // d S33 / d C33
  double energy;
  double Q;      // Fung exponent
};
inline void fung_ps(const Material& m, const double C2[3], double x, FungPS& r, bool want_energy) {
  const double C11 = C2[0], C22 = C2[1], C12 = C2[2];
  const double a = C11 * C22 - C12 * C12;
  const double ax = a * x;
  const double J = std::sqrt(ax);
  const double j = 1.0 / std::cbrt(ax);  // J^(-2/3)
  // deviatoric Green strain, components 11, 22, 33, 12 (Voigt slots 0, 1, 2, 3)
  const double E[4] = {0.5 * (j * C11 - 1.0), 0.5 * (j * C22 - 1.0), 0.5 * (j * x - 1.0), 0.5 * j * C12};
  const double dE[4] = {-j * C11 / (6 * x), -j * C22 / (6 * x), j / 3.0, -j * C12 / (6 * x)};
  static const int slot[4] = {0, 1, 2, 3};
  const double w[4] = {1, 1, 1, 2};
  double bE[4], bdE[4], Q = 0, dQ = 0;
  for (int I = 0; I < 4; ++I) {
    double s = 0, sd = 0;
    for (int K = 0; K < 4; ++K) {
      double b = m.fung_b[slot[I]][slot[K]];
      s += b * w[K] * E[K];
      sd += b * w[K] * dE[K];
    }
    bE[I] = s;
    bdE[I] = sd;
    Q += w[I] * E[I] * s;
    dQ += 2 * w[I] * dE[I] * s;
  }
  Q = std::min(Q, 600.0);
  const double ceQ = m.fung_c * std::exp(Q);
  r.Q = Q;
  double Sb[4], dSb[4];
  for (int I = 0; I < 4; ++I) { Sb[I] = ceQ * bE[I]; dSb[I] = ceQ * (bE[I] * dQ + bdE[I]); }
  const double D = m.fung_D > 0 ? m.fung_D : 1.0 / (10.0 * std::max(m.shear_modulus(), 1e-30));
  // S = j (Sb - T/3 Cinv) + (J^2 - 1)/D Cinv, T = Sb:C ; in-plane Cinv = adj/a, Cinv33 = 1/x
  const double T = Sb[0] * C11 + Sb[1] * C22 + Sb[2] * x + 2 * Sb[3] * C12;
  const double Ci11 = C22 / a, Ci22 = C11 / a, Ci12 = -C12 / a, Ci33 = 1.0 / x;
  const double pv = (ax - 1.0) / D;
  r.S[0] = j * (Sb[0] - T / 3 * Ci11) + pv * Ci11;
  r.S[1] = j * (Sb[1] - T / 3 * Ci22) + pv * Ci22;
  r.S[2] = j * (Sb[3] - T / 3 * Ci12) + pv * Ci12;
  r.S[3] = j * (Sb[2] - T / 3 * Ci33) + pv * Ci33;
  const double dj = -j / (3 * x);
  const double dT = dSb[0] * C11 + dSb[1] * C22 + dSb[2] * x + Sb[2] + 2 * dSb[3] * C12;
  // S33 = j (Sb33 - T/(3x)) + a/D - 1/(D x)
  r.dS33 = dj * (Sb[2] - T / (3 * x)) + j * (dSb[2] - dT / (3 * x) + T / (3 * x * x)) + 1.0 / (D * x * x);
  r.energy = want_energy ? 0.5 * m.fung_c * (std::exp(Q) - 1) + ((ax - 1) * 0.5 - std::log(J)) / D : 0.0;
}

}  // namespace

void hyper_cauchy_3d(const Material& m, const Mat3& F, Sym6& sigma, double* energy) {
  Sym6 C = Sym6::from(F.T() * F);
  Sym6 S;
  hyper_pk2(m, C, S, energy);
  double J = F.det();
  sigma = Sym6::from(F * S.full() * F.T() * (1.0 / J));
}

bool hyper_ps(const Material& m, const double C2[3], double& lam3, double S[3], double* energy, double* qexp) {
  if (qexp) *qexp = -1;
  double detC2 = C2[0] * C2[1] - C2[2] * C2[2];
  if (!(detC2 > 1e-12) || !(C2[0] > 0) || !(C2[1] > 0)) {
    S[0] = S[1] = S[2] = 0;
    lam3 = 1;
    if (energy) *energy = 0;
    return false;
  }
  if (m.fung) {  // fast path: Newton on x = lambda3^2 with the exact derivative, bracketed
    double x0 = 1.0 / detC2;  // incompressible estimate
    double x = (lam3 > 0 && std::fabs(lam3 * lam3 / x0 - 1) < 0.4) ? lam3 * lam3 : x0;
    FungPS r;
    fung_ps(m, C2, x, r, false);
    double scale = std::fabs(r.S[0]) + std::fabs(r.S[1]) + std::fabs(r.S[2]) + 1e-9 * m.shear_modulus() + 1e-300;
    // absolute floor: round-off of the volumetric term (J^2 - 1)/D near J = 1 (unloaded points)
    const double Dv = m.fung_D > 0 ? m.fung_D : 1.0 / (10.0 * std::max(m.shear_modulus(), 1e-30));
    const double tol = std::max(1e-10 * scale, 1e-14 / Dv), tol_loose = std::max(1e-6 * scale, 1e-12 / Dv);
    double lo = 0, hi = 1e300;
    bool ok = false;
    for (int it = 0; it < 60; ++it) {
      double f = r.S[3];
      if (std::fabs(f) <= tol) { ok = true; break; }
      if (f > 0) hi = std::min(hi, x); else lo = std::max(lo, x);
      double xn = r.dS33 > 0 ? x - f / r.dS33 : -1;
      if (!(xn > lo && xn < hi)) xn = hi < 1e299 ? 0.5 * (lo + hi) : (f > 0 ? 0.5 * (x + lo) : 2 * x);
      bool tiny = std::fabs(xn - x) <= 1e-15 * x;
      x = xn;
      fung_ps(m, C2, x, r, false);
      if (tiny) { ok = std::fabs(r.S[3]) <= tol_loose; break; }
    }
    if (energy) fung_ps(m, C2, x, r, true), *energy = r.energy;
    if (qexp) *qexp = r.Q;
    lam3 = std::sqrt(x);
    S[0] = r.S[0]; S[1] = r.S[1]; S[2] = r.S[2];
    return ok;
  }
  auto s33 = [&](double l3, Sym6& Sout, double* en) {
    Sym6 C;
    C[0] = C2[0]; C[1] = C2[1]; C[2] = l3 * l3; C[3] = C2[2];
    hyper_pk2(m, C, Sout, en);
    return Sout[2];
  };
  Sym6 Sx;
  double guess = 1.0 / std::sqrt(std::max(detC2, 1e-30));  // incompressible estimate
  if (lam3 > 0 && std::fabs(lam3 / guess - 1) < 0.2) guess = lam3;
  double x = guess, fx = s33(x, Sx, nullptr);
  double scale = std::fabs(Sx[0]) + std::fabs(Sx[1]) + std::fabs(Sx[3]) + 1e-9 * m.shear_modulus() + 1e-300;
  double lo = 0, hi = 1e300;
  bool ok = false;
  // safeguarded Newton (S33 increases with the thickness stretch)
  for (int it = 0; it < 60; ++it) {
    if (std::fabs(fx) <= 1e-10 * scale) { ok = true; break; }
    if (fx > 0) hi = std::min(hi, x); else lo = std::max(lo, x);
    double h = 1e-7 * x;
    double fp = s33(x + h, Sx, nullptr);
    double d = (fp - fx) / h;
    double xn = (d > 0) ? x - fx / d : (fx > 0 ? 0.5 * (x + lo) : 2 * x);
    if (xn <= lo || xn >= hi) xn = hi < 1e299 ? 0.5 * (lo + hi) : (fx > 0 ? 0.5 * (x + lo) : 2 * x);
    if (std::fabs(xn - x) <= 1e-14 * x) { x = xn; fx = s33(x, Sx, nullptr); ok = std::fabs(fx) <= 1e-6 * scale; break; }
    x = xn;
    fx = s33(x, Sx, nullptr);
  }
  lam3 = x;
  s33(lam3, Sx, energy);
  S[0] = Sx[0]; S[1] = Sx[1]; S[2] = Sx[3];
  return ok;
}

}  // namespace abq
