#include "material/material.hpp"

#include <cmath>

namespace abq {

namespace {
inline double ddot(const Sym6& a, const Sym6& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]);
}
}  // namespace

double hgo_fibre_modulus(const Material& m, const Mat3& F, const double* fib) {
  if (!m.hgo || !fib) return 0;
  const Mat3 C = F.T() * F;
  const double J = std::max(F.det(), 1e-12), Jm23 = std::pow(J, -2.0 / 3.0);
  const double I1b = Jm23 * (C.m[0][0] + C.m[1][1] + C.m[2][2]), ka = m.hgo_kappa;
  double kmax = 0;
  for (int f = 0; f < m.hgo_nfib; ++f) {
    const Vec3 a(fib[3 * f], fib[3 * f + 1], fib[3 * f + 2]);
    const double I4 = Jm23 * dot(a, C * a);
    const double E = std::max(0.0, ka * (I1b - 3) + (1 - 3 * ka) * (I4 - 1));
    const double ex = std::exp(std::min(m.hgo_k2 * E * E, 600.0));
    // d2U/dE2 (dE/deps)^2 + dU/dE d2E/deps2 along the fibre, dE/deps = 2 (1 - 3 kappa) I4
    const double c = 2 * (1 - 3 * ka) * I4;
    kmax = std::max(kmax, m.hgo_k1 * ex * (1 + 2 * m.hgo_k2 * E * E) * c * c + 2 * m.hgo_k1 * E * ex * c);
  }
  return kmax;
}

double dilatational_modulus(const Material& m) {
  if (m.sma) {
    const double E = std::max(m.sma_EA, m.sma_EM), nu = m.sma_nuA;
    return E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
  }
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

// Linear kinematic hardening (Prager): yield on the shifted deviator xi = dev(sigma) - alpha, backstress increment
// (2/3) C deps_p, radial return with dp = f / (3 G + C).
double j2_kinematic_3d(const Material& m, const Sym6& d, Sym6& s, Sym6& alpha, double& ep) {
  const double E = m.E, nu = m.nu;
  const double G = E / (2 * (1 + nu)), lam = E * nu / ((1 + nu) * (1 - 2 * nu));
  const double tr = d[0] + d[1] + d[2];
  Sym6 t = s;
  for (int i = 0; i < 3; ++i) t[i] += lam * tr + 2 * G * d[i];
  for (int i = 3; i < 6; ++i) t[i] += 2 * G * d[i];
  const double p = (t[0] + t[1] + t[2]) / 3.0;
  Sym6 xi = t;
  for (int i = 0; i < 3; ++i) xi[i] -= p;
  for (int i = 0; i < 6; ++i) xi[i] -= alpha[i];
  const double q = std::sqrt(1.5 * ddot(xi, xi)), sy = m.plastic_sy.front();
  if (q <= sy) { s = t; return 0; }
  const double dp = (q - sy) / (3 * G + m.kin_C);
  for (int i = 0; i < 6; ++i) {
    const double n = xi[i] / q;  // 3/2 n = flow direction per unit dp
    t[i] -= 3 * G * dp * n;
    alpha[i] += m.kin_C * dp * n;
  }
  s = t;
  ep += dp;
  return sy * dp;
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

// ------------------------------------------------------------------ superelasticity
// Generalized-plasticity model of the austenite-martensite transformation (Auricchio and Taylor, 1997):
//  - elasticity with moduli mixed linearly between austenite and martensite, sigma = C(xi) (eps - etr);
//  - loading function of Drucker-Prager type F = q + (tan b / 3) tr(sigma), tan b = 3 (s_cs - s_ls) / (s_cs + s_ls)
//    from the start stresses in compression and tension; thresholds are the uniaxial stresses times (1 + tan b / 3)
//    shifted by the stress-temperature slopes;
//  - transformation strain along dF/dsigma, scaled to eps_L in uniaxial tension; the reverse transformation retraces
//    the transformation strain;
//  - linear kinetics integrated exactly over the increment: forward dxi / (1 - xi) = dF / (R_f - F), reverse
//    dxi / xi = dF / (F - R_f).
namespace {

struct SmaState {
  double lam, mu;
};
inline SmaState sma_moduli(const Material& m, double xi) {
  const double E = m.sma_EA + xi * (m.sma_EM - m.sma_EA), nu = m.sma_nuA + xi * (m.sma_nuM - m.sma_nuA);
  return {E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))};
}
inline void sma_stress(const Material& m, double xi, const Sym6& ee, Sym6& sig) {
  const SmaState k = sma_moduli(m, xi);
  const double tr = ee[0] + ee[1] + ee[2];
  for (int i = 0; i < 6; ++i) sig[i] = 2 * k.mu * ee[i];
  for (int i = 0; i < 3; ++i) sig[i] += k.lam * tr;
}
inline double sma_F(const Sym6& s, double tb) {
  const double p = (s[0] + s[1] + s[2]) / 3;
  double j2 = 0;
  for (int i = 0; i < 3; ++i) j2 += (s[i] - p) * (s[i] - p);
  for (int i = 3; i < 6; ++i) j2 += 2 * s[i] * s[i];
  return std::sqrt(1.5 * j2) + tb * p;
}

}  // namespace

void sma_update(const Material& m, const Sym6& eps, Sym6& etr, double& xi, double T, const Sym6& sig_old, Sym6& sig) {
  const double tb = 3 * (m.sma_sCLs - m.sma_sLs) / std::max(m.sma_sCLs + m.sma_sLs, 1e-30);
  const double k = 1 + tb / 3;
  const double dT = T - m.sma_T0;
  const double RsAS = (m.sma_sLs + m.sma_CL * dT) * k, RfAS = (m.sma_sLe + m.sma_CL * dT) * k;
  const double RsSA = (m.sma_sUs + m.sma_CU * dT) * k, RfSA = (m.sma_sUe + m.sma_CU * dT) * k;
  const double Fn = sma_F(sig_old, tb);
  const double xin = xi;
  const Sym6 etr_n = etr;
  // stress for a fraction x with the transformation strain moved along dir; returns F with the deviator projected on
  // the unit deviatoric direction nd (signed, so that F decreases monotonically as the transformation proceeds), or
  // the plain F when nd is null
  auto state = [&](double x, const Sym6& dir, double scale, Sym6& e_out, Sym6& s_out, const Sym6* nd) {
    for (int i = 0; i < 6; ++i) e_out[i] = etr_n[i] + scale * (x - xin) * dir[i];
    Sym6 ee;
    for (int i = 0; i < 6; ++i) ee[i] = eps[i] - e_out[i];
    sma_stress(m, x, ee, s_out);
    if (!nd) return sma_F(s_out, tb);
    const double p = (s_out[0] + s_out[1] + s_out[2]) / 3;
    double sn = 0;
    for (int i = 0; i < 6; ++i) sn += (i < 3 ? 1.0 : 2.0) * (s_out[i] - (i < 3 ? p : 0.0)) * (*nd)[i];
    return std::sqrt(1.5) * sn + tb * p;
  };
  auto unit_dev = [](const Sym6& a, Sym6& nd) {
    const double p = (a[0] + a[1] + a[2]) / 3;
    double n2 = 0;
    for (int i = 0; i < 6; ++i) { nd[i] = a[i] - (i < 3 ? p : 0.0); n2 += (i < 3 ? 1.0 : 2.0) * nd[i] * nd[i]; }
    if (!(n2 > 0)) return false;
    for (int i = 0; i < 6; ++i) nd[i] /= std::sqrt(n2);
    return true;
  };
  Sym6 e_tr, s_tr, zero;
  const double Ftr = state(xin, zero, 0.0, e_tr, s_tr, nullptr);
  sig = s_tr;
  // forward transformation
  if (Ftr > Fn && Ftr > RsAS && xin < 1) {
    Sym6 N;  // dF/dsigma at the trial stress, scaled to eps_L in uniaxial tension
    const double p = (s_tr[0] + s_tr[1] + s_tr[2]) / 3, q = sma_F(s_tr, 0.0);
    if (q <= 0) return;
    for (int i = 0; i < 6; ++i) N[i] = 1.5 * (s_tr[i] - (i < 3 ? p : 0.0)) / q / k;
    for (int i = 0; i < 3; ++i) N[i] += tb / 3 / k;
    Sym6 nd;
    if (!unit_dev(s_tr, nd)) return;
    const double F0 = std::max(Fn, RsAS);
    auto g = [&](double x) {
      Sym6 e, s;
      const double F = std::min(state(x, N, m.sma_epsL, e, s, &nd), RfAS);
      return (1 - x) * (RfAS - F0) - (1 - xin) * (RfAS - F);
    };
    double lo = xin, hi = 1.0;  // g decreases from g(xi_n) > 0
    if (g(hi) >= 0) lo = hi;
    else
      for (int it = 0; it < 80 && hi - lo > 1e-14; ++it) {
        const double mid = 0.5 * (lo + hi);
        (g(mid) > 0 ? lo : hi) = mid;
      }
    xi = lo;
    state(xi, N, m.sma_epsL, etr, sig, nullptr);
    return;
  }
  // reverse transformation: the transformation strain shrinks in proportion to xi
  if (Ftr < Fn && Ftr < RsSA && xin > 0) {
    Sym6 D;
    for (int i = 0; i < 6; ++i) D[i] = etr_n[i] / xin;
    Sym6 nd;
    if (!unit_dev(etr_n, nd)) return;
    const double F0 = std::min(Fn, RsSA);
    auto g = [&](double x) {
      Sym6 e, s;
      const double F = std::max(state(x, D, 1.0, e, s, &nd), RfSA);
      return x * (F0 - RfSA) - xin * (F - RfSA);
    };
    double lo = 0.0, hi = xin;  // g increases with xi up to g(xi_n) > 0
    if (g(lo) >= 0) hi = lo;
    else
      for (int it = 0; it < 80 && hi - lo > 1e-14; ++it) {
        const double mid = 0.5 * (lo + hi);
        (g(mid) > 0 ? hi : lo) = mid;
      }
    xi = hi;
    state(xi, D, 1.0, etr, sig, nullptr);
  }
}

// ------------------------------------------------------------------ hyperelasticity
namespace {

// 2nd Piola-Kirchhoff stress (tensor components, material basis) from C (Sym6).
void hyper_pk2(const Material& m, const Sym6& C, Sym6& S, double* energy, const double* fib = nullptr) {
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
  } else if (m.hgo) {
    // U = C10 (I1b - 3) + k1 / (2 k2) sum (exp(k2 <E_a>^2) - 1) + ((J^2 - 1) / 2 - ln J) / D,
    // E_a = kappa (I1b - 3) + (1 - 3 kappa)(I4b_a - 1); a fibre family acts in tension only (E_a > 0)
    const double I1b = Jm23 * (C[0] + C[1] + C[2]);
    const double ka = m.hgo_kappa;
    double dI1 = m.hgo_C10;
    double Ufib = 0;
    double dI4[3] = {0, 0, 0};
    for (int f = 0; f < m.hgo_nfib && fib; ++f) {
      const double* a = fib + 3 * f;
      const double aCa = C[0] * a[0] * a[0] + C[1] * a[1] * a[1] + C[2] * a[2] * a[2] +
                         2 * (C[3] * a[0] * a[1] + C[4] * a[0] * a[2] + C[5] * a[1] * a[2]);
      const double E = ka * (I1b - 3) + (1 - 3 * ka) * (Jm23 * aCa - 1);
      if (E <= 0) continue;
      const double ex = std::exp(std::min(m.hgo_k2 * E * E, 600.0));
      Ufib += m.hgo_k1 / (2 * m.hgo_k2) * (ex - 1);
      const double dUdE = m.hgo_k1 * E * ex;
      dI1 += dUdE * ka;
      dI4[f] = dUdE * (1 - 3 * ka);
    }
    for (int i = 0; i < 3; ++i) Sb[i] = 2 * dI1;
    for (int i = 3; i < 6; ++i) Sb[i] = 0;
    for (int f = 0; f < m.hgo_nfib && fib; ++f) {
      const double* a = fib + 3 * f;
      Sb[0] += 2 * dI4[f] * a[0] * a[0]; Sb[1] += 2 * dI4[f] * a[1] * a[1]; Sb[2] += 2 * dI4[f] * a[2] * a[2];
      Sb[3] += 2 * dI4[f] * a[0] * a[1]; Sb[4] += 2 * dI4[f] * a[0] * a[2]; Sb[5] += 2 * dI4[f] * a[1] * a[2];
    }
    Udev = m.hgo_C10 * (I1b - 3) + Ufib;
    const double D = m.hgo_D > 0 ? m.hgo_D : 1.0 / (10.0 * std::max(m.hgo_C10, 1e-30));
    dUdJ = (J - 1.0 / J) / D;
    Uvol = ((J * J - 1) * 0.5 - std::log(J)) / D;
  } else if (m.mooney || m.yeoh) {
    const double I1b = Jm23 * (C[0] + C[1] + C[2]);
    double dU1 = m.C10, dU2 = 0;
    if (m.mooney) {
      dU2 = m.C01;
      for (int i = 0; i < 3; ++i) Sb[i] = 2 * (dU1 + dU2 * I1b);
      for (int i = 0; i < 6; ++i) Sb[i] -= 2 * dU2 * Jm23 * C[i];
      double trC2 = 0;
      for (int i = 0; i < 3; ++i) trC2 += C[i] * C[i];
      for (int i = 3; i < 6; ++i) trC2 += 2 * C[i] * C[i];
      const double I2b = 0.5 * (I1b * I1b - Jm23 * Jm23 * trC2);
      Udev = m.C10 * (I1b - 3) + m.C01 * (I2b - 3);
    } else {
      const double x = I1b - 3;
      dU1 = m.C10 + 2 * m.C20 * x + 3 * m.C30 * x * x;
      for (int i = 0; i < 3; ++i) Sb[i] = 2 * dU1;
      Udev = m.C10 * x + m.C20 * x * x + m.C30 * x * x * x;
    }
    const double Ds[3] = {m.D1 > 0 ? m.D1 : 1.0 / (10.0 * std::max(m.C10, 1e-30)), m.D2, m.D3};
    dUdJ = 0; Uvol = 0;
    for (int i = 0; i < (m.yeoh ? 3 : 1); ++i) {
      if (!(Ds[i] > 0)) continue;
      const int n = 2 * (i + 1);
      dUdJ += n * std::pow(J - 1, n - 1) / Ds[i];
      Uvol += std::pow(J - 1, n) / Ds[i];
    }
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

// eigenvalues and eigenvectors (columns of V) of a symmetric 3 x 3 matrix (cyclic Jacobi)
static void sym_eigen(Mat3 A, double* w, Mat3& V) {
  V = Mat3::identity();
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = std::fabs(A.m[0][1]) + std::fabs(A.m[0][2]) + std::fabs(A.m[1][2]);
    if (off < 1e-15 * (std::fabs(A.m[0][0]) + std::fabs(A.m[1][1]) + std::fabs(A.m[2][2]) + 1e-300)) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (A.m[p][q] == 0) continue;
        const double th = 0.5 * (A.m[q][q] - A.m[p][p]) / A.m[p][q];
        const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double c = 1 / std::sqrt(t * t + 1), sn = t * c;
        for (int k = 0; k < 3; ++k) {  // A <- A J
          const double akp = A.m[k][p], akq = A.m[k][q];
          A.m[k][p] = c * akp - sn * akq; A.m[k][q] = sn * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {  // A <- J^T A
          const double apk = A.m[p][k], aqk = A.m[q][k];
          A.m[p][k] = c * apk - sn * aqk; A.m[q][k] = sn * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
          const double vkp = V.m[k][p], vkq = V.m[k][q];
          V.m[k][p] = c * vkp - sn * vkq; V.m[k][q] = sn * vkp + c * vkq;
        }
      }
  }
  for (int i = 0; i < 3; ++i) w[i] = A.m[i][i];
}

void hyper_cauchy_3d(const Material& m, const Mat3& F, Sym6& sigma, double* energy, const double* fib) {
  if (m.ogden) {  // principal stretches of b = F F^T: tau_i = sum_p 2 mu_p / alpha_p (lb_i^a_p - mean_j lb_j^a_p)
    const Mat3 b = F * F.T();
    double w[3];
    Mat3 V;
    sym_eigen(b, w, V);
    const double J = F.det(), Jm13 = std::pow(J, -1.0 / 3.0);
    double lb[3];
    for (int i = 0; i < 3; ++i) lb[i] = Jm13 * std::sqrt(std::max(w[i], 1e-300));
    double tau[3] = {0, 0, 0}, U = 0;
    for (size_t p = 0; p < m.og_mu.size(); ++p) {
      const double mu = m.og_mu[p], a = m.og_alpha[p];
      double la[3], mean = 0;
      for (int i = 0; i < 3; ++i) { la[i] = std::pow(lb[i], a); mean += la[i] / 3; }
      for (int i = 0; i < 3; ++i) tau[i] += 2 * mu / a * (la[i] - mean);
      U += 2 * mu / (a * a) * (3 * mean - 3);
    }
    double dUdJ = 0;
    for (size_t i = 0; i < m.og_D.size(); ++i) {
      double D = m.og_D[i];
      if (i == 0 && !(D > 0)) D = 2.0 / (20 * m.shear_modulus());
      if (!(D > 0)) continue;
      const int n = 2 * (int)(i + 1);
      dUdJ += n * std::pow(J - 1, n - 1) / D;
      U += std::pow(J - 1, n) / D;
    }
    Mat3 s;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double v = 0;
        for (int k = 0; k < 3; ++k) v += tau[k] / J * V.m[i][k] * V.m[j][k];
        s.m[i][j] = v + (i == j ? dUdJ : 0.0);
      }
    sigma = Sym6::from(s);
    if (energy) *energy = U;
    return;
  }
  Sym6 C = Sym6::from(F.T() * F);
  Sym6 S;
  hyper_pk2(m, C, S, energy, fib);
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
