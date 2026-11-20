// S4R / S3R / S3RS shells: Belytschko-Lin-Tsay style corotational one-point element (C0 triangle for S3R, S3RS),
// Simpson integration through the thickness, stiffness hourglass control derived from exact element
// modes (in-plane bending of the membrane, curvature variation of the rotations, linear transverse
// shear sampled at the edge midpoints) and a small drilling stiffness.
//  - elastic / J2 plastic materials: hypoelastic update of the corotational stresses at each
//    section point (plane stress, thickness strain integrated);
//  - hyperelastic (Fung, neo-Hooke): total Lagrangian membrane stretch C relative to the reference
//    material basis + incremental curvature, plane stress solved for the thickness stretch.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "element/groups.hpp"
#include "material/material.hpp"

namespace abq {

namespace {

constexpr int MAXPT = 15;

struct Geo {
  Vec3 e1, e2, e3;   // local frame
  double xl[4][3];   // local coordinates (relative to centroid)
  double B1[4], B2[4];
  double A;
  double gam[4];     // hourglass projection (quads)
};

// local frame and one-point gradients for a quad (n=4) or triangle (n=3)
void geometry(const double* x, int n, Geo& g) {
  Vec3 p[4];
  for (int I = 0; I < n; ++I) p[I] = load3(x + 3 * I);
  Vec3 c;
  for (int I = 0; I < n; ++I) c += p[I];
  c = c / n;
  Vec3 s;
  if (n == 4) {
    g.e3 = normalized(cross(p[2] - p[0], p[3] - p[1]));
    s = (p[1] + p[2] - p[0] - p[3]);
  } else {
    g.e3 = normalized(cross(p[1] - p[0], p[2] - p[0]));
    s = p[1] - p[0];
  }
  s -= g.e3 * dot(s, g.e3);
  g.e1 = normalized(s);
  g.e2 = cross(g.e3, g.e1);
  for (int I = 0; I < n; ++I) {
    Vec3 d = p[I] - c;
    g.xl[I][0] = dot(d, g.e1); g.xl[I][1] = dot(d, g.e2); g.xl[I][2] = dot(d, g.e3);
  }
  if (n == 4) {
    double x13 = g.xl[0][0] - g.xl[2][0], x24 = g.xl[1][0] - g.xl[3][0];
    double y13 = g.xl[0][1] - g.xl[2][1], y24 = g.xl[1][1] - g.xl[3][1];
    g.A = 0.5 * (x13 * y24 - x24 * y13);
    double i2A = 1.0 / (2 * g.A);
    g.B1[0] = (g.xl[1][1] - g.xl[3][1]) * i2A; g.B1[1] = (g.xl[2][1] - g.xl[0][1]) * i2A;
    g.B1[2] = (g.xl[3][1] - g.xl[1][1]) * i2A; g.B1[3] = (g.xl[0][1] - g.xl[2][1]) * i2A;
    g.B2[0] = (g.xl[3][0] - g.xl[1][0]) * i2A; g.B2[1] = (g.xl[0][0] - g.xl[2][0]) * i2A;
    g.B2[2] = (g.xl[1][0] - g.xl[3][0]) * i2A; g.B2[3] = (g.xl[2][0] - g.xl[0][0]) * i2A;
    const double h[4] = {1, -1, 1, -1};
    double hx = 0, hy = 0;
    for (int I = 0; I < 4; ++I) { hx += h[I] * g.xl[I][0]; hy += h[I] * g.xl[I][1]; }
    for (int I = 0; I < 4; ++I) g.gam[I] = 0.25 * (h[I] - hx * g.B1[I] - hy * g.B2[I]);
  } else {
    double x1 = g.xl[0][0], x2 = g.xl[1][0], x3 = g.xl[2][0], y1 = g.xl[0][1], y2 = g.xl[1][1], y3 = g.xl[2][1];
    g.A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    double i2A = 1.0 / (2 * g.A);
    g.B1[0] = (y2 - y3) * i2A; g.B1[1] = (y3 - y1) * i2A; g.B1[2] = (y1 - y2) * i2A;
    g.B2[0] = (x3 - x2) * i2A; g.B2[1] = (x1 - x3) * i2A; g.B2[2] = (x2 - x1) * i2A;
    g.B1[3] = g.B2[3] = 0;
    for (int I = 0; I < 4; ++I) g.gam[I] = 0;
  }
}

// max eigenvalue of D*G, D symmetric 3x3 (Voigt, engineering shear), G = sum_I B_I B_I^T
// largest eigenvalue of a symmetric 3x3 matrix (closed form, trigonometric solution of the cubic)
inline double max_eig_sym3(double a11, double a22, double a33, double a12, double a13, double a23) {
  double p1 = a12 * a12 + a13 * a13 + a23 * a23;
  double q = (a11 + a22 + a33) / 3.0;
  double b11 = a11 - q, b22 = a22 - q, b33 = a33 - q;
  double p2 = b11 * b11 + b22 * b22 + b33 * b33 + 2 * p1;
  if (p2 <= 1e-28 * (q * q) || p2 <= 1e-300) return std::max(a11, std::max(a22, a33));
  double p = std::sqrt(p2 / 6.0);
  double detB = b11 * (b22 * b33 - a23 * a23) - a12 * (a12 * b33 - a23 * a13) + a13 * (a12 * a23 - b22 * a13);
  double r = std::clamp(detB / (2 * p * p * p), -1.0, 1.0);
  return q + 2 * p * std::cos(std::acos(r) / 3.0);
}

double max_eig_DG(const double D[3][3], const Geo& g, int n) {
  double G[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int I = 0; I < n; ++I) {
    double bx = g.B1[I], by = g.B2[I];
    double B[3][2] = {{bx, 0}, {0, by}, {by, bx}};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) G[i][j] += B[i][0] * B[j][0] + B[i][1] * B[j][1];
  }
  double L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  double tr = G[0][0] + G[1][1] + G[2][2];
  for (int i = 0; i < 3; ++i) G[i][i] += 1e-12 * tr;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = G[i][j];
      for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
      L[i][j] = (i == j) ? std::sqrt(std::max(s, 1e-300)) : s / L[j][j];
    }
  double S[3][3];  // L^T D L
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) s += L[a][i] * D[a][b] * L[b][j];
      S[i][j] = s;
    }
  return max_eig_sym3(S[0][0], S[1][1], S[2][2], 0.5 * (S[0][1] + S[1][0]), 0.5 * (S[0][2] + S[2][0]), 0.5 * (S[1][2] + S[2][1]));
}

// gradients of the bilinear shape functions at corner J of a quad (local coordinates g.xl): a = dN/dx, b = dN/dy
void corner_gradients(const Geo& g, int J, double a[4], double b[4]) {
  static const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
  double dNx[4], dNe[4], j11 = 0, j12 = 0, j21 = 0, j22 = 0;
  for (int I = 0; I < 4; ++I) {
    dNx[I] = 0.25 * XI[I] * (1 + ET[I] * ET[J]);
    dNe[I] = 0.25 * ET[I] * (1 + XI[I] * XI[J]);
    j11 += dNx[I] * g.xl[I][0]; j12 += dNx[I] * g.xl[I][1];
    j21 += dNe[I] * g.xl[I][0]; j22 += dNe[I] * g.xl[I][1];
  }
  const double det = j11 * j22 - j12 * j21;
  for (int I = 0; I < 4; ++I) {
    a[I] = (j22 * dNx[I] - j12 * dNe[I]) / det;
    b[I] = (-j21 * dNx[I] + j11 * dNe[I]) / det;
  }
}

void iso_plane_stress_D(double Eps, double nu, double D[3][3]) {  // Eps = E/(1-nu^2)
  D[0][0] = D[1][1] = Eps; D[0][1] = D[1][0] = nu * Eps;
  D[0][2] = D[1][2] = D[2][0] = D[2][1] = 0;
  D[2][2] = Eps * (1 - nu) / 2;
}

struct Pt {      // section point state
  double s[3];   // corotational Cauchy stress 11, 22, 12
  double e33;    // integrated thickness strain (hypoelastic) / ln(lambda3) (hyperelastic)
  double ep;
  double lam3;
  double jv;     // hyperelastic: volume ratio J of the last solution (predictor of the thickness stretch)
};

// Hourglass stiffness of a quadrilateral from the incompatible-mode element (QM6: Wilson modes 1 - xi^2,
// 1 - eta^2 with Taylor's correction, so that the patch test is passed), statically condensed and projected
// on the two hourglass displacement modes (u_x and u_y following h = (1, -1, 1, -1) without its linear
// part). Plane stress, unit Young's modulus and thickness; k = {k_xx, k_yy, k_xy}. For a rectangle this is
// the exact energy of pure in-plane bending of an element.
void qm6_hourglass(const double xl[4][3], double nu, double k[3]) {
  const double c = 1.0 / (1 - nu * nu);
  const double D[3][3] = {{c, c * nu, 0}, {c * nu, c, 0}, {0, 0, c * (1 - nu) / 2}};
  const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
  auto jac = [&](double xi, double eta, double J[2][2], double dNx[4], double dNy[4]) {
    double dxi[4], deta[4];
    for (int I = 0; I < 4; ++I) { dxi[I] = 0.25 * XI[I] * (1 + eta * ET[I]); deta[I] = 0.25 * ET[I] * (1 + xi * XI[I]); }
    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0;
    for (int I = 0; I < 4; ++I) {
      J[0][0] += dxi[I] * xl[I][0]; J[0][1] += dxi[I] * xl[I][1];
      J[1][0] += deta[I] * xl[I][0]; J[1][1] += deta[I] * xl[I][1];
    }
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    for (int I = 0; I < 4; ++I) {
      dNx[I] = (J[1][1] * dxi[I] - J[0][1] * deta[I]) / det;
      dNy[I] = (-J[1][0] * dxi[I] + J[0][0] * deta[I]) / det;
    }
    return det;
  };
  double J0[2][2], n0x[4], n0y[4];
  const double det0 = jac(0, 0, J0, n0x, n0y);
  double K[12][12] = {};
  const double g = 1.0 / std::sqrt(3.0);
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b) {
      double xi = a ? g : -g, eta = b ? g : -g, J[2][2], dNx[4], dNy[4];
      double det = jac(xi, eta, J, dNx, dNy);
      // incompatible modes: derivatives with the centre Jacobian, scaled by det0 / det (Taylor)
      double p1x = (J0[1][1] * (-2 * xi)) / det0 * det0 / det, p1y = (-J0[1][0] * (-2 * xi)) / det0 * det0 / det;
      double p2x = (-J0[0][1] * (-2 * eta)) / det0 * det0 / det, p2y = (J0[0][0] * (-2 * eta)) / det0 * det0 / det;
      double B[3][12] = {};
      for (int I = 0; I < 4; ++I) {
        B[0][2 * I] = dNx[I]; B[1][2 * I + 1] = dNy[I]; B[2][2 * I] = dNy[I]; B[2][2 * I + 1] = dNx[I];
      }
      B[0][8] = p1x; B[2][8] = p1y;   // u = a1 (1 - xi^2)
      B[0][9] = p2x; B[2][9] = p2y;   // u = a2 (1 - eta^2)
      B[1][10] = p1y; B[2][10] = p1x; // v = a3 (1 - xi^2)
      B[1][11] = p2y; B[2][11] = p2x; // v = a4 (1 - eta^2)
      for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j) {
          double s = 0;
          for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q) s += B[p][i] * D[p][q] * B[q][j];
          K[i][j] += s * det;
        }
    }
  // static condensation of the 4 internal modes
  double Kii[4][4], inv[4][4];
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) Kii[i][j] = K[8 + i][8 + j];
  {  // Gauss-Jordan inverse
    double M[4][8];
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 8; ++j) M[i][j] = j < 4 ? Kii[i][j] : (j - 4 == i);
    for (int c0 = 0; c0 < 4; ++c0) {
      int piv = c0;
      for (int r = c0 + 1; r < 4; ++r) if (std::fabs(M[r][c0]) > std::fabs(M[piv][c0])) piv = r;
      for (int j = 0; j < 8; ++j) std::swap(M[c0][j], M[piv][j]);
      double d = M[c0][c0];
      for (int j = 0; j < 8; ++j) M[c0][j] /= d;
      for (int r = 0; r < 4; ++r)
        if (r != c0) { double f = M[r][c0]; for (int j = 0; j < 8; ++j) M[r][j] -= f * M[c0][j]; }
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) inv[i][j] = M[i][4 + j];
  }
  double Kc[8][8];
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) {
      double s = K[i][j];
      for (int p = 0; p < 4; ++p)
        for (int q = 0; q < 4; ++q) s -= K[i][8 + p] * inv[p][q] * K[8 + q][j];
      Kc[i][j] = s;
    }
  // hourglass displacement mode: h without its linear part, normalised so that gamma . d = 1
  const double h[4] = {1, -1, 1, -1};
  double bh1 = 0, bh2 = 0;
  for (int I = 0; I < 4; ++I) { bh1 += n0x[I] * h[I]; bh2 += n0y[I] * h[I]; }
  double d[4], hx = 0, hy = 0, B1[4], B2[4];
  for (int I = 0; I < 4; ++I) { d[I] = h[I] - bh1 * xl[I][0] - bh2 * xl[I][1]; B1[I] = n0x[I]; B2[I] = n0y[I]; }
  for (int I = 0; I < 4; ++I) { hx += h[I] * xl[I][0]; hy += h[I] * xl[I][1]; }
  double gd = 0;
  for (int I = 0; I < 4; ++I) gd += 0.25 * (h[I] - hx * B1[I] - hy * B2[I]) * d[I];
  for (int I = 0; I < 4; ++I) d[I] /= gd;
  double kxx = 0, kyy = 0, kxy = 0;
  for (int I = 0; I < 4; ++I)
    for (int J = 0; J < 4; ++J) {
      kxx += d[I] * Kc[2 * I][2 * J] * d[J];
      kyy += d[I] * Kc[2 * I + 1][2 * J + 1] * d[J];
      kxy += d[I] * Kc[2 * I][2 * J + 1] * d[J];
    }
  k[0] = kxx; k[1] = kyy; k[2] = kxy;
}

// Transverse shear of the quadrilateral sampled along its edges (assumed natural strain): the tangential shear
// of an edge, w_J - w_I + (mean rotation of the end nodes) x (edge vector), vanishes for every field with linear
// rotations and quadratic deflection, so that constant-curvature bending of a distorted element is free of
// transverse shear. Edges 0->1 and 3->2 carry the xi strain, 0->3 and 1->2 the eta strain; the centre values give
// the uniform shear (Cartesian, through the centre Jacobian), their differences the linear variations
// (stabilisation). G[c][I][k]: derivative of component c (0: gamma_xz, 1: gamma_yz at the centre, 2: variation
// of the xi strain along eta, 3: variation of the eta strain along xi) with respect to dof k of node I
// (0: w, 1: theta_x, 2: theta_y). inv[4]: inverse centre Jacobian (gamma_x = inv0 e_xi + inv1 e_eta,
// gamma_y = inv2 e_xi + inv3 e_eta).
void edge_shear_ops(const double xl[4][3], double G[4][4][3], double inv[4]) {
  static const int ED[4][2] = {{0, 1}, {3, 2}, {0, 3}, {1, 2}};
  static const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
  double xx = 0, xe = 0, yx = 0, ye = 0;
  for (int I = 0; I < 4; ++I) {
    xx += 0.25 * XI[I] * xl[I][0]; xe += 0.25 * ET[I] * xl[I][0];
    yx += 0.25 * XI[I] * xl[I][1]; ye += 0.25 * ET[I] * xl[I][1];
  }
  const double det = xx * ye - yx * xe;
  inv[0] = ye / det; inv[1] = -yx / det; inv[2] = -xe / det; inv[3] = xx / det;
  double Ed[4][4][3] = {};  // derivative of each edge value
  for (int k = 0; k < 4; ++k) {
    const int I = ED[k][0], J = ED[k][1];
    const double dx = xl[J][0] - xl[I][0], dy = xl[J][1] - xl[I][1];
    Ed[k][J][0] += 1; Ed[k][I][0] -= 1;
    Ed[k][I][1] -= 0.5 * dy; Ed[k][J][1] -= 0.5 * dy;
    Ed[k][I][2] += 0.5 * dx; Ed[k][J][2] += 0.5 * dx;
  }
  for (int I = 0; I < 4; ++I)
    for (int q = 0; q < 3; ++q) {
      const double exi = 0.25 * (Ed[0][I][q] + Ed[1][I][q]), eeta = 0.25 * (Ed[2][I][q] + Ed[3][I][q]);
      G[0][I][q] = inv[0] * exi + inv[1] * eeta;
      G[1][I][q] = inv[2] * exi + inv[3] * eeta;
      G[2][I][q] = 0.25 * (Ed[1][I][q] - Ed[0][I][q]);
      G[3][I][q] = 0.25 * (Ed[3][I][q] - Ed[2][I][q]);
    }
}

// transverse shear of the reduced-integration quad, from the covariant shear strains of the four edges
// (e = w_K - w_I + 1/2 (t_x (theta_y,I + theta_y,K) - t_y (theta_x,I + theta_x,K)), t = x_K - x_I; zero for
// every Kirchhoff state). G[0], G[1]: shear strains gamma_xz, gamma_yz at the centre from the mean edge strains
// u_xi = (e_bottom + e_top) / 2, u_eta = (e_left + e_right) / 2 (gamma = J0^-1 (u_xi, u_eta) / 2, rows of J0 =
// dx/dxi, dx/deta). G[2]: one stabilised mode, half the sum of the edge-strain variations v_xi = e_top -
// e_bottom, v_eta = e_right - e_left, less its part along the centre gradient (projection of the bilinear mode,
// h.x J0^-1 (u_xi, u_eta) / 2): for w alone it is sum gamma_I w_I with the projected hourglass vector gamma.
// G[3] unused. Returns gamma.gamma (4 for a parallelogram).
double quad_shear_ops(const Geo& g, double G[4][4][3]) {
  static const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
  static const int EDGE[4][2] = {{0, 1}, {3, 2}, {0, 3}, {1, 2}};  // bottom, top, left, right
  double E[4][4][3];  // edge strain operators
  for (int k = 0; k < 4; ++k) {
    for (int I = 0; I < 4; ++I) E[k][I][0] = E[k][I][1] = E[k][I][2] = 0.0;
    const int a = EDGE[k][0], b = EDGE[k][1];
    const double tx = g.xl[b][0] - g.xl[a][0], ty = g.xl[b][1] - g.xl[a][1];
    E[k][b][0] += 1.0; E[k][a][0] -= 1.0;
    for (int N : {a, b}) { E[k][N][1] += -0.5 * ty; E[k][N][2] += 0.5 * tx; }
  }
  double jxi[2] = {0, 0}, jet[2] = {0, 0}, hx = 0, hy = 0;
  for (int I = 0; I < 4; ++I) {
    for (int c = 0; c < 2; ++c) { jxi[c] += 0.25 * XI[I] * g.xl[I][c]; jet[c] += 0.25 * ET[I] * g.xl[I][c]; }
    hx += XI[I] * ET[I] * g.xl[I][0]; hy += XI[I] * ET[I] * g.xl[I][1];
  }
  const double det = jxi[0] * jet[1] - jxi[1] * jet[0];
  // J0^-1 with J0 = [jxi; jet]: inverse = [jet1 -jxi1; -jet0 jxi0] / det
  const double i00 = jet[1] / det, i01 = -jxi[1] / det, i10 = -jet[0] / det, i11 = jxi[0] / det;
  // projection: q -= (hx, hy) . J0^-1 (u_xi, u_eta) / 2
  const double cxi = 0.5 * (hx * i00 + hy * i10), ceta = 0.5 * (hx * i01 + hy * i11);
  double gg = 0;
  for (int I = 0; I < 4; ++I) {
    for (int q = 0; q < 3; ++q) {
      const double uxi = 0.5 * (E[0][I][q] + E[1][I][q]), ueta = 0.5 * (E[2][I][q] + E[3][I][q]);
      const double vxi = E[1][I][q] - E[0][I][q], veta = E[3][I][q] - E[2][I][q];
      G[0][I][q] = 0.5 * (i00 * uxi + i01 * ueta);
      G[1][I][q] = 0.5 * (i10 * uxi + i11 * ueta);
      G[2][I][q] = 0.5 * (vxi + veta) - cxi * uxi - ceta * ueta;
      G[3][I][q] = 0.0;
    }
    gg += G[2][I][0] * G[2][I][0];
  }
  return gg;
}

// S3RS transverse shear from the three edge strains: the shear field a + c (-y, x) (coordinates from the centroid)
// whose tangential component at the middle of each edge reproduces the edge strain; the element keeps the uniform
// part a only (one point at the centroid, no stiffness on the circulation c). G[c][I][k]: derivative of a_x (c = 0),
// a_y (c = 1) and c (c = 2) with respect to dof k of node I (0: w, 1: theta_x, 2: theta_y); G[3] = 0. Reproduces
// the single-element Abaqus/Explicit S3RS stiffness to 1e-7 (5 shapes, thickness 0.01 to 1)
void tri_edge_shear_ops(const Geo& g, double G[4][4][3]) {
  double M[3][3], R[3][3][3];  // rows: edges; R[edge][node][dof]
  for (int k = 0; k < 3; ++k) {
    const int a = k, b = (k + 1) % 3;
    const double tx = g.xl[b][0] - g.xl[a][0], ty = g.xl[b][1] - g.xl[a][1];
    const double mx = 0.5 * (g.xl[a][0] + g.xl[b][0]), my = 0.5 * (g.xl[a][1] + g.xl[b][1]);
    M[k][0] = tx; M[k][1] = ty; M[k][2] = mx * ty - my * tx;
    for (int I = 0; I < 3; ++I) R[k][I][0] = R[k][I][1] = R[k][I][2] = 0.0;
    R[k][b][0] += 1.0; R[k][a][0] -= 1.0;
    for (int N : {a, b}) { R[k][N][1] += -0.5 * ty; R[k][N][2] += 0.5 * tx; }
  }
  // inverse of M by cofactors
  const double c00 = M[1][1] * M[2][2] - M[1][2] * M[2][1], c01 = M[1][2] * M[2][0] - M[1][0] * M[2][2];
  const double c02 = M[1][0] * M[2][1] - M[1][1] * M[2][0];
  const double det = M[0][0] * c00 + M[0][1] * c01 + M[0][2] * c02;
  double Mi[3][3];
  Mi[0][0] = c00 / det; Mi[1][0] = c01 / det; Mi[2][0] = c02 / det;
  Mi[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / det;
  Mi[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det;
  Mi[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) / det;
  Mi[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
  Mi[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / det;
  Mi[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;
  for (int c = 0; c < 3; ++c)
    for (int I = 0; I < 3; ++I)
      for (int q = 0; q < 3; ++q) G[c][I][q] = Mi[c][0] * R[0][I][q] + Mi[c][1] * R[1][I][q] + Mi[c][2] * R[2][I][q];
  for (int I = 0; I < 3; ++I) G[3][I][0] = G[3][I][1] = G[3][I][2] = 0.0;
  for (int q = 0; q < 4; ++q) for (int a = 0; a < 3; ++a) G[q][3][a] = 0.0;
}

// transverse shear of the triangle as a degenerate quad: corners (k+1, k+2, k, k), the degenerate corner k being
// the apex of the most isosceles view of the triangle (the other two angles closest). Homogeneous shear from the
// covariant edge strains (no shear for constant curvature), with the circulation of the edge strains ("crop
// circle" mode) added with the coefficients c (cxi, ceta), plus a weak stiffness on that mode (0.001 of the
// homogeneous one, weight returned in cw: K_cc = kts cw). G[0], G[1]: gamma_xz, gamma_yz; G[2]: circulation.
// Together with the thin-shell factor this gives the Abaqus S3R shear stiffness (exact for isosceles
// triangles, within a few % otherwise; single-element probes)
void tri_degenerate_shear_ops(const Geo& g, int k, const double c[2], double G[4][4][3], double& cw) {
  static const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
  const int o[4] = {(k + 1) % 3, (k + 2) % 3, k, k};
  auto edge = [&](int a, int b, double e[3][3]) {  // e[node][q] for the edge a -> b of the triangle
    for (int I = 0; I < 3; ++I) e[I][0] = e[I][1] = e[I][2] = 0.0;
    if (a == b) return;
    const double tx = g.xl[b][0] - g.xl[a][0], ty = g.xl[b][1] - g.xl[a][1];
    e[b][0] += 1.0; e[a][0] -= 1.0;
    for (int N : {a, b}) { e[N][1] += -0.5 * ty; e[N][2] += 0.5 * tx; }
  };
  double eA[3][3], eB[3][3], eC[3][3];
  edge(o[0], o[3], eA); edge(o[0], o[1], eB); edge(o[1], o[2], eC);  // the top edge (o[3] -> o[2]) is degenerate
  double jxi[2] = {0, 0}, jet[2] = {0, 0};
  for (int I = 0; I < 4; ++I)
    for (int d = 0; d < 2; ++d) { jxi[d] += 0.25 * XI[I] * g.xl[o[I]][d]; jet[d] += 0.25 * ET[I] * g.xl[o[I]][d]; }
  const double det = jxi[0] * jet[1] - jxi[1] * jet[0];
  const double i00 = jet[1] / det, i01 = -jxi[1] / det, i10 = -jet[0] / det, i11 = jxi[0] / det;
  for (int I = 0; I < 3; ++I)
    for (int q = 0; q < 3; ++q) {
      const double gcc = 0.5 * (-eA[I][q] + eB[I][q] + eC[I][q]);
      const double u1 = 0.25 * eB[I][q] + c[0] * gcc, u2 = 0.25 * (eA[I][q] + eC[I][q]) + c[1] * gcc;
      G[0][I][q] = i00 * u1 + i01 * u2;
      G[1][I][q] = i10 * u1 + i11 * u2;
      G[2][I][q] = gcc;
      G[3][I][q] = 0.0;
    }
  for (int q = 0; q < 4; ++q) for (int a = 0; a < 3; ++a) G[q][3][a] = 0.0;
  // weak stiffness of the circulation: 0.001 Cs_eff / 12 with Cs_eff the mean of the inverse metric diagonal
  const double m00 = jxi[0] * jxi[0] + jxi[1] * jxi[1], m11 = jet[0] * jet[0] + jet[1] * jet[1];
  const double m01 = jxi[0] * jet[0] + jxi[1] * jet[1], dm = m00 * m11 - m01 * m01;
  cw = 4.0 * det * 0.001 / 12.0 * 0.5 * (m11 + m00) / dm;
}

// apex and crop-circle coefficients of the degenerate view of a triangle (reference geometry): apex k where the
// other two angles are closest; c from the skew s and aspect q of the degenerate quad (Jacobian rows j_xi, j_eta:
// s = j_xi.j_eta / |j_xi|^2, q = |j_xi x j_eta| / |j_xi|^2), fitted on single-element Abaqus probes of 40
// triangles: ceta = 1.63 s / (2.963 - q^2 + 1.413 s^2) (bounded by 0.3), cxi = -1.44 ceta^2; zero for isosceles
void tri_collapse_setup(const Geo& g, int& kout, double c[2]) {
  double ang[3];
  for (int k = 0; k < 3; ++k) {
    const double ax = g.xl[(k + 1) % 3][0] - g.xl[k][0], ay = g.xl[(k + 1) % 3][1] - g.xl[k][1];
    const double bx = g.xl[(k + 2) % 3][0] - g.xl[k][0], by = g.xl[(k + 2) % 3][1] - g.xl[k][1];
    ang[k] = std::atan2(std::fabs(ax * by - ay * bx), ax * bx + ay * by);
  }
  int k = 0;
  double best = 1e300;
  for (int j = 0; j < 3; ++j) {
    const double d = std::fabs(ang[(j + 1) % 3] - ang[(j + 2) % 3]);
    if (d < best - 1e-12) { best = d; k = j; }
  }
  static const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
  const int o[4] = {(k + 1) % 3, (k + 2) % 3, k, k};
  double jxi[2] = {0, 0}, jet[2] = {0, 0};
  for (int I = 0; I < 4; ++I)
    for (int d = 0; d < 2; ++d) { jxi[d] += 0.25 * XI[I] * g.xl[o[I]][d]; jet[d] += 0.25 * ET[I] * g.xl[o[I]][d]; }
  const double n2 = jxi[0] * jxi[0] + jxi[1] * jxi[1];
  const double sk = (jxi[0] * jet[0] + jxi[1] * jet[1]) / n2, qa = std::fabs(jxi[0] * jet[1] - jxi[1] * jet[0]) / n2;
  double ce = 1.63 * sk / (2.963 - qa * qa + 1.413 * sk * sk);
  if (!(std::fabs(ce) <= 0.3)) ce = sk >= 0 ? 0.3 : -0.3;
  c[0] = -1.44 * ce * ce;
  c[1] = ce;
  kout = k;
}


struct ElState {
  Pt pt[MAXPT];
  double Q[2];                 // transverse shear resultants
  double hgm[2], hgb[2], hgs[2];  // hourglass generalized forces (membrane, rotations, transverse shear)
  double hgmv[2], hgbv[2];        // relaxing parts of the membrane and rotation hourglass forces (explicit)
  double md[4];                // drilling moments
  double em[3];                // integrated membrane strain (e11, e22, 2e12)
  double kap[3];               // integrated curvature (k11, k22, 2k12)
  double Meff;                 // effective in-plane modulus (stable increment, hourglass)
  int age;                     // increments since the effective modulus update
  double Cref[2][3];           // C of the two outer section points at the last update
  double Qref[2];              // Fung exponent of the two outer points at the last update
};

class ShellGroup : public ElementGroup {
 public:
  std::vector<ElState> st;
  std::vector<const Section*> sec;
  std::vector<const Material*> mat;
  std::vector<double> ref;   // per element: dN/dX1 (4), dN/dX2 (4) in the reference material basis, A0
  std::vector<double> khg;   // quads: hourglass stiffness (unit modulus and thickness) xx, yy, xy
  std::vector<double> irot;  // rotary inertia per element node (unscaled)
  std::vector<double> mfr;   // share of the element mass carried by each node
  std::vector<double> oax;   // per element: output axis 1 in the initial element frame (cos, sin)

  static void simpson(int n, double* z, double* w) {
    if (n == 1) { z[0] = 0; w[0] = 1; return; }
    for (int k = 0; k < n; ++k) {
      z[k] = -1 + 2.0 * k / (n - 1);
      w[k] = (k == 0 || k == n - 1) ? 1 : (k % 2 ? 4 : 2);
    }
    double s = 0;
    for (int k = 0; k < n; ++k) s += w[k];
    for (int k = 0; k < n; ++k) w[k] /= s;
  }

  bool hyper(int e) const { return mat[e]->fung || mat[e]->neohooke; }
  double nu_eff(int e) const { return mat[e]->elastic ? mat[e]->nu : (hyper(e) ? 0.49 : 0.3); }

  std::vector<double> ts_factor;  // transverse shear stiffness factor (thin elements)
  std::vector<int> tri_k;         // triangles: degenerate corner of the shear formulation
  std::vector<double> tri_c;      // triangles: crop-circle coefficients (2 per element)
  double ts_stiff(int e, int dir) const {
    const Section& s = *sec[e];
    double k = s.has_ts ? (dir == 0 ? s.ts_k11 : s.ts_k22) : 5.0 / 6.0 * mat[e]->shear_modulus() * s.thickness;
    return ts_factor.empty() ? k : k * ts_factor[e];
  }

  void setup(const Model& m) override {
    model = &m;
    rot = true;
    sec.resize(size());
    mat.resize(size());
    st.assign(size(), ElState());
    ref.assign(size() * 9, 0.0);
    irot.assign(size() * nen, 0.0);
    mfr.assign(size() * nen, 1.0 / nen);
    oax.assign(size() * 2, 0.0);
    if (nen == 4) khg.assign(size() * 3, 0.0);
    if (nen == 4) {  // S4R: transverse shear stiffness kappa G t / (1 + 0.25e-4 A / t^2) (default stiffness)
      ts_factor.assign(size(), 1.0);
      for (int e = 0; e < size(); ++e) {
        const Section& s = m.sections[m.elems[ids[e]].section];
        if (type != ElType::S4R || s.has_ts) continue;
        const Vec3 p0 = m.X[conn[e * 4]], p1 = m.X[conn[e * 4 + 1]], p2 = m.X[conn[e * 4 + 2]], p3 = m.X[conn[e * 4 + 3]];
        const double area = 0.5 * norm(cross(p2 - p0, p3 - p1));
        ts_factor[e] = 1.0 / (1.0 + 0.25e-4 * area / (s.thickness * s.thickness));
      }
    }
    if (nen == 3) {  // S3R: thin-shell factor as S4R, degenerate corner and crop-circle coefficients; S3RS factor
      ts_factor.assign(size(), 1.0);
      tri_k.assign(size(), 0);
      tri_c.assign(size() * 2, 0.0);
      for (int e = 0; e < size(); ++e) {
        const Section& s = m.sections[m.elems[ids[e]].section];
        double x[12] = {};
        for (int I = 0; I < 3; ++I) store3(x + 3 * I, m.X[conn[e * 3 + I]]);
        Geo g0;
        geometry(x, 3, g0);
        tri_collapse_setup(g0, tri_k[e], &tri_c[e * 2]);
        if (s.has_ts) continue;
        const double t = s.thickness, fp = 1.0 / (1.0 + 0.25e-4 * g0.A / (t * t));
        if (type != ElType::S3RS) { ts_factor[e] = fp; continue; }
        // S3RS: full shear stiffness on thick elements (t > L / 5, L the shortest edge; no thin-shell factor); on
        // thin ones the shear acts as a penalty scaled by t / (2.4 L) and by the thin-shell factor, independent of
        // the elastic constants (single-element probes: 6 shapes, thickness swept across the switch, three moduli
        // and Poisson ratios)
        double L2 = 1e300;
        for (int I = 0; I < 3; ++I) {
          const double dx = g0.xl[(I + 1) % 3][0] - g0.xl[I][0], dy = g0.xl[(I + 1) % 3][1] - g0.xl[I][1];
          L2 = std::min(L2, dx * dx + dy * dy);
        }
        const double L = std::sqrt(L2);
        ts_factor[e] = 5.0 * t > L ? 1.0 : fp * t / (2.4 * L);
      }
    }
    for (int e = 0; e < size(); ++e) {
      const Section& s = m.sections[m.elems[ids[e]].section];
      sec[e] = &s;
      mat[e] = &m.materials[s.material];
      if (s.nip > MAXPT) throw std::runtime_error("shell section: at most 15 integration points");
      if (s.nip % 2 == 0 && s.nip > 1) throw std::runtime_error("shell section: Simpson integration needs an odd number of points");
      if (s.thickness <= 0) throw std::runtime_error("shell section: thickness must be positive");
      ElState& S = st[e];
      std::memset(&S, 0, sizeof(S));
      for (int k = 0; k < MAXPT; ++k) S.pt[k].lam3 = 1.0;
      S.Meff = plane_stress_modulus(*mat[e]);
      double X[12];
      for (int I = 0; I < nen; ++I) store3(X + 3 * I, m.X[conn[e * nen + I]]);
      Geo g;
      geometry(X, nen, g);
      if (g.A <= 0) throw std::runtime_error("element " + m.elem_name(ids[e]) + ": zero or negative area");
      if (nen == 4 && (s.hg_enhanced || type == ElType::S4)) qm6_hourglass(g.xl, nu_eff(e), &khg[e * 3]);
      Vec3 m1 = g.e1;
      if (s.orientation >= 0) {
        Vec3 a1 = m.orients[s.orientation].R.row(0);
        Vec3 p = a1 - g.e3 * dot(a1, g.e3);
        if (norm(p) < 1e-6 * norm(a1)) {  // axis 1 normal to the shell: use axis 2 rotated
          Vec3 a2 = m.orients[s.orientation].R.row(1);
          p = cross(a2 - g.e3 * dot(a2, g.e3), g.e3);
        }
        m1 = normalized(p);
      }
      double c = dot(g.e1, m1), sn = dot(g.e2, m1);  // m1 = c e1 + sn e2
      // stress output axes: the orientation if any, otherwise the projection of the global x axis on the surface
      // (of the global z axis when x is within 0.1 degree of the normal); they follow the element rotation
      Vec3 o1 = m1;
      if (s.orientation < 0) {
        o1 = Vec3(1, 0, 0) - g.e3 * g.e3.x;
        if (norm(o1) < std::sin(0.1 * M_PI / 180.0)) o1 = Vec3(0, 0, 1) - g.e3 * g.e3.z;
        o1 = normalized(o1);
      }
      oax[e * 2] = dot(g.e1, o1); oax[e * 2 + 1] = dot(g.e2, o1);
      for (int I = 0; I < nen; ++I) {
        ref[e * 9 + I] = c * g.B1[I] + sn * g.B2[I];
        ref[e * 9 + 4 + I] = -sn * g.B1[I] + c * g.B2[I];
      }
      ref[e * 9 + 8] = g.A;
      // rotary inertia of a node (the same about the three axes): its translational mass times
      // A/12 + ((L^2/16)^4 + (t^2/6)^4)^(1/4), L the shortest edge of a triangle, the area over the longest median
      // of a quadrilateral (the distance between opposite edges of a parallelogram); fitted on single-element
      // Abaqus/Explicit runs (initial angular acceleration under nodal moments, 11 shapes, thickness 0.01 to 8,
      // several moduli and Poisson ratios: 1e-7). Rotations then do not control the stable increment
      const double t = s.thickness, mtr = mat[e]->density * t * g.A / nen;
      double L2;
      if (nen == 3) {
        L2 = 1e300;
        for (int I = 0; I < 3; ++I) {
          const double dx = g.xl[(I + 1) % 3][0] - g.xl[I][0], dy = g.xl[(I + 1) % 3][1] - g.xl[I][1];
          L2 = std::min(L2, dx * dx + dy * dy);
        }
      } else {
        double med2 = 0;
        for (int I = 0; I < 2; ++I) {  // median joining the middles of edges I-(I+1) and (I+2)-(I+3)
          const double dx = 0.5 * (g.xl[(I + 2) % 4][0] + g.xl[(I + 3) % 4][0] - g.xl[I][0] - g.xl[I + 1][0]);
          const double dy = 0.5 * (g.xl[(I + 2) % 4][1] + g.xl[(I + 3) % 4][1] - g.xl[I][1] - g.xl[I + 1][1]);
          med2 = std::max(med2, dx * dx + dy * dy);
        }
        L2 = g.A * g.A / std::max(med2, 1e-300);
      }
      const double a4 = L2 * L2 / 256.0, b4 = t * t * t * t / 36.0;
      // S3RS: nodal mass times max(A / 4, t^2 / 6) (single-element Abaqus/Explicit runs, 5 shapes, 3 thicknesses: 1e-5)
      const double r2 = type == ElType::S3RS ? std::max(g.A / 4.0, t * t / 6.0)
                                             : g.A / 12.0 + std::pow(a4 * a4 + b4 * b4, 0.25);
      if (nen == 4) {
        // nodal masses of a quadrilateral: integral of the shape function over the element (row sums of the
        // consistent mass matrix, as Abaqus lumps them: 0.278 / 0.222 of the mass on a trapezoid of sides 3 and
        // 1.5); dA = (j0 + j1 xi + j2 eta) dxi deta
        static const double XI[4] = {-1, 1, 1, -1}, ET[4] = {-1, -1, 1, 1};
        double ax = 0, bx = 0, cx = 0, ay = 0, by = 0, cy = 0;
        for (int I = 0; I < 4; ++I) {
          ax += 0.25 * XI[I] * g.xl[I][0]; bx += 0.25 * ET[I] * g.xl[I][0]; cx += 0.25 * XI[I] * ET[I] * g.xl[I][0];
          ay += 0.25 * XI[I] * g.xl[I][1]; by += 0.25 * ET[I] * g.xl[I][1]; cy += 0.25 * XI[I] * ET[I] * g.xl[I][1];
        }
        const double j0 = ax * by - bx * ay, j1 = ax * cy - cx * ay, j2 = cx * by - bx * cy;
        if (std::fabs(j0) > 0)
          for (int I = 0; I < 4; ++I) mfr[e * 4 + I] = (j0 + (j1 * XI[I] + j2 * ET[I]) / 3.0) / (4.0 * j0);
      }
      for (int I = 0; I < nen; ++I) irot[e * nen + I] = mtr * nen * mfr[e * nen + I] * r2;
    }
  }

  void save(BinOut& o) const override { ElementGroup::save(o); o.vec(st); o.vec(irot); o.vec(mfr); }
  void load(BinIn& in) override { ElementGroup::load(in); in.vec(st); in.vec(irot); in.vec(mfr); }

  double alpha_damping(int e) const override { return mat[e]->damp_alpha; }
  double contact_thickness(int e) const override { return sec[e]->thickness; }

  void mass(int e, const double* X, double* mt, double* mr) const override {
    double A = ref[e * 9 + 8];
    double m = mat[e]->density * sec[e]->thickness * A;
    for (int I = 0; I < nen; ++I) { mt[I] = m * mfr[e * nen + I]; mr[I] = irot[e * nen + I]; }
  }

  // largest eigenvalue of the plane-stress hyperelastic tangent dS/dE at C (finite differences)
  static double hyper_tangent_max(const Material& m, const double C[3], double lam3) {
    double S0[3], l = lam3;
    hyper_ps(m, C, l, S0);
    double D[3][3];
    for (int j = 0; j < 3; ++j) {
      double h = 1e-6;
      double Cp[3] = {C[0], C[1], C[2]};
      if (j < 2) Cp[j] += 2 * h; else Cp[2] += h;
      double S1[3], lp = l;
      hyper_ps(m, Cp, lp, S1);
      for (int i = 0; i < 3; ++i) D[i][j] = (S1[i] - S0[i]) / h;
    }
    return max_eig_sym3(D[0][0], D[1][1], D[2][2], 0.5 * (D[0][1] + D[1][0]), 0.5 * (D[0][2] + D[2][0]), 0.5 * (D[1][2] + D[2][1]));
  }

  void compute(int e, const ElemIn& in, ElemOut& out, bool commit) override {
    const int n = nen;
    const Section& sc = *sec[e];
    const Material& m = *mat[e];
    ElState S = st[e];
    const double t0 = sc.thickness;
    const double rho = m.density * msf[e];
    double xm[12], x1[12];
    for (int I = 0; I < n; ++I)
      for (int k = 0; k < 3; ++k) {
        double base = in.nlgeom ? in.xo[3 * I + k] : in.X[3 * I + k];
        xm[3 * I + k] = base + (in.nlgeom ? 0.5 * in.du[6 * I + k] : 0.0);
        x1[3 * I + k] = base + (in.nlgeom ? in.du[6 * I + k] : 0.0);
      }
    Geo gm, g1;
    geometry(xm, n, gm);
    geometry(x1, n, g1);
    double du[4][3], dth[4][3];
    for (int I = 0; I < n; ++I) {
      Vec3 u = load3(in.du + 6 * I), r = load3(in.du + 6 * I + 3);
      du[I][0] = dot(u, gm.e1); du[I][1] = dot(u, gm.e2); du[I][2] = dot(u, gm.e3);
      dth[I][0] = dot(r, gm.e1); dth[I][1] = dot(r, gm.e2); dth[I][2] = dot(r, gm.e3);
    }
    double de[3] = {0, 0, 0}, dk[3] = {0, 0, 0}, dg[2] = {0, 0};
    double thx = 0, thy = 0, wz = 0;
    for (int I = 0; I < n; ++I) {
      de[0] += gm.B1[I] * du[I][0];
      de[1] += gm.B2[I] * du[I][1];
      de[2] += gm.B2[I] * du[I][0] + gm.B1[I] * du[I][1];
      dk[0] += gm.B1[I] * dth[I][1];
      dk[1] += -gm.B2[I] * dth[I][0];
      dk[2] += gm.B2[I] * dth[I][1] - gm.B1[I] * dth[I][0];
      dg[0] += gm.B1[I] * du[I][2];
      dg[1] += gm.B2[I] * du[I][2];
      thx += dth[I][0]; thy += dth[I][1];
      wz += 0.5 * (gm.B1[I] * du[I][1] - gm.B2[I] * du[I][0]);
    }
    dg[0] += thy / n;
    dg[1] += -thx / n;
    // uniform transverse shear from the edge strains (edge_shear_ops, tri_shear_ops)
    double Gm[4][4][3], invm[4] = {0, 0, 0, 0}, cwm = 0;
    double ggm = 4.0;
    if (n == 4 && type == ElType::S4) edge_shear_ops(gm.xl, Gm, invm);
    else if (n == 4) ggm = quad_shear_ops(gm, Gm);
    else if (type == ElType::S3RS) tri_edge_shear_ops(gm, Gm);
    else tri_degenerate_shear_ops(gm, tri_k[e], &tri_c[e * 2], Gm, cwm);
    dg[0] = dg[1] = 0;
    for (int I = 0; I < n; ++I)
      for (int c = 0; c < 2; ++c) dg[c] += Gm[c][I][0] * du[I][2] + Gm[c][I][1] * dth[I][0] + Gm[c][I][2] * dth[I][1];
    for (int i = 0; i < 3; ++i) { S.em[i] += de[i]; S.kap[i] += dk[i]; }

    double dtemp = 0;
    if (in.temp0 && in.temp1) { for (int I = 0; I < n; ++I) dtemp += in.temp1[I] - in.temp0[I]; dtemp /= n; }

    double zk[MAXPT], wk[MAXPT];
    const int np = sc.nip;
    simpson(np, zk, wk);
    double N[3] = {0, 0, 0}, M[3] = {0, 0, 0};
    double tcur = t0;
    double stretch_fac = 1.0;  // push-forward of the material tangent (hyperelastic)
    double meff_now = -1;      // current effective modulus (hyperelastic), S.Meff otherwise
    double mis = 0, epm = 0, lemax = -1e30;
    const bool hyp = hyper(e);
    if (hyp) {
      const double* B0 = &ref[e * 9];
      double F[2][2] = {{0, 0}, {0, 0}};
      for (int I = 0; I < n; ++I)
        for (int a = 0; a < 2; ++a) { F[a][0] += g1.xl[I][a] * B0[I]; F[a][1] += g1.xl[I][a] * B0[4 + I]; }
      double detF = F[0][0] * F[1][1] - F[0][1] * F[1][0];
      // fibre deformation gradient F(z) = (I + z kappa) F_m (kappa in the current element frame)
      double k11 = S.kap[0], k22 = S.kap[1], k12 = 0.5 * S.kap[2];
      double Cout[2][3], lout[2], Jout[2];  // outer section points (stable increment)
      double c1max = 0, l3max = 0, Qout[2] = {-1, -1};
      for (int k = 0; k < np; ++k) {
        double z = zk[k] * 0.5 * t0;
        double Mz[2][2] = {{1 + z * k11, z * k12}, {z * k12, 1 + z * k22}};
        double Fz[2][2];
        for (int a2 = 0; a2 < 2; ++a2)
          for (int b2 = 0; b2 < 2; ++b2) Fz[a2][b2] = Mz[a2][0] * F[0][b2] + Mz[a2][1] * F[1][b2];
        double detz = Fz[0][0] * Fz[1][1] - Fz[0][1] * Fz[1][0];
        if (detz <= 1e-3 * std::fabs(detF)) {  // fibre folded through itself: keep the mid-surface state
          for (int a2 = 0; a2 < 2; ++a2) for (int b2 = 0; b2 < 2; ++b2) Fz[a2][b2] = F[a2][b2];
          detz = detF;
          out.failed = true;
        }
        double C[3] = {Fz[0][0] * Fz[0][0] + Fz[1][0] * Fz[1][0], Fz[0][1] * Fz[0][1] + Fz[1][1] * Fz[1][1],
                       Fz[0][0] * Fz[0][1] + Fz[1][0] * Fz[1][1]};
        double l3 = S.pt[k].lam3, Sp[3], qe = -1;
        const double a2 = C[0] * C[1] - C[2] * C[2];
        if (S.pt[k].jv > 0 && a2 > 0) l3 = S.pt[k].jv / std::sqrt(a2);  // J nearly constant (incompressibility)
        if (!hyper_ps(m, C, l3, Sp, nullptr, &qe)) out.failed = true;
        if (a2 > 0) S.pt[k].jv = l3 * std::sqrt(a2);
        S.pt[k].lam3 = l3;
        S.pt[k].e33 = std::log(l3);
        if (k == 0 || k == np - 1) {
          int o = k == 0 ? 0 : 1;
          Qout[o] = qe;
          for (int i = 0; i < 3; ++i) Cout[o][i] = C[i];
          lout[o] = l3;
          Jout[o] = std::max(detz * l3, 1e-12);
        }
        double Sm[2][2] = {{Sp[0], Sp[2]}, {Sp[2], Sp[1]}};
        double FS[2][2];
        for (int a2 = 0; a2 < 2; ++a2) for (int b2 = 0; b2 < 2; ++b2) FS[a2][b2] = Fz[a2][0] * Sm[0][b2] + Fz[a2][1] * Sm[1][b2];
        double sg[2][2];
        for (int a2 = 0; a2 < 2; ++a2)
          for (int b2 = 0; b2 < 2; ++b2) sg[a2][b2] = (FS[a2][0] * Fz[b2][0] + FS[a2][1] * Fz[b2][1]) / (detz * l3);
        S.pt[k].s[0] = sg[0][0]; S.pt[k].s[1] = sg[1][1]; S.pt[k].s[2] = sg[0][1];
        double l1, l2, cs, sn;
        eig2sym(C[0], C[2], C[1], l1, l2, cs, sn);
        c1max = std::max(c1max, l1);
        l3max = std::max(l3max, l3);
      }
      lemax = std::max(0.5 * std::log(std::max(c1max, 1e-30)), std::log(l3max));
      tcur = t0 * S.pt[np / 2].lam3;
      // effective modulus for the stable increment: largest pushed-forward plane-stress tangent of the
      // two outer section points (with bending, the exponential stiffening is strongest there; the
      // mid-surface value underestimates the bending frequency by orders of magnitude)
      double dC = 0;
      for (int o = 0; o < 2; ++o)
        for (int i = 0; i < 3; ++i) dC = std::max(dC, std::fabs(Cout[o][i] - S.Cref[o][i]));
      // between finite-difference updates the modulus follows the exponential stiffening exp(Q - Qref)
      double dQ = 0;
      bool qjump = false;
      if (Qout[0] >= 0 && Qout[1] >= 0)
        for (int o = 0; o < 2; ++o) { dQ = std::max(dQ, Qout[o] - S.Qref[o]); qjump |= std::fabs(Qout[o] - S.Qref[o]) > 0.4; }
      if (S.age == 0 || S.age >= 50 || dC > 2e-3 || qjump) {
        double Mx = 1e-30;
        for (int o = 0; o < 2; ++o) {
          double l1, l2, cs, sn;
          eig2sym(Cout[o][0], Cout[o][2], Cout[o][1], l1, l2, cs, sn);
          Mx = std::max(Mx, hyper_tangent_max(m, Cout[o], lout[o]) * std::max(l1, 1e-12) * std::max(l1, 1e-12) / Jout[o]);
          for (int i = 0; i < 3; ++i) S.Cref[o][i] = Cout[o][i];
          S.Qref[o] = Qout[o];
        }
        S.Meff = Mx;
        S.age = 0;
        dQ = 0;
      }
      meff_now = S.Meff * std::exp(std::min(dQ, 50.0));
      ++S.age;
    } else {
      double tc = t0 * std::exp(S.pt[np / 2].e33);
      for (int k = 0; k < np; ++k) {
        Pt& p = S.pt[k];
        double z = zk[k] * 0.5 * tc;
        double d11 = de[0] + z * dk[0] - m.alpha * dtemp;
        double d22 = de[1] + z * dk[1] - m.alpha * dtemp;
        double d12 = 0.5 * (de[2] + z * dk[2]);
        double d33;
        double w = j2_update_ps(m, d11, d22, d12, p.s, p.ep, d33);
        p.e33 += d33 + m.alpha * dtemp;
        out.e_plastic += w * wk[k] * tc * g1.A;
        if (k == 0 || k == np - 1) epm = std::max(epm, p.ep);  // output at the bottom and top section points
        double e11 = S.em[0] + z * S.kap[0], e22 = S.em[1] + z * S.kap[1], e12 = 0.5 * (S.em[2] + z * S.kap[2]);
        double l1, l2, cs, sn;
        eig2sym(e11, e12, e22, l1, l2, cs, sn);
        lemax = std::max(lemax, std::max(l1, p.e33));
      }
      tcur = t0 * std::exp(S.pt[np / 2].e33);
      S.Meff = plane_stress_modulus(m);
    }
    for (int k = 0; k < np; ++k) {
      const Pt& p = S.pt[k];
      double z = zk[k] * 0.5 * tcur;
      for (int i = 0; i < 3; ++i) { N[i] += wk[k] * tcur * p.s[i]; M[i] += wk[k] * tcur * p.s[i] * z; }
      double vm = std::sqrt(std::max(0.0, p.s[0] * p.s[0] - p.s[0] * p.s[1] + p.s[1] * p.s[1] + 3 * p.s[2] * p.s[2]));
      if (k == 0 || k == np - 1) mis = std::max(mis, vm);  // as PEEQ: bottom and top section points
    }
    // transverse shear (elastic resultants)
    S.Q[0] += ts_stiff(e, 0) * dg[0];
    S.Q[1] += ts_stiff(e, 1) * dg[1];
    const double Q[2] = {S.Q[0], S.Q[1]};

    // stable increment: element eigenvalue of membrane + transverse shear, rotary check
    double smax_p = 0;  // largest principal membrane stress (geometric stiffness)
    {
      double l1, l2, cs, sn;
      eig2sym(N[0] / tcur, N[2] / tcur, N[1] / tcur, l1, l2, cs, sn);
      smax_p = std::max(0.0, l1);
    }
    const double Me = (meff_now > 0 ? meff_now : S.Meff) * stretch_fac + smax_p, nu = nu_eff(e);
    double D[3][3];
    iso_plane_stress_D(Me, nu, D);
    double lam = max_eig_DG(D, g1, n);
    double BB = 0;
    for (int I = 0; I < n; ++I) BB += g1.B1[I] * g1.B1[I] + g1.B2[I] * g1.B2[I];
    double kts = std::max(ts_stiff(e, 0), ts_stiff(e, 1));
    const double mtr = rho * t0 * g1.A / n;
    double Imin = 1e300;
    for (int I = 0; I < n; ++I) Imin = std::min(Imin, irot[e * n + I] * msf[e]);
    // hourglass (stiffness form) and drilling stiffnesses, evaluated here so that they enter the bound
    //  - membrane modes u_x, u_y and rotation modes theta_y, -theta_x (k, kb: 2x2 matrices per unit modulus,
    //    scaled by E t and E t^3/12):
    //    * default: Flanagan-Belytschko form, G/E A sum(B1^2 + B2^2) times 0.005 for the membrane modes (the
    //      documented default factor) and 0.18 for the rotation modes, same for x and y. The membrane value is
    //      much softer than the physical stiffness of in-plane bending: coarse meshes are more flexible, and
    //      in-plane bending of a single row of elements is mostly carried by the drilling stiffness; these
    //      factors reproduce the Abaqus results of the validation suite (plates, roofs, cylinders, hemisphere,
    //      in-plane and twisted cantilevers, thick plate);
    //    * enhanced hourglass control (*SECTION CONTROLS, HOURGLASS=ENHANCED): stiffness of the hourglass
    //      modes from the incompatible-mode element (qm6_hourglass, initial geometry) with the modulus matrix
    //      E [1 nu 0; nu 1 0; 0 0 (1-nu)/2] (that is (1 - nu^2) times the plane-stress matrix), which gives the
    //      Abaqus results for rectangular and trapezoidal elements;
    //  - transverse shear: the edge strains (edge_shear_ops) vary linearly across the element; their
    //    variations get the stiffness of the exact integral of the linear field (kts A / 12 per variation
    //    for a rectangle; the uniform twisting mode w = c x y with consistent rotations has no edge strain
    //    and stays free).
    double kxx = 0, kyy = 0, kxy = 0, klam = 0, kbx = 0, kby = 0, kbxy = 0, kblam = 0, gg = 0;
    // transverse shear of quads (edge strains): kv = stiffness of the two variations, hv = sum over the dofs of
    // (derivative^2 / nodal mass), hq = the same for the uniform shear (bounds on the element frequency)
    double kv[2] = {0, 0}, hv[2] = {0, 0}, hq = 0;
    const double Eh = Me * (1 - nu * nu), Em = Eh * tcur, Db = Eh * tcur * tcur * tcur / 12.0;
    if (n == 4) {
      const double r = sc.hg_scale;
      for (int I = 0; I < 4; ++I) gg += gm.gam[I] * gm.gam[I];
      if (sc.hg_enhanced || type == ElType::S4) {  // S4 (full integration): stiffness of a physical element
        const double f = r * (1 - nu * nu);
        kxx = kbx = f * khg[e * 3]; kyy = kby = f * khg[e * 3 + 1]; kxy = kbxy = f * khg[e * 3 + 2];
      } else {
        double bb = 0;
        for (int I = 0; I < 4; ++I) bb += gm.B1[I] * gm.B1[I] + gm.B2[I] * gm.B2[I];
        const double fb = r / (2 * (1 + nu)) * gm.A * bb;
        // static procedures: membrane 0.005 G (documented default), rotations 0.18 (fitted on Abaqus/Standard).
        // Explicit: instantaneous stiffnesses 0.2 G t A sum|B|^2 (membrane) and 0.2 G t^3 A sum|B|^2 (rotations)
        // for the vectors gamma, which relax (below); single-element Abaqus/Explicit runs, HOURGLASS=STIFFNESS
        kxx = kyy = (in.dynamic ? 0.2 : 0.005) * fb;
        kbx = kby = (in.dynamic ? 2.4 : 0.18) * fb;
      }
      klam = 0.5 * (kxx + kyy) + std::sqrt(0.25 * (kxx - kyy) * (kxx - kyy) + kxy * kxy);
      kblam = 0.5 * (kbx + kby) + std::sqrt(0.25 * (kbx - kby) * (kbx - kby) + kbxy * kbxy);
    }
    {
      // variation of the edge strains. Quad: linear across the element, integral of eta^2 (xi^2) = A / 3;
      // triangle: rotational part of the shear field, polar moment of the area
      if (n == 4 && type == ElType::S4) {
        kv[0] = sc.hg_scale * kts * gm.A / 3.0 * (invm[0] * invm[0] + invm[2] * invm[2]);
        kv[1] = sc.hg_scale * kts * gm.A / 3.0 * (invm[1] * invm[1] + invm[3] * invm[3]);
      } else if (n == 4) {
        // S4R: stiffness of the stabilised mode kts A sum|grad N_I|^2 / (12 (gamma.gamma - 2)) with the centre
        // gradients (kts A sum / 24 on a parallelogram); with the uniform part, the stiffness matrix of the Abaqus
        // element (single-element probes on rectangles, parallelograms, trapezoids, kites, distorted quads)
        double bb = 0;
        for (int I = 0; I < 4; ++I) bb += Gm[0][I][0] * Gm[0][I][0] + Gm[1][I][0] * Gm[1][I][0];
        kv[0] = sc.hg_scale * kts * gm.A * bb / (12.0 * std::max(ggm - 2.0, 1e-3));
        kv[1] = 0.0;
      } else {
        kv[0] = kts * cwm;
      }
      double hrot[2] = {0, 0};
      for (int I = 0; I < n; ++I) {
        for (int c = 0; c < 2; ++c) {
          hq += Gm[c][I][0] * Gm[c][I][0] / mtr;
          hrot[c] += (Gm[c][I][1] * Gm[c][I][1] + Gm[c][I][2] * Gm[c][I][2]) / Imin;
          hv[c] += Gm[2 + c][I][0] * Gm[2 + c][I][0] / mtr +
                   (Gm[2 + c][I][1] * Gm[2 + c][I][1] + Gm[2 + c][I][2] * Gm[2 + c][I][2]) / Imin;
        }
      }
      hq += std::max(hrot[0], hrot[1]);  // rotations: each component acts mostly on one rotation
    }
    // drilling stiffness: nodal rotation about the normal against the in-plane spin of the displacement field.
    // Quads: spin of the bilinear field at each corner, kd = G t A t^2 / (24 (t^2 + 0.1 A)) (about 0.42 G t^3 for
    // thin shells, G t A / 24 for thick ones); with the membrane and rotation hourglass stiffnesses above this
    // gives the stiffness matrix of the Abaqus element (single-element probes: hourglass and uniform rotation
    // patterns, rectangles 1x1 to 10x2, t 0.1 to 1). Triangles: spin of the linear field (constant), a tenth of
    // the bending stiffness E t^3/12 A sum(B1^2 + B2^2)
    const double Gd = Eh / (2.0 * (1.0 + nu));
    const double kd = Gd * tcur * g1.A * tcur * tcur / (24.0 * (tcur * tcur + 0.1 * g1.A));
    // bound on the element frequency: membrane | bending + transverse shear (coupled w / rotations);
    // hourglass modes: nodal stiffness k gamma gamma^T, largest eigenvalue k_max sum(gamma^2)
    // hourglass modes are orthogonal to the constant-strain modes (exactly for parallelograms): the largest
    // membrane frequency is the larger of the two families, not their sum
    // (drilling: its coupling to the in-plane displacements, corner gradients at most twice the mean ones)
    double om_m2 = std::max(n * lam / rho, Em * klam * gg / mtr) + (n == 3 ? 4.2 : 1.0) * kd * BB / mtr;
    // bending side kept additive (conservative): with the max of the shear families the increment grows and
    // the margin left for penalty contact is too small (energy created in the plate impact test)
    double om_b2 = g1.A * tcur * tcur * tcur * lam / (12.0 * Imin) +
                   kts * g1.A * hq +
                   Db * kblam * gg / Imin + 2.0 * std::max(kv[0] * hv[0], kv[1] * hv[1]) + (n == 3 ? 4.2 : 2.0) * kd / Imin;
    // Linear bulk viscosity of shells (Abaqus/Explicit doc, "Bulk viscosity"): membrane pressure
    // p = b1 rho c_d L_e eps_vol_dot and "pressure moment" m = b1 h0^2/12 rho c_d L kappa_dot
    // (kappa_dot = k11_dot + k22_dot), both with the dilatational wave speed c_d of the material
    // (bulk modulus included: for nearly incompressible hyperelastic leaflets this is the dominant
    // damping of curvature changes). The damping enters the stable increment through
    // dt = 2/omega (sqrt(1 + xi^2) - xi), xi = nu / (2 omega).
    double dtm, dtb;
    {
      double om_m = std::sqrt(std::max(om_m2, 1e-300)), om_b = std::sqrt(std::max(om_b2, 1e-300));
      double xim = 0, xib = 0;
      if (in.dynamic && in.bv1 > 0) {
        const double cd = std::sqrt(dilatational_modulus(m) / rho);
        const double Le = std::sqrt(2.0 / std::max(BB, 1e-300));
        const double kv = in.bv1 * rho * cd * Le;
        if (in.dt > 0) {
          double dvol = de[0] + de[1];  // in-plane dilatation rate
          double pbv = kv * dvol / in.dt;
          N[0] += pbv * tcur; N[1] += pbv * tcur;
          out.e_visc += pbv * dvol * tcur * g1.A;
          double mbv = kv * t0 * t0 / 12.0 * (dk[0] + dk[1]) / in.dt;
          M[0] += mbv * tcur; M[1] += mbv * tcur;
          out.e_visc += mbv * tcur * (dk[0] + dk[1]) * g1.A;
        }
        double Dv[3][3] = {{kv, kv, 0}, {kv, kv, 0}, {0, 0, 0}};
        double nu_m = n * max_eig_DG(Dv, g1, n) / rho;
        double nu_b = g1.A * tcur * kv * t0 * t0 / 12.0 * BB / Imin;
        xim = std::max(in.bv1, 0.5 * nu_m / om_m);
        xib = 0.5 * nu_b / om_b;
      }
      dtm = 2.0 / om_m * (std::sqrt(1 + xim * xim) - xim);
      dtb = 2.0 / om_b * (std::sqrt(1 + xib * xib) - xib);
      static int trace_dt = std::getenv("ABQ_TRACE_DT") ? std::atoi(std::getenv("ABQ_TRACE_DT")) : -1;
      if (commit && model->elems[ids[e]].label == trace_dt)
        std::fprintf(stderr, "DTTRACE om_m %.3e xim %.3f dtm %.3e | om_b %.3e xib %.3f dtb %.3e | Imin %.3e mtr %.3e\n", om_m, xim, dtm, om_b, xib, dtb, Imin, mtr);
    }
    double dt_el = std::min(dtm, dtb);
    out.dt_stable = dt_el;

    // local nodal forces (end configuration)
    double fl[4][6];
    std::memset(fl, 0, sizeof(fl));
    const double A = g1.A;
    double G1[4][4][3], inv1[4], cw1;
    if (n == 4 && type == ElType::S4) edge_shear_ops(g1.xl, G1, inv1);
    else if (n == 4) quad_shear_ops(g1, G1);
    else if (type == ElType::S3RS) tri_edge_shear_ops(g1, G1);
    else tri_degenerate_shear_ops(g1, tri_k[e], &tri_c[e * 2], G1, cw1);
    for (int I = 0; I < n; ++I) {
      double b1 = g1.B1[I], b2 = g1.B2[I];
      fl[I][0] += A * (b1 * N[0] + b2 * N[2]);
      fl[I][1] += A * (b2 * N[1] + b1 * N[2]);
      fl[I][3] += A * (-b2 * M[1] - b1 * M[2]);
      fl[I][4] += A * (b1 * M[0] + b2 * M[2]);
      for (int q = 0; q < 3; ++q) fl[I][2 + q] += A * (G1[0][I][q] * Q[0] + G1[1][I][q] * Q[1]);
    }
    if (n == 3 && type != ElType::S3RS) {  // rotational part of the shear field of the triangle (physical shear energy)
      double dgc = 0;
      for (int I = 0; I < 3; ++I) dgc += Gm[2][I][0] * du[I][2] + Gm[2][I][1] * dth[I][0] + Gm[2][I][2] * dth[I][1];
      S.hgs[0] += kv[0] * dgc;
      for (int I = 0; I < 3; ++I)
        for (int q = 0; q < 3; ++q) fl[I][2 + q] += G1[2][I][q] * S.hgs[0];
    }
    if (n == 4) {  // hourglass control (stiffness form, incremental)
      // warping correction: on a warped element (nodes off the mean plane by z_I), a rigid rotation moves the
      // nodes in the plane by -/+ z_I theta; the mean nodal rotation increment is removed before the membrane
      // hourglass modes are measured
      double thx = 0, thy = 0;
      for (int I = 0; I < 4; ++I) { thx += 0.25 * dth[I][0]; thy += 0.25 * dth[I][1]; }
      double dqx = 0, dqy = 0, dqtx = 0, dqty = 0, dgx = 0, dgy = 0;
      for (int I = 0; I < 4; ++I) {
        const double zI = gm.xl[I][2];
        dqx += gm.gam[I] * (du[I][0] - zI * thy); dqy += gm.gam[I] * (du[I][1] + zI * thx);
        dqtx += gm.gam[I] * dth[I][0]; dqty += gm.gam[I] * dth[I][1];
        // variations of the edge shear strains across the element
        dgx += Gm[2][I][0] * du[I][2] + Gm[2][I][1] * dth[I][0] + Gm[2][I][2] * dth[I][1];
        dgy += Gm[3][I][0] * du[I][2] + Gm[3][I][1] * dth[I][0] + Gm[3][I][2] * dth[I][1];
      }
      double o[6] = {S.hgm[0], S.hgm[1], S.hgb[0], S.hgb[1], S.hgs[0], S.hgs[1]};
      // membrane modes (u_x, u_y); rotation modes: theta_y acts as u_x and -theta_x as u_y
      S.hgm[0] += Em * (kxx * dqx + kxy * dqy); S.hgm[1] += Em * (kxy * dqx + kyy * dqy);
      S.hgb[0] += Db * (kby * dqtx - kbxy * dqty); S.hgb[1] += Db * (kbx * dqty - kbxy * dqtx);
      S.hgs[0] += kv[0] * dgx; S.hgs[1] += kv[1] * dgy;
      // relaxed stiffness of the membrane and rotation modes (explicit default, HOURGLASS=RELAX STIFFNESS), as for
      // C3D8R: force = Q_e - 0.9 Q_v, Q_v following Q_e with tau = 500 times the undamped stable increment of the
      // element, L / c_p with L^2 = max(1 / sum|B|^2, 0.92 L_m^2), L_m from the one-point membrane eigenvalue (the
      // Abaqus element increment); transverse shear does not relax (single-element Abaqus/Explicit runs)
      double hf[4] = {S.hgm[0], S.hgm[1], S.hgb[0], S.hgb[1]}, hfo[4] = {o[0], o[1], o[2], o[3]};
      if (in.dynamic && sc.hg_relax && !sc.hg_enhanced && type != ElType::S4) {
        double rb = 0, ra = 0;
        if (in.dt > 0) {
          const double L2 = std::max(1.0 / std::max(BB, 1e-300), 0.92 * 4.0 * Me / std::max(n * lam, 1e-300));
          const double tau = 500.0 * std::sqrt(L2 * rho / Me), x = in.dt / tau;
          rb = x < 1e-3 ? x * (1.0 - x * (0.5 - x / 6.0)) : 1.0 - std::exp(-x);
          ra = 1.0 - rb / x;
        }
        double* v[4] = {&S.hgmv[0], &S.hgmv[1], &S.hgbv[0], &S.hgbv[1]};
        for (int k = 0; k < 4; ++k) {
          const double vo = *v[k];
          *v[k] = vo * (1.0 - rb) + hfo[k] * rb + (hf[k] - hfo[k]) * ra;
          hfo[k] -= 0.9 * vo;
          hf[k] -= 0.9 * *v[k];
        }
      }
      out.e_art += 0.5 * ((hfo[0] + hf[0]) * dqx + (hfo[1] + hf[1]) * dqy + (hfo[2] + hf[2]) * dqtx +
                          (hfo[3] + hf[3]) * dqty + (o[4] + S.hgs[0]) * dgx + (o[5] + S.hgs[1]) * dgy);
      // no viscous term on the hourglass modes: the force of a single S4R launched at 1000 mm/s is the stiffness
      // force K u from the first increment (Abaqus/Explicit, free elements and a clamped plate kicked at one node:
      // same displacements to 1e-5, same spectrum up to the highest plate modes)
      const double Gx = S.hgs[0], Gy = S.hgs[1];
      // moments conjugate to the warping correction: d(q_x)/d(theta_y) = -sum(gamma z)/4 at every node,
      // d(q_y)/d(theta_x) = +sum(gamma z)/4 (without them the stiffness is not symmetric and the rigid
      // rotations of a warped element have a negative energy)
      double gz = 0;
      for (int I = 0; I < 4; ++I) gz += g1.gam[I] * g1.xl[I][2];
      for (int I = 0; I < 4; ++I) {
        double gI = g1.gam[I];
        fl[I][0] += gI * hf[0]; fl[I][1] += gI * hf[1];
        fl[I][4] -= 0.25 * gz * hf[0]; fl[I][3] += 0.25 * gz * hf[1];
        fl[I][3] += gI * hf[2];
        fl[I][4] += gI * hf[3];
        for (int q = 0; q < 3; ++q) fl[I][2 + q] += G1[2][I][q] * Gx + G1[3][I][q] * Gy;
      }
    }
    // drilling. Static procedures: corner rotation against the spin of the bilinear field at that corner (quads),
    // weighted constraints of the triangle (Abaqus/Standard stiffness matrices). Explicit: every corner rotation
    // against the mean spin of the element (centre gradients: the hourglass modes do not load the drilling); the
    // triangle as a degenerate quadrilateral whose doubled corner is the first node of the connectivity (Abaqus/
    // Explicit runs: kd per corner, 2 kd at the first node, no coupling between corners; strip of triangles: the
    // reaction moments of the nodes match only with the first node doubled)
    if (n == 3 && in.dynamic) {
      const int kc = 0;
      double dsum = 0, msum = 0;
      for (int I = 0; I < 3; ++I) {
        const double wI = I == kc ? 2.0 : 1.0, dphi = dth[I][2] - wz, old = S.md[I];
        S.md[I] += wI * kd * dphi;
        out.e_art += 0.5 * (old + S.md[I]) * dphi;
        fl[I][5] += S.md[I];
        msum += S.md[I];
        dsum += dphi;
      }
      (void)dsum;
      for (int I = 0; I < 3; ++I) {
        fl[I][0] += 0.5 * g1.B2[I] * msum;
        fl[I][1] += -0.5 * g1.B1[I] * msum;
      }
    } else if (n == 4) {
      // warping correction as for the membrane hourglass modes: the in-plane motion -/+ z_I theta of the warped
      // nodes under a rigid rotation about an in-plane axis is removed before the spins are measured
      double thxm = 0, thym = 0;
      for (int I = 0; I < 4; ++I) { thxm += 0.25 * dth[I][0]; thym += 0.25 * dth[I][1]; }
      double mx = 0, my = 0;  // conjugate moments of the correction (spread over the nodes)
      for (int J = 0; J < 4; ++J) {
        double a[4], b[4];
        if (in.dynamic) for (int I = 0; I < 4; ++I) { a[I] = gm.B1[I]; b[I] = gm.B2[I]; }
        else corner_gradients(gm, J, a, b);
        double wJ = 0;
        for (int I = 0; I < 4; ++I) {
          const double zI = gm.xl[I][2];
          wJ += 0.5 * (a[I] * (du[I][1] + zI * thxm) - b[I] * (du[I][0] - zI * thym));
        }
        const double dphi = dth[J][2] - wJ;
        const double old = S.md[J];
        S.md[J] += kd * dphi;
        out.e_art += 0.5 * (old + S.md[J]) * dphi;
        fl[J][5] += S.md[J];
        if (in.dynamic) for (int I = 0; I < 4; ++I) { a[I] = g1.B1[I]; b[I] = g1.B2[I]; }
        else corner_gradients(g1, J, a, b);
        for (int I = 0; I < 4; ++I) {
          fl[I][0] += 0.5 * b[I] * S.md[J]; fl[I][1] += -0.5 * a[I] * S.md[J];
          mx -= 0.5 * a[I] * g1.xl[I][2] * S.md[J]; my -= 0.5 * b[I] * g1.xl[I][2] * S.md[J];
        }
      }
      for (int I = 0; I < 4; ++I) { fl[I][3] += 0.25 * mx; fl[I][4] += 0.25 * my; }
    } else {  // triangles: rotations of the second and third corners against the element spin (2 kd each), plus
      // a weak constraint on (2 theta_1 + theta_2 + theta_3) / 4 (0.0125 x 16 kd): the stiffness matrix of the
      // Abaqus triangle (single-element probes; it depends on which corner comes first in the connectivity)
      const double d1 = dth[1][2] - wz, d2 = dth[2][2] - wz, d3 = 2 * dth[0][2] + dth[1][2] + dth[2][2] - 4 * wz;
      const double o1 = S.md[0], o2 = S.md[1], o3 = S.md[2];
      S.md[0] += 2 * kd * d1; S.md[1] += 2 * kd * d2; S.md[2] += 0.0125 * kd * d3;
      out.e_art += 0.5 * ((o1 + S.md[0]) * d1 + (o2 + S.md[1]) * d2 + (o3 + S.md[2]) * d3);
      const double m[3] = {2 * S.md[2], S.md[0] + S.md[2], S.md[1] + S.md[2]};
      double sum = 0;
      for (int I = 0; I < 3; ++I) { fl[I][5] += m[I]; sum += m[I]; }
      for (int I = 0; I < 3; ++I) {
        fl[I][0] += 0.5 * g1.B2[I] * sum;
        fl[I][1] += -0.5 * g1.B1[I] * sum;
      }
    }
    for (int I = 0; I < n; ++I) {
      Vec3 f = g1.e1 * fl[I][0] + g1.e2 * fl[I][1] + g1.e3 * fl[I][2];
      Vec3 mo = g1.e1 * fl[I][3] + g1.e2 * fl[I][4] + g1.e3 * fl[I][5];
      store3(out.f + 6 * I, f);
      store3(out.f + 6 * I + 3, mo);
    }
    static int trace_label = std::getenv("ABQ_TRACE_ELEM") ? std::atoi(std::getenv("ABQ_TRACE_ELEM")) : -1;
    if (commit && model->elems[ids[e]].label == trace_label) {
      std::fprintf(stderr, "ETRACE dt=%.3e Meff=%.3e Me=%.3e kap=(%.3e %.3e %.3e) em=(%.3e %.3e %.3e) N=(%.3e %.3e %.3e) M=(%.3e %.3e %.3e) Q=(%.3e %.3e) hg=(%.2e %.2e %.2e %.2e %.2e) md=(%.2e %.2e %.2e %.2e) A=%.4e fail=%d lam3=(%.3f %.3f %.3f)\n",
                   out.dt_stable, S.Meff, Me, S.kap[0], S.kap[1], S.kap[2], S.em[0], S.em[1], S.em[2], N[0], N[1], N[2], M[0], M[1], M[2], Q[0], Q[1],
                   S.hgm[0], S.hgm[1], S.hgb[0], S.hgb[1], S.hgs[0], S.md[0], S.md[1], S.md[2], S.md[3], g1.A, (int)out.failed, S.pt[0].lam3, S.pt[np/2].lam3, S.pt[np-1].lam3);
    }
    if (commit && type == ElType::S4) {
      // S4 output at its four integration points: the linear variation of the strains across the element is the
      // field of the hourglass modes (total amplitudes from the stored forces, elastic), added to the stresses of
      // the bottom and top section points; largest von Mises stress over the 4 x 2 points
      const double dm = kxx * kyy - kxy * kxy, db = kbx * kby - kbxy * kbxy;
      const double qx = dm > 0 ? (kyy * S.hgm[0] - kxy * S.hgm[1]) / (Em * dm) : 0;
      const double qy = dm > 0 ? (kxx * S.hgm[1] - kxy * S.hgm[0]) / (Em * dm) : 0;
      const double qtx = db > 0 ? (kbx * S.hgb[0] + kbxy * S.hgb[1]) / (Db * db) : 0;
      const double qty = db > 0 ? (kbxy * S.hgb[0] + kby * S.hgb[1]) / (Db * db) : 0;
      static const double XI4[4] = {-1, 1, 1, -1}, ET4[4] = {-1, -1, 1, 1};
      const double h[4] = {1, -1, 1, -1};
      double bh1 = 0, bh2 = 0, d[4], gd = 0;
      for (int I = 0; I < 4; ++I) { bh1 += g1.B1[I] * h[I]; bh2 += g1.B2[I] * h[I]; }
      for (int I = 0; I < 4; ++I) d[I] = h[I] - bh1 * g1.xl[I][0] - bh2 * g1.xl[I][1];
      for (int I = 0; I < 4; ++I) gd += g1.gam[I] * d[I];
      if (std::fabs(gd) > 1e-300) for (int I = 0; I < 4; ++I) d[I] /= gd;
      double Dp[3][3];
      iso_plane_stress_D(Me, nu, Dp);
      const double gp = 1.0 / std::sqrt(3.0);
      for (int ig = 0; ig < 4; ++ig) {
        const double xi = gp * XI4[ig], eta = gp * ET4[ig];
        double J[2][2] = {{0, 0}, {0, 0}}, dxi[4], deta[4];
        for (int I = 0; I < 4; ++I) {
          dxi[I] = 0.25 * XI4[I] * (1 + eta * ET4[I]); deta[I] = 0.25 * ET4[I] * (1 + xi * XI4[I]);
          J[0][0] += dxi[I] * g1.xl[I][0]; J[0][1] += dxi[I] * g1.xl[I][1];
          J[1][0] += deta[I] * g1.xl[I][0]; J[1][1] += deta[I] * g1.xl[I][1];
        }
        const double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (!(std::fabs(det) > 1e-300)) continue;
        double gx = 0, gy = 0;
        for (int I = 0; I < 4; ++I) {
          gx += (J[1][1] * dxi[I] - J[0][1] * deta[I]) / det * d[I];
          gy += (-J[1][0] * dxi[I] + J[0][0] * deta[I]) / det * d[I];
        }
        const double dem[3] = {qx * gx, qy * gy, qx * gy + qy * gx};
        const double dka[3] = {qty * gx, -qtx * gy, qty * gy - qtx * gx};
        for (int k : {0, np - 1}) {
          const double z = zk[k] * 0.5 * tcur;
          double sg[3];
          for (int i = 0; i < 3; ++i) {
            sg[i] = S.pt[k].s[i];
            for (int j = 0; j < 3; ++j) sg[i] += Dp[i][j] * (dem[j] + z * dka[j]);
          }
          mis = std::max(mis, std::sqrt(std::max(0.0, sg[0] * sg[0] - sg[0] * sg[1] + sg[1] * sg[1] + 3 * sg[2] * sg[2])));
        }
      }
    }
    if (commit) {
      st[e] = S;
      float* o = &this->out[e * N_OUT];
      o[OUT_MISES] = (float)mis;
      const double ca = oax[e * 2], sa = oax[e * 2 + 1], cc = ca * ca, ss = sa * sa, cs = ca * sa;
      o[OUT_S11] = (float)((cc * N[0] + ss * N[1] + 2 * cs * N[2]) / tcur);
      o[OUT_S22] = (float)((ss * N[0] + cc * N[1] - 2 * cs * N[2]) / tcur);
      o[OUT_S12] = (float)((cs * (N[1] - N[0]) + (cc - ss) * N[2]) / tcur);
      o[OUT_S33] = 0;
      o[OUT_S13] = (float)((ca * Q[0] + sa * Q[1]) / tcur); o[OUT_S23] = (float)((ca * Q[1] - sa * Q[0]) / tcur);
      o[OUT_PEEQ] = (float)epm;
      o[OUT_LEMAX] = (float)lemax;
      o[OUT_STH] = (float)tcur;
    }
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_shell_group() { return std::make_unique<ShellGroup>(); }

}  // namespace abq
