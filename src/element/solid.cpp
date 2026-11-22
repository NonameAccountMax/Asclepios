// Continuum elements C3D8R (uniform strain, Flanagan-Belytschko hourglass control), C3D8I (incompatible modes),
// C3D8 (selective
// reduced / B-bar, 2x2x2) and C3D4 (constant strain). Incremental corotational formulation:
// strain increment from the mid-increment configuration, stresses rotated with Hughes-Winget.
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "element/groups.hpp"
#include "material/material.hpp"

namespace abq {

namespace {

const double HEX_XI[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                             {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
const double HG[4][8] = {{1, 1, -1, -1, -1, -1, 1, 1},
                         {1, -1, -1, 1, -1, 1, 1, -1},
                         {1, -1, 1, -1, 1, -1, 1, -1},
                         {-1, 1, -1, 1, 1, -1, 1, -1}};
const double G8 = 0.5773502691896258;

// shape function derivatives of the trilinear hex at (r,s,t)
void hex_dN(double r, double s, double t, double dN[8][3]) {
  for (int I = 0; I < 8; ++I) {
    double a = HEX_XI[I][0], b = HEX_XI[I][1], c = HEX_XI[I][2];
    dN[I][0] = 0.125 * a * (1 + b * s) * (1 + c * t);
    dN[I][1] = 0.125 * b * (1 + a * r) * (1 + c * t);
    dN[I][2] = 0.125 * c * (1 + a * r) * (1 + b * s);
  }
}

// physical gradients b[I] = dN_I/dx at a point; returns det J
double hex_grad(const double* x, double r, double s, double t, Vec3 b[8]) {
  double dN[8][3];
  hex_dN(r, s, t, dN);
  Mat3 J;  // J[i][j] = dx_j / dxi_i
  for (int I = 0; I < 8; ++I)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) J.m[i][j] += dN[I][i] * x[3 * I + j];
  double det = J.det();
  Mat3 Ji = J.inverse();
  for (int I = 0; I < 8; ++I) {
    Vec3 g(dN[I][0], dN[I][1], dN[I][2]);
    b[I] = Ji * g;
  }
  return det;
}

// uniform (volume averaged) gradient and volume of a hex
double hex_mean_grad(const double* x, Vec3 b[8]) {
  Vec3 bb[8];
  double V = 0;
  for (int I = 0; I < 8; ++I) b[I] = Vec3();
  for (int g = 0; g < 8; ++g) {
    double d = hex_grad(x, HEX_XI[g][0] * G8, HEX_XI[g][1] * G8, HEX_XI[g][2] * G8, bb);
    V += d;
    for (int I = 0; I < 8; ++I) b[I] += bb[I] * d;
  }
  for (int I = 0; I < 8; ++I) b[I] = b[I] / V;
  return V;
}

double tet_grad(const double* x, Vec3 b[4]) {
  Vec3 p0 = load3(x), p1 = load3(x + 3), p2 = load3(x + 6), p3 = load3(x + 9);
  Mat3 J = Mat3::rows(p1 - p0, p2 - p0, p3 - p0);  // rows: dx/dxi_k
  double det = J.det();
  Mat3 Ji = J.inverse();
  Vec3 g[4] = {Vec3(-1, -1, -1), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
  for (int I = 0; I < 4; ++I) b[I] = Ji * g[I];
  return det / 6.0;
}

// element length of the explicit stable increment built on the largest face of a hex (area from the diagonals):
// L^2 = c(s) V^2 / A_max^2, with s = A_max^1.5 / V - 1 >= 0 measuring the elongation of the element (0 for a
// cube) and c rising from 0.504 to 0.670 with s. Fitted on the element increments of single-element Abaqus runs
// (boxes from 1x1x1 to 10x1x1, sheared, tapered and randomly distorted bricks, nu from 0 to 0.45); the stable
// length is the larger of this one and the uniform-gradient bound
double face_length(const double* x, double V) {
  static const int F[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  double amax = 0;
  for (const auto& f : F) {
    Vec3 d1 = load3(x + 3 * f[2]) - load3(x + 3 * f[0]), d2 = load3(x + 3 * f[3]) - load3(x + 3 * f[1]);
    amax = std::max(amax, 0.5 * norm(cross(d1, d2)));
  }
  if (amax <= 0 || V <= 0) return 0.0;
  double s = std::max(0.0, amax * std::sqrt(amax) / V - 1.0);
  double c = 0.504 + 0.166 * (1.0 - std::exp(-2.5 * s));
  return std::sqrt(c) * V / amax;
}


// solve A X = B in place (Gauss elimination with partial pivoting), A n x n, B n x m (row-major)
void dense_solve(std::vector<double>& A, std::vector<double>& B, int n, int m) {
  for (int c = 0; c < n; ++c) {
    int pv = c;
    for (int r = c + 1; r < n; ++r) if (std::fabs(A[r * n + c]) > std::fabs(A[pv * n + c])) pv = r;
    if (pv != c) {
      for (int j = 0; j < n; ++j) std::swap(A[c * n + j], A[pv * n + j]);
      for (int j = 0; j < m; ++j) std::swap(B[c * m + j], B[pv * m + j]);
    }
    const double d = A[c * n + c];
    if (std::fabs(d) < 1e-300) continue;
    for (int r = 0; r < n; ++r) {
      if (r == c) continue;
      const double f = A[r * n + c] / d;
      if (f == 0) continue;
      for (int j = c; j < n; ++j) A[r * n + j] -= f * A[c * n + j];
      for (int j = 0; j < m; ++j) B[r * m + j] -= f * B[c * m + j];
    }
  }
  for (int r = 0; r < n; ++r) {
    const double d = A[r * n + r];
    for (int j = 0; j < m; ++j) B[r * m + j] = std::fabs(d) > 1e-300 ? B[r * m + j] / d : 0.0;
  }
}

// Enhanced hourglass control of the reduced-integration brick: stiffness of its 12 hourglass modes taken from
// the incompatible-mode brick (2x2x2 points, 9 incompatible modes condensed, linear elastic lam, mu), written for the generalized hourglass displacements q_ai = sum_I u_Ii g_aI used by the element
// (g_a = (h_a - (h_a . x_I) b_I) / 8, b mean gradients): Kq = W^T Kc W with W the dual basis of the g_a
void enhanced_hourglass(const double* xe, double lam, double mu, double Kq[144]) {
  // the incompatible-mode stiffness and its projection are evaluated on the parallelepiped spanned by the centre
  // Jacobian (constant Jacobian: no distortion stiffening; reproduces the Abaqus results on trapezoidal elements)
  double x[24];
  {
    double dN[8][3];
    hex_dN(0, 0, 0, dN);
    Mat3 J;
    Vec3 c0;
    for (int I = 0; I < 8; ++I) {
      c0 += load3(xe + 3 * I) * 0.125;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) J.m[i][j] += dN[I][i] * xe[3 * I + j];
    }
    for (int I = 0; I < 8; ++I)
      for (int j = 0; j < 3; ++j) {
        double v = c0[j];
        for (int i = 0; i < 3; ++i) v += HEX_XI[I][i] * J.m[i][j];
        x[3 * I + j] = v;
      }
  }
  Vec3 b[8][8], gi[8][3];
  double d[8];
  for (int g = 0; g < 8; ++g) d[g] = hex_grad(x, HEX_XI[g][0] * G8, HEX_XI[g][1] * G8, HEX_XI[g][2] * G8, b[g]);
  Mat3 J0;
  {
    double dN[8][3];
    hex_dN(0, 0, 0, dN);
    for (int I = 0; I < 8; ++I)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) J0.m[i][j] += dN[I][i] * x[3 * I + j];
  }
  const double det0 = J0.det();
  const Mat3 Ji0 = J0.inverse();
  for (int g = 0; g < 8; ++g)
    for (int k = 0; k < 3; ++k) {
      Vec3 v;
      v[k] = -2.0 * HEX_XI[g][k] * G8;
      gi[g][k] = (Ji0 * v) * (det0 / d[g]);
    }
  // stiffness of 24 nodal + 9 internal dofs
  const int N = 33;
  std::vector<double> K(N * N, 0.0);
  auto strain = [](const Vec3& gr, int c, double v[6]) {
    for (int i = 0; i < 6; ++i) v[i] = 0;
    v[c] = gr[c];
    if (c == 0) { v[3] = gr.y; v[4] = gr.z; }
    if (c == 1) { v[3] = gr.x; v[5] = gr.z; }
    if (c == 2) { v[4] = gr.x; v[5] = gr.y; }
  };
  for (int g = 0; g < 8; ++g) {
    double Bm[33][6], DB[33][6];
    for (int I = 0; I < 8; ++I) for (int c = 0; c < 3; ++c) strain(b[g][I], c, Bm[3 * I + c]);
    for (int k = 0; k < 3; ++k) for (int c = 0; c < 3; ++c) strain(gi[g][k], c, Bm[24 + 3 * k + c]);
    for (int r = 0; r < N; ++r) {
      const double tr = Bm[r][0] + Bm[r][1] + Bm[r][2];
      for (int i = 0; i < 3; ++i) DB[r][i] = lam * tr + 2 * mu * Bm[r][i];
      for (int i = 3; i < 6; ++i) DB[r][i] = mu * Bm[r][i];
    }
    for (int r = 0; r < N; ++r)
      for (int c = 0; c < N; ++c) {
        double v = 0;
        for (int i = 0; i < 6; ++i) v += Bm[r][i] * DB[c][i];
        K[r * N + c] += v * d[g];
      }
  }
  // condensation: Kc = Kuu - Kua Kaa^-1 Kau
  std::vector<double> Kaa(81), Kau(9 * 24);
  for (int a = 0; a < 9; ++a) {
    for (int b2 = 0; b2 < 9; ++b2) Kaa[a * 9 + b2] = K[(24 + a) * N + 24 + b2];
    for (int u = 0; u < 24; ++u) Kau[a * 24 + u] = K[(24 + a) * N + u];
  }
  std::vector<double> X = Kau;
  dense_solve(Kaa, X, 9, 24);
  double Kc[24][24];
  for (int r = 0; r < 24; ++r)
    for (int c = 0; c < 24; ++c) {
      double v = K[r * N + c];
      for (int a = 0; a < 9; ++a) v -= K[r * N + 24 + a] * X[a * 24 + c];
      Kc[r][c] = v;
    }
  // hourglass vectors g_a (mean gradients) and their dual basis W (8 x 4): g^T W = I
  Vec3 bm[8];
  hex_mean_grad(x, bm);
  double gv[4][8];
  for (int a = 0; a < 4; ++a) {
    double hx[3] = {0, 0, 0};
    for (int I = 0; I < 8; ++I) for (int j = 0; j < 3; ++j) hx[j] += HG[a][I] * x[3 * I + j];
    for (int I = 0; I < 8; ++I) gv[a][I] = 0.125 * (HG[a][I] - (hx[0] * bm[I].x + hx[1] * bm[I].y + hx[2] * bm[I].z));
  }
  std::vector<double> GtG(16), W(32);
  for (int a = 0; a < 4; ++a)
    for (int c = 0; c < 4; ++c) {
      double v = 0;
      for (int I = 0; I < 8; ++I) v += gv[a][I] * gv[c][I];
      GtG[a * 4 + c] = v;
    }
  std::vector<double> Bm(16, 0.0);
  for (int a = 0; a < 4; ++a) Bm[a * 4 + a] = 1.0;
  dense_solve(GtG, Bm, 4, 4);  // (G^T G)^-1
  for (int I = 0; I < 8; ++I)
    for (int a = 0; a < 4; ++a) {
      double v = 0;
      for (int c = 0; c < 4; ++c) v += gv[c][I] * Bm[c * 4 + a];
      W[I * 4 + a] = v;
    }
  // Kq[(a,i),(c,j)] = sum_IJ W_Ia W_Jc Kc[(I,i),(J,j)]
  for (int a = 0; a < 4; ++a)
    for (int i = 0; i < 3; ++i)
      for (int c = 0; c < 4; ++c)
        for (int j = 0; j < 3; ++j) {
          double v = 0;
          for (int I = 0; I < 8; ++I)
            for (int J = 0; J < 8; ++J) v += W[I * 4 + a] * W[J * 4 + c] * Kc[3 * I + i][3 * J + j];
          Kq[(3 * a + i) * 12 + 3 * c + j] = v;
        }
}


struct PointState {  // per integration point
  Sym6 sig, eps;
  double ep = 0;
  double F[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double sdv = 0;
};

class SolidGroup : public ElementGroup {
 public:
  int npt = 1;
  std::vector<PointState> st;    // size()*npt
  std::vector<double> hgq;       // C3D8R: 12 per element (4 modes x 3 comps): generalized hourglass strains
                                 // (default control) or total forces (enhanced control)
  std::vector<double> hgl;       // C3D8R, explicit: element length in the initial configuration (relaxation time)
  std::vector<double> hgv;       // C3D8R, explicit: relaxed part of the hourglass strains (12 per element)
  std::vector<double> hgf;       // C3D8R, explicit: applied hourglass forces (12 per element)
  std::vector<char> hgr;         // C3D8R, explicit: relaxed hourglass stiffness (default) or constant
  std::vector<double> hgs;       // hourglass scale of the section (*SECTION CONTROLS)
  std::vector<char> fdt;         // stable increment of C3D8R with the face-based length (Section::face_dt)
  std::vector<double> hgK;       // C3D8R with enhanced hourglass control: 12 x 12 per element, empty otherwise
  std::vector<int> hgKi;         // index of the element's matrix in hgK, -1 for the default control
  std::vector<const Material*> mat;
  std::vector<int> orient;
  std::vector<double> fib;       // Holzapfel-Gasser-Ogden: 3 fibre directions per element (global, reference)

  void setup(const Model& m) override {
    model = &m;
    nen = type == ElType::C3D4 ? 4 : 8;
    npt = (type == ElType::C3D8 || type == ElType::C3D8I) ? 8 : 1;
    rot = false;
    st.assign(size() * npt, PointState());
    if (type == ElType::C3D8R) {
      hgq.assign(size() * 12, 0.0); hgv.assign(size() * 12, 0.0); hgf.assign(size() * 12, 0.0); hgr.assign(size(), 1); hgl.assign(size(), 0.0); hgs.assign(size(), 1.0); fdt.assign(size(), 1);
      hgKi.assign(size(), -1); hgK.clear();
    }
    mat.resize(size());
    orient.resize(size());
    for (int e = 0; e < size(); ++e) {
      const Section& s = m.sections[m.elems[ids[e]].section];
      mat[e] = &m.materials[s.material];
      orient[e] = s.orientation;
      if (mat[e]->hgo) {
        // fibre families: the LOCAL DIRECTIONS of the orientation (default: local axes 1, 2, 3), in global axes
        if (fib.empty()) fib.assign(size() * 9, 0.0);
        const Orientation* o = s.orientation >= 0 ? &m.orients[s.orientation] : nullptr;
        const Mat3 R = o ? o->R : Mat3::identity();
        for (int f = 0; f < 3; ++f) {
          Vec3 l = o && f < (int)o->local_dirs.size() ? o->local_dirs[f] : Vec3(f == 0, f == 1, f == 2);
          Vec3 a = normalized(R.row(0) * l.x + R.row(1) * l.y + R.row(2) * l.z);
          store3(&fib[e * 9 + 3 * f], a);
        }
      }
      if (type == ElType::C3D8R) {
        double x[24];
        for (int I = 0; I < 8; ++I) store3(x + 3 * I, m.X[conn[e * 8 + I]]);
        Vec3 b[8];
        double V = hex_mean_grad(x, b);
        double bb = 0;
        for (int I = 0; I < 8; ++I) bb += norm2(b[I]);
        // explicit procedures: the relaxation time of the hourglass forces is set by the initial element (the
        // undamped stable-increment length of the undeformed brick), the stiffness follows the current one
        hgl[e] = std::max(2.0 / std::sqrt(8 * bb), s.face_dt ? face_length(x, V) : 0.0);
        hgr[e] = s.hg_relax;
        hgs[e] = s.hg_scale;
        fdt[e] = s.face_dt;
        // enhanced hourglass control: requested (*SECTION CONTROLS, HOURGLASS=ENHANCED) or, as the Abaqus default,
        // for hyperelastic materials
        if (s.hg_enhanced || mat[e]->hyper()) {
          const double mu0 = mat[e]->shear_modulus(), lam0 = std::max(0.0, dilatational_modulus(*mat[e]) - 2 * mu0);
          hgKi[e] = (int)(hgK.size() / 144);
          hgK.resize(hgK.size() + 144);
          enhanced_hourglass(x, lam0, mu0, &hgK[hgKi[e] * 144]);
          for (int k = 0; k < 144; ++k) hgK[hgKi[e] * 144 + k] *= s.hg_scale;
        }
      }
    }
  }

  void save(BinOut& o) const override { ElementGroup::save(o); o.vec(st); o.vec(hgq); o.vec(hgv); o.vec(hgf); o.vec(hgl); o.vec(hgs); }
  void load(BinIn& in) override { ElementGroup::load(in); in.vec(st); in.vec(hgq); in.vec(hgv); in.vec(hgf); in.vec(hgl); in.vec(hgs); }
  double alpha_damping(int e) const override { return mat[e]->damp_alpha; }

  void mass(int e, const double* X, double* mt, double* mr) const override {
    double V = type == ElType::C3D4 ? tet_vol(X) : hex_vol(X);
    double m = mat[e]->density * V / nen;
    for (int I = 0; I < nen; ++I) { mt[I] = m; mr[I] = 0; }
  }
  static double hex_vol(const double* X) { Vec3 b[8]; return hex_mean_grad(X, b); }
  static double tet_vol(const double* X) { Vec3 b[4]; return tet_grad(X, b); }

  // material update at one point; returns plastic dissipation density
  double update_point(const Material& m, PointState& p, const Sym6& deps, const Vec3& w, const Mat3& dF, bool rotate,
                      double dtemp, const double* fibres = nullptr) {
    Mat3 Q = rotate ? hughes_winget(w) : Mat3::identity();
    if (rotate) { p.sig = abq::rotate(Q, p.sig); p.eps = abq::rotate(Q, p.eps); }
    for (int i = 0; i < 6; ++i) p.eps[i] += deps[i];
    if (m.hyper()) {
      Mat3 F;
      std::memcpy(&F.m[0][0], p.F, sizeof(p.F));
      F = dF * F;
      std::memcpy(p.F, &F.m[0][0], sizeof(p.F));
      hyper_cauchy_3d(m, F, p.sig, nullptr, fibres);
      return 0;
    }
    Sym6 d = deps;
    if (m.alpha != 0) for (int i = 0; i < 3; ++i) d[i] -= m.alpha * dtemp;
    return j2_update_3d(m, d, p.sig, p.ep);
  }

  void compute(int e, const ElemIn& in, ElemOut& out, bool commit) override {
    const Material& m = *mat[e];
    const int n = nen;
    double xm[24], x1[24];
    const double* xref = in.nlgeom ? nullptr : in.X;
    for (int I = 0; I < n; ++I)
      for (int k = 0; k < 3; ++k) {
        xm[3 * I + k] = xref ? xref[3 * I + k] : in.xo[3 * I + k] + 0.5 * in.du[6 * I + k];
        x1[3 * I + k] = xref ? xref[3 * I + k] : in.xo[3 * I + k] + in.du[6 * I + k];
      }
    double dtemp = 0;
    if (in.temp0 && in.temp1) { for (int I = 0; I < n; ++I) dtemp += in.temp1[I] - in.temp0[I]; dtemp /= n; }
    for (int i = 0; i < 6 * n; ++i) out.f[i] = 0;
    const bool rotate = in.nlgeom;
    const double rho = m.density * msf[e];
    double M = dilatational_modulus(m);
    if (m.hgo && !fib.empty()) {  // fibres stiffen with the stretch: current modulus along the fibres
      Mat3 F0;
      std::memcpy(&F0.m[0][0], st[e * npt].F, sizeof(st[e * npt].F));
      M += hgo_fibre_modulus(m, F0, &fib[e * 9]);
    }
    double Le, V1, Ldt = 0;  // Le: bulk viscosity and stable increment; Ldt: stable increment when larger
    double trd = 0;
    double ep_max = 0, mis_max = 0;
    Sym6 sig_avg, eps_avg;

    if (type == ElType::C3D8R || type == ElType::C3D4) {
      Vec3 bm[8], b1[8];
      if (type == ElType::C3D8R) { hex_mean_grad(xm, bm); V1 = hex_mean_grad(x1, b1); }
      else { tet_grad(xm, bm); V1 = tet_grad(x1, b1); }
      {
        double bb = 0;
        for (int I = 0; I < n; ++I) bb += norm2(b1[I]);
        Le = 2.0 / std::sqrt(n * bb);  // Flanagan-Belytschko bound on the element frequency
        // the face-based length lengthens the stable increment only: in the bulk viscosity it would raise the
        // damping of the highest modes and make the run unstable at the global increment
        if (type == ElType::C3D8R && fdt[e]) Ldt = face_length(x1, V1);
      }
      // displacement-increment gradient at mid configuration
      Mat3 L;
      for (int I = 0; I < n; ++I)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) L.m[i][j] += in.du[6 * I + i] * bm[I][j];
      Sym6 deps = Sym6::from(L);
      Vec3 w(L.m[2][1] - L.m[1][2], L.m[0][2] - L.m[2][0], L.m[1][0] - L.m[0][1]);
      w *= 0.5;
      // deformation-gradient increment relative to start configuration: dF = I + L (1 + L/2)^-1 approx
      Mat3 dF = Mat3::identity();
      if (m.hyper()) {
        Vec3 bo[8];
        double xo[24];
        for (int i = 0; i < 3 * n; ++i) xo[i] = xref ? xref[i] : in.xo[i];
        if (type == ElType::C3D8R) hex_mean_grad(xo, bo); else tet_grad(xo, bo);
        for (int I = 0; I < n; ++I)
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) dF.m[i][j] += in.du[6 * I + i] * bo[I][j];
      }
      PointState p = st[e];
      double wp = update_point(m, p, deps, w, dF, rotate, dtemp, fib.empty() ? nullptr : &fib[e * 9]);
      // explicit: while the point flows plastically, the increments of the applied hourglass forces are scaled by
      // max(0.1, H / (3 G + H)) (single-element Abaqus/Explicit runs: 0.1 for perfect plasticity and for a slope H
      // of 0.08 E, with HOURGLASS=STIFFNESS and RELAX STIFFNESS alike: an hourglass deformation added during the
      // flow gets 0.1 of the stiffness, the forces built while elastic are kept and relax 10 times slower)
      double hg_flow = 1.0;
      if (in.dynamic && p.ep > st[e].ep) {
        double H = 0;
        m.yield(p.ep, &H);
        const double G = m.shear_modulus();
        hg_flow = std::max(0.1, std::min(1.0, H / (3 * G + H)));
      }
      out.e_plastic += wp * V1;
      trd = deps.trace();
      Sym6 s = p.sig;
      // bulk viscosity
      if (in.dynamic && in.dt > 0 && (in.bv1 > 0 || in.bv2 > 0)) {
        double rate = trd / in.dt;
        double c = std::sqrt(M / rho);
        const double Lb = in.initial ? std::max(Le, Ldt) : Le;  // increment 0: stable-increment length, linear term
        double pbv = in.bv1 * rho * c * Lb * rate +
                     (in.initial ? 0.0 : rho * in.bv2 * in.bv2 * Le * Le * rate * std::fabs(std::min(0.0, rate)));
        for (int i = 0; i < 3; ++i) s[i] += pbv;
        out.e_visc += -(-pbv) * trd * V1;  // pbv * dvol (always >= 0)
      }
      for (int I = 0; I < n; ++I) {
        Vec3 f = Vec3(s[0] * b1[I].x + s[3] * b1[I].y + s[4] * b1[I].z, s[3] * b1[I].x + s[1] * b1[I].y + s[5] * b1[I].z,
                      s[4] * b1[I].x + s[5] * b1[I].y + s[2] * b1[I].z) * V1;
        add3(out.f + 6 * I, f);
      }
      // hourglass (C3D8R)
      if (type == ElType::C3D8R) {
        double* q = &hgq[e * 12];
        double* qv = &hgv[e * 12];
        double* qf = &hgf[e * 12];
        double qn[12], qvn[12], qfn[12];
        std::memcpy(qn, q, sizeof(qn));
        std::memcpy(qvn, qv, sizeof(qvn));
        std::memcpy(qfn, qf, sizeof(qfn));
        if (rotate) {
          Mat3 Q = hughes_winget(w);
          for (int a = 0; a < 4; ++a) {
            store3(qn + 3 * a, Q * load3(qn + 3 * a)); store3(qvn + 3 * a, Q * load3(qvn + 3 * a));
            store3(qfn + 3 * a, Q * load3(qfn + 3 * a));
          }
        }
        // relaxed stiffness (explicit default, HOURGLASS=RELAX STIFFNESS): standard linear solid on each mode,
        // relaxed strain q_v following q with the relaxation time tau (exponential kernel, integrated exactly for
        // a linear q over the increment), r = 0.9 (long-term stiffness 0.1 times the instantaneous one); tau = 1000
        // times the undamped stable increment of the initial element (independent of later stretching). Fitted on
        // single-element Abaqus/Explicit runs (hourglass mode prescribed, then held: relaxation of the reaction,
        // independent of the step time, scaling with size, density and modulus; bricks stretched beforehand)
        const bool relax = in.dynamic && hgr[e] && hgKi[e] < 0;
        const double rr = relax ? 0.9 : 0.0;
        double ra = 0, rb = 0;
        if (relax && in.dt > 0) {
          const double tau = 1000.0 * hgl[e] / std::sqrt(M / rho), x = in.dt / tau;
          rb = x < 1e-3 ? x * (1.0 - x * (0.5 - x / 6.0)) : 1.0 - std::exp(-x);  // 1 - exp(-x)
          ra = 1.0 - rb / x;      // weight of the increment of Q_e
        }
        // stiffness per spatial direction i. Explicit: k (instantaneous value, isotropic).
        // Implicit procedures, materials other than hyperelastic: documented default 0.005 G, for hourglass
        // vectors with unit nodal values, with a factor per direction proportional to sum_I (dN_I/dx_i)^2 (current
        // mean gradients): k_i = c G V sum_I (dN_I/dx_i)^2, c = 0.0754 for the vectors gamma below (1/8 included),
        // which reproduces the Abaqus results of bending, torsion, pressure and the NAFEMS LE10 plate
        // explicit: instantaneous stiffness 8/15 G V sum_I |grad N_I|^2 on the current configuration (isotropic;
        // single-element Abaqus/Explicit runs, HOURGLASS=STIFFNESS and RELAX STIFFNESS, bricks stretched from
        // lambda = 0.5 to 1.5: exact to 1e-4), G the elastic shear modulus also during plastic flow
        double bbc = 0;
        for (int I = 0; I < 8; ++I) bbc += norm2(b1[I]);
        const double kc = 8.0 / 15.0 * hgs[e] * mat[e]->shear_modulus() * V1 * bbc;
        double kd3[3] = {kc, kc, kc};
        if (!in.dynamic && !(mat[e]->hyper())) {
          double Si[3] = {0, 0, 0};
          for (int I = 0; I < 8; ++I) { Si[0] += b1[I].x * b1[I].x; Si[1] += b1[I].y * b1[I].y; Si[2] += b1[I].z * b1[I].z; }
          const double c = 0.0754 * hgs[e] * mat[e]->shear_modulus() * V1;
          for (int i = 0; i < 3; ++i) kd3[i] = c * Si[i];
        }
        double gall[4][8], dqall[12];
        for (int a = 0; a < 4; ++a) {
          // gamma_a = (h_a - (h_a . x_I) b_I) / 8 at end configuration
          double hx[3] = {0, 0, 0};
          for (int I = 0; I < 8; ++I) for (int j = 0; j < 3; ++j) hx[j] += HG[a][I] * x1[3 * I + j];
          double* g = gall[a];
          for (int I = 0; I < 8; ++I) g[I] = 0.125 * (HG[a][I] - (hx[0] * b1[I].x + hx[1] * b1[I].y + hx[2] * b1[I].z));
          Vec3 dq;
          for (int I = 0; I < 8; ++I) dq += load3(in.du + 6 * I) * g[I];
          store3(dqall + 3 * a, dq);
        }
        const double* Kq = hgKi[e] >= 0 ? &hgK[hgKi[e] * 144] : nullptr;
        for (int a = 0; a < 4; ++a) {
          const Vec3 dq = load3(dqall + 3 * a);
          Vec3 dQ(dq.x * kd3[0], dq.y * kd3[1], dq.z * kd3[2]);
          if (Kq) {  // enhanced control: coupled stiffness of the 12 hourglass modes
            double v[3] = {0, 0, 0};
            for (int i = 0; i < 3; ++i)
              for (int c = 0; c < 12; ++c) v[i] += Kq[(3 * a + i) * 12 + c] * dqall[c];
            dQ = Vec3(v[0], v[1], v[2]);
          }
          Vec3 qold = load3(qn + 3 * a), vold = load3(qvn + 3 * a), fold = load3(qfn + 3 * a);
          Vec3 qnew, vnew, fnew;
          if (Kq) {  // enhanced control: total forces of the coupled stiffness
            qnew = qold + dQ;
            vnew = vold;
            fnew = qnew;
          } else {
            // default: the relaxation acts on the generalized hourglass strains (q: total, q_v: relaxed part
            // following q with the time tau); force increments k (dq - r dq_v) at the current stiffness, scaled
            // during plastic flow. Same as relaxing the forces while k is constant; differs when the element
            // shape changes (single-element Abaqus/Explicit runs: compression to lambda = 0.5 together with
            // an hourglass mode, elastic and plastic, reproduced to 4 digits)
            qnew = qold + dq;
            vnew = vold * (1.0 - rb) + qold * rb + dq * ra;
            const Vec3 dv = vnew - vold;
            fnew = fold + (dQ - Vec3(dv.x * kd3[0], dv.y * kd3[1], dv.z * kd3[2]) * rr) * hg_flow;
          }
          out.e_art += 0.5 * dot(fold + fnew, dq);
          store3(qn + 3 * a, qnew);
          store3(qvn + 3 * a, vnew);
          store3(qfn + 3 * a, fnew);
          for (int I = 0; I < 8; ++I) add3(out.f + 6 * I, fnew * gall[a][I]);
        }
        if (commit) { std::memcpy(q, qn, sizeof(qn)); std::memcpy(qv, qvn, sizeof(qvn)); std::memcpy(qf, qfn, sizeof(qfn)); }
      }
      if (commit) st[e] = p;
      sig_avg = p.sig; eps_avg = p.eps; ep_max = p.ep; mis_max = mises(p.sig);
    } else if (type == ElType::C3D8I) {  // incompatible modes (Wilson-Taylor), full integration
      // gradients at the 2x2x2 points (middle and end configurations); incompatible displacement modes
      // alpha_k (1 - xi_k^2) per component, gradients with the centre Jacobian scaled by det0/det (patch test)
      Vec3 bm[8][8], b1[8][8], gi[8][3];
      double dm[8], d1[8];
      V1 = 0;
      for (int g = 0; g < 8; ++g) {
        double r = HEX_XI[g][0] * G8, s = HEX_XI[g][1] * G8, t = HEX_XI[g][2] * G8;
        dm[g] = hex_grad(xm, r, s, t, bm[g]);
        d1[g] = hex_grad(x1, r, s, t, b1[g]);
        V1 += d1[g];
      }
      Mat3 J0;
      {
        double dN[8][3];
        hex_dN(0, 0, 0, dN);
        for (int I = 0; I < 8; ++I)
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) J0.m[i][j] += dN[I][i] * xm[3 * I + j];
      }
      const double det0 = J0.det();
      const Mat3 Ji0 = J0.inverse();
      for (int g = 0; g < 8; ++g)
        for (int k = 0; k < 3; ++k) {
          Vec3 d;
          d[k] = -2.0 * HEX_XI[g][k] * G8;
          gi[g][k] = (Ji0 * d) * (det0 / dm[g]);
        }
      // elastic modulus matrix (engineering shear) for the condensation of the incompatible modes
      const double Gm = m.shear_modulus(), lam = std::max(0.0, M - 2 * Gm);
      auto eps_col = [](const Vec3& b, int c, double v[6]) {  // engineering strain of a unit displacement c
        for (int i = 0; i < 6; ++i) v[i] = 0;
        v[c] = b[c];
        if (c == 0) { v[3] = b.y; v[4] = b.z; }
        if (c == 1) { v[3] = b.x; v[5] = b.z; }
        if (c == 2) { v[4] = b.x; v[5] = b.y; }
      };
      auto Dmul = [&](const double v[6], double o[6]) {
        const double tr = v[0] + v[1] + v[2];
        for (int i = 0; i < 3; ++i) o[i] = lam * tr + 2 * Gm * v[i];
        for (int i = 3; i < 6; ++i) o[i] = Gm * v[i];
      };
      double Kaa[9][9] = {}, rhs[9] = {};
      for (int g = 0; g < 8; ++g) {
        double ga[9][6], dga[9][6];
        for (int k = 0; k < 3; ++k)
          for (int c = 0; c < 3; ++c) { eps_col(gi[g][k], c, ga[3 * k + c]); Dmul(ga[3 * k + c], dga[3 * k + c]); }
        double deps_u[6] = {0, 0, 0, 0, 0, 0}, v[6];
        for (int I = 0; I < 8; ++I)
          for (int c = 0; c < 3; ++c) {
            eps_col(bm[g][I], c, v);
            for (int i = 0; i < 6; ++i) deps_u[i] += v[i] * in.du[6 * I + c];
          }
        const Sym6& so = st[e * 8 + g].sig;
        const double sv[6] = {so[0], so[1], so[2], so[3], so[4], so[5]};
        for (int a = 0; a < 9; ++a) {
          double ra = 0, ku = 0;
          for (int i = 0; i < 6; ++i) { ra += ga[a][i] * sv[i]; ku += dga[a][i] * deps_u[i]; }
          rhs[a] -= (ra + ku) * dm[g];
          for (int b2 = 0; b2 < 9; ++b2) {
            double kab = 0;
            for (int i = 0; i < 6; ++i) kab += dga[a][i] * ga[b2][i];
            Kaa[a][b2] += kab * dm[g];
          }
        }
      }
      double dal[9];
      {  // Gauss elimination with partial pivoting
        double Mx[9][10];
        for (int i = 0; i < 9; ++i) { for (int j = 0; j < 9; ++j) Mx[i][j] = Kaa[i][j]; Mx[i][9] = rhs[i]; }
        for (int c0 = 0; c0 < 9; ++c0) {
          int pv = c0;
          for (int r = c0 + 1; r < 9; ++r) if (std::fabs(Mx[r][c0]) > std::fabs(Mx[pv][c0])) pv = r;
          for (int j = 0; j < 10; ++j) std::swap(Mx[c0][j], Mx[pv][j]);
          const double dd = Mx[c0][c0];
          if (std::fabs(dd) < 1e-300) continue;
          for (int r = 0; r < 9; ++r) {
            if (r == c0) continue;
            const double f = Mx[r][c0] / dd;
            for (int j = c0; j < 10; ++j) Mx[r][j] -= f * Mx[c0][j];
          }
        }
        for (int i = 0; i < 9; ++i) dal[i] = std::fabs(Mx[i][i]) > 1e-300 ? Mx[i][9] / Mx[i][i] : 0.0;
      }
      for (int g = 0; g < 8; ++g) {
        Mat3 L;
        for (int I = 0; I < 8; ++I)
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) L.m[i][j] += in.du[6 * I + i] * bm[g][I][j];
        for (int k = 0; k < 3; ++k)
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) L.m[i][j] += dal[3 * k + i] * gi[g][k][j];
        Sym6 deps = Sym6::from(L);
        Vec3 w(0.5 * (L.m[2][1] - L.m[1][2]), 0.5 * (L.m[0][2] - L.m[2][0]), 0.5 * (L.m[1][0] - L.m[0][1]));
        Mat3 dF = Mat3::identity();
        if (m.hyper())
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) dF.m[i][j] += L.m[i][j];
        PointState p = st[e * 8 + g];
        double wp = update_point(m, p, deps, w, dF, rotate, dtemp, fib.empty() ? nullptr : &fib[e * 9]);
        out.e_plastic += wp * d1[g];
        const Sym6& s = p.sig;
        for (int I = 0; I < 8; ++I) {
          const Vec3& b = b1[g][I];
          Vec3 f(s[0] * b.x + s[3] * b.y + s[4] * b.z, s[3] * b.x + s[1] * b.y + s[5] * b.z, s[4] * b.x + s[5] * b.y + s[2] * b.z);
          add3(out.f + 6 * I, f * d1[g]);
        }
        if (commit) st[e * 8 + g] = p;
        for (int i = 0; i < 6; ++i) { sig_avg[i] += p.sig[i] / 8; eps_avg[i] += p.eps[i] / 8; }
        ep_max = std::max(ep_max, p.ep);
        mis_max = std::max(mis_max, mises(p.sig));
        trd += deps.trace() * d1[g];
      }
      trd /= V1;
      Vec3 bbar[8];
      for (int g = 0; g < 8; ++g) for (int I = 0; I < 8; ++I) bbar[I] += b1[g][I] * (d1[g] / V1);
      {
        double bb = 0;
        for (int I = 0; I < 8; ++I) bb += norm2(bbar[I]);
        Le = 2.0 / std::sqrt(8 * bb);
      }
      if (in.dynamic && in.dt > 0 && (in.bv1 > 0 || in.bv2 > 0)) {
        double rate = trd / in.dt;
        double c = std::sqrt(M / rho);
        const double Lb = in.initial ? std::max(Le, Ldt) : Le;  // increment 0: stable-increment length, linear term
        double pbv = in.bv1 * rho * c * Lb * rate +
                     (in.initial ? 0.0 : rho * in.bv2 * in.bv2 * Le * Le * rate * std::fabs(std::min(0.0, rate)));
        for (int I = 0; I < 8; ++I) add3(out.f + 6 * I, bbar[I] * (pbv * V1));
        out.e_visc += pbv * trd * V1;
      }
    } else {  // C3D8 selective reduced (B-bar)
      Vec3 bm[8][8], b1[8][8];
      double dm[8], d1[8];
      Vec3 bbar_m[8], bbar_1[8];
      double Vm = 0;
      V1 = 0;
      for (int I = 0; I < 8; ++I) { bbar_m[I] = Vec3(); bbar_1[I] = Vec3(); }
      for (int g = 0; g < 8; ++g) {
        double r = HEX_XI[g][0] * G8, s = HEX_XI[g][1] * G8, t = HEX_XI[g][2] * G8;
        dm[g] = hex_grad(xm, r, s, t, bm[g]);
        d1[g] = hex_grad(x1, r, s, t, b1[g]);
        Vm += dm[g]; V1 += d1[g];
        for (int I = 0; I < 8; ++I) { bbar_m[I] += bm[g][I] * dm[g]; bbar_1[I] += b1[g][I] * d1[g]; }
      }
      for (int I = 0; I < 8; ++I) { bbar_m[I] = bbar_m[I] / Vm; bbar_1[I] = bbar_1[I] / V1; }
      {
        double bb = 0;
        for (int I = 0; I < 8; ++I) bb += norm2(bbar_1[I]);
        Le = 2.0 / std::sqrt(8 * bb);
      }
      double vol_bar = 0;
      for (int I = 0; I < 8; ++I) vol_bar += dot(bbar_m[I], load3(in.du + 6 * I));
      for (int g = 0; g < 8; ++g) {
        Mat3 L;
        for (int I = 0; I < 8; ++I)
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) L.m[i][j] += in.du[6 * I + i] * bm[g][I][j];
        Sym6 deps = Sym6::from(L);
        double corr = (vol_bar - deps.trace()) / 3.0;
        for (int i = 0; i < 3; ++i) deps[i] += corr;
        Vec3 w(0.5 * (L.m[2][1] - L.m[1][2]), 0.5 * (L.m[0][2] - L.m[2][0]), 0.5 * (L.m[1][0] - L.m[0][1]));
        Mat3 dF = Mat3::identity();
        if (m.hyper()) {
          double xo[24];
          for (int i = 0; i < 24; ++i) xo[i] = xref ? xref[i] : in.xo[i];
          Vec3 bo[8];
          hex_grad(xo, HEX_XI[g][0] * G8, HEX_XI[g][1] * G8, HEX_XI[g][2] * G8, bo);
          for (int I = 0; I < 8; ++I)
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) dF.m[i][j] += in.du[6 * I + i] * bo[I][j];
        }
        PointState p = st[e * 8 + g];
        double wp = update_point(m, p, deps, w, dF, rotate, dtemp, fib.empty() ? nullptr : &fib[e * 9]);
        out.e_plastic += wp * d1[g];
        const Sym6& s = p.sig;
        double trs = s.trace() / 3.0;
        for (int I = 0; I < 8; ++I) {
          const Vec3& b = b1[g][I];
          Vec3 f(s[0] * b.x + s[3] * b.y + s[4] * b.z, s[3] * b.x + s[1] * b.y + s[5] * b.z, s[4] * b.x + s[5] * b.y + s[2] * b.z);
          f += (bbar_1[I] - b) * trs;
          add3(out.f + 6 * I, f * d1[g]);
        }
        if (commit) st[e * 8 + g] = p;
        for (int i = 0; i < 6; ++i) { sig_avg[i] += p.sig[i] / 8; eps_avg[i] += p.eps[i] / 8; }
        ep_max = std::max(ep_max, p.ep);
        mis_max = std::max(mis_max, mises(p.sig));
      }
      trd = vol_bar;
      if (in.dynamic && in.dt > 0 && (in.bv1 > 0 || in.bv2 > 0)) {
        double rate = trd / in.dt;
        double c = std::sqrt(M / rho);
        const double Lb = in.initial ? std::max(Le, Ldt) : Le;  // increment 0: stable-increment length, linear term
        double pbv = in.bv1 * rho * c * Lb * rate +
                     (in.initial ? 0.0 : rho * in.bv2 * in.bv2 * Le * Le * rate * std::fabs(std::min(0.0, rate)));
        for (int I = 0; I < 8; ++I) add3(out.f + 6 * I, bbar_1[I] * (pbv * V1));
        out.e_visc += pbv * trd * V1;
      }
    }
    // stable increment
    Le = std::max(Le, Ldt);
    double c = std::sqrt(M / rho);
    double xi = 0;
    // linear term also at a zero increment (initial estimate, mass scaling): the damping is that of the element at
    // rest; the quadratic term needs a compression rate
    if (in.dynamic) xi = in.bv1 - (in.dt > 0 ? in.bv2 * in.bv2 * (Le / c) * std::min(0.0, trd / in.dt) : 0.0);
    if (xi < 0) xi = 0;
    out.dt_stable = Le / c * (std::sqrt(1 + xi * xi) - xi) * std::sqrt(1.0);
    if (commit) {
      float* o = &this->out[e * N_OUT];
      o[OUT_MISES] = (float)mis_max;
      for (int i = 0; i < 6; ++i) o[OUT_S11 + i] = (float)sig_avg[i];
      o[OUT_PEEQ] = (float)ep_max;
      double ev[3];
      eig3sym(eps_avg.full(), ev);
      o[OUT_LEMAX] = (float)std::max(ev[0], std::max(ev[1], ev[2]));
    }
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_solid_group() { return std::make_unique<SolidGroup>(); }

}  // namespace abq
