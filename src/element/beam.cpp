// B31: two-node Timoshenko beam, corotational incremental formulation, one integration point
// along the axis (uniform shear, no locking), cross-section integration points for CIRC / PIPE /
// RECT (J2 plasticity in sigma11 + torsion shear, thermal expansion). Purely elastic sections use
// resultants. Section frame rotated with the nodes (shear along the section directions), transverse shear stiffness
// with Abaqus' slenderness compensation factor (0.25) of each bending plane.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

#include "element/groups.hpp"
#include "material/material.hpp"

namespace abq {

namespace {

struct SecPt { double y1, y2, w; };
constexpr int MAXSP = 32;

struct SectionProps {
  double A = 0, I11 = 0, I22 = 0, J = 0, Ip = 0;  // I11: about n1 (fibres at y2), I22: about n2
  double kappa = 0.85;  // RECT (free RECT elements in Abaqus/Explicit: transverse motions exact with 0.85)
  double rmax = 0;
  bool torsion_points = true;
  std::vector<SecPt> pts;
};

SectionProps section_props(const Section& s) {
  SectionProps p;
  if (s.shape == "CIRC" || s.shape == "PIPE") {
    double R = s.dims.at(0);
    const bool solid = s.shape == "CIRC";
    // PIPE is thin-walled (outer radius R, wall t): A = 2 pi r t, I = pi r^3 t, J = 2 I on the mean radius r,
    // section points on the mean radius (Abaqus/Standard cantilevers: bending and torsion stiffness, stresses)
    const double t = solid ? 0.0 : s.dims.at(1), rm = R - 0.5 * t;
    p.A = solid ? M_PI * R * R : 2 * M_PI * rm * t;
    p.I11 = p.I22 = solid ? M_PI / 4 * std::pow(R, 4) : M_PI * rm * rm * rm * t;
    p.J = p.Ip = 2 * p.I11;
    p.kappa = solid ? 0.89 : 0.53;
    p.rmax = R;
    // Abaqus default section points (beam in space): 8 circumferential points starting on the section axes;
    // radially, CIRC: the centre (output only, no weight, not stored), R/2 and R, with the weights that make the
    // area and the second moment exact (R^2/3 and R^2/6 per radian); PIPE: the mean radius
    const int nth = 8;
    const double rr[3] = {0.0, solid ? 0.5 * R : rm, R};
    for (int i = 0; i < 3; ++i) {
      const double r = rr[i];
      const double w = solid ? (i == 1 ? R * R / 3.0 : (i == 2 ? R * R / 6.0 : 0.0)) : (i == 1 ? rm * t : 0.0);
      if (w <= 0) continue;
      for (int k = 0; k < nth; ++k) {
        const double th = 2 * M_PI * k / nth;
        p.pts.push_back({r * std::cos(th), r * std::sin(th), w * 2 * M_PI / nth});
      }
    }
  } else {  // RECT: a along n1, b along n2
    double a = s.dims.at(0), b = s.dims.at(1);
    p.A = a * b;
    p.I11 = a * b * b * b / 12.0;
    p.I22 = b * a * a * a / 12.0;
    double lo = std::max(a, b), sh = std::min(a, b), r = sh / lo;
    p.J = lo * sh * sh * sh * (1.0 / 3.0 - 0.21 * r * (1 - r * r * r * r / 12.0));
    p.Ip = p.I11 + p.I22;
    p.rmax = 0.5 * std::sqrt(a * a + b * b);
    p.torsion_points = false;
    const int n = 5;
    const double w1[n] = {1, 4, 2, 4, 1};
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        p.pts.push_back({-0.5 * a + a * i / (n - 1), -0.5 * b + b * j / (n - 1),
                         (a / (3.0 * (n - 1))) * w1[i] * (b / (3.0 * (n - 1))) * w1[j]});
  }
  return p;
}

struct BeamState {
  double n1[3];                  // section frame (axis, n1, n2 = axis x n1), rotated by the nodal rotations: the
  double ax[3];                  // section axis departs from the chord by the transverse shear strains
  double e0;                     // integrated axial strain
  double N, M1, M2, T, V1, V2;   // resultants
};

class BeamGroup : public ElementGroup {
 public:
  std::vector<SectionProps> props;
  std::vector<int> pidx;
  std::vector<BeamState> st;
  std::vector<double> pts;   // per element: npt x (s11, t12, t13, ep)
  std::vector<int> poff;
  std::vector<const Material*> mat;
  std::vector<const Section*> sec;
  std::vector<double> irot, L0;

  double poisson(int e) const { return sec[e]->poisson > 0 ? sec[e]->poisson : mat[e]->nu; }

  // transverse shear stiffness with the slenderness compensation of the bending plane it belongs to: I is the
  // second moment of the coupled bending (I22 for shear along n1, I11 for shear along n2; free RECT elements in
  // Abaqus/Explicit: transverse motions of each plane follow their own factor)
  double shear_stiff(int e, double L, double G, double I) const {
    const SectionProps& P = props[pidx[e]];
    double fp = 1.0 / (1.0 + 0.25 * L * L * P.A / (12.0 * I));
    return fp * P.kappa * G * P.A;
  }

  void setup(const Model& m) override {
    model = &m;
    rot = true;
    std::vector<int> sec_to_prop(m.sections.size(), -1);
    pidx.resize(size()); mat.resize(size()); sec.resize(size()); st.assign(size(), BeamState());
    irot.assign(size(), 0); L0.assign(size(), 0); poff.assign(size() + 1, 0);
    for (int e = 0; e < size(); ++e) {
      int si = m.elems[ids[e]].section;
      const Section& s = m.sections[si];
      sec[e] = &s;
      if (sec_to_prop[si] < 0) { sec_to_prop[si] = (int)props.size(); props.push_back(section_props(s)); }
      pidx[e] = sec_to_prop[si];
      mat[e] = &m.materials[s.material];
      if (mat[e]->hgo) throw std::runtime_error("Holzapfel-Gasser-Ogden material " + mat[e]->name + ": solid elements only");
      if (!mat[e]->elastic) throw std::runtime_error("beam elements need an *ELASTIC material (" + mat[e]->name + ")");
      const SectionProps& P = props[pidx[e]];
      Vec3 a = m.X[conn[2 * e]], b = m.X[conn[2 * e + 1]];
      double L = norm(b - a);
      if (L <= 0) throw std::runtime_error("beam element " + m.elem_name(ids[e]) + " has zero length");
      L0[e] = L;
      Vec3 t = (b - a) / L;
      Vec3 n1 = s.n1 - t * dot(s.n1, t);
      if (norm(n1) < 1e-6) {
        Vec3 g = std::fabs(t.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
        n1 = g - t * dot(g, t);
      }
      n1 = normalized(n1);
      BeamState& S = st[e];
      std::memset(&S, 0, sizeof(S));
      store3(S.n1, n1);
      store3(S.ax, t);
      bool points = mat[e]->plastic() || mat[e]->alpha != 0;
      poff[e + 1] = poff[e] + (points ? (int)P.pts.size() : 0);
      const double rho = mat[e]->density;
      // nodal rotary inertia. Default (Abaqus EXACT rotary inertia of Timoshenko beams): the lumped inertia of the
      // half element, anisotropic and carried by the beam frame, rho L / 2 (I11 + I22, I11, I22) about (axis, n1,
      // n2) (rotary_tensor; free single-element Abaqus/Explicit runs: torsional rotations exact to 2e-4 for CIRC
      // and PIPE). ROTARY INERTIA=ISOTROPIC: rho L / 2 Imax max((2 + 3u) / (2 + u), 6u / (1 + u)) about every
      // axis, u = A L^2 / (24 Imin); independent of the elastic constants and of the density (single-element
      // Abaqus/Explicit runs: 3 sections, L / r from 0.5 to 50, three moduli: 1e-5)
      const double Imin = std::min(P.I11, P.I22), Imax = std::max(P.I11, P.I22);
      const double u = P.A * L * L / (24.0 * Imin);
      const double Iiso = rho * L / 2 * Imax * std::max((2 + 3 * u) / (2 + u), 6 * u / (1 + u));
      irot[e] = s.rotary_isotropic ? Iiso : 0.0;
    }
    pts.assign((size_t)poff[size()] * 4, 0.0);
  }

  void save(BinOut& o) const override { ElementGroup::save(o); o.vec(st); o.vec(pts); o.vec(irot); }
  void load(BinIn& in) override { ElementGroup::load(in); in.vec(st); in.vec(pts); in.vec(irot); }

  double alpha_damping(int e) const override { return mat[e]->damp_alpha; }
  double contact_thickness(int e) const override { return 2 * props[pidx[e]].rmax; }

  void mass(int e, const double* X, double* mt, double* mr) const override {
    double m = mat[e]->density * props[pidx[e]].A * L0[e] / 2;
    mt[0] = mt[1] = m;
    mr[0] = mr[1] = irot[e];
  }

  bool rotary_tensor(int e, const double* x, double* J6) const override {
    if (sec[e]->rotary_isotropic) return false;
    const SectionProps& P = props[pidx[e]];
    (void)x;
    const Vec3 t = load3(st[e].ax), n1 = load3(st[e].n1), n2 = cross(t, n1);
    const double c = mat[e]->density * L0[e] / 2, a = c * (P.I11 + P.I22), b1 = c * P.I11, b2 = c * P.I22;
    auto add = [&](const Vec3& d, double s) {
      J6[0] += s * d.x * d.x; J6[1] += s * d.y * d.y; J6[2] += s * d.z * d.z;
      J6[3] += s * d.x * d.y; J6[4] += s * d.x * d.z; J6[5] += s * d.y * d.z;
    };
    for (int k = 0; k < 6; ++k) J6[k] = 0;
    add(t, a); add(n1, b1); add(n2, b2);
    return true;
  }

  void compute(int e, const ElemIn& in, ElemOut& out, bool commit) override {
    const Material& m = *mat[e];
    const SectionProps& P = props[pidx[e]];
    BeamState S = st[e];
    const double E = m.E, G = E / (2 * (1 + poisson(e)));
    const double rho = m.density * msf[e];
    const bool nl = in.nlgeom;
    Vec3 xo1 = load3(nl ? in.xo : in.X), xo2 = load3(nl ? in.xo + 3 : in.X + 3);
    Vec3 du1 = load3(in.du), du2 = load3(in.du + 6), dr1 = load3(in.du + 3), dr2 = load3(in.du + 9);
    Vec3 x1a = nl ? xo1 + du1 : xo1, x1b = nl ? xo2 + du2 : xo2;
    Vec3 xma = nl ? xo1 + du1 * 0.5 : xo1, xmb = nl ? xo2 + du2 * 0.5 : xo2;
    double Lo = norm(xo2 - xo1), Lm = norm(xmb - xma), L1 = norm(x1b - x1a);
    Vec3 tm = (xmb - xma) / Lm, t1 = (x1b - x1a) / L1;
    Vec3 n1o = load3(S.n1), axo = load3(S.ax);
    Vec3 rbar = (dr1 + dr2) * 0.5;
    // section frame rotated by the mean nodal rotation (at mid-increment and at the end)
    auto frame = [&](const Vec3& rv, Vec3& a, Vec3& n1, Vec3& n2) {
      a = axo; n1 = n1o;
      if (nl) { const Mat3 R = rotation_from_vector(rv); a = R * axo; n1 = R * n1o; }
      a = normalized(a);
      n1 = normalized(n1 - a * dot(n1, a));
      n2 = cross(a, n1);
    };
    Vec3 axm, n1m, n2m, ax1, n11, n21;
    frame(rbar * 0.5, axm, n1m, n2m);
    frame(rbar, ax1, n11, n21);
    Vec3 dd = du2 - du1, dr = dr2 - dr1;
    double de0 = nl ? std::log(L1 / Lo) : dot(dd, tm) / Lm;
    double dk1 = dot(dr, n1m) / Lm, dk2 = dot(dr, n2m) / Lm, dphi = dot(dr, axm) / Lm;
    // transverse shear strains gamma_a = n_a . (x2 - x1) / L on the rotated section directions
    double dg1 = dot(dd, n1m) / Lm + dot(rbar, cross(n1m, tm)), dg2 = dot(dd, n2m) / Lm + dot(rbar, cross(n2m, tm));
    double dtemp = 0;
    if (in.temp0 && in.temp1) dtemp = 0.5 * (in.temp1[0] + in.temp1[1] - in.temp0[0] - in.temp0[1]);
    S.e0 += de0;
    const double Ks1 = shear_stiff(e, L1, G, P.I22), Ks2 = shear_stiff(e, L1, G, P.I11);
    S.V1 += Ks1 * dg1;
    S.V2 += Ks2 * dg2;
    double N, M1, M2, T;
    double mis = 0, epm = 0, smax = 0;
    const int np = poff[e + 1] - poff[e];
    double loc[MAXSP * 4];
    if (np == 0) {  // elastic resultants
      S.N += E * P.A * de0;
      S.M1 += E * P.I11 * dk1;
      S.M2 += E * P.I22 * dk2;
      S.T += G * P.J * dphi;
      N = S.N; M1 = S.M1; M2 = S.M2; T = S.T;
      double r = P.rmax;
      smax = std::fabs(N) / P.A + std::sqrt(M1 * M1 / (P.I11 * P.I11) + M2 * M2 / (P.I22 * P.I22)) * r;
      mis = std::sqrt(smax * smax + 3 * std::pow(T * r / std::max(P.J, 1e-300), 2));
    } else {
      std::memcpy(loc, &pts[(size_t)poff[e] * 4], sizeof(double) * 4 * np);
      N = M1 = M2 = T = 0;
      for (int k = 0; k < np; ++k) {
        const SecPt& q = P.pts[k];
        double* s = loc + 4 * k;
        double d11 = de0 - q.y1 * dk2 + q.y2 * dk1 - m.alpha * dtemp;
        double g12 = P.torsion_points ? -q.y2 * dphi : 0.0, g13 = P.torsion_points ? q.y1 * dphi : 0.0;
        const double so[3] = {s[0], s[1], s[2]};
        double w = j2_update_beam(m, G, d11, g12, g13, s[0], s[1], s[2], s[3]);
        out.e_plastic += w * q.w * L1;
        {  // internal energy at the point: elastic energy change plus dissipation (thermal strain work removed)
          auto We = [&](double a, double b, double c) { return a * a / (2 * E) + (b * b + c * c) / (2 * G); };
          const double trap = 0.5 * ((so[0] + s[0]) * (d11 + m.alpha * dtemp) + (so[1] + s[1]) * g12 + (so[2] + s[2]) * g13);
          out.e_corr += (We(s[0], s[1], s[2]) - We(so[0], so[1], so[2]) + w - trap) * q.w * L1;
        }
        N += q.w * s[0];
        M1 += q.w * s[0] * q.y2;
        M2 -= q.w * s[0] * q.y1;
        T += q.w * (q.y1 * s[2] - q.y2 * s[1]);
        mis = std::max(mis, std::sqrt(s[0] * s[0] + 3 * (s[1] * s[1] + s[2] * s[2])));
        epm = std::max(epm, s[3]);
        smax = std::max(smax, std::fabs(s[0]));
      }
      if (!P.torsion_points) { S.T += G * P.J * dphi; T = S.T; }
      S.N = N; S.M1 = M1; S.M2 = M2;
    }
    const double c = std::sqrt(E / rho);
    const double N_out = N;  // section force output without the bulk viscosity (as Abaqus/Explicit SF)
    if (in.dynamic && in.dt > 0 && in.bv1 > 0) {
      // linear bulk viscosity: b1 rho c L on the axial strain rate, half of it on the curvature rates (moments
      // 0.5 b1 rho c L I k'), none on transverse shear and torsion (single-element Abaqus/Explicit runs, viscous
      // dissipation of each mode: linear in b1, independent of the Poisson ratio)
      const double bv = in.bv1 * rho * c * L1 / in.dt;
      const double pbv = bv * de0, mb1 = 0.5 * bv * P.I11 * dk1, mb2 = 0.5 * bv * P.I22 * dk2;
      N += pbv * P.A;
      M1 += mb1;
      M2 += mb2;
      out.e_visc += (pbv * P.A * de0 + mb1 * dk1 + mb2 * dk2) * L1;
    }
    // axial force along the chord, shear forces along the section directions (their couple balanced by the
    // nodal moments)
    const Vec3 Fs = n11 * S.V1 + n21 * S.V2;
    Vec3 F = t1 * N + Fs;
    Vec3 Mm = n11 * M1 + n21 * M2 + ax1 * T;
    Vec3 Ms = cross(Fs, t1) * (0.5 * L1);
    store3(out.f, -F);
    store3(out.f + 6, F);
    store3(out.f + 3, -Mm + Ms);
    store3(out.f + 9, Mm + Ms);
    // stable increment, as estimated by Abaqus/Explicit (free single elements: CIRC, PIPE, RECT, L / r from 0.5 to
    // 50, two Poisson ratios, prescribed transverse shear stiffness, with and without bulk viscosity: 2e-3 on the
    // increment). For the rotary inertia Ir = rho L / 2 Ib about each bending axis, the larger of two modes, scaled
    // by 1.118:
    //  - opposite nodal rotations: bending stiffness of the OTHER axis plus the shear energy of the linear rotation
    //    field (integrated exactly) with the shear stiffness compensated for axis b, (2 E Ia / L + Ks(Ib) L / 6) / Ir;
    //  - equal rotations with opposite deflections, shear stiffness compensated for the other axis,
    //    Ks(Ia) L (2 / (L^2 m) + 1 / (2 Ir)), m the half element mass.
    // The cross pairing (not the coupling of the element's own modes) is what the estimate does: RECT sections with
    // the physical pairing predict increments up to 7 times too large. The axial wave (0.99 L / c) only governs with
    // ISOTROPIC rotary inertia, where the rotations never limit. The bulk viscosity factor applies to every mode.
    const double xi = in.dynamic ? in.bv1 : 0.0;
    double dt_el = 0.99 * L1 / c;
    if (!sec[e]->rotary_isotropic) {
      const double mt = rho * P.A * L0[e] / 2, cI = rho * L0[e] / 2;
      auto rot_mode = [&](double Ia, double Ib, double Ksa, double Ksb) {
        const double Ir = cI * Ib;
        return std::max((2 * E * Ia / L1 + Ksb * L1 / 6) / Ir, Ksa * L1 * (2.0 / (L1 * L1 * mt) + 1.0 / (2 * Ir)));
      };
      const double om_t = 2 * G * P.J / (L1 * cI * (P.I11 + P.I22));  // torsion (never governs for these sections)
      // Ks1 is compensated with I22, Ks2 with I11
      const double om2 = 1.118 * std::max({rot_mode(P.I22, P.I11, Ks1, Ks2), rot_mode(P.I11, P.I22, Ks2, Ks1), om_t});
      dt_el = std::min(dt_el, 2.0 / std::sqrt(om2));
    }
    out.dt_stable = dt_el * (std::sqrt(1 + xi * xi) - xi);
    if (commit) {
      store3(S.n1, n11);
      store3(S.ax, ax1);
      st[e] = S;
      if (np) std::memcpy(&pts[(size_t)poff[e] * 4], loc, sizeof(double) * 4 * np);
      float* o = &this->out[e * N_OUT];
      o[OUT_MISES] = (float)mis;
      o[OUT_S11] = (float)smax;
      o[OUT_PEEQ] = (float)epm;
      o[OUT_LEMAX] = (float)S.e0;
      o[OUT_SF] = (float)N_out;
      o[OUT_STH] = (float)(2 * P.rmax);
    }
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_beam_group() { return std::make_unique<BeamGroup>(); }

}  // namespace abq
