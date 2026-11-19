// B31: two-node Timoshenko beam, corotational incremental formulation, one integration point
// along the axis (uniform shear, no locking), cross-section integration points for CIRC / PIPE /
// RECT (J2 plasticity in sigma11 + torsion shear, thermal expansion). Purely elastic sections use
// resultants. Transverse shear stiffness with Abaqus' slenderness compensation factor (0.25).
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
  double kappa = 5.0 / 6.0;
  double rmax = 0;
  bool torsion_points = true;
  std::vector<SecPt> pts;
};

SectionProps section_props(const Section& s) {
  SectionProps p;
  if (s.shape == "CIRC" || s.shape == "PIPE") {
    double R = s.dims.at(0);
    double Ri = s.shape == "PIPE" ? R - s.dims.at(1) : 0.0;
    p.A = M_PI * (R * R - Ri * Ri);
    p.I11 = p.I22 = M_PI / 4 * (std::pow(R, 4) - std::pow(Ri, 4));
    p.J = p.Ip = 2 * p.I11;
    p.kappa = s.shape == "CIRC" ? 0.89 : 0.53;
    p.rmax = R;
    // Abaqus default section points (beam in space): 8 circumferential points starting on the section axes;
    // radially, CIRC: the centre (output only, no weight, not stored), R/2 and R, with the weights that make the
    // area and the second moment exact (R^2/3 and R^2/6 per radian); PIPE: Simpson's rule at the inner, middle
    // and outer radii (exact for both)
    const int nth = 8;
    const bool solid = s.shape == "CIRC";
    const double h = R - Ri;
    const double rr[3] = {Ri, 0.5 * (R + Ri), R};
    for (int i = 0; i < 3; ++i) {
      const double r = rr[i];
      const double w = solid ? (i == 1 ? R * R / 3.0 : (i == 2 ? R * R / 6.0 : 0.0)) : h / 6.0 * (i == 1 ? 4.0 : 1.0) * r;
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
  double n1[3];
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

  double shear_stiff(int e, double L, double G) const {
    const SectionProps& P = props[pidx[e]];
    double I = std::min(P.I11, P.I22);
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
      bool points = mat[e]->plastic() || mat[e]->alpha != 0;
      poff[e + 1] = poff[e] + (points ? (int)P.pts.size() : 0);
      const double E = mat[e]->E, G = E / (2 * (1 + poisson(e))), rho = mat[e]->density;
      double Ks = shear_stiff(e, L, G);
      double om_ax2 = 4 * E / (rho * L * L);
      double Ineed = 1.2 * (2 * std::max(E * std::max(P.I11, P.I22), G * P.J) / L + Ks * L / 2) / om_ax2;
      // nodal rotary inertia: by default the physical flexural inertia rho I L / 2 (Abaqus default EXACT rotary
      // inertia of Timoshenko beams, here isotropic with the larger flexural moment); with ROTARY INERTIA=ISOTROPIC
      // a scaled inertia chosen so that the rotations do not limit the stable increment
      // (isotropic: the smallest scaling of the flexural inertia that keeps the rotational bound at the axial one)
      const double Iflex = rho * L / 2 * std::max(P.I11, P.I22);
      irot[e] = s.rotary_isotropic ? std::max(Iflex, Ineed / 1.2) : Iflex;
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
    Vec3 n1o = load3(S.n1);
    Vec3 rbar = (dr1 + dr2) * 0.5;
    auto frame = [&](const Vec3& t, const Vec3& rv, Vec3& n1, Vec3& n2) {
      Vec3 v = nl ? rotation_from_vector(rv) * n1o : n1o;
      v -= t * dot(v, t);
      n1 = normalized(v);
      n2 = cross(t, n1);
    };
    Vec3 n1m, n2m, n11, n21;
    frame(tm, rbar * 0.5, n1m, n2m);
    frame(t1, rbar, n11, n21);
    Vec3 dd = du2 - du1, dr = dr2 - dr1;
    double de0 = nl ? std::log(L1 / Lo) : dot(dd, tm) / Lm;
    double dk1 = dot(dr, n1m) / Lm, dk2 = dot(dr, n2m) / Lm, dphi = dot(dr, tm) / Lm;
    double dg1 = dot(dd, n1m) / Lm - dot(rbar, n2m), dg2 = dot(dd, n2m) / Lm + dot(rbar, n1m);
    double dtemp = 0;
    if (in.temp0 && in.temp1) dtemp = 0.5 * (in.temp1[0] + in.temp1[1] - in.temp0[0] - in.temp0[1]);
    S.e0 += de0;
    const double Ks = shear_stiff(e, L1, G);
    S.V1 += Ks * dg1;
    S.V2 += Ks * dg2;
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
        double w = j2_update_beam(m, G, d11, g12, g13, s[0], s[1], s[2], s[3]);
        out.e_plastic += w * q.w * L1;
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
    if (in.dynamic && in.dt > 0 && in.bv1 > 0) {  // linear bulk viscosity on the axial rate
      double pbv = in.bv1 * rho * c * L1 * de0 / in.dt;
      N += pbv * P.A;
      out.e_visc += pbv * P.A * de0 * L1;
    }
    Vec3 F = t1 * N + n11 * S.V1 + n21 * S.V2;
    Vec3 Mm = n11 * M1 + n21 * M2 + t1 * T;
    Vec3 Ms = (n21 * (-0.5 * S.V1) + n11 * (0.5 * S.V2)) * L1;
    store3(out.f, -F);
    store3(out.f + 6, F);
    store3(out.f + 3, -Mm + Ms);
    store3(out.f + 9, Mm + Ms);
    // stable increment: axial wave and rotational (bending + coupled shear) bound
    double xi = in.dynamic ? in.bv1 : 0.0;
    double dt_ax = L1 / c * (std::sqrt(1 + xi * xi) - xi);
    double In = irot[e] * msf[e], mt = rho * P.A * L0[e] / 2;
    double om_r2 = 2 * std::max(E * std::max(P.I11, P.I22), G * P.J) / (L1 * In) + Ks * L1 * (2.0 / (L1 * L1 * mt) + 1.0 / (2 * In));
    out.dt_stable = std::min(dt_ax, 2.0 / std::sqrt(om_r2));
    if (commit) {
      store3(S.n1, n11);
      st[e] = S;
      if (np) std::memcpy(&pts[(size_t)poff[e] * 4], loc, sizeof(double) * 4 * np);
      float* o = &this->out[e * N_OUT];
      o[OUT_MISES] = (float)mis;
      o[OUT_S11] = (float)smax;
      o[OUT_PEEQ] = (float)epm;
      o[OUT_LEMAX] = (float)S.e0;
      o[OUT_SF] = (float)N;
      o[OUT_STH] = (float)(2 * P.rmax);
    }
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_beam_group() { return std::make_unique<BeamGroup>(); }

}  // namespace abq
