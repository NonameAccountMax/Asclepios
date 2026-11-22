// T3D2 trusses: axial force only, three translations per node. Hypoelastic axial stress on the logarithmic
// strain increment of the chord (current length; reference length without NLGEOM), uniaxial J2 plasticity with
// the isotropic hardening table, constant cross-section. Explicit: linear bulk viscosity b1 rho c L on the axial
// strain rate; stable increment L / c.
#include <cmath>

#include "element/groups.hpp"

namespace abq {

namespace {

class TrussGroup : public ElementGroup {
 public:
  std::vector<const Material*> mat;
  std::vector<double> area, L0, sig, ep;

  void setup(const Model& m) override {
    model = &m;
    nen = 2;
    rot = false;
    mat.resize(size()); area.resize(size()); L0.resize(size());
    sig.assign(size(), 0.0); ep.assign(size(), 0.0);
    for (int e = 0; e < size(); ++e) {
      const Section& s = m.sections[m.elems[ids[e]].section];
      mat[e] = &m.materials[s.material];
      area[e] = s.area;
      L0[e] = norm(m.X[conn[2 * e + 1]] - m.X[conn[2 * e]]);
    }
  }
  void save(BinOut& o) const override { ElementGroup::save(o); o.vec(sig); o.vec(ep); }
  void load(BinIn& in) override { ElementGroup::load(in); in.vec(sig); in.vec(ep); }
  void mass(int e, const double* X, double* mt, double* mr) const override {
    const double m = mat[e]->density * area[e] * norm(load3(X + 3) - load3(X)) * 0.5;
    mt[0] = mt[1] = m;
    mr[0] = mr[1] = 0;
  }
  double alpha_damping(int e) const override { return mat[e]->damp_alpha; }

  void compute(int e, const ElemIn& in, ElemOut& out, bool commit) override {
    for (int i = 0; i < 12; ++i) out.f[i] = 0;
    const Material& m = *mat[e];
    const Vec3 a0 = load3(in.xo), b0 = load3(in.xo + 3);
    const Vec3 a1 = a0 + load3(in.du), b1 = b0 + load3(in.du + 6);
    double L, dl;
    Vec3 n;
    if (in.nlgeom) {
      const double Lo = norm(b0 - a0);
      L = norm(b1 - a1);
      dl = std::log(L / std::max(Lo, 1e-300));
      n = (b1 - a1) / std::max(L, 1e-300);
    } else {
      const Vec3 X0 = load3(in.X), X1 = load3(in.X + 3);
      L = L0[e];
      n = (X1 - X0) / std::max(L, 1e-300);
      dl = dot(load3(in.du + 6) - load3(in.du), n) / L;
    }
    if (in.temp0 && in.temp1) dl -= m.alpha * 0.5 * (in.temp1[0] + in.temp1[1] - in.temp0[0] - in.temp0[1]);
    // elastic predictor, uniaxial return mapping
    double s = sig[e] + m.E * dl, p = ep[e];
    if (m.plastic()) {
      const double sy = m.yield(p);
      if (std::fabs(s) > sy) {
        double dp = 0;
        for (int it = 0; it < 50; ++it) {
          double H;
          const double syk = m.yield(p + dp, &H);
          const double r = std::fabs(s) - m.E * dp - syk;
          const double ddp = r / (m.E + H);
          dp += ddp;
          if (dp < 0) dp = 0;
          if (std::fabs(ddp) <= 1e-12 * (p + dp + 1e-12)) break;
        }
        s -= (s > 0 ? 1.0 : -1.0) * m.E * dp;
        out.e_plastic = m.yield(p + 0.5 * dp) * dp * area[e] * L;
        p += dp;
      }
    }
    const double c = std::sqrt(m.E / std::max(m.density * msf[e], 1e-300));
    double sv = 0;  // linear bulk viscosity on the axial strain rate
    if (in.dynamic && in.dt > 0 && in.bv1 > 0) {
      sv = in.bv1 * m.density * msf[e] * c * L * dl / in.dt;
      out.e_visc = sv * dl * area[e] * L;
    }
    const double N = (s + sv) * area[e];
    add3(out.f, n * (-N));
    add3(out.f + 3 * 2, n * N);
    out.dt_stable = L / c;
    if (commit) {
      sig[e] = s; ep[e] = p;
      float* o = &this->out[e * N_OUT];
      o[OUT_S11] = (float)s; o[OUT_MISES] = (float)std::fabs(s); o[OUT_PEEQ] = (float)p; o[OUT_SF] = (float)(s * area[e]);
    }
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_truss_group() { return std::make_unique<TrussGroup>(); }

}  // namespace abq
