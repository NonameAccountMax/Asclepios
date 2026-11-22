// Discrete elements: SPRINGA (axial spring along the current line of its nodes, linear or tabulated force against
// the change of length, linear beyond the table), DASHPOTA (axial damper on the relative velocity) and MASS (point
// mass). The springs report their tangent stiffness to the driver for the nodal frequency bound.
#include <cmath>

#include "element/groups.hpp"

namespace abq {

namespace {

class DiscreteGroup : public ElementGroup {
 public:
  std::vector<const Section*> sec;
  std::vector<double> L0;

  void setup(const Model& m) override {
    model = &m;
    rot = false;
    sec.resize(size());
    L0.assign(size(), 0.0);
    for (int e = 0; e < size(); ++e) {
      sec[e] = &m.sections[m.elems[ids[e]].section];
      if (nen == 2) L0[e] = norm(m.X[conn[2 * e + 1]] - m.X[conn[2 * e]]);
    }
  }
  void mass(int e, const double* X, double* mt, double* mr) const override {
    for (int I = 0; I < nen; ++I) { mt[I] = 0; mr[I] = 0; }
    if (type == ElType::MASS) mt[0] = sec[e]->point_mass;
  }

  static double table_force(const Section& s, double u, double* k) {
    const auto& U = s.spring_u;
    const auto& F = s.spring_F;
    if (U.size() < 2) { *k = 0; return 0; }
    size_t i = 1;
    while (i + 1 < U.size() && u > U[i]) ++i;  // segment [i-1, i], extended linearly beyond the ends
    *k = (F[i] - F[i - 1]) / (U[i] - U[i - 1]);
    return F[i - 1] + *k * (u - U[i - 1]);
  }

  void compute(int e, const ElemIn& in, ElemOut& out, bool commit) override {
    for (int i = 0; i < 6 * nen; ++i) out.f[i] = 0;
    out.dt_stable = 1e30;
    if (type == ElType::MASS) return;
    const Section& s = *sec[e];
    const Vec3 a = load3(in.xo) + load3(in.du), b = load3(in.xo + 3) + load3(in.du + 6);
    const Vec3 d = in.nlgeom ? b - a : load3(in.X + 3) - load3(in.X);
    const double L = norm(d);
    const Vec3 n = L > 0 ? d / L : Vec3(1, 0, 0);
    double F = 0, k = 0;
    if (type == ElType::SPRINGA) {
      const double u = in.nlgeom ? L - L0[e] : dot(b - a - (load3(in.X + 3) - load3(in.X)), n);
      if (!s.spring_F.empty()) F = table_force(s, u, &k);
      else { k = s.spring_k; F = k * u; }
      out.dt_stable = k > 0 ? -k : 1e30;  // spring-like: frequency bound from the nodal masses
    } else {  // DASHPOTA
      const double v = in.dt > 0 && !in.initial ? dot(load3(in.du + 6) - load3(in.du), n) / in.dt : 0.0;
      F = s.dashpot_c * v;
      out.e_visc = F * v * in.dt;
    }
    add3(out.f, n * (-F));
    add3(out.f + 6, n * F);
    if (commit) { float* o = &this->out[e * N_OUT]; o[OUT_SF] = (float)F; o[OUT_S11] = (float)F; }
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_discrete_group() { return std::make_unique<DiscreteGroup>(); }

}  // namespace abq
