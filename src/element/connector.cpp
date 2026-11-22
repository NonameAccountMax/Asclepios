// CONN3D2 connectors: AXIAL (force along the current line between the nodes, relative to a
// constitutive reference length; linear damping on the rate of change of length), CARTESIAN (relative displacement in a local system attached to node a) and BUSHING (CARTESIAN plus
// the relative rotation in the same axes, components 4-6). Elastic (linear or tabular, field-variable dependent)
// behaviour. JOIN, HINGE and BEAM are constraints (links.cpp).
#include <cmath>

#include "element/groups.hpp"

namespace abq {

namespace {

class ConnectorGroup : public ElementGroup {
 public:
  std::vector<const Section*> sec;
  std::vector<double> L0;      // initial length
  std::vector<double> frame;   // CARTESIAN: 9 per element (local axes rows), rotated with node a
  std::vector<double> force;   // last force components (3)
  std::vector<double> relrot;  // BUSHING: relative rotation in the local axes (3)
  std::vector<double> rota;    // CARTESIAN / BUSHING without NLGEOM: accumulated rotation of node a (3)

  void setup(const Model& m) override {
    model = &m;
    rot = true;  // CARTESIAN frames follow node a rotation
    sec.resize(size());
    L0.resize(size());
    frame.assign(size() * 9, 0.0);
    force.assign(size() * 3, 0.0);
    relrot.assign(size() * 3, 0.0);
    rota.assign(size() * 3, 0.0);
    for (int e = 0; e < size(); ++e) {
      const Section& s = m.sections[m.elems[ids[e]].section];
      sec[e] = &s;
      Vec3 a = m.X[conn[2 * e]], b = m.X[conn[2 * e + 1]];
      L0[e] = norm(b - a);
      Mat3 R = s.orientation >= 0 ? m.orients[s.orientation].R : Mat3::identity();
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) frame[e * 9 + 3 * i + j] = R.m[i][j];
    }
  }

  void save(BinOut& o) const override { ElementGroup::save(o); o.vec(frame); o.vec(force); o.vec(relrot); o.vec(rota); }
  void load(BinIn& in) override { ElementGroup::load(in); in.vec(frame); in.vec(force); in.vec(relrot); in.vec(rota); }
  void mass(int e, const double* X, double* mt, double* mr) const override {
    mt[0] = mt[1] = 0;
    mr[0] = mr[1] = 0;
  }

  void compute(int e, const ElemIn& in, ElemOut& out, bool commit) override {
    for (int i = 0; i < 12; ++i) out.f[i] = 0;
    const Section& s = *sec[e];
    if (s.behavior < 0) return;
    const ConnectorBehavior& B = model->conn_behaviors[s.behavior];
    Vec3 xa = load3(in.xo) + load3(in.du), xb = load3(in.xo + 3) + load3(in.du + 6);
    double fv = in.fv1 ? 0.5 * (in.fv1[0] + in.fv1[1]) : 0.0;
    double kmax = 0;
    float* o = &this->out[e * N_OUT];
    if (s.ctype == Section::AXIAL) {
      Vec3 d = xb - xa;
      double L = norm(d);
      Vec3 n = L > 0 ? d / L : Vec3(1, 0, 0);
      const ConnectorComponent& c = B.comp[0];
      double ref = c.has_ref ? c.ref : L0[e];
      double k = 0;
      double F = c.elastic ? c.force(L - ref, fv, &k) : 0.0;
      kmax = std::fabs(k);
      if (c.damping != 0 && in.dt > 0 && !in.initial) {  // viscous force on the rate of change of length
        const double Lo = norm(load3(in.xo + 3) - load3(in.xo));
        const double Fd = c.damping * (L - Lo) / in.dt;
        F += Fd;
        out.e_visc = Fd * (L - Lo);
      }

      add3(out.f, n * (-F));
      add3(out.f + 6, n * F);
      if (commit) { o[OUT_SF] = (float)F; o[OUT_S11] = (float)(L - ref); }
    } else if (s.ctype == Section::CARTESIAN || s.ctype == Section::BUSHING) {
      double* R = &frame[e * 9];
      Mat3 Rn;
      for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Rn.m[i][j] = R[3 * i + j];
      if (in.nlgeom) {
        Mat3 Q = rotation_from_vector(load3(in.du + 3));
        Rn = Rn * Q.T();  // rows rotated: e_i' = Q e_i
      }
      // position of node b in the axes that follow node a, less its initial position: with NLGEOM the axes rotate
      // with node a, otherwise u_b - u_a - theta_a x d (small rotations, theta_a accumulated)
      const Vec3 d0 = load3(in.X + 3) - load3(in.X);
      Vec3 ul;
      if (in.nlgeom) {
        const Mat3 R0 = s.orientation >= 0 ? model->orients[s.orientation].R : Mat3::identity();
        ul = Rn * (xb - xa) - R0 * d0;
      } else {
        const Vec3 tha = load3(&rota[3 * e]) + load3(in.du + 3);
        ul = Rn * ((xb - xa) - d0 - cross(tha, d0));
      }
      Vec3 fl;
      for (int i = 0; i < 3; ++i) {
        const ConnectorComponent& c = B.comp[i];
        if (!c.elastic) continue;
        double k = 0;
        double ref = c.has_ref ? c.ref : 0.0;
        fl[i] = c.force(ul[i] - ref, fv, &k);
        kmax = std::max(kmax, std::fabs(k));
      }
      Vec3 fg = Rn.Tmul(fl);
      add3(out.f, -fg);
      add3(out.f + 6, fg);
      add3(out.f + 3, -cross(xb - xa, fg));  // moment of the force pair across the offset, on node a (equilibrium)
      if (commit && !in.nlgeom) add3(&rota[3 * e], load3(in.du + 3));
      if (s.ctype == Section::BUSHING) {  // rotational components 4-6 on the relative rotation in the local axes
        double* rr = &relrot[e * 3];
        const Vec3 dphi = Rn * (load3(in.du + 9) - load3(in.du + 3));
        Vec3 ml;
        for (int i = 0; i < 3; ++i) {
          const ConnectorComponent& c = B.comp[3 + i];
          if (!c.elastic) continue;
          double k = 0;
          ml[i] = c.force(rr[i] + dphi[i], fv, &k);
          kmax = std::max(kmax, std::fabs(k));
        }
        const Vec3 mg = Rn.Tmul(ml);
        add3(out.f + 3, -mg);
        add3(out.f + 9, mg);
        if (commit) for (int i = 0; i < 3; ++i) rr[i] += dphi[i];
      }
      if (commit) {
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) R[3 * i + j] = Rn.m[i][j];
        o[OUT_SF] = (float)norm(fl);
        o[OUT_S11] = (float)ul.x; o[OUT_S22] = (float)ul.y; o[OUT_S33] = (float)ul.z;
      }
    }
    // stable increment of a spring between the two nodal masses is handled by the driver
    // (needs nodal masses): store the stiffness in dt_stable as a negative marker
    out.dt_stable = kmax > 0 ? -kmax : 1e30;
  }
};

}  // namespace

std::unique_ptr<ElementGroup> make_connector_group() { return std::make_unique<ConnectorGroup>(); }

}  // namespace abq
