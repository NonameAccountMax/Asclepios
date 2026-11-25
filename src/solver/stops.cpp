// Connector stops in explicit dynamics (*CONNECTOR STOP), kinematic: a stop is an inequality constraint on the
// relative position of the connector (AXIAL: its length). At the end of each increment, after the prediction of the
// next one, a connector whose predicted length leaves the admissible range is brought back to the limit by
// corrections of the accelerations of its two nodes along its line, shared by their masses so that momentum is
// conserved (a prescribed or massless node does not move); the arrival on a stop is a plastic impact, as in
// Abaqus/Explicit. The corrections are returned as contact-like forces (mass x correction; the force taken by a
// prescribed node as its reaction).
#include <cmath>

#include "core/log.hpp"
#include "solver/solvers.hpp"

namespace abq {

bool has_connector_stops(const Model& m) {
  for (const auto& el : m.elems) {
    if (el.type != ElType::CONN3D2 || el.section < 0 || el.nodes.size() != 2) continue;
    const Section& s = m.sections[el.section];
    if (s.behavior < 0) continue;
    const ConnectorBehavior& B = m.conn_behaviors[s.behavior];
    for (int k = 0; k < 6; ++k)
      if (B.comp[k].stop) {
        if (s.ctype == Section::AXIAL && k == 0) return true;
        static bool once = false;
        if (!once) warn("connector stops are supported on the length of AXIAL connectors only: other stops ignored");
        once = true;
      }
  }
  return false;
}

void connector_stops(const Analysis& A, double c, std::vector<double>& dpred, std::vector<double>& a,
                     std::vector<double>& fk) {
  const Model& m = A.m;
  auto X = [&](int n) {
    return Vec3(A.X[3 * n] + A.u[6 * n] + dpred[6 * n], A.X[3 * n + 1] + A.u[6 * n + 1] + dpred[6 * n + 1],
                A.X[3 * n + 2] + A.u[6 * n + 2] + dpred[6 * n + 2]);
  };
  for (const auto& el : m.elems) {
    if (el.type != ElType::CONN3D2 || el.section < 0 || el.nodes.size() != 2) continue;
    const Section& s = m.sections[el.section];
    if (s.ctype != Section::AXIAL || s.behavior < 0) continue;
    const ConnectorComponent& cc = m.conn_behaviors[s.behavior].comp[0];
    if (!cc.stop) continue;
    const int na = el.nodes[0], nb = el.nodes[1];
    const Vec3 d = X(nb) - X(na);
    const double l = norm(d);
    if (!(l > 0)) continue;
    const double over = l > cc.stop_hi ? l - cc.stop_hi : l < cc.stop_lo ? l - cc.stop_lo : 0.0;
    if (over == 0) continue;
    const Vec3 n = d / l;
    // mobility of each node along the line: inverse mass of its free translational dofs projected on n
    // (prescribed dofs and massless nodes do not move)
    auto mobility = [&](int nd) {
      if (!(A.mt[nd] > 0)) return 0.0;
      double w = 0;
      for (int k = 0; k < 3; ++k) if (!A.fixed[6 * nd + k]) w += n[k] * n[k];
      return w / A.mt[nd];
    };
    const double ia = mobility(na), ib = mobility(nb);
    if (ia + ib == 0) continue;
    // force pair along n restoring the limit in the predicted configuration (relative motion -over along n)
    const double lam = over / (ia + ib) / c;
    const Vec3 fa = n * lam, fb = n * (-lam);
    for (int k = 0; k < 3; ++k) {
      const int ia6 = 6 * na + k, ib6 = 6 * nb + k;
      if (A.mt[na] > 0 && !A.fixed[ia6]) { a[ia6] += fa[k] / A.mt[na]; dpred[ia6] += c * fa[k] / A.mt[na]; }
      if (A.mt[nb] > 0 && !A.fixed[ib6]) { a[ib6] += fb[k] / A.mt[nb]; dpred[ib6] += c * fb[k] / A.mt[nb]; }
      fk[ia6] += fa[k];
      fk[ib6] += fb[k];
    }
  }
}

}  // namespace abq
