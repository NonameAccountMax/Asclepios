// Heat conduction in the elements with temperature degrees of freedom (DC3D8, C3D8T): trilinear temperature,
// conduction matrix k grad N . grad N integrated with 2 x 2 x 2 points on the reference geometry, lumped heat
// capacity rho c integral of N. Prescribed temperatures: *BOUNDARY on degree of freedom 11 (kept with the
// prescribed temperature fields). Steady state: K T = 0 with the prescribed values; transient: implicit Euler,
// (C / dt + K) T1 = C / dt T0. *HEAT TRANSFER steps, and the thermal half of *COUPLED TEMPERATURE-DISPLACEMENT
// increments (the temperature does not depend on the deformation here: no heat generation, no gap conductance).
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cmath>

#include "core/log.hpp"
#include "output/output.hpp"
#include "solver/solvers.hpp"

namespace abq {

namespace {

// conduction matrix (8 x 8) and nodal capacities of a trilinear brick
void brick_heat(const double* X, double k, double rc, double* Ke, double* ce) {
  static const int nat[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}};
  const double g = 1.0 / std::sqrt(3.0);
  for (int i = 0; i < 64; ++i) Ke[i] = 0;
  for (int i = 0; i < 8; ++i) ce[i] = 0;
  for (int p = 0; p < 8; ++p) {
    const double r = nat[p][0] * g, s = nat[p][1] * g, t = nat[p][2] * g;
    double dN[8][3], N[8];
    for (int a = 0; a < 8; ++a) {
      const double xa = nat[a][0], ya = nat[a][1], za = nat[a][2];
      N[a] = 0.125 * (1 + r * xa) * (1 + s * ya) * (1 + t * za);
      dN[a][0] = 0.125 * xa * (1 + s * ya) * (1 + t * za);
      dN[a][1] = 0.125 * ya * (1 + r * xa) * (1 + t * za);
      dN[a][2] = 0.125 * za * (1 + r * xa) * (1 + s * ya);
    }
    Mat3 J;
    for (int a = 0; a < 8; ++a)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) J.m[i][j] += dN[a][i] * X[3 * a + j];
    const double det = J.det();
    const Mat3 Ji = J.inverse();
    Vec3 b[8];
    for (int a = 0; a < 8; ++a) b[a] = Ji * Vec3(dN[a][0], dN[a][1], dN[a][2]);
    for (int a = 0; a < 8; ++a) {
      ce[a] += rc * N[a] * det;
      for (int c = 0; c < 8; ++c) Ke[8 * a + c] += k * dot(b[a], b[c]) * det;
    }
  }
}

}  // namespace

struct HeatSystem {
  int nn = 0;
  Eigen::SparseMatrix<double> K;
  std::vector<double> C;
  std::vector<char> node;  // node with a temperature degree of freedom
};

namespace {

HeatSystem heat_system(const Analysis& A) {
  const Model& m = A.m;
  HeatSystem H;
  H.nn = A.nn;
  H.C.assign(A.nn, 0.0);
  H.node.assign(A.nn, 0);
  std::vector<Eigen::Triplet<double>> trip;
  for (const auto& el : m.elems) {
    if (!el.heat || el.section < 0) continue;
    if (el.nodes.size() != 8) throw std::runtime_error("heat transfer: 8-node bricks only");
    const int mi = m.sections[el.section].material;
    if (mi < 0) continue;
    const Material& mat = m.materials[mi];
    if (mat.conductivity <= 0) throw std::runtime_error("material " + mat.name + ": *CONDUCTIVITY required for heat transfer elements");
    double X[24], Ke[64], ce[8];
    for (int a = 0; a < 8; ++a)
      for (int i = 0; i < 3; ++i) X[3 * a + i] = A.X[3 * el.nodes[a] + i];
    brick_heat(X, mat.conductivity, mat.density * mat.specific_heat, Ke, ce);
    for (int a = 0; a < 8; ++a) {
      H.node[el.nodes[a]] = 1;
      H.C[el.nodes[a]] += ce[a];
      for (int c = 0; c < 8; ++c) trip.push_back({el.nodes[a], el.nodes[c], Ke[8 * a + c]});
    }
  }
  H.K.resize(A.nn, A.nn);
  H.K.setFromTriplets(trip.begin(), trip.end());
  return H;
}

}  // namespace

bool has_heat_elements(const Model& m) {
  for (const auto& el : m.elems) if (el.heat) return true;
  return false;
}

// Temperatures at the end of an increment (step time st, total time tt) into A.temp from A.temp_prev; dt = 0:
// steady state. Returns the largest temperature change of a free node.
double heat_increment(Analysis& A, double dt, double st, double tt) {
  if (!A.heat) A.heat = std::make_shared<HeatSystem>(heat_system(A));
  const HeatSystem& H = *A.heat;
  const int nn = A.nn;
  std::vector<char> pres(nn, 0);
  std::vector<double> T = A.temp_prev;
  for (auto& p : A.loads.temps) { T[p.node] = A.loads.value(p, st, tt, A.m); pres[p.node] = 1; }
  std::vector<int> eq(nn, -1);
  int neq = 0;
  for (int n = 0; n < nn; ++n) if (H.node[n] && !pres[n]) eq[n] = neq++;
  double dTmax = 0;
  if (neq > 0) {
    std::vector<Eigen::Triplet<double>> trip;
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(neq);
    for (int c = 0; c < H.K.outerSize(); ++c)
      for (Eigen::SparseMatrix<double>::InnerIterator it(H.K, c); it; ++it) {
        const int i = (int)it.row(), j = (int)it.col();
        if (eq[i] < 0) continue;
        if (eq[j] >= 0) trip.push_back({eq[i], eq[j], it.value()});
        else rhs[eq[i]] -= it.value() * T[j];
      }
    if (dt > 0)
      for (int n = 0; n < nn; ++n)
        if (eq[n] >= 0) { trip.push_back({eq[n], eq[n], H.C[n] / dt}); rhs[eq[n]] += H.C[n] / dt * A.temp_prev[n]; }
    Eigen::SparseMatrix<double> Kf(neq, neq);
    Kf.setFromTriplets(trip.begin(), trip.end());
    Kf.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.analyzePattern(Kf);
    lu.factorize(Kf);
    if (lu.info() != Eigen::Success) throw std::runtime_error("heat transfer: singular conduction matrix (a region without prescribed temperature?)");
    const Eigen::VectorXd x = lu.solve(rhs);
    for (int n = 0; n < nn; ++n)
      if (eq[n] >= 0) { T[n] = x[eq[n]]; dTmax = std::max(dTmax, std::fabs(T[n] - A.temp_prev[n])); }
  }
  A.temp = T;
  return dTmax;
}

void run_heat_step(Analysis& A, const Step& st, int si) {
  A.begin_step(st, si, st.amp_ramp == -1 ? st.heat_steady : st.amp_ramp == 1);  // default: STEP when transient
  OutputWriter& W = *A.writer;
  W.sta_header(si, st.name, false);
  info("\n STEP %d  %s  (heat transfer, %s)", si + 1, st.name.c_str(), st.heat_steady ? "steady state" : "transient");
  if (A.total_increments == 0 && si == 0) W.frame(si, st.name, 0.0, A.total_time);
  double t = 0, dt = std::min(st.dt0, st.period);
  int inc = 0;
  const int nfr = st.output.field_intervals, ffreq = st.output.field_frequency;
  const double fr_dt = nfr > 0 ? st.period / nfr : 0;
  double next_frame = fr_dt;
  while (t < st.period * (1 - 1e-12)) {
    if (inc >= st.maxinc) throw std::runtime_error("too many increments in step " + st.name);
    dt = std::min(dt, st.period - t);
    if (st.period - t - dt < 1e-9 * st.period) dt = st.period - t;
    const double dt_auto = dt;
    if (nfr > 0 && t + dt > next_frame + 1e-9 * st.period && next_frame - t > 1e-9 * st.period) dt = next_frame - t;
    const double tn = t + dt, ttn = A.total_time + dt;
    A.temp_prev = A.temp;
    const double dT = heat_increment(A, st.heat_steady ? 0.0 : dt, tn, ttn);
    if (!st.heat_steady && st.deltmx > 0 && dT > st.deltmx * (1 + 1e-6)) {  // repeated with a shorter increment
      A.temp = A.temp_prev;
      dt *= std::max(0.25, 0.85 * st.deltmx / dT);
      if (dt < st.dtmin) throw std::runtime_error("time increment smaller than the minimum in step " + st.name);
      continue;
    }
    ++inc;
    t = tn;
    A.total_time = ttn;
    ++A.total_increments;
    W.sta_line("%5d %5d %3d %6d %5d %10.3E %10.3E %10.3E", si + 1, inc, 1, 0, 1, A.total_time, t, dt);
    W.history(si, inc, t, A.total_time, dt, 0.0);
    const bool last = t >= st.period * (1 - 1e-12);
    if (ffreq > 0 ? (inc % ffreq == 0 || last) : (nfr <= 0 || t >= next_frame - 1e-9 * st.period || last)) {
      W.frame(si, st.name, t, A.total_time);
      if (nfr > 0) while (next_frame <= t + 1e-9 * st.period) next_frame += fr_dt;
    }
    dt = dt_auto;
    if (!st.heat_steady && st.deltmx > 0 && dT > 0) dt = std::min(dt * 1.5, dt * 0.85 * st.deltmx / dT);
    dt = std::min(dt, st.dtmax);
  }
  W.sta_line("\n THE STEP %d HAS COMPLETED: %d increments", si + 1, inc);
}

}  // namespace abq
