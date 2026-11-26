// Arc-length procedure (*STATIC, RIKS): the load proportionality factor lambda is an unknown; external loads
// f = f_base + lambda (f_step - f_base). Arc length in the space of lambda and of the translations scaled by u_s,
// the norm of the translations of the first linear solution for lambda = 1: dl^2 = dlambda^2 + |du|^2 / u_s^2.
// Predictor along the tangent, in the direction of the previous increment; corrections in the plane normal to the
// current increment (Riks-Wempner with an updated normal); full Newton tangent (finite differences).
// The step ends at the maximum load factor, when the given displacement is reached, or at the maximum number of
// increments.
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cmath>

#include "core/log.hpp"
#include "output/output.hpp"
#include "solver/solvers.hpp"

namespace abq {

void run_riks_step(Analysis& A, const Step& st, int si) {
  const int nn = A.nn, n6 = 6 * nn;
  A.nlgeom = st.nlgeom;
  A.mat_time = 0.0;
  OutputWriter& W = *A.writer;
  W.sta_header(si, st.name, false);
  info("\n STEP %d  %s  (arc length, %s)", si + 1, st.name.c_str(), st.nlgeom ? "nlgeom" : "small displacements");
  std::vector<double> fbase(n6), fstep(n6);
  std::fill(A.du.begin(), A.du.end(), 0.0);
  A.external_forces(A.loads.period, A.total_time, fbase);
  A.begin_step(st, si, true);
  if (A.total_increments == 0 && si == 0) W.frame(si, st.name, 0.0, A.total_time);
  for (auto& b : A.loads.bcs)
    if (std::fabs(A.loads.value(b, st.period, A.total_time + st.period, A.m) - A.u[6 * b.node + b.dof]) > 1e-12)
      warn("*STATIC, RIKS: prescribed displacement changes ignored (node %d)", A.m.node_label[b.node]);
  A.external_forces(st.period, A.total_time + st.period, fstep);
  std::vector<double> P(n6);
  for (int i = 0; i < n6; ++i) P[i] = fstep[i] - fbase[i];
  const DofMap dm = make_dofmap(A);
  const double hfd = fd_step(A);
  info("  %d equations", dm.neq);
  if (dm.neq == 0) throw std::runtime_error("*STATIC, RIKS: no free degree of freedom");
  std::vector<char> trans(dm.neq, 0);  // translational equations (arc length metric)
  for (int i = 0; i < n6; ++i) if (dm.eq[i] >= 0 && i % 6 < 3) trans[dm.eq[i]] = 1;
  auto to_eq = [&](std::vector<double> f) {
    A.gather_links(f);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(dm.neq);
    for (int i = 0; i < n6; ++i) if (dm.eq[i] >= 0) r[dm.eq[i]] += f[i];
    return r;
  };
  const Eigen::VectorXd Pq = to_eq(P);
  auto tdot = [&](const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    double s = 0;
    for (int i = 0; i < dm.neq; ++i) if (trans[i]) s += a[i] * b[i];
    return s;
  };
  auto factor = [&](Eigen::SparseLU<Eigen::SparseMatrix<double>>& lu) {
    std::vector<Eigen::Triplet<double>> trip;
    assemble_tangent(A, dm, 0.0, hfd, trip);
    Eigen::SparseMatrix<double> K(dm.neq, dm.neq);
    K.setFromTriplets(trip.begin(), trip.end());
    for (int i = 0; i < dm.neq; ++i) K.coeffRef(i, i) += 1e-12 * (1 + std::fabs(K.coeff(i, i)));
    K.makeCompressed();
    lu.analyzePattern(K);
    lu.factorize(K);
    return lu.info() == Eigen::Success;
  };
  auto apply_eq = [&](const Eigen::VectorXd& x, std::vector<double>& v) {
    for (int i = 0; i < n6; ++i) if (dm.eq[i] >= 0) v[i] += x[dm.eq[i]];
    A.spread_links(v);
  };

  double lam = 0, us = 0, dl = st.dt0, arc = 0;
  const double dlmax = st.dtmax, dlmin = st.dtmin;
  Eigen::VectorXd du_prev = Eigen::VectorXd::Zero(dm.neq);
  double dlam_prev = 0;
  int inc = 0, fast = 0;
  std::vector<double> fint_c = A.fint, fext_c = A.fext, R(n6);
  const int ffreq = st.output.field_frequency > 0 ? st.output.field_frequency : 1;
  const double sign_end = st.riks_u < 0 ? -1.0 : 1.0;
  bool done = false;
  while (!done) {
    if (inc >= st.maxinc) { warn("*STATIC, RIKS: maximum number of increments (%d) reached", st.maxinc); break; }
    ++inc;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    std::fill(A.du.begin(), A.du.end(), 0.0);
    if (!factor(lu)) throw std::runtime_error("*STATIC, RIKS: singular tangent");
    const Eigen::VectorXd ut = lu.solve(Pq);
    if (us == 0) {
      us = std::sqrt(tdot(ut, ut));
      if (us == 0) us = 1;
    }
    const double us2 = us * us;
    double sgn = 1.0;
    if (inc > 1 && tdot(du_prev, ut) / us2 + dlam_prev < 0) sgn = -1.0;
    double dlam = sgn * dl / std::sqrt(1 + tdot(ut, ut) / us2);
    Eigen::VectorXd dU = dlam * ut;
    std::fill(A.du.begin(), A.du.end(), 0.0);
    apply_eq(dU, A.du);
    bool converged = false;
    int it;
    double qt = 0;
    for (it = 1; it <= 16; ++it) {
      A.eval_elements(0.0, 0, 0, false, false, A.fint);
      double fsum = 0, rmax = 0;
      int fcnt = 0;
      for (int i = 0; i < n6; ++i) {
        R[i] = fbase[i] + (lam + dlam) * P[i] - A.fint[i];
        if (i % 6 < 3 || A.node_rot[i / 6]) { const double fa = std::fabs(A.fint[i]); if (fa > 0) { fsum += fa; ++fcnt; } }
      }
      const Eigen::VectorXd Rq = to_eq(R);
      for (int i = 0; i < dm.neq; ++i) rmax = std::max(rmax, std::fabs(Rq[i]));
      if (!std::isfinite(rmax)) break;
      qt = std::max(qt, fcnt ? fsum / fcnt : 0.0);
      if (it > 1 && rmax <= 0.005 * std::max(qt, 1e-300)) { converged = true; break; }
      if (!factor(lu)) break;
      const Eigen::VectorXd dR = lu.solve(Rq), dP = lu.solve(Pq);
      const double dln = -(tdot(dU, dR) / us2) / (tdot(dU, dP) / us2 + dlam);
      const Eigen::VectorXd c = dR + dln * dP;
      dU += c;
      dlam += dln;
      std::fill(A.du.begin(), A.du.end(), 0.0);
      apply_eq(dU, A.du);
      double cmax = 0, umax = 0;
      for (int i = 0; i < dm.neq; ++i) if (trans[i]) { cmax = std::max(cmax, std::fabs(c[i])); umax = std::max(umax, std::fabs(dU[i])); }
      if (!std::isfinite(cmax)) break;
      if (it > 1 && rmax <= 0.005 * std::max(qt, 1e-300) && cmax <= 0.01 * std::max(umax, 1e-300)) { converged = true; break; }
    }
    if (!converged) {
      --inc;
      dl *= 0.25;
      W.sta_line("%5d %5d  cutback of the arc length to %.3e", si + 1, inc + 1, dl);
      if (dl < dlmin) throw std::runtime_error("*STATIC, RIKS: arc length increment below the minimum");
      continue;
    }
    // commit
    Analysis::ElemSum S = A.eval_elements(0.0, 0, 0, false, true, A.fint);
    lam += dlam;
    std::vector<double> fe(n6), fr(n6);
    for (int i = 0; i < n6; ++i) { fe[i] = fbase[i] + lam * P[i]; fr[i] = A.fint[i] - fe[i]; }
    A.fext = fe;
    A.gather_links(fr);
    double dWi = 0, dWe = 0;
    for (int i = 0; i < n6; ++i) {
      const double d = A.du[i];
      dWi += 0.5 * (fint_c[i] + A.fint[i]) * d;
      dWe += 0.5 * (fext_c[i] + A.fext[i]) * d;
      A.rf[i] = A.fixed[i] ? fr[i] : (A.report_rf[i] ? A.fint[i] - A.fext[i] : 0.0);
      A.u[i] += A.du[i];
    }
    A.commit_rotations();
    A.en.plastic += S.e_plastic;
    A.en.internal += dWi + S.e_corr;
    A.en.external += dWe;
    fint_c = A.fint;
    fext_c = A.fext;
    du_prev = dU;
    dlam_prev = dlam;
    const double dla = std::sqrt(dlam * dlam + tdot(dU, dU) / us2);
    arc += dla;
    A.total_time += dla;
    ++A.total_increments;
    W.sta_line("%5d %5d %3d %6d %5d %10.3E %10.3E LPF %10.4E", si + 1, inc, 1, 0, it - 1, A.total_time, arc, lam);
    info("  increment %4d  arc length %.4g  LPF %.6g  iterations %d", inc, arc, lam, it - 1);
    W.history(si, inc, arc, A.total_time, dla, 0.0);
    if (st.riks_lpf_max > 0 && lam >= st.riks_lpf_max) done = true;
    if (st.riks_node >= 0 && sign_end * (A.u[6 * st.riks_node + st.riks_dof] - st.riks_u) >= 0) done = true;
    if (inc % ffreq == 0 || done) W.frame(si, st.name, arc, A.total_time);
    // arc length control: longer after quick convergence, shorter after slow
    if (it - 1 <= 4) ++fast; else fast = 0;
    if (fast >= 2) { dl = std::min(dl * 1.5, dlmax); fast = 0; }
    else if (it - 1 > 8) dl = std::max(dl * 0.75, dlmin);
  }
  std::fill(A.du.begin(), A.du.end(), 0.0);
  W.sta_line("\n THE STEP %d HAS COMPLETED: %d increments, LPF %.6g", si + 1, inc, lam);
}

}  // namespace abq
