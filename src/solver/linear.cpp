// Linear perturbation procedures about the base state (the state reached at the end of the last general step):
// *FREQUENCY (natural modes), *BUCKLE (eigenvalue buckling), *STATIC in a PERTURBATION step, *MODAL DYNAMIC and
// *STEADY STATE DYNAMICS (mode superposition on the modes of the last *FREQUENCY step).
// The loads and boundary conditions of a perturbation step act on that step only; its response is not added to
// the base state. Stiffness: element tangents by finite differences (as the static procedure); masses: lumped, with
// the physical rotary inertia (as implicit dynamics).
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <random>

#include "core/binio.hpp"
#include "core/log.hpp"
#include "output/output.hpp"
#include "solver/solvers.hpp"

namespace abq {

namespace {

using SpMat = Eigen::SparseMatrix<double>;

SpMat tangent_matrix(Analysis& A, const DofMap& dm, double hfd) {
  std::fill(A.du.begin(), A.du.end(), 0.0);
  std::vector<Eigen::Triplet<double>> trip;
  assemble_tangent(A, dm, 0.0, hfd, trip);
  SpMat K(dm.neq, dm.neq);
  K.setFromTriplets(trip.begin(), trip.end());
  K.makeCompressed();
  return K;
}

// Generalised symmetric eigenproblem B x = mu K x with K positive definite on the equations: the p eigenpairs of
// largest |mu| (lowest |lambda| = 1 / |mu| of K x = lambda B x), by subspace iteration with Rayleigh-Ritz projection
// (q = max(2p, p + 8) vectors). Columns of X are K-orthonormal.
struct EigenPairs {
  std::vector<double> mu;
  Eigen::MatrixXd X;
  int iterations = 0;
};
EigenPairs subspace_iteration(const SpMat& K, const std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)>& Bmul,
                              const Eigen::VectorXd& bdiag, int p, int max_iter) {
  const int n = (int)K.rows();
  p = std::min(p, n);
  const int q = std::min(n, std::max(2 * p, p + 8));
  Eigen::SparseLU<SpMat> lu;
  lu.analyzePattern(K);
  lu.factorize(K);
  if (lu.info() != Eigen::Success) throw std::runtime_error("eigenvalue extraction: singular stiffness (free body or mechanism)");
  // starting vectors: the diagonal of B, unit vectors at the dofs of largest b_ii / k_ii, then random ones
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, q);
  X.col(0) = bdiag.cwiseAbs();
  if (X.col(0).norm() == 0) X.col(0).setOnes();
  std::vector<std::pair<double, int>> ratio;
  for (int i = 0; i < n; ++i) ratio.push_back({-std::fabs(bdiag[i]) / std::max(std::fabs(K.coeff(i, i)), 1e-300), i});
  std::sort(ratio.begin(), ratio.end());
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);
  for (int c = 1; c < q; ++c) {
    if (c < q - 2 && c - 1 < n && ratio[c - 1].first < 0) X(ratio[c - 1].second, c) = 1.0;
    else for (int i = 0; i < n; ++i) X(i, c) = uni(rng);
  }
  EigenPairs r;
  std::vector<double> prev(q, 0.0);
  for (int it = 1; it <= max_iter; ++it) {
    const Eigen::MatrixXd Y = Bmul(X);
    Eigen::MatrixXd Xb(n, q);
    for (int c = 0; c < q; ++c) Xb.col(c) = lu.solve(Y.col(c));
    Eigen::MatrixXd Kp = Xb.transpose() * Y;   // Xb^T K Xb
    Eigen::MatrixXd Bp = Xb.transpose() * Bmul(Xb);
    Kp = 0.5 * (Kp + Kp.transpose());
    Bp = 0.5 * (Bp + Bp.transpose());
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Bp, Kp);
    if (ges.info() != Eigen::Success) throw std::runtime_error("eigenvalue extraction: projected problem failed");
    std::vector<int> order(q);
    for (int c = 0; c < q; ++c) order[c] = c;
    const Eigen::VectorXd ev = ges.eigenvalues();
    std::sort(order.begin(), order.end(), [&](int a, int b) { return std::fabs(ev[a]) > std::fabs(ev[b]); });
    Eigen::MatrixXd V(q, q);
    std::vector<double> mu(q);
    for (int c = 0; c < q; ++c) { V.col(c) = ges.eigenvectors().col(order[c]); mu[c] = ev[order[c]]; }
    X = Xb * V;
    bool conv = it > 1;
    for (int c = 0; c < p && conv; ++c) conv = std::fabs(mu[c] - prev[c]) <= 1e-12 * std::fabs(mu[c]);
    prev = mu;
    r.iterations = it;
    if (conv) break;
    if (it == max_iter) warn("eigenvalue extraction: no convergence in %d iterations", max_iter);
  }
  r.mu.assign(prev.begin(), prev.begin() + p);
  r.X = X.leftCols(p);
  return r;
}

// external forces of the loads of the step at step time t, less those of the base state
void perturbation_load(Analysis& A, double t, const std::vector<double>& fbase, std::vector<double>& f) {
  A.external_forces(t, A.total_time + t, f);
  for (size_t i = 0; i < f.size(); ++i) f[i] -= fbase[i];
}

Eigen::VectorXd to_eq(const Analysis& A, const DofMap& dm, std::vector<double> f) {
  A.gather_links(f);
  Eigen::VectorXd r = Eigen::VectorXd::Zero(dm.neq);
  for (size_t i = 0; i < f.size(); ++i) if (dm.eq[i] >= 0) r[dm.eq[i]] += f[i];
  return r;
}

std::vector<double> from_eq(const Analysis& A, const DofMap& dm, const Eigen::VectorXd& x) {
  std::vector<double> u(6 * (size_t)A.nn, 0.0);
  for (size_t i = 0; i < u.size(); ++i) if (dm.eq[i] >= 0) u[i] = x[dm.eq[i]];
  A.spread_links(u);
  return u;
}

// frame of a perturbation field u (6nn), reactions rf; the base state is restored afterwards
void perturbation_frame(Analysis& A, int si, const std::string& name, double t, const std::vector<double>& u,
                        const std::vector<double>* rf, const std::vector<double>* v = nullptr) {
  std::vector<double> u0 = A.u, rf0 = A.rf, v0 = A.v;
  A.u = u;
  if (rf) A.rf = *rf; else std::fill(A.rf.begin(), A.rf.end(), 0.0);
  if (v) A.v = *v; else std::fill(A.v.begin(), A.v.end(), 0.0);
  A.writer->frame(si, name, t, A.total_time + t);
  A.u = u0; A.rf = rf0; A.v = v0;
}

// reactions of a perturbation displacement u under the perturbation load f: element forces of u (small
// displacements) less the load, on the prescribed dofs
std::vector<double> perturbation_reactions(Analysis& A, const std::vector<double>& u, const std::vector<double>& f) {
  const bool nl = A.nlgeom;
  A.du = u;
  A.nlgeom = false;
  std::vector<double> fi(u.size());
  A.eval_elements(0.0, 0, 0, false, false, fi);
  A.nlgeom = nl;
  std::fill(A.du.begin(), A.du.end(), 0.0);
  std::vector<double> r(u.size());
  for (size_t i = 0; i < u.size(); ++i) r[i] = fi[i] - f[i];
  A.gather_links(r);
  for (size_t i = 0; i < u.size(); ++i) r[i] = A.fixed[i] ? r[i] : 0.0;
  return r;
}

double damping_ratio(const Step& st, int mode1) {
  double z = 0;
  for (const auto& d : st.modal_damping)
    if (mode1 >= d.first && mode1 <= d.last) z = d.ratio;
  return z;
}

void frequency_step(Analysis& A, const Step& st, int si, const DofMap& dm, double hfd) {
  OutputWriter& W = *A.writer;
  const SpMat K = tangent_matrix(A, dm, hfd);
  std::vector<double> mdof;
  dof_masses(A, mdof);
  Eigen::VectorXd md = Eigen::VectorXd::Zero(dm.neq);
  for (int i = 0; i < 6 * A.nn; ++i) if (dm.eq[i] >= 0) md[dm.eq[i]] += mdof[i];
  auto Mmul = [&](const Eigen::MatrixXd& X) { return Eigen::MatrixXd(md.asDiagonal() * X); };
  const EigenPairs ep = subspace_iteration(K, Mmul, md, st.neig, 300);
  A.modes.clear();
  A.mode_omega.clear();
  info("  %d eigenvalues (%d subspace iterations)", (int)ep.mu.size(), ep.iterations);
  info("  MODE   EIGENVALUE       FREQUENCY (RAD/TIME)   FREQUENCY (CYCLES/TIME)");
  W.sta_line("  MODE   EIGENVALUE   FREQUENCY (CYCLES/TIME)");
  int kept = 0;
  for (size_t k = 0; k < ep.mu.size(); ++k) {
    const double mu = ep.mu[k];
    if (!(mu > 0)) continue;
    const double w2 = 1.0 / mu, w = std::sqrt(w2), f = w / (2 * M_PI);
    if (st.eig_min > 0 && f < st.eig_min) continue;
    if (st.eig_max > 0 && f > st.eig_max) continue;
    ++kept;
    // mass normalisation: x^T M x = mu for K-orthonormal x
    Eigen::VectorXd phi = ep.X.col(k) / std::sqrt(mu);
    std::vector<double> u = from_eq(A, dm, phi);
    A.modes.push_back(u);
    A.mode_omega.push_back(w);
    info("  %4d  %14.6e  %14.6e  %14.6e", kept, w2, w, f);
    W.sta_line("  %4d  %12.5e  %12.5e", kept, w2, f);
    W.next_frame_mode(kept, f);
    perturbation_frame(A, si, st.name, (double)kept, u, nullptr);
  }
}

void buckle_step(Analysis& A, const Step& st, int si, const DofMap& dm, double hfd, const std::vector<double>& fbase) {
  OutputWriter& W = *A.writer;
  const bool nl = A.nlgeom;
  A.nlgeom = true;  // stress stiffness of the base state and of the perturbation load
  const SpMat KA = tangent_matrix(A, dm, hfd);
  std::vector<double> f(6 * A.nn);
  perturbation_load(A, st.period, fbase, f);
  Eigen::SparseLU<SpMat> lu;
  lu.analyzePattern(KA);
  lu.factorize(KA);
  if (lu.info() != Eigen::Success) throw std::runtime_error("*BUCKLE: singular base state stiffness");
  std::vector<double> uP = from_eq(A, dm, lu.solve(to_eq(A, dm, f)));
  double umax = 0;
  for (size_t i = 0; i < uP.size(); ++i) if (i % 6 < 3) umax = std::max(umax, std::fabs(uP[i]));
  if (umax == 0) throw std::runtime_error("*BUCKLE: the perturbation load produces no displacement");
  // stress stiffness of the perturbation load: the element states take the stresses of the (scaled) linear
  // response, the tangent is evaluated again at the base geometry; states restored afterwards
  const double s = 0.01 * (hfd / 1e-7) / umax;
  const std::string tmp = A.job + ".bkl.tmp";
  {
    BinOut o(tmp);
    for (auto& g : A.groups) g->save(o);
  }
  for (size_t i = 0; i < uP.size(); ++i) A.du[i] = s * uP[i];
  A.nlgeom = false;
  std::vector<double> scratch(6 * A.nn);
  A.eval_elements(0.0, 0, 0, false, true, scratch);
  A.nlgeom = true;
  const SpMat KB = tangent_matrix(A, dm, hfd);
  {
    BinIn in(tmp);
    for (auto& g : A.groups) g->load(in);
  }
  std::remove(tmp.c_str());
  A.nlgeom = nl;
  const SpMat KG = (KB - KA) / s;
  // K x = lambda (-KG) x
  const SpMat Bm = -KG;
  auto Bmul = [&](const Eigen::MatrixXd& X) { return Eigen::MatrixXd(Bm * X); };
  const EigenPairs ep = subspace_iteration(KA, Bmul, Eigen::VectorXd(Bm.diagonal()), st.neig, 300);
  info("  MODE   EIGENVALUE (buckling load factor), %d subspace iterations", ep.iterations);
  W.sta_line("  MODE   EIGENVALUE");
  int kept = 0;
  for (size_t k = 0; k < ep.mu.size(); ++k) {
    const double lam = 1.0 / ep.mu[k];
    if (st.eig_max > 0 && std::fabs(lam) > st.eig_max) continue;
    ++kept;
    std::vector<double> u = from_eq(A, dm, ep.X.col(k));
    double m = 0;  // largest displacement component scaled to 1
    for (size_t i = 0; i < u.size(); ++i) if (i % 6 < 3 && std::fabs(u[i]) > std::fabs(m)) m = u[i];
    if (m != 0) for (double& x : u) x /= m;
    info("  %4d  %14.6e", kept, lam);
    W.sta_line("  %4d  %12.5e", kept, lam);
    W.next_frame_mode(kept, lam);
    perturbation_frame(A, si, st.name, (double)kept, u, nullptr);
  }
}

void static_perturbation_step(Analysis& A, const Step& st, int si, const DofMap& dm, double hfd,
                              const std::vector<double>& fbase) {
  const SpMat K = tangent_matrix(A, dm, hfd);
  std::vector<double> f(6 * A.nn);
  perturbation_load(A, st.period, fbase, f);
  Eigen::SparseLU<SpMat> lu;
  lu.analyzePattern(K);
  lu.factorize(K);
  if (lu.info() != Eigen::Success) throw std::runtime_error("static perturbation: singular stiffness");
  const std::vector<double> u = from_eq(A, dm, lu.solve(to_eq(A, dm, f)));
  const std::vector<double> rf = perturbation_reactions(A, u, f);
  // a perturbation step takes no time (frame at the end of the base state, as the reference solver writes it)
  perturbation_frame(A, si, st.name, 2.220446049250313e-16, u, &rf);
}

// mode superposition: modal force of mode k for a nodal force vector f
double modal_force(const std::vector<double>& phi, const std::vector<double>& f) {
  double s = 0;
  for (size_t i = 0; i < f.size(); ++i) s += phi[i] * f[i];
  return s;
}

void modal_dynamic_step(Analysis& A, const Step& st, int si, const std::vector<double>& fbase) {
  if (A.modes.empty()) throw std::runtime_error("*MODAL DYNAMIC: no modes (a *FREQUENCY step must come first)");
  const int nm = (int)A.modes.size(), n6 = 6 * A.nn;
  const double dt = st.dt0;
  const int nstep = std::max(1, (int)std::lround(st.period / dt));
  std::vector<double> q(nm, 0.0), qd(nm, 0.0), p0(nm), p1(nm), f(n6);
  auto modal_loads = [&](double t, std::vector<double>& p) {
    perturbation_load(A, t, fbase, f);
    std::vector<double> g = f;
    A.gather_links(g);
    for (int k = 0; k < nm; ++k) p[k] = modal_force(A.modes[k], g);
  };
  // exact response of q'' + 2 z w q' + w^2 q = p(t) for p linear over the increment (unit modal mass)
  struct Coef { double a, b, c, d, ad, bd, cd, dd; };
  std::vector<Coef> cf(nm);
  for (int k = 0; k < nm; ++k) {
    const double w = A.mode_omega[k], z = std::min(damping_ratio(st, k + 1), 0.999999), kk = w * w;
    const double sq = std::sqrt(1 - z * z), wd = w * sq, e = std::exp(-z * w * dt), sn = std::sin(wd * dt),
                 cs = std::cos(wd * dt);
    Coef& C = cf[k];
    C.a = e * (z / sq * sn + cs);
    C.b = e * sn / wd;
    C.c = (2 * z / (w * dt) + e * (((1 - 2 * z * z) / (wd * dt) - z / sq) * sn - (1 + 2 * z / (w * dt)) * cs)) / kk;
    C.d = (1 - 2 * z / (w * dt) + e * ((2 * z * z - 1) / (wd * dt) * sn + 2 * z / (w * dt) * cs)) / kk;
    C.ad = -e * w / sq * sn;
    C.bd = e * (cs - z / sq * sn);
    C.cd = (-1 / dt + e * ((w / sq + z / (dt * sq)) * sn + cs / dt)) / kk;
    C.dd = (1 - e * (z / sq * sn + cs)) / (kk * dt);
  }
  OutputWriter& W = *A.writer;
  const int nfr = st.output.field_intervals, ffreq = st.output.field_frequency;
  const double fr_dt = nfr > 0 ? st.period / nfr : 0;
  double next_frame = fr_dt;
  auto output = [&](double t) {
    std::vector<double> u(n6, 0.0), v(n6, 0.0);
    for (int k = 0; k < nm; ++k)
      for (int i = 0; i < n6; ++i) { u[i] += A.modes[k][i] * q[k]; v[i] += A.modes[k][i] * qd[k]; }
    perturbation_frame(A, si, st.name, t, u, nullptr, &v);
  };
  modal_loads(0.0, p0);
  for (int n = 1; n <= nstep; ++n) {
    const double t = n * dt;
    modal_loads(t, p1);
    for (int k = 0; k < nm; ++k) {
      const Coef& C = cf[k];
      const double qn = C.a * q[k] + C.b * qd[k] + C.c * p0[k] + C.d * p1[k];
      const double vn = C.ad * q[k] + C.bd * qd[k] + C.cd * p0[k] + C.dd * p1[k];
      q[k] = qn; qd[k] = vn;
    }
    p0 = p1;
    const bool last = n == nstep;
    if (ffreq > 0 ? (n % ffreq == 0 || last) : (nfr <= 0 || t >= next_frame - 1e-9 * st.period || last)) {
      output(t);
      if (nfr > 0) while (next_frame <= t + 1e-9 * st.period) next_frame += fr_dt;
    }
    W.history(si, n, t, A.total_time + t, dt, 0.0);
  }
  info("  %d increments of %g on %d modes", nstep, dt, nm);
}

void steady_state_step(Analysis& A, const Step& st, int si, const std::vector<double>& fbase) {
  if (A.modes.empty()) throw std::runtime_error("*STEADY STATE DYNAMICS: no modes (a *FREQUENCY step must come first)");
  const int nm = (int)A.modes.size(), n6 = 6 * A.nn;
  std::vector<double> f(n6);
  perturbation_load(A, st.period, fbase, f);
  A.gather_links(f);
  std::vector<double> pk(nm);
  for (int k = 0; k < nm; ++k) pk[k] = modal_force(A.modes[k], f);
  // frequency points: the range divided at the eigenfrequencies inside it, N points per interval (ends included),
  // logarithmic scale, clustered towards the ends by the bias b: s = 1/2 + 1/2 |y|^(1/b) sign(y), y in [-1, 1]
  const double f0 = st.ssd_f0, f1 = std::max(st.ssd_f1, st.ssd_f0);
  std::vector<double> anchors = {f0};
  std::vector<double> fe;
  for (double w : A.mode_omega) fe.push_back(w / (2 * M_PI));
  std::sort(fe.begin(), fe.end());
  for (double x : fe) if (x > f0 * (1 + 1e-9) && x < f1 * (1 - 1e-9)) anchors.push_back(x);
  if (f1 > f0) anchors.push_back(f1);
  std::vector<double> freqs = {f0};
  const int N = std::max(2, st.ssd_points);
  const double b = st.ssd_bias > 0 ? st.ssd_bias : 1.0;
  for (size_t j = 0; j + 1 < anchors.size(); ++j) {
    const double la = std::log(std::max(anchors[j], 1e-300)), lb = std::log(std::max(anchors[j + 1], 1e-300));
    for (int k = 1; k < N; ++k) {
      const double y = -1 + 2.0 * k / (N - 1);
      const double sgn = y < 0 ? -1.0 : 1.0, sfrac = 0.5 + 0.5 * sgn * std::pow(std::fabs(y), 1.0 / b);
      freqs.push_back(k == N - 1 ? anchors[j + 1] : std::exp(la + sfrac * (lb - la)));
    }
  }
  if (f0 <= 0) warn("*STEADY STATE DYNAMICS: logarithmic frequency scale with a zero lower frequency");
  // response u = sum phi_k p_k / (w_k^2 - W^2 + 2 i z_k w_k W): real part as U, imaginary part not written
  for (double fr : freqs) {
    const double W2 = std::pow(2 * M_PI * fr, 2), Wc = 2 * M_PI * fr;
    std::vector<double> u(n6, 0.0);
    for (int k = 0; k < nm; ++k) {
      const double w = A.mode_omega[k], z = damping_ratio(st, k + 1);
      const double re = w * w - W2, im = 2 * z * w * Wc, den = re * re + im * im;
      const double hr = re / den;
      for (int i = 0; i < n6; ++i) u[i] += A.modes[k][i] * pk[k] * hr;
    }
    perturbation_frame(A, si, st.name, fr, u, nullptr);
  }
  info("  %d frequencies from %g to %g", (int)freqs.size(), f0, f1);
}

}  // namespace

void run_linear_step(Analysis& A, const Step& st, int si) {
  const bool nl0 = A.nlgeom;
  A.nlgeom = st.nlgeom;
  A.mat_time = 0.0;
  OutputWriter& W = *A.writer;
  const char* what = st.proc == Step::FREQUENCY ? "frequency" : st.proc == Step::BUCKLE ? "buckle"
                   : st.proc == Step::MODAL_DYNAMIC ? "modal dynamic" : st.proc == Step::STEADY_STATE ? "steady state dynamics"
                   : "static perturbation";
  info("\n STEP %d  %s  (%s, linear perturbation)", si + 1, st.name.c_str(), what);
  W.sta_header(si, st.name, false);
  // base state: loads and boundary conditions restored after the step
  std::vector<double> fbase(6 * A.nn);
  std::fill(A.du.begin(), A.du.end(), 0.0);
  A.external_forces(A.loads.period, A.total_time, fbase);
  const LoadState loads0 = A.loads;
  A.begin_step(st, si, st.proc != Step::MODAL_DYNAMIC);
  // frame of the base state at the start of the step (the reference writes one per step)
  W.frame(si, st.name, 0.0, A.total_time);
  const DofMap dm = make_dofmap(A);
  info("  %d equations", dm.neq);
  const double hfd = fd_step(A);
  if (dm.neq == 0) throw std::runtime_error("step " + st.name + ": no free degree of freedom");
  switch (st.proc) {
    case Step::FREQUENCY: frequency_step(A, st, si, dm, hfd); break;
    case Step::BUCKLE: buckle_step(A, st, si, dm, hfd, fbase); break;
    case Step::MODAL_DYNAMIC: modal_dynamic_step(A, st, si, fbase); break;
    case Step::STEADY_STATE: steady_state_step(A, st, si, fbase); break;
    default: static_perturbation_step(A, st, si, dm, hfd, fbase); break;
  }
  std::fill(A.du.begin(), A.du.end(), 0.0);
  A.loads = loads0;
  A.update_fixed(st);
  A.nlgeom = nl0;
  W.sta_line("\n THE STEP %d HAS COMPLETED", si + 1);
}

}  // namespace abq
