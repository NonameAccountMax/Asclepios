// Static implicit procedure (*STATIC): full Newton-Raphson on the incremental (corotational)
// element formulation shared with the explicit driver. Element tangents by central finite
// differences of the internal forces; sparse direct solve (Eigen SparseLU).
// Default Abaqus/Standard controls: R_n = 0.005, C_n = 0.01, I_0 = 4, I_R = 8, I_C = 16,
// I_G = 4 (two consecutive increments -> x1.5), cutback D_f = 0.25.
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <omp.h>

#include <chrono>
#include <cmath>

#include "contact/contact.hpp"
#include "core/log.hpp"
#include "output/output.hpp"
#include "solver/solvers.hpp"

namespace abq {

namespace {

double wall() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct DofMap {
  std::vector<int> eq;  // 6nn -> equation or -1
  int neq = 0;
};

// element tangent by central differences of compute(commit=false)
void element_tangent(ElementGroup& G, int e, ElemIn in, double* du, int nen, double hscale, double* K) {
  const int nd = 6 * nen;
  std::vector<double> fp(nd), fm(nd);
  ElemOut out;
  for (int j = 0; j < nd; ++j) {
    double h = (j % 6) < 3 ? hscale : 1e-6;
    double save = du[j];
    du[j] = save + h;
    out = ElemOut(); out.f = fp.data();
    G.compute(e, in, out, false);
    du[j] = save - h;
    out = ElemOut(); out.f = fm.data();
    G.compute(e, in, out, false);
    du[j] = save;
    for (int i = 0; i < nd; ++i) K[i * nd + j] = (fp[i] - fm[i]) / (2 * h);
  }
}

}  // namespace

void run_static_step(Analysis& A, const Step& st, int si) {
  const int nn = A.nn;
  A.nlgeom = st.nlgeom;
  A.begin_step(st, si, st.amp_ramp != 0);
  if (st.contact_defined) warn("general contact is not supported in *STATIC steps: ignored");
  OutputWriter& W = *A.writer;
  W.sta_header(si, st.name, false);
  info("\n STEP %d  %s  (static, %s)  period %g", si + 1, st.name.c_str(), st.nlgeom ? "nlgeom" : "small displacements", st.period);
  if (A.total_increments == 0 && si == 0) W.frame(si, st.name, 0.0, A.total_time);

  // active dofs: translations of nodes carrying elements, rotations of shell/beam nodes
  std::vector<char> has_el(nn, 0);
  for (auto& gp : A.groups) for (int n : gp->conn) has_el[n] = 1;
  DofMap dm;
  dm.eq.assign(6 * nn, -1);
  for (int n = 0; n < nn; ++n)
    for (int k = 0; k < 6; ++k) {
      if (!has_el[n] || (k >= 3 && !A.node_rot[n])) continue;
      if (A.fixed[6 * n + k]) continue;
      dm.eq[6 * n + k] = dm.neq++;
    }
  // tied translations (CONN3D2 JOIN) share the equations of the first node
  std::vector<char> slave(6 * nn, 0);
  for (auto [na, nb] : A.ties)
    for (int k = 0; k < 3; ++k) {
      int ia = 6 * na + k, ib = 6 * nb + k;
      if (dm.eq[ia] < 0 || dm.eq[ib] < 0) continue;
      slave[ib] = 1;
    }
  if (!A.ties.empty()) {  // renumber without the slave equations, then point the slaves to their master
    int neq = 0;
    std::vector<int> old = dm.eq;
    for (int i = 0; i < 6 * nn; ++i) dm.eq[i] = (old[i] >= 0 && !slave[i]) ? neq++ : -1;
    dm.neq = neq;
    for (auto [na, nb] : A.ties)
      for (int k = 0; k < 3; ++k)
        if (slave[6 * nb + k]) dm.eq[6 * nb + k] = dm.eq[6 * na + k];
  }
  info("  %d equations", dm.neq);
  // characteristic length for finite-difference steps
  double lmin = 1e30;
  for (auto& gp : A.groups)
    for (int e = 0; e < gp->size(); ++e) {
      int a = gp->conn[e * gp->nen], b = gp->conn[e * gp->nen + 1];
      double l = std::sqrt(std::pow(A.X[3 * a] - A.X[3 * b], 2) + std::pow(A.X[3 * a + 1] - A.X[3 * b + 1], 2) +
                           std::pow(A.X[3 * a + 2] - A.X[3 * b + 2], 2));
      if (l > 0) lmin = std::min(lmin, l);
    }
  const double hfd = 1e-7 * lmin;

  double t = 0, dt = std::min(st.dt0, st.period);
  int inc = 0, fast = 0;
  double qtilde = 0;
  int nq = 0;
  const double t0wall = wall();
  std::vector<double> R(6 * nn), du_save;
  // last converged nodal forces, for the work of internal forces, external loads and reactions (trapezoidal rule)
  std::vector<double> fint_c = A.fint, fext_c = A.fext;
  const int nfr = st.output.field_intervals;
  const double fr_dt = nfr > 0 ? st.period / nfr : 0;
  double next_frame = fr_dt;

  while (t < st.period * (1 - 1e-12)) {
    if (inc >= st.maxinc) throw std::runtime_error("too many increments in step " + st.name + " (INC=" + std::to_string(st.maxinc) + ")");
    dt = std::min(dt, st.period - t);
    if (st.period - t - dt < 1e-9 * st.period) dt = st.period - t;
    // output times are reached exactly (as with TIME MARKS=YES): the increment is shortened, not the next ones
    const double dt_auto = dt;
    if (nfr > 0 && t + dt > next_frame + 1e-9 * st.period && next_frame - t > 1e-9 * st.period) dt = next_frame - t;
    const bool cut_to_frame = dt < dt_auto;
    const double tn = t + dt, ttn = A.total_time + dt;
    ++inc;
    // predictor: prescribed increments only
    std::fill(A.du.begin(), A.du.end(), 0.0);
    for (auto& b : A.loads.bcs) {
      double target = A.loads.value(b, tn, ttn, A.m);
      A.du[6 * b.node + b.dof] = b.velocity ? target * dt : target - A.u[6 * b.node + b.dof];
    }
    A.temp_prev = A.temp;
    A.temperatures(tn, ttn);
    bool converged = false;
    int it = 0;
    double rmax = 0, cmax = 0, dumax = 0;
    for (it = 1; it <= 16; ++it) {
      // residual
      A.eval_elements(dt, 0, 0, false, false, A.fint);
      A.external_forces(tn, ttn, A.fext);
      double fsum = 0;
      int fcnt = 0;
      rmax = 0;
      std::vector<double> req(dm.neq, 0.0);  // residual per equation (tied dofs summed)
      bool finite = true;
      for (int i = 0; i < 6 * nn; ++i) {
        R[i] = A.fext[i] - A.fint[i];
        finite = finite && std::isfinite(R[i]);
        if (dm.eq[i] >= 0) req[dm.eq[i]] += R[i];
        if (i % 6 < 3 || A.node_rot[i / 6]) {  // time-averaged force: all active dofs, reactions included
          double fa = std::fabs(A.fint[i]);
          if (fa > 0) { fsum += fa; ++fcnt; }
        }
      }
      for (double r : req) rmax = std::max(rmax, std::fabs(r));
      double qcur = fcnt ? fsum / fcnt : 0;
      double qt = nq ? std::max(qtilde, qcur) : qcur;
      if (qt <= 0) qt = 1e-300;
      // a diverged iterate (a force not finite, reactions included) is never accepted: NaN would be skipped by
      // the maxima below, and an infinite reaction would make the relative criterion trivially satisfied
      if (!finite || !std::isfinite(rmax) || !std::isfinite(qcur)) break;
      static const bool trace_newton = std::getenv("ABQ_TRACE_NEWTON") != nullptr;
      if (trace_newton) std::fprintf(stderr, "NEWTON inc %d it %d dt %.3e rmax %.3e qt %.3e cmax %.3e dumax %.3e\n", inc, it, dt, rmax, qt, cmax, dumax);
      if (it > 1 && rmax <= 0.005 * qt && cmax <= 0.01 * std::max(dumax, 1e-300)) { converged = true; break; }
      if (rmax <= 1e-8 * qt && it > 1) { converged = true; break; }
      if (dm.neq == 0) { converged = true; break; }  // fully prescribed
      // tangent
      std::vector<Eigen::Triplet<double>> trip;
      for (auto& gp : A.groups) {
        ElementGroup& G = *gp;
        const int nen = G.nen, nd = 6 * nen;
        std::vector<std::vector<Eigen::Triplet<double>>> local(A.nthreads);
#pragma omp parallel num_threads(A.nthreads)
        {
          std::vector<double> Xe(3 * nen), xo(3 * nen), due(nd), t0(nen), t1(nen), fv(nen), K(nd * nd);
          auto& loc = local[omp_get_thread_num()];
#pragma omp for schedule(dynamic, 32)
          for (int e = 0; e < G.size(); ++e) {
            const int* cn = &G.conn[e * nen];
            for (int I = 0; I < nen; ++I) {
              int n = cn[I];
              for (int k = 0; k < 3; ++k) { Xe[3 * I + k] = A.X[3 * n + k]; xo[3 * I + k] = A.X[3 * n + k] + A.u[6 * n + k]; }
              for (int k = 0; k < 6; ++k) due[6 * I + k] = A.du[6 * n + k];
              t0[I] = A.temp_prev[n]; t1[I] = A.temp[n]; fv[I] = A.fv1[n];
            }
            ElemIn in;
            in.X = Xe.data(); in.xo = xo.data(); in.du = due.data();
            if (A.has_temp) { in.temp0 = t0.data(); in.temp1 = t1.data(); }
            in.fv1 = fv.data(); in.dt = dt; in.dynamic = false; in.nlgeom = A.nlgeom;
            element_tangent(G, e, in, due.data(), nen, hfd, K.data());
            for (int i = 0; i < nd; ++i) {
              int gi = dm.eq[6 * cn[i / 6] + i % 6];
              if (gi < 0) continue;
              for (int j = 0; j < nd; ++j) {
                int gj = dm.eq[6 * cn[j / 6] + j % 6];
                if (gj < 0 || K[i * nd + j] == 0) continue;
                loc.push_back({gi, gj, K[i * nd + j]});
              }
            }
          }
        }
        for (auto& l : local) trip.insert(trip.end(), l.begin(), l.end());
      }
      Eigen::SparseMatrix<double> Kg(dm.neq, dm.neq);
      Kg.setFromTriplets(trip.begin(), trip.end());
      // tiny regularisation of dofs without stiffness (e.g. drilling of isolated nodes)
      for (int i = 0; i < dm.neq; ++i) Kg.coeffRef(i, i) += 1e-12 * (1 + std::fabs(Kg.coeff(i, i)));
      Kg.makeCompressed();
      Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
      lu.analyzePattern(Kg);
      lu.factorize(Kg);
      if (lu.info() != Eigen::Success) { warn("increment %d: singular tangent", inc); break; }
      Eigen::VectorXd rhs(dm.neq);
      rhs.setZero();
      for (int i = 0; i < 6 * nn; ++i) if (dm.eq[i] >= 0) rhs[dm.eq[i]] += R[i];
      Eigen::VectorXd c = lu.solve(rhs);
      cmax = 0;
      for (int i = 0; i < 6 * nn; ++i)
        if (dm.eq[i] >= 0) { A.du[i] += c[dm.eq[i]]; cmax = std::max(cmax, std::fabs(c[dm.eq[i]])); }
      dumax = 0;
      for (int i = 0; i < 6 * nn; ++i) if (i % 6 < 3) dumax = std::max(dumax, std::fabs(A.du[i]));
      if (!std::isfinite(cmax)) break;
    }
    if (!converged) {
      dt *= 0.25;
      --inc;
      W.sta_line("%5d %5d  cutback (residual %.3e)", si + 1, inc + 1, rmax);
      if (dt < st.dtmin) throw std::runtime_error("time increment smaller than the minimum in step " + st.name);
      A.temp = A.temp_prev;
      continue;
    }
    // commit
    Analysis::ElemSum S = A.eval_elements(dt, 0, 0, false, true, A.fint);
    A.en.plastic += S.e_plastic;
    double dWi = 0, dWe = 0;
    for (int i = 0; i < 6 * nn; ++i) {
      const double d = A.du[i];
      dWi += 0.5 * (fint_c[i] + A.fint[i]) * d;
      dWe += 0.5 * (fext_c[i] + A.fext[i]) * d;
      if (A.fixed[i]) dWe += 0.5 * (fint_c[i] - fext_c[i] + A.fint[i] - A.fext[i]) * d;  // work of the reactions
      A.rf[i] = A.fixed[i] ? A.fint[i] - A.fext[i] : 0.0;
      A.u[i] += A.du[i];
    }
    A.en.internal += dWi;
    A.en.external += dWe;
    fint_c = A.fint;
    fext_c = A.fext;
    double fsum = 0;
    int fcnt = 0;
    for (int i = 0; i < 6 * nn; ++i) if (A.fint[i] != 0) { fsum += std::fabs(A.fint[i]); ++fcnt; }
    if (fcnt) { qtilde = (qtilde * nq + fsum / fcnt) / (nq + 1); ++nq; }
    t = tn;
    A.total_time = ttn;
    ++A.total_increments;
    W.sta_line("%5d %5d %3d %6d %5d %10.3E %10.3E %10.3E", si + 1, inc, 1, 0, it - 1, A.total_time, t, dt);
    info("  increment %4d  t=%.4g  dt=%.4g  iterations %d", inc, t, dt, it - 1);
    W.history(si, inc, t, A.total_time, dt, 0.0);
    bool last = t >= st.period * (1 - 1e-12);
    const int ffreq = st.output.field_frequency;
    if (ffreq > 0 ? (inc % ffreq == 0 || last) : (nfr <= 0 || t >= next_frame - 1e-9 * st.period || last)) {
      W.frame(si, st.name, t, A.total_time);
      if (nfr > 0) while (next_frame <= t + 1e-9 * st.period) next_frame += fr_dt;
    }
    // automatic incrementation
    if (cut_to_frame) dt = dt_auto;
    if (it - 1 <= 4) ++fast; else fast = 0;
    if (fast >= 2) { dt = std::min(dt * 1.5, st.dtmax); fast = 0; }
    else if (it - 1 > 8) dt = std::max(dt * 0.75, st.dtmin);
  }
  std::fill(A.du.begin(), A.du.end(), 0.0);
  W.sta_line("\n THE STEP %d HAS COMPLETED: %d increments, wallclock %.1f s", si + 1, inc, wall() - t0wall);
}

}  // namespace abq
