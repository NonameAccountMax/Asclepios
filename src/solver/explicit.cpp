// Explicit dynamics: central differences with half-step velocities (leapfrog), lumped masses,
// element-by-element stable increment, variable mass scaling, Rayleigh mass damping.
#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "contact/contact.hpp"
#include "core/binio.hpp"
#include "core/log.hpp"
#include "output/output.hpp"
#include "solver/analysis.hpp"
#include "solver/solvers.hpp"

namespace abq {

namespace {

double wall() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// variable mass scaling: element density multipliers so that dt_e >= target (BELOW MIN) or all = target (UNIFORM)
bool apply_mass_scaling(Analysis& A, const Step::MassScaling& ms, double dt_target) {
  if (!ms.on) return false;
  bool changed = false;
  for (auto& gp : A.groups) {
    auto& G = *gp;
    if (G.type == ElType::CONN3D2) continue;
    for (int e = 0; e < G.size(); ++e) {
      double dtu = G.dt_last[e];  // unscaled estimate
      double f = 1.0;
      if (ms.factor > 0 && ms.dt <= 0) f = ms.factor;
      else if (dtu < dt_target || !ms.below_min) f = std::max(1.0, (dt_target / dtu) * (dt_target / dtu));
      if (!ms.below_min && ms.dt > 0) f = (dt_target / dtu) * (dt_target / dtu);
      if (std::fabs(f - G.msf[e]) > 1e-12 * f) { G.msf[e] = f; changed = true; }
    }
  }
  if (changed) A.lump_masses();
  return changed;
}

double total_mass(const Analysis& A) {
  double s = 0;
  for (double x : A.mt) s += x;
  return s;
}

}  // namespace

// loop state of an explicit step, stored in restart files
struct Cursor {
  double t, dt_prev, dt, dt_el, next_frame, next_hist, next_print, mass0, dt_last_inc;
  int inc, ms_counter, nfail_tot, nfail_inc, crit;
};

// Global estimate of the stable increment: largest eigenvalue of M^-1 K by power iteration, K x from central
// differences of the internal forces at the current state (no viscous terms). The element-by-element bound is
// conservative; the global value (with the Abaqus safety factor, see global_update) is used once two successive
// estimates agree within 0.5 %, as a ratio to the element bound so that later changes of the element bound are
// followed.
struct GlobalDt {
  std::vector<double> x;  // iteration vector (6 nn)
  double lam = 0, ratio = 1;
  int ok = 0, countdown = 0;
};

double power_iterations(Analysis& A, std::vector<double>& x, int iters, double dt) {
  const int nd = 6 * A.nn;
  std::vector<double> fp(nd), fm(nd), du_save = A.du;
  double lam = 0, h = 0;
  for (int i = 0; i < 3 * A.nn; ++i) h = std::max(h, std::fabs(A.X[i]));
  h = 1e-7 * std::max(h, 1e-30);
  for (int it = 0; it < iters; ++it) {
    double xmax = 0;
    for (int i = 0; i < nd; ++i) {
      const int n = i / 6, k = i % 6;
      const double mm = k < 3 ? A.mt[n] : A.mr[n];
      if (A.fixed[i] || !(mm > 0)) x[i] = 0;
      xmax = std::max(xmax, std::fabs(x[i]));
    }
    if (!(xmax > 0)) return 0;
    for (int i = 0; i < nd; ++i) x[i] /= xmax;
    for (int i = 0; i < nd; ++i) A.du[i] = h * x[i];
    A.eval_elements(dt, 0, 0, false, false, fp);
    for (int i = 0; i < nd; ++i) A.du[i] = -h * x[i];
    A.eval_elements(dt, 0, 0, false, false, fm);
    double xkx = 0, xmx = 0;
    for (int i = 0; i < nd; ++i) {
      const int n = i / 6, k = i % 6;
      const double mm = k < 3 ? A.mt[n] : A.mr[n];
      if (A.fixed[i] || !(mm > 0)) { fp[i] = 0; continue; }
      const double kx = (fp[i] - fm[i]) / (2 * h);
      xkx += x[i] * kx;
      xmx += x[i] * mm * x[i];
      fp[i] = kx / mm;  // next iterate M^-1 K x
    }
    lam = xmx > 0 ? xkx / xmx : 0;
    x.swap(fp);
  }
  A.du = du_save;
  return lam;
}

void run_explicit_step(Analysis& A, const Step& st, int si) {
  const int nn = A.nn;
  const bool resume = A.resume_step == si && !A.restart_in.empty();
  const bool first = !resume && A.total_increments == 0 && si == 0;
  A.nlgeom = st.nlgeom;
  if (!resume) A.begin_step(st, si, st.amp_ramp == 1);
  if (!A.contact) A.contact = std::make_unique<GeneralContact>(A);
  if (resume) {  // contact definition in force at this step
    for (int k = si; k >= 0; --k)
      if (A.m.steps[k].contact_defined) { A.contact->define(A.m.steps[k]); break; }
  } else if (st.contact_defined) {
    A.contact->define(st);
  }
  if (!resume) {
    A.temperatures(0.0, A.total_time);
    A.temp_prev = A.temp;
  }

  OutputWriter& W = *A.writer;
  W.sta_header(si, st.name, true);
  info("\n STEP %d  %s  (dynamic, explicit)  period %g", si + 1, st.name.c_str(), st.period);

  static double mass0 = 0;
  static double dt_last_inc = 0;
  // FREQUENCY applies to Abaqus/Standard field output only: ignored here, as Abaqus/Explicit does (20 intervals)
  if (st.output.field_frequency > 0)
    warn("step %s: *OUTPUT, FIELD, FREQUENCY ignored in an explicit step (NUMBER INTERVAL or TIME INTERVAL)", st.name.c_str());
  const int nfr = st.output.field_intervals > 0 || st.output.field_time_interval > 0 ? st.output.field_intervals : 20;
  const double fr_dt = nfr > 0 ? st.period / nfr : (st.output.field_time_interval > 0 ? st.output.field_time_interval : st.period);
  const int nhist = std::max(1, st.output.history_intervals);
  const double hist_dt = st.period / nhist;
  std::vector<double> fint_old, fext_old, fcon_old, fd(6 * nn, 0.0);
  double next_frame, next_hist, next_print, dt_prev, t, dt_el, dt;
  int crit = -1, inc, ms_counter = 0, nfail_tot = 0, nfail_inc = 0;
  double ke_rot_last = 0;
  GlobalDt gdt;
  // models with rotational dofs (shells, beams) keep the element bound: their rotational modes are limited by the
  // damping of the bulk viscosity on curvature rates, which the element bound includes and the global one does not
  bool has_rot = false;
  for (int n = 0; n < nn; ++n) has_rot = has_rot || A.node_rot[n];
  const bool use_global = !st.element_by_element && !(st.fixed_dt && st.user_dt > 0) && !has_rot;
  // damping of the highest mode by the linear bulk viscosity, as in the element bound
  const double gdamp = std::sqrt(1 + st.bv1 * st.bv1) - st.bv1;
  auto global_update = [&](int iters) {
    if ((int)gdt.x.size() != 6 * nn) {  // start vector with signs alternating pseudo-randomly
      gdt.x.assign(6 * nn, 0.0);
      unsigned s0 = 4096;
      for (int i = 0; i < 6 * nn; ++i) { s0 = (s0 * 150889u + 1u) % 714025u; gdt.x[i] = (s0 % 2 ? 1.0 : -1.0) * (0.5 + s0 / 1428050.0); }
    }
    double lam = power_iterations(A, gdt.x, iters, dt);
    bool agree = gdt.lam > 0 && lam > 0 && std::fabs(lam - gdt.lam) <= 0.005 * lam;
    gdt.lam = lam;
    if (agree && lam > 0 && dt_el > 0) {
      // 0.917 of the converged estimate, and the element bound unless the global one is at least 5 % larger:
      // the global increments of Abaqus/Explicit (exact largest eigenvalue of clamped C3D8R blocks computed
      // separately: ratios 0.914 to 0.919; global estimate dropped below 1.05 times the element bound)
      const double r = 0.917 * 2.0 / std::sqrt(lam) * gdamp / dt_el;
      gdt.ratio = r < 1.05 ? 1.0 : std::min(1.5, r);
      gdt.ok = 1;
    }
  };
  // penalty contact adds up to ~0.6/dt^2 to the nodal frequencies: keep a margin
  double cfac = 1.0;

  if (resume) {
    BinIn in(A.restart_in);
    in.tag("asclepios-restart");
    int sfile;
    in.pod(sfile);
    A.load_state(in);
    A.contact->load(in);
    in.tag("explicit");
    Cursor c;
    in.pod(c);
    in.vec(fint_old); in.vec(fext_old); in.vec(fcon_old);
    in.vec(gdt.x); in.pod(gdt.lam); in.pod(gdt.ratio); in.pod(gdt.ok); in.pod(gdt.countdown);
    t = c.t; dt_prev = c.dt_prev; dt = c.dt; dt_el = c.dt_el; next_frame = c.next_frame; next_hist = c.next_hist;
    next_print = c.next_print; mass0 = c.mass0; dt_last_inc = c.dt_last_inc; inc = c.inc; ms_counter = c.ms_counter;
    nfail_tot = c.nfail_tot; nfail_inc = c.nfail_inc; crit = c.crit;
    cfac = A.contact->active() ? 0.93 : 1.0;
    A.next_save = A.save_interval > 0 ? (std::floor(A.total_time / A.save_interval + 1e-9) + 1) * A.save_interval : 0;
    info("  resumed from %s at step time %.6e (total %.6e, increment %d)", A.restart_in.c_str(), t, A.total_time, inc);
    W.frame(si, st.name, t, A.total_time);
  } else {
    if (first) {
      std::fill(A.du.begin(), A.du.end(), 0.0);
      A.lump_masses();
      // initial element pass (zero increment) for stable increments and initial forces; the element increments
      // include the damping of the linear bulk viscosity (they drive the variable mass scaling)
      A.eval_elements(0.0, st.bv1, st.bv2, true, true, A.fint);
      mass0 = total_mass(A);
    }
    // mass scaling at the start of the step
    if (st.mscale.on) {
      double target = st.mscale.dt;
      apply_mass_scaling(A, st.mscale, target);
      info("  mass scaling: target dt %g, mass change %+.3f %%", target, 100 * (total_mass(A) / mass0 - 1));
    } else if (si > 0) {
      // mass scaling of a previous step is not propagated
      bool any = false;
      for (auto& gp : A.groups) for (auto& f : gp->msf) if (f != 1.0) { f = 1.0; any = true; }
      if (any) A.lump_masses();
    }
    // stable increment for the current masses (zero increment pass, nothing committed)
    std::fill(A.du.begin(), A.du.end(), 0.0);
    std::vector<double> scratch(6 * nn, 0.0);
    auto fresh = A.eval_elements(0.0, st.bv1, st.bv2, true, false, scratch);
    fint_old = A.fint;
    A.external_forces(0.0, A.total_time, A.fext);
    fext_old = A.fext;
    if (first) {
      if (A.contact->active()) A.contact->compute(1e-30, A.fcon, A.en, 0);
      fcon_old = A.fcon;
      // viscous forces of the initial velocities: an element pass on the increment v h over a vanishing time h,
      // less the forces of the initial state. Abaqus/Explicit at increment 0 (Taylor impact: reaction and first
      // increment exact): linear bulk viscosity on the stable-increment length in the reactions, and the damping
      // of the quadratic term at the initial compression rate in the first stable increment
      std::vector<double> fvis(6 * nn, 0.0);
      if ((st.bv1 > 0 || st.bv2 > 0) && std::any_of(A.v.begin(), A.v.end(), [](double x) { return x != 0.0; })) {
        const double h = 1e-6 * fresh.dt_min;
        for (int i = 0; i < 6 * nn; ++i) A.du[i] = h * A.v[i];
        A.initial_pass = true;
        auto vis = A.eval_elements(h, st.bv1, st.bv2, true, false, fvis);
        A.initial_pass = false;
        std::fill(A.du.begin(), A.du.end(), 0.0);
        for (int i = 0; i < 6 * nn; ++i) fvis[i] -= A.fint[i];
        if (vis.dt_min < fresh.dt_min) { fresh.dt_min = vis.dt_min; fresh.dt_elem = vis.dt_elem; }
      }
      for (int n = 0; n < nn; ++n)
        for (int k = 0; k < 6; ++k) {
          double mm = k < 3 ? A.mt[n] : A.mr[n];
          int i = 6 * n + k;
          A.a[i] = (mm > 0 && !A.fixed[i]) ? (A.fext[i] + A.fcon[i] - A.fint[i]) / mm : 0.0;
          // reactions of the initial state (as Abaqus reports them at increment 0: e.g. the share of a suddenly
          // applied load carried by the prescribed nodes)
          A.rf[i] = A.fixed[i] ? -(A.fext[i] + A.fcon[i] - A.fint[i] - fvis[i]) : 0.0;
        }
      W.frame(si, st.name, 0.0, A.total_time);
    } else {
      fcon_old = A.fcon;
    }
    next_frame = fr_dt;
    next_hist = 0;
    next_print = 0.05 * st.period;
    dt_prev = first ? 0.0 : dt_last_inc;
    t = 0;
    inc = 0;
    dt_el = fresh.dt_min;
    crit = fresh.dt_elem;
    cfac = A.contact->active() ? 0.93 : 1.0;
    dt = st.fixed_dt && st.user_dt > 0 ? st.user_dt : dt_el * std::min(1.0, st.dt_scale) * cfac;
    if (A.save_interval > 0 && A.next_save <= A.total_time) A.next_save = (std::floor(A.total_time / A.save_interval + 1e-9) + 1) * A.save_interval;
  }
  auto save_restart = [&]() {
    std::string path = A.job + ".rst";
    {
      BinOut o(path + ".tmp");
      o.tag("asclepios-restart");
      o.pod(si);
      A.save_state(o);
      A.contact->save(o);
      o.tag("explicit");
      Cursor c{t, dt_prev, dt, dt_el, next_frame, next_hist, next_print, mass0, dt_last_inc, inc, ms_counter, nfail_tot, nfail_inc, crit};
      o.pod(c);
      o.vec(fint_old); o.vec(fext_old); o.vec(fcon_old);
      o.vec(gdt.x); o.pod(gdt.lam); o.pod(gdt.ratio); o.pod(gdt.ok); o.pod(gdt.countdown);
    }
    std::rename((path + ".tmp").c_str(), path.c_str());
    // keep a copy per save time as well (job_t<time>.rst) for debugging restarts
    if (std::getenv("ABQ_KEEP_RESTARTS")) {
      char buf[64];
      std::snprintf(buf, sizeof buf, "_t%.6f.rst", A.total_time);
      std::FILE* src = std::fopen(path.c_str(), "rb");
      std::FILE* dst = std::fopen((A.job + buf).c_str(), "wb");
      if (src && dst) { char b[1 << 16]; size_t k; while ((k = std::fread(b, 1, sizeof b, src)) > 0) std::fwrite(b, 1, k, dst); }
      if (src) std::fclose(src);
      if (dst) std::fclose(dst);
    }
  };
  const double t0wall = wall();
  int trace = -1;
  if (const char* tr = std::getenv("ABQ_TRACE_NODE")) {
    for (int n = 0; n < nn; ++n) if (A.m.node_name(n) == tr) trace = n;
  }

  while (t < st.period * (1 - 1e-12)) {
    // increment size (avoid a tiny last increment: spread the remainder over up to 10 increments)
    double rem = st.period - t;
    if (dt >= rem) dt = rem;
    else if (rem < 10 * dt) { double nrem = std::ceil(rem / dt); dt = rem / nrem; }
    // with TIME MARKS=YES the field output times are reached exactly: the increments before an output time are made
    // equal (remainder spread over up to 10 increments; shortening only the last one, periodically, made an
    // undamped beam model unstable); otherwise the frame is written at the end of the increment that reaches or
    // passes the output time
    const double dt_full = dt;
    if (st.output.time_marks && next_frame < st.period * (1 - 1e-12) && t + 10 * dt > next_frame &&
        next_frame - t > 1e-6 * dt) {
      const double r = next_frame - t;
      dt = r / std::ceil(r / dt * (1 - 1e-9));
    }
    const double tn = t + dt, ttn = A.total_time + dt;

    // velocities at n+1/2
    double hdt = 0.5 * (dt_prev + dt);
#pragma omp parallel for num_threads(A.nthreads) schedule(static)
    for (int i = 0; i < 6 * nn; ++i) A.v[i] += hdt * A.a[i];
    // prescribed motion
    for (auto& b : A.loads.bcs) {
      int i = 6 * b.node + b.dof;
      double target = A.loads.value(b, tn, ttn, A.m);
      if (b.velocity) A.v[i] = target;
      // with an amplitude, the velocity follows the amplitude (magnitude x rate of the amplitude): a node whose
      // condition is redefined keeps its current position as the reference, as in Abaqus/Explicit
      else if (b.amp >= 0) A.v[i] = (target - A.loads.value(b, t, A.total_time, A.m)) / dt;
      else A.v[i] = (target - A.u[i]) / dt;
    }
    for (auto& b : A.m.model_bcs) { int i = 6 * b.node + b.dof; A.v[i] = (b.value - A.u[i]) / dt; }
#pragma omp parallel for num_threads(A.nthreads) schedule(static)
    for (int i = 0; i < 6 * nn; ++i) A.du[i] = dt * A.v[i];

    // predefined fields at the end of the increment
    A.temp_prev = A.temp;
    A.temperatures(tn, ttn);

    // internal forces
    Analysis::ElemSum S = A.eval_elements(dt, st.bv1, st.bv2, true, true, A.fint);
    nfail_inc += S.failed > 0;
    nfail_tot += S.failed;
    A.external_forces(tn, ttn, A.fext);
    if (A.contact->active()) A.contact->compute(dt, A.fcon, A.en, A.total_increments);
    // damping forces (mass proportional, half-step velocity)
#pragma omp parallel for num_threads(A.nthreads) schedule(static)
    for (int n = 0; n < nn; ++n)
      for (int k = 0; k < 3; ++k) fd[6 * n + k] = A.cdamp[n] * A.v[6 * n + k];

    // energies (trapezoidal rule on nodal forces)
    double dWi = 0, dWe = 0, dWc = 0, dWd = 0;
#pragma omp parallel for num_threads(A.nthreads) reduction(+ : dWi, dWe, dWc, dWd) schedule(static)
    for (int i = 0; i < 6 * nn; ++i) {
      double d = A.du[i];
      dWi += 0.5 * (fint_old[i] + A.fint[i]) * d;
      if (A.fixed[i]) dWe += 0.5 * (fint_old[i] - fext_old[i] - fcon_old[i] + A.fint[i] - A.fext[i] - A.fcon[i]) * d;  // reactions
      dWe += 0.5 * (fext_old[i] + A.fext[i]) * d;
      dWc += 0.5 * (fcon_old[i] + A.fcon[i]) * d;
      dWd += fd[i] * d;
    }
    A.en.viscous += S.e_visc;
    A.en.plastic += S.e_plastic;
    A.en.artificial += S.e_art;
    A.en.internal += dWi - S.e_visc;
    A.en.external += dWe;
    A.en.contact += -dWc;  // work absorbed by contact (penalty energy + friction)
    A.en.damping += dWd;

    if (trace >= 0) {
      int n = trace;
      std::fprintf(stderr, "TRACE t=%.6e dt=%.3e v=(%.3e %.3e %.3e) w=(%.3e %.3e %.3e) fint=(%.3e %.3e %.3e) m=(%.3e %.3e %.3e) fext=(%.3e %.3e %.3e) fcon=(%.3e %.3e %.3e) m=%.3e I=%.3e\n",
                   tn, dt, A.v[6*n], A.v[6*n+1], A.v[6*n+2], A.v[6*n+3], A.v[6*n+4], A.v[6*n+5], A.fint[6*n], A.fint[6*n+1], A.fint[6*n+2],
                   A.fint[6*n+3], A.fint[6*n+4], A.fint[6*n+5], A.fext[6*n], A.fext[6*n+1], A.fext[6*n+2], A.fcon[6*n], A.fcon[6*n+1], A.fcon[6*n+2], A.mt[n], A.mr[n]);
    }
    // accelerations; reactions on prescribed dofs
    double ke = 0, ker = 0;
#pragma omp parallel for num_threads(A.nthreads) schedule(static)
    for (int n = 0; n < nn; ++n)
      for (int k = 0; k < 6; ++k) {
        int i = 6 * n + k;
        double r = A.fext[i] + A.fcon[i] - A.fint[i] - fd[i];
        double mm = k < 3 ? A.mt[n] : A.mr[n];
        if (A.fixed[i]) { A.rf[i] = -r; A.a[i] = 0; }
        else { A.rf[i] = 0; A.a[i] = mm > 0 ? r / mm : 0.0; }
        if (std::isnan(A.a[i])) A.a[i] = 0;
      }
    if (!A.ties.empty()) A.tie_accelerations();
#pragma omp parallel for num_threads(A.nthreads) reduction(+ : ke, ker) schedule(static)
    for (int n = 0; n < nn; ++n)
      for (int k = 0; k < 6; ++k) {
        int i = 6 * n + k;
        double mm = k < 3 ? A.mt[n] : A.mr[n];
        double vf = A.v[i] + 0.5 * dt * A.a[i];
        ke += 0.5 * mm * vf * vf;
        if (k >= 3) ker += 0.5 * mm * vf * vf;
      }
#pragma omp parallel for num_threads(A.nthreads) schedule(static)
    for (int i = 0; i < 6 * nn; ++i) A.u[i] += A.du[i];
    fint_old = A.fint;
    fext_old = A.fext;
    fcon_old = A.fcon;
    A.en.kinetic = ke;
    ke_rot_last = ker;
    t = tn;
    A.total_time = ttn;
    ++inc;
    ++A.total_increments;
    dt_prev = dt;
    dt_last_inc = dt;

    if (!std::isfinite(ke) || ke > 1e30) {
      int worst = -1; double vmax = -1;
      for (int n = 0; n < nn; ++n) {
        double vv = std::fabs(A.v[6 * n]) + std::fabs(A.v[6 * n + 1]) + std::fabs(A.v[6 * n + 2]);
        if (!std::isfinite(vv)) vv = 1e300;
        if (vv > vmax) { vmax = vv; worst = n; }
      }
      std::string els;
      for (auto& gp : A.groups)
        for (int e = 0; e < gp->size(); ++e)
          for (int I = 0; I < gp->nen; ++I)
            if (gp->conn[e * gp->nen + I] == worst) els += std::string(" ") + eltype_name(gp->type) + ":" + A.m.elem_name(gp->ids[e]);
      throw std::runtime_error("the analysis is unstable at time " + std::to_string(A.total_time) + " (increment " +
                               std::to_string(inc) + "); fastest node " + A.m.node_name(worst) + ", elements" + els);
    }

    // mass scaling update
    if (st.mscale.on && st.mscale.frequency > 0 && ++ms_counter >= st.mscale.frequency) {
      ms_counter = 0;
      apply_mass_scaling(A, st.mscale, st.mscale.dt);
    }
    // next increment
    crit = S.dt_elem;
    dt_el = S.dt_min;
    if (use_global && --gdt.countdown <= 0) {  // a few iterations at the start, then every 20 increments
      global_update(gdt.x.empty() ? 12 : 2);
      gdt.countdown = gdt.ok ? 20 : 1;
    }
    double dt_next = st.fixed_dt && st.user_dt > 0 ? st.user_dt
                                                   : dt_el * (gdt.ok ? gdt.ratio : 1.0) * std::min(1.0, st.dt_scale) * cfac;
    // limit the growth of the increment (a sudden x10 jump after a stiff transient destabilises the
    // viscous terms, which are evaluated with the rate du/dt of the new increment)
    if (!(st.fixed_dt && st.user_dt > 0)) dt_next = std::min(dt_next, 1.25 * dt_full);
    dt = dt_next;

    bool last = t >= st.period * (1 - 1e-12);
    if (t >= next_hist - 1e-12 * st.period || last) {
      W.history(si, inc, t, A.total_time, dt_prev, ke);
      double pm = 100 * (total_mass(A) / mass0 - 1);
      W.sta_line("%6d %10d %10.3E %10.3E %10.1f %10.3E %10s %11.3E %11.3E %9.2E %8d", si + 1, inc, t, A.total_time, wall() - t0wall,
                 dt_prev, crit >= 0 ? A.m.elem_name(crit).c_str() : "-", ke, A.en.total(ke), pm, A.contact->ncontacts());
      while (next_hist <= t + 1e-12 * st.period) next_hist += hist_dt;
    }
    if (t >= next_frame - 1e-9 * st.period || last) {
      W.frame(si, st.name, t, A.total_time);
      while (next_frame <= t + 1e-9 * st.period) next_frame += fr_dt;
    }
    if (t >= next_print || last) {
      info("  t=%10.4e  inc %9d  dt %10.3e  KE %10.3e (rot %8.2e)  ETOTAL %10.3e  wall %7.1fs", t, inc, dt_prev, ke, ke_rot_last, A.en.total(ke), wall() - t0wall);
      while (next_print <= t) next_print += 0.05 * st.period;
    }
    if (A.save_interval > 0 && A.total_time >= A.next_save - 1e-12 && !last) {
      save_restart();
      while (A.next_save <= A.total_time + 1e-12) A.next_save += A.save_interval;
    }
  }
  if (nfail_tot) warn("step %d: material constitutive routines failed to converge at %d points over %d increments", si + 1, nfail_tot, nfail_inc);
  W.sta_line("\n THE STEP %d HAS COMPLETED: %d increments, wallclock %.1f s", si + 1, inc, wall() - t0wall);
}

}  // namespace abq
