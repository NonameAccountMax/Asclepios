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
#include <cstdlib>

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

// equations of a dof with their coefficients: itself, or the terms of its link (prescribed terms excluded)
int eq_terms(const Analysis& A, const DofMap& dm, int i, int* eqs, double* cs) {
  const int l = A.link_of[i];
  if (l < 0) { if (dm.eq[i] < 0) return 0; eqs[0] = dm.eq[i]; cs[0] = 1.0; return 1; }
  int n = 0;
  const auto& L = A.links[l];
  for (int t = L.first; t < L.first + L.count && n < 64; ++t) {
    const int g = dm.eq[A.link_terms[t].dof];
    if (g < 0) continue;
    eqs[n] = g; cs[n] = A.link_terms[t].c; ++n;
  }
  return n;
}

// Contact in the implicit procedures: penalty on the contact points of the general contact (facets only), normal
// stiffness 7 times the representative stiffness of the underlying elements (largest translational diagonal term
// of the element tangents at the slave node and at the master nodes, the softer side; the factor reproduces the
// reaction and the contact energy of a punch pressed on a block in the reference solver, whose default is 10 times
// its own measure of the element stiffness), starting quadratically over 1e-6 of the mean facet length; Coulomb
// friction with a penalty limited to an elastic slip of 0.5 % of the mean facet length. Contact forces enter the
// internal forces (residual and reactions); their tangent is k g g^T in the normal direction (and in the tangent
// plane while sticking); points within 1e-3 of the mean facet length of contact keep a stiffness in the tangent
// (stable active set during the iterations).
struct ImplicitContact {
  GeneralContact* gc = nullptr;
  std::vector<GeneralContact::Point> pts;
  std::vector<double> kn, kt;      // per point: normal and tangential tangent stiffnesses
  std::vector<Vec3> ftan;          // per point: tangential force
  std::vector<Vec3> ft_c;          // per node: committed tangential force (friction)
  std::vector<double> kd;          // per node: representative translational stiffness of its elements
  double mu = 0, slip = 0;
  double epot = 0;                 // elastic energy of the normal penalty (last evaluation)
  double d0 = 0;                   // overclosure over which the penalty force starts quadratically
  const SurfaceInteraction* soft = nullptr;
  double clear = 0;
  bool active() const { return gc != nullptr; }
};

void contact_forces(Analysis& A, ImplicitContact& C, std::vector<double>& fcon, bool commit,
                    std::vector<double>* ftan = nullptr) {
  std::fill(fcon.begin(), fcon.end(), 0.0);
  if (ftan) std::fill(ftan->begin(), ftan->end(), 0.0);
  // points near contact keep their stiffness in the tangent (no force): stable active set during the iterations
  C.d0 = 1e-6 * C.gc->mean_facet_length();
  C.gc->implicit_points(C.pts, 1e-3 * C.gc->mean_facet_length());
  const size_t np = C.pts.size();
  C.kn.assign(np, 0.0); C.kt.assign(np, 0.0); C.ftan.assign(np, Vec3());
  C.epot = 0;
  std::vector<Vec3> ft_new;
  if (commit) ft_new.assign(A.nn, Vec3());
  for (size_t i = 0; i < np; ++i) {
    const auto& p = C.pts[i];
    const double ks = C.kd[p.slave];
    double km = 0;
    for (int a = 0; a < p.nm; ++a) km += p.w[a] * C.kd[p.mn[a]];
    const double krep = ks > 0 && km > 0 ? std::min(ks, km) : std::max(ks, km);
    double f, kk;
    if (C.soft) {
      double dp = 0;
      f = p.area * C.soft->pressure(p.pen - C.clear, &dp);
      kk = p.area * dp;
    } else {
      // linear penalty with a quadratic start over a small overclosure d0 (continuous slope: Newton iterations
      // do not chatter on the nodes at the edge of the contact zone)
      const double k = 7.0 * krep, d0 = C.d0;
      if (p.pen <= 0) { f = 0; kk = k; }
      else if (p.pen < d0) { f = 0.5 * k * p.pen * p.pen / d0; kk = k * p.pen / d0; }
      else { f = k * (p.pen - 0.5 * d0); kk = k; }
    }
    if (!(f > 0)) {  // near contact: stiffness only
      if (!C.soft) C.kn[i] = 0.5 * kk;
      continue;
    }
    C.epot += 0.5 * f * (C.soft ? p.pen - C.clear : p.pen);
    Vec3 F = p.n * f;
    if (C.mu > 0 && C.slip > 0) {  // Coulomb friction, elastic slip up to C.slip
      Vec3 dus = load3(&A.du[6 * p.slave]);
      for (int a = 0; a < p.nm; ++a) dus = dus - load3(&A.du[6 * p.mn[a]]) * p.w[a];
      dus = dus - p.n * dot(dus, p.n);
      Vec3 f0 = C.ft_c[p.slave];
      f0 = f0 - p.n * dot(f0, p.n);
      const double kt = C.mu * f / C.slip;
      Vec3 ft = f0 - dus * kt;
      const double lim = C.mu * f, mag = norm(ft);
      if (mag > lim && mag > 0) ft = ft * (lim / mag);
      else C.kt[i] = kt;
      C.ftan[i] = ft;
      F = F + ft;
      if (commit) ft_new[p.slave] = ft;
    }
    C.kn[i] = kk;
    // internal force convention: the contact force F acts on the slave, -w F on the master nodes
    for (int k = 0; k < 3; ++k) {
      fcon[6 * p.slave + k] -= F[k];
      for (int a = 0; a < p.nm; ++a) fcon[6 * p.mn[a] + k] += p.w[a] * F[k];
      if (ftan) {
        (*ftan)[6 * p.slave + k] -= C.ftan[i][k];
        for (int a = 0; a < p.nm; ++a) (*ftan)[6 * p.mn[a] + k] += p.w[a] * C.ftan[i][k];
      }
    }
  }
  if (commit) C.ft_c = ft_new;
}

void contact_tangent(const Analysis& A, const DofMap& dm, const ImplicitContact& C,
                     std::vector<Eigen::Triplet<double>>& trip) {
  for (size_t i = 0; i < C.pts.size(); ++i) {
    const auto& p = C.pts[i];
    if (C.kn[i] <= 0 && C.kt[i] <= 0) continue;
    // directions: the normal, and the tangent plane while sticking
    Vec3 dirs[3];
    double kdir[3];
    int nd = 0;
    dirs[nd] = p.n; kdir[nd++] = C.kn[i];
    if (C.kt[i] > 0) {
      Vec3 t1 = std::fabs(p.n.x) < 0.9 ? cross(p.n, Vec3(1, 0, 0)) : cross(p.n, Vec3(0, 1, 0));
      t1 = normalized(t1);
      dirs[nd] = t1; kdir[nd++] = C.kt[i];
      dirs[nd] = cross(p.n, t1); kdir[nd++] = C.kt[i];
    }
    int nodes[5];
    double wts[5];
    nodes[0] = p.slave; wts[0] = 1.0;
    for (int a = 0; a < p.nm; ++a) { nodes[1 + a] = p.mn[a]; wts[1 + a] = -p.w[a]; }
    const int nnod = 1 + p.nm;
    for (int d = 0; d < nd; ++d) {
      int eq[15];
      double gv[15];
      int ng = 0;
      for (int a = 0; a < nnod; ++a)
        for (int k = 0; k < 3; ++k) {
          int e2[64];
          double c2[64];
          const int nt = eq_terms(A, dm, 6 * nodes[a] + k, e2, c2);
          for (int t = 0; t < nt && ng < 15; ++t) { eq[ng] = e2[t]; gv[ng] = wts[a] * dirs[d][k] * c2[t]; ++ng; }
        }
      for (int x = 0; x < ng; ++x)
        for (int y = 0; y < ng; ++y) trip.push_back({eq[x], eq[y], kdir[d] * gv[x] * gv[y]});
    }
  }
}

// representative stiffness of each node: largest translational diagonal term of the element tangents
void node_stiffness(const Analysis& A, const DofMap& dm, const std::vector<Eigen::Triplet<double>>& trip,
                    std::vector<double>& kd) {
  std::vector<double> diag(dm.neq, 0.0);
  for (const auto& t : trip) if (t.row() == t.col()) diag[t.row()] += t.value();
  kd.assign(A.nn, 0.0);
  for (int n = 0; n < A.nn; ++n)
    for (int k = 0; k < 3; ++k) {
      const int e = dm.eq[6 * n + k];
      if (e >= 0) kd[n] = std::max(kd[n], diag[e]);
    }
}

double dt0_for_contact(const Step& st) { return std::min(st.dt0, st.period); }

}  // namespace

DofMap make_dofmap(const Analysis& A) {
  const int nn = A.nn;
  std::vector<char> has_el(nn, 0), is_term(6 * nn, 0);
  for (auto& gp : A.groups) for (int n : gp->conn) has_el[n] = 1;
  for (int j : A.link_term_dofs) is_term[j] = 1;  // e.g. the reference node of a kinematic coupling
  DofMap dm;
  dm.eq.assign(6 * nn, -1);
  for (int n = 0; n < nn; ++n)
    for (int k = 0; k < 6; ++k) {
      const bool pres = k == 3 && A.any_pres && A.node_pres[n];  // C3D4H pressure unknown
      if (!is_term[6 * n + k] && (!has_el[n] || (k >= 3 && !A.node_rot[n] && !pres))) continue;
      if (A.fixed[6 * n + k]) continue;
      dm.eq[6 * n + k] = dm.neq++;
    }
  // dependent dofs (surface ties, JOIN connectors) are eliminated: their rows and columns go to their terms
  if (!A.links.empty()) {
    int neq = 0;
    for (int i = 0; i < 6 * nn; ++i) dm.eq[i] = (dm.eq[i] >= 0 && A.link_of[i] < 0) ? neq++ : -1;
    dm.neq = neq;
  }
  return dm;
}

double fd_step(const Analysis& A) {
  double lmin = 1e30;
  for (auto& gp : A.groups)
    for (int e = 0; e < gp->size(); ++e) {
      int a = gp->conn[e * gp->nen], b = gp->conn[e * gp->nen + 1];
      double l = std::sqrt(std::pow(A.X[3 * a] - A.X[3 * b], 2) + std::pow(A.X[3 * a + 1] - A.X[3 * b + 1], 2) +
                           std::pow(A.X[3 * a + 2] - A.X[3 * b + 2], 2));
      if (l > 0) lmin = std::min(lmin, l);
    }
  return 1e-7 * lmin;
}

void assemble_tangent(Analysis& A, const DofMap& dm, double dt, double hfd, std::vector<Eigen::Triplet<double>>& trip) {
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
        in.fv1 = fv.data(); in.dt = dt; in.dynamic = false; in.nlgeom = A.nlgeom; in.dt_mat = dt * A.mat_time;
        element_tangent(G, e, in, due.data(), nen, hfd, K.data());
        int ei[64], ej[64];
        double ci[64], cj[64];
        for (int i = 0; i < nd; ++i) {
          const int ni = eq_terms(A, dm, 6 * cn[i / 6] + i % 6, ei, ci);
          if (!ni) continue;
          for (int j = 0; j < nd; ++j) {
            const double kij = K[i * nd + j];
            if (kij == 0) continue;
            const int nj = eq_terms(A, dm, 6 * cn[j / 6] + j % 6, ej, cj);
            for (int a = 0; a < ni; ++a)
              for (int b = 0; b < nj; ++b) loc.push_back({ei[a], ej[b], ci[a] * cj[b] * kij});
          }
        }
      }
    }
    for (auto& l : local) trip.insert(trip.end(), l.begin(), l.end());
  }
}

void dof_masses(Analysis& A, std::vector<double>& mdof) {
  const int nn = A.nn;
  for (auto& gp : A.groups) gp->physical_rotary = true;
  A.lump_masses();
  for (auto& gp : A.groups) gp->physical_rotary = false;
  mdof.assign(6 * nn, 0.0);
  for (int n = 0; n < nn; ++n)
    for (int k = 0; k < 6; ++k)
      mdof[6 * n + k] = k < 3 ? A.mt[n] : A.mr[n] + (A.node_aniso[n] ? A.jrot[6 * (size_t)n + k - 3] : 0.0);
  for (const auto& L : A.links)
    for (int t = L.first; t < L.first + L.count; ++t) mdof[A.link_terms[t].dof] += A.link_terms[t].mc * mdof[L.dof];
  for (const auto& L : A.links) mdof[L.dof] = 0;
}

void run_static_step(Analysis& A, const Step& st, int si) {
  const int nn = A.nn;
  A.nlgeom = st.nlgeom;
  // time-dependent materials: real time in implicit dynamics and *VISCO, instantaneous response in *STATIC (long-term
  // with LONG TERM)
  A.mat_time = st.implicit_dynamic || st.visco ? 1.0 : st.long_term ? 1e30 : 0.0;
  // default amplitude: RAMP in static steps, STEP in dynamic ones
  // (transient coupled temperature-displacement: STEP, as transient heat transfer)
  const bool step_default = st.implicit_dynamic || (st.coupled && !st.heat_steady);
  A.begin_step(st, si, step_default ? st.amp_ramp == 1 : st.amp_ramp != 0);
  // contact: general contact of the step, or the contact pairs treated as general contact between their surfaces
  ImplicitContact IC;
  if (st.contact_defined || !st.contact_pairs.empty()) {
    Step cs = st;
    if (!st.contact_defined) {
      cs.contact_defined = true;
      cs.contact_all_exterior = false;
      cs.contact_surfaces.clear();
      for (const auto& cp : st.contact_pairs) {
        cs.contact_surfaces.push_back(cp.first);
        if (cp.second != cp.first) cs.contact_surfaces.push_back(cp.second);
        if (cp.interaction >= 0) cs.contact_property = cp.interaction;
      }
      warn("*CONTACT PAIR in an implicit step: penalty contact between the surfaces of the pairs");
    }
    if (!A.contact) A.contact = std::make_unique<GeneralContact>(A);
    if (A.contact->define(cs)) {
      IC.gc = A.contact.get();
      IC.mu = IC.gc->friction();
      IC.soft = IC.gc->softened();
      IC.clear = IC.gc->clearance();
      IC.ft_c.assign(nn, Vec3());
    }
  }
  OutputWriter& W = *A.writer;
  W.sta_header(si, st.name, false);
  info("\n STEP %d  %s  (%s, %s)  period %g", si + 1, st.name.c_str(), st.implicit_dynamic ? "implicit dynamic" : st.visco ? "visco" : "static",
       st.nlgeom ? "nlgeom" : "small displacements", st.period);
  if (A.total_increments == 0 && si == 0) W.frame(si, st.name, 0.0, A.total_time);

  const DofMap dm = make_dofmap(A);
  info("  %d equations", dm.neq);
  const double hfd = fd_step(A);
  std::vector<double> fcon(6 * nn, 0.0), fcon_c(6 * nn, 0.0), ftan(6 * nn, 0.0), ftan_c(6 * nn, 0.0);
  if (IC.active()) {  // representative stiffness of the nodes from the element tangents of the start state
    std::fill(A.du.begin(), A.du.end(), 0.0);
    std::vector<Eigen::Triplet<double>> trip;
    assemble_tangent(A, dm, dt0_for_contact(st), hfd, trip);
    node_stiffness(A, dm, trip, IC.kd);
    IC.slip = 0.005 * A.contact->mean_facet_length();
  }

  // implicit dynamics (HHT-alpha): lumped masses, initial accelerations from equilibrium
  const bool dyn = st.implicit_dynamic;
  const double ha = st.hht_alpha, hb = 0.25 * (1 - ha) * (1 - ha), hg = 0.5 - ha;
  std::vector<double> mdof(6 * nn, 0.0), fint0, fext0;
  if (dyn) {
    for (auto& gp : A.groups) gp->physical_rotary = true;
    A.lump_masses();
    for (auto& gp : A.groups) gp->physical_rotary = false;
    for (int n = 0; n < nn; ++n)
      for (int k = 0; k < 6; ++k)
        mdof[6 * n + k] = k < 3 ? A.mt[n] : A.mr[n] + (A.node_aniso[n] ? A.jrot[6 * (size_t)n + k - 3] : 0.0);
    std::fill(A.du.begin(), A.du.end(), 0.0);
    A.spread_links(A.v);
    A.eval_elements(0.0, 0, 0, false, false, A.fint);
    A.external_forces(0.0, A.total_time, A.fext);
    std::vector<double> r0(6 * nn), m0 = mdof;
    for (int i = 0; i < 6 * nn; ++i) r0[i] = A.fext[i] - A.fint[i] - A.damp_coef(i, true) * A.v[i];
    A.gather_links(r0);
    for (const auto& L : A.links)
      for (int t = L.first; t < L.first + L.count; ++t) m0[A.link_terms[t].dof] += A.link_terms[t].mc * mdof[L.dof];
    for (int i = 0; i < 6 * nn; ++i) A.a[i] = (!A.fixed[i] && A.link_of[i] < 0 && m0[i] > 0) ? r0[i] / m0[i] : 0.0;
    A.spread_links(A.a);
    fint0 = A.fint; fext0 = A.fext;
    // mass of each equation (lumped onto the independent dofs)
    for (const auto& L : A.links)
      for (int t = L.first; t < L.first + L.count; ++t) mdof[A.link_terms[t].dof] += A.link_terms[t].mc * mdof[L.dof];
    for (const auto& L : A.links) mdof[L.dof] = 0;
  }

  double t = 0, dt = std::min(st.dt0, st.period);
  int inc = 0, fast = 0;
  double dT_inc = 0;  // coupled temperature-displacement: largest temperature change of the increment
  double kscale = 0;  // largest diagonal term of the first tangent of the increment
  double qtilde = 0;
  int nq = 0;
  const double t0wall = wall();
  std::vector<double> R(6 * nn), du_save, inert(6 * nn, 0.0);
  // converged increment of the free dofs and its length: linear extrapolation in time of the previous increment
  // starts the Newton iterations (static procedure, as the default EXTRAPOLATION=LINEAR of Abaqus/Standard)
  std::vector<double> du_last;
  double dt_last = 0;
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
    if (!dyn && !du_last.empty() && dt_last > 0) {
      const double r = dt / dt_last;
      for (int i = 0; i < 6 * nn; ++i) if (!A.fixed[i] && dm.eq[i] >= 0) A.du[i] = r * du_last[i];
    }
    A.spread_links(A.du);
    A.temp_prev = A.temp;
    if (st.coupled) {  // temperatures of the increment first (independent of the deformation), then the mechanics
      dT_inc = heat_increment(A, st.heat_steady ? 0.0 : dt, tn, ttn);
      if (!st.heat_steady && st.deltmx > 0 && dT_inc > st.deltmx * (1 + 1e-6)) {
        A.temp = A.temp_prev;
        --inc;
        dt *= std::max(0.25, 0.85 * st.deltmx / dT_inc);
        if (dt < st.dtmin) throw std::runtime_error("time increment smaller than the minimum in step " + st.name);
        continue;
      }
    } else A.temperatures(tn, ttn);
    bool converged = false;
    int it = 0;
    double qmax_inc = 0;
    double rmax = 0, cmax = 0, dumax = 0;
    for (it = 1; it <= 16; ++it) {
      // residual
      A.eval_elements(dt, 0, 0, false, false, A.fint);
      if (IC.active()) {
        contact_forces(A, IC, fcon, false);
        for (int i = 0; i < 6 * nn; ++i) A.fint[i] += fcon[i];
      }
      A.external_forces(tn, ttn, A.fext);
      double fsum = 0;
      int fcnt = 0;
      rmax = 0;
      std::vector<double> req(dm.neq, 0.0);  // residual per equation (tied dofs summed)
      bool finite = true;
      if (!dyn) for (int i = 0; i < 6 * nn; ++i) R[i] = A.fext[i] - A.fint[i];
      else {
        // (1 + a)(fext - fint)_{n+1} - a (fext - fint)_n - M a_{n+1} - C ((1 + a) v_{n+1} - a v_n), Newmark a, v
        for (int i = 0; i < 6 * nn; ++i) {
          const double an = (A.du[i] - dt * A.v[i] - dt * dt * (0.5 - hb) * A.a[i]) / (hb * dt * dt);
          const double vn = A.v[i] + dt * ((1 - hg) * A.a[i] + hg * an);
          const double c = A.damp_coef(i, true);
          R[i] = (1 + ha) * (A.fext[i] - A.fint[i]) - ha * (fext0[i] - fint0[i]) - c * ((1 + ha) * vn - ha * A.v[i]);
          inert[i] = (i % 6 < 3 ? A.mt[i / 6] : A.mr[i / 6] + (A.node_aniso[i / 6] ? A.jrot[6 * (size_t)(i / 6) + i % 6 - 3] : 0.0)) * an;
        }
        A.gather_links(inert);
        for (int i = 0; i < 6 * nn; ++i) R[i] -= inert[i];
      }
      A.gather_links(R);
      for (int i = 0; i < 6 * nn; ++i) {
        finite = finite && std::isfinite(R[i]);
        if (dm.eq[i] >= 0) req[dm.eq[i]] += R[i];
        if (i % 6 < 3 || A.node_rot[i / 6]) {  // time-averaged force: all active dofs, reactions included
          double fa = std::fabs(A.fint[i]);
          if (fa > 0) { fsum += fa; ++fcnt; }
        }
      }
      int rarg = -1;
      for (int q = 0; q < dm.neq; ++q) if (std::fabs(req[q]) > rmax) { rmax = std::fabs(req[q]); rarg = q; }
      double qcur = fcnt ? fsum / fcnt : 0;
      if (it == 1) qmax_inc = qcur;
      double qt = nq ? std::max(qtilde, qcur) : qcur;
      // a state released to zero force (e.g. initial stresses relaxed) keeps a reference: the force at the start of
      // the increment
      qt = std::max(qt, 1e-3 * qmax_inc);
      // no force at all (free flight before an impact): the residual is compared with a force scale of the model,
      // 1e-12 x largest stiffness of the last tangent x shortest element edge
      qt = std::max(qt, 1e-12 * kscale * (hfd / 1e-7));
      if (qt <= 0) qt = 1e-300;
      // a diverged iterate (a force not finite, reactions included) is never accepted: NaN would be skipped by
      // the maxima below, and an infinite reaction would make the relative criterion trivially satisfied
      if (!finite || !std::isfinite(rmax) || !std::isfinite(qcur)) break;
      static const bool trace_newton = std::getenv("ABQ_TRACE_NEWTON") != nullptr;
      if (trace_newton) {
        int dof = -1;
        for (int i = 0; i < 6 * nn && rarg >= 0; ++i) if (dm.eq[i] == rarg) { dof = i; break; }
        std::fprintf(stderr, "NEWTON inc %d it %d dt %.3e rmax %.3e (node %d dof %d) qt %.3e cmax %.3e dumax %.3e\n", inc, it, dt, rmax,
                     dof >= 0 ? A.m.node_label[dof / 6] : -1, dof % 6 + 1, qt, cmax, dumax);
      }
      if (it > 1 && rmax <= 0.005 * qt && cmax <= 0.01 * std::max(dumax, 1e-300)) { converged = true; break; }
      if (rmax <= 1e-8 * qt && it > 1) { converged = true; break; }
      if (dm.neq == 0) { converged = true; break; }  // fully prescribed
      // tangent
      std::vector<Eigen::Triplet<double>> trip;
      assemble_tangent(A, dm, dt, hfd, trip);
      if (IC.active()) {
        node_stiffness(A, dm, trip, IC.kd);
        contact_tangent(A, dm, IC, trip);
      }
      if (dyn) {
        for (auto& tr : trip) tr = Eigen::Triplet<double>(tr.row(), tr.col(), (1 + ha) * tr.value());
        for (int i = 0; i < 6 * nn; ++i) {
          if (dm.eq[i] < 0) continue;
          const double c = A.damp_coef(i, true);
          trip.push_back({dm.eq[i], dm.eq[i], mdof[i] / (hb * dt * dt) + (1 + ha) * hg / (hb * dt) * c});
        }
      }
      Eigen::SparseMatrix<double> Kg(dm.neq, dm.neq);
      Kg.setFromTriplets(trip.begin(), trip.end());
      // tiny regularisation of dofs without stiffness (e.g. drilling of isolated nodes)
      for (int i = 0; i < dm.neq; ++i) Kg.coeffRef(i, i) += 1e-12 * (1 + std::fabs(Kg.coeff(i, i)));
      Kg.makeCompressed();
      if (it == 1) {  // scale of the increment: its first tangent (not a diverging iterate)
        kscale = 0;
        for (int i = 0; i < dm.neq; ++i) kscale = std::max(kscale, std::fabs(Kg.coeff(i, i)));
      }
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
        if (dm.eq[i] >= 0) {
          A.du[i] += c[dm.eq[i]];
          if (!(A.any_pres && i % 6 == 3 && A.node_pres[i / 6])) cmax = std::max(cmax, std::fabs(c[dm.eq[i]]));
        }
      A.spread_links(A.du);
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
    double dWc = 0, dWf = 0;  // work of the contact forces: normal (penalty) and friction
    if (IC.active()) {
      contact_forces(A, IC, fcon, true, &ftan);
      for (int i = 0; i < 6 * nn; ++i) {
        A.fint[i] += fcon[i];
        dWf += 0.5 * (ftan_c[i] + ftan[i]) * A.du[i];
        dWc += 0.5 * (fcon_c[i] + fcon[i]) * A.du[i];
      }
      dWc -= dWf;
      fcon_c = fcon;
      ftan_c = ftan;
    }
    A.en.plastic += S.e_plastic;
    A.en.viscous += S.e_visc;
    double dWi = 0, dWe = 0;
    std::vector<double> fr(6 * nn);  // nodal reactions, the forces of the dependent dofs moved onto their terms
    for (int i = 0; i < 6 * nn; ++i) fr[i] = A.fint[i] - A.fext[i];
    // implicit dynamics: reactions of the HHT-weighted equilibrium
    if (dyn)
      for (int i = 0; i < 6 * nn; ++i) fr[i] = (1 + ha) * (A.fint[i] - A.fext[i]) - ha * (fint0[i] - fext0[i]);
    A.gather_links(fr);
    for (int i = 0; i < 6 * nn; ++i) {
      const double d = A.du[i];
      dWi += 0.5 * (fint_c[i] + A.fint[i]) * d;
      dWe += 0.5 * (fext_c[i] + A.fext[i]) * d;
      if (A.fixed[i]) dWe += 0.5 * (fint_c[i] - fext_c[i] + A.fint[i] - A.fext[i]) * d;  // work of the reactions
      A.rf[i] = A.fixed[i] ? fr[i] : (A.report_rf[i] ? A.fint[i] - A.fext[i] : 0.0);
      A.u[i] += A.du[i];
    }
    du_last = A.du;
    dt_last = dt;
    A.commit_rotations();
    if (dyn) {
      double ke = 0, dWd = 0;
      for (int i = 0; i < 6 * nn; ++i) {
        const double an = (A.du[i] - dt * A.v[i] - dt * dt * (0.5 - hb) * A.a[i]) / (hb * dt * dt);
        const double vn = A.v[i] + dt * ((1 - hg) * A.a[i] + hg * an);
        const double c = A.damp_coef(i, true);
        dWd += c * 0.5 * (vn + A.v[i]) * A.du[i];
        A.a[i] = an; A.v[i] = vn;
      }
      A.spread_links(A.v);
      for (int i = 0; i < 6 * nn; ++i) {
        const double mm = i % 6 < 3 ? A.mt[i / 6] : A.mr[i / 6] + (A.node_aniso[i / 6] ? A.jrot[6 * (size_t)(i / 6) + i % 6 - 3] : 0.0);
        ke += 0.5 * mm * A.v[i] * A.v[i];
      }
      A.en.kinetic = ke;
      A.en.damping += dWd;
      fint0 = A.fint; fext0 = A.fext;
    }
    A.en.internal += dWi - dWc - dWf - S.e_visc + S.e_corr;
    // contact energy: elastic energy of the penalty springs (a state function: an impact within one increment makes
    // the trapezoidal work of the contact forces meaningless) plus the friction dissipation, also reported as ALLFD
    A.en.friction += dWf;
    if (IC.active()) A.en.contact = IC.epot + A.en.friction;
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
    W.history(si, inc, t, A.total_time, dt, dyn ? A.en.kinetic : 0.0);
    bool last = t >= st.period * (1 - 1e-12);
    const int ffreq = st.output.field_frequency;
    if (ffreq > 0 ? (inc % ffreq == 0 || last) : (nfr <= 0 || t >= next_frame - 1e-9 * st.period || last)) {
      W.frame(si, st.name, t, A.total_time);
      if (nfr > 0) while (next_frame <= t + 1e-9 * st.period) next_frame += fr_dt;
    }
    // automatic incrementation
    if (cut_to_frame) dt = dt_auto;
    if (it - 1 <= 4) ++fast; else fast = 0;
    if (fast >= (dyn ? 1 : 2)) { dt = std::min(dt * 1.5, st.dtmax); fast = 0; }
    else if (it - 1 > 8) dt = std::max(dt * 0.75, st.dtmin);
    if (st.coupled && !st.heat_steady && st.deltmx > 0 && dT_inc > 0) dt = std::min(dt, dt_last * 0.85 * st.deltmx / dT_inc);
  }
  std::fill(A.du.begin(), A.du.end(), 0.0);
  W.sta_line("\n THE STEP %d HAS COMPLETED: %d increments, wallclock %.1f s", si + 1, inc, wall() - t0wall);
}

}  // namespace abq
