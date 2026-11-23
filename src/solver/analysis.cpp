#include "solver/analysis.hpp"

#include <stdexcept>

#include <omp.h>

#include <algorithm>
#include <functional>
#include <cmath>
#include <map>

#include "contact/contact.hpp"
#include "core/log.hpp"
#include "output/output.hpp"

namespace abq {

double LoadState::value(double v0, double v1, int amp, double st, double tt, const Model& m) const {
  if (amp == -2) return v1;
  if (amp >= 0) return v1 * m.amplitudes[amp].eval(st, tt);
  if (ramp) return v0 + (v1 - v0) * std::min(1.0, period > 0 ? st / period : 1.0);
  return v1;
}

Analysis::Analysis(Model& model) : m(model) {}

Analysis::~Analysis() = default;

void Analysis::setup() {
  nn = m.nn();
  tie_surfaces();  // may move slave nodes (ADJUST=YES) before the reference geometry is stored
  X.resize(3 * nn);
  for (int i = 0; i < nn; ++i) store3(&X[3 * i], m.X[i]);
  u.assign(6 * nn, 0); v.assign(6 * nn, 0); a.assign(6 * nn, 0); du.assign(6 * nn, 0);
  fint.assign(6 * nn, 0); fext.assign(6 * nn, 0); fcon.assign(6 * nn, 0); rf.assign(6 * nn, 0);
  mt.assign(nn, 0); mr.assign(nn, 0); cdamp.assign(nn, 0);
  temp.assign(nn, 0); temp_prev.assign(nn, 0); fv1.assign(nn, 0);
  fixed.assign(6 * nn, 0);
  node_rot.assign(nn, 0);
  node_aniso.assign(nn, 0);
  jrot.assign(6 * (size_t)nn, 0.0);
  groups = make_element_groups(m);
  joins.clear();
  for (const auto& el : m.elems)
    if (el.type == ElType::CONN3D2 && el.nodes.size() == 2 && el.section >= 0 &&
        m.sections[el.section].ctype == Section::JOIN)
      joins.push_back({el.nodes[0], el.nodes[1]});
  efor.resize(groups.size());
  for (size_t g = 0; g < groups.size(); ++g) {
    auto& G = *groups[g];
    efor[g].assign((size_t)G.size() * G.nen * 6, 0.0);
    if (G.rot && G.type != ElType::CONN3D2)
      for (int n : G.conn) node_rot[n] = 1;
    if (G.rot) {
      std::vector<double> xe(3 * G.nen);
      double J6[6];
      for (int e = 0; e < G.size(); ++e) {
        for (int I = 0; I < G.nen; ++I) for (int k = 0; k < 3; ++k) xe[3 * I + k] = X[3 * G.conn[e * G.nen + I] + k];
        if (G.rotary_tensor(e, xe.data(), J6))
          for (int I = 0; I < G.nen; ++I) node_aniso[G.conn[e * G.nen + I]] = 1;
      }
    }
    int cnt = G.size();
    info("  %-8s %8d elements", eltype_name(G.type), cnt);
  }
  // lumped nodal masses of every element (body forces), without mass scaling
  elem_mass_off.assign(m.ne(), -1);
  elem_mass.clear();
  for (auto& gp : groups) {
    auto& G = *gp;
    std::vector<double> Xe(3 * G.nen), me(G.nen), mre(G.nen);
    for (int e = 0; e < G.size(); ++e) {
      for (int I = 0; I < G.nen; ++I) for (int k = 0; k < 3; ++k) Xe[3 * I + k] = X[3 * G.conn[e * G.nen + I] + k];
      G.mass(e, Xe.data(), me.data(), mre.data());
      elem_mass_off[G.ids[e]] = (int)elem_mass.size();
      for (int I = 0; I < G.nen; ++I) elem_mass.push_back(me[I] / std::max(G.msf[e], 1e-300));
    }
  }
  elem_rigid.assign(m.ne(), 0);
  build_constraints();
  group_rigid.assign(groups.size(), {});
  for (size_t g = 0; g < groups.size(); ++g) {
    const auto& G = *groups[g];
    bool any = false;
    std::vector<char> r(G.size(), 0);
    for (int e = 0; e < G.size(); ++e) any |= (r[e] = elem_rigid[G.ids[e]]) != 0;
    if (any) group_rigid[g] = std::move(r);
  }
  for (auto& t : m.init_temperature) { temp[t.node] = t.value; has_temp = true; }
  temp_prev = temp;
  for (auto& iv : m.init_velocity) v[6 * iv.node + iv.dof] = iv.value;
  for (auto& s : m.steps) for (auto& l : s.loads) if (l.kind == LoadBlock::TEMPERATURE) has_temp = true;
  for (auto& mat : m.materials) if (mat.damp_beta != 0) warn("material %s: stiffness proportional damping (BETA) ignored", mat.name.c_str());
}

void Analysis::lump_masses() {
  std::fill(mt.begin(), mt.end(), 0.0);
  std::fill(mr.begin(), mr.end(), 0.0);
  std::fill(cdamp.begin(), cdamp.end(), 0.0);
  for (auto& gp : groups) {
    auto& G = *gp;
    std::vector<double> Xe(3 * G.nen), me(G.nen), mre(G.nen);
    for (int e = 0; e < G.size(); ++e) {
      for (int I = 0; I < G.nen; ++I) for (int k = 0; k < 3; ++k) Xe[3 * I + k] = X[3 * G.conn[e * G.nen + I] + k];
      G.mass(e, Xe.data(), me.data(), mre.data());
      double al = G.alpha_damping(e);
      for (int I = 0; I < G.nen; ++I) {
        int n = G.conn[e * G.nen + I];
        mt[n] += me[I] * G.msf[e];
        mr[n] += mre[I] * G.msf[e];
        cdamp[n] += al * me[I] * G.msf[e];
      }
    }
  }
  rotary_tensors(false);
}

void Analysis::rotary_tensors(bool with_increment) {
  any_aniso = false;
  for (int n = 0; n < nn; ++n)
    if (node_aniso[n]) { any_aniso = true; for (int k = 0; k < 6; ++k) jrot[6 * n + k] = 0.0; }
  if (!any_aniso) return;
  for (auto& gp : groups) {
    auto& G = *gp;
    if (!G.rot) continue;
    std::vector<double> xe(3 * G.nen);
    double J6[6];
    for (int e = 0; e < G.size(); ++e) {
      for (int I = 0; I < G.nen; ++I) {
        const int n = G.conn[e * G.nen + I];
        for (int k = 0; k < 3; ++k) xe[3 * I + k] = X[3 * n + k] + u[6 * n + k] + (with_increment ? du[6 * n + k] : 0.0);
      }
      if (!G.rotary_tensor(e, xe.data(), J6)) continue;
      for (int I = 0; I < G.nen; ++I) {
        const int n = G.conn[e * G.nen + I];
        for (int k = 0; k < 6; ++k) jrot[6 * n + k] += J6[k] * G.msf[e];
      }
    }
  }
}

namespace {
Mat3 node_inertia(const std::vector<double>& jrot, double mr, int n) {
  const double* j = &jrot[6 * (size_t)n];
  Mat3 J;
  J.m[0][0] = j[0] + mr; J.m[1][1] = j[1] + mr; J.m[2][2] = j[2] + mr;
  J.m[0][1] = J.m[1][0] = j[3]; J.m[0][2] = J.m[2][0] = j[4]; J.m[1][2] = J.m[2][1] = j[5];
  return J;
}
}  // namespace

void Analysis::rotary_accelerations(int n, const double* r, const double* w, double* acc, double* react) const {
  const Mat3 J = node_inertia(jrot, mr[n], n);
  const Vec3 om(w[0], w[1], w[2]);
  const Vec3 g = cross(om, J * om);  // gyroscopic moment
  double rhs[3] = {r[0] - g.x, r[1] - g.y, r[2] - g.z};
  // fixed dofs: zero acceleration, the reduced system is solved on the free ones
  Mat3 K = J;
  double b[3];
  for (int i = 0; i < 3; ++i) {
    const bool fx = fixed[6 * n + 3 + i];
    b[i] = fx ? 0.0 : rhs[i];
    if (fx)
      for (int j = 0; j < 3; ++j) { K.m[i][j] = K.m[j][i] = 0.0; }
  }
  for (int i = 0; i < 3; ++i) if (fixed[6 * n + 3 + i]) K.m[i][i] = 1.0;
  const double d = K.det();
  Vec3 x = std::fabs(d) > 0 ? K.inverse() * Vec3(b[0], b[1], b[2]) : Vec3();
  const double xa[3] = {x.x, x.y, x.z};
  for (int i = 0; i < 3; ++i) {
    const bool fx = fixed[6 * n + 3 + i];
    acc[i] = fx ? 0.0 : xa[i];
    double ja = 0;
    for (int j = 0; j < 3; ++j) ja += J.m[i][j] * (fixed[6 * n + 3 + j] ? 0.0 : xa[j]);
    react[i] = fx ? -(rhs[i] - ja) : 0.0;
  }
}

double Analysis::rotary_energy(int n, const double* w) const {
  const Mat3 J = node_inertia(jrot, mr[n], n);
  const Vec3 om(w[0], w[1], w[2]);
  return 0.5 * dot(om, J * om);
}

void Analysis::temperatures(double st, double tt) {
  if (!has_temp && loads.fields.empty()) return;
  for (auto& p : loads.temps) temp[p.node] = loads.value(p, st, tt, m);
  for (auto& p : loads.fields) fv1[p.node] = loads.value(p, st, tt, m);
}

void Analysis::begin_step(const Step& st, int si, bool ramp_default) {
  // freeze the values reached at the end of the previous step
  double T = loads.period;
  auto freeze = [&](std::vector<Prescribed>& v) {
    for (auto& p : v) { double x = loads.value(p, T, total_time, m); p.v0 = p.v1 = x; p.amp = -2; }
  };
  if (si > 0) {
    // prescribed displacements continue from the position actually reached (in explicit steps, a displacement
    // condition with an amplitude prescribes the velocity, so the node may not sit at magnitude x amplitude)
    for (auto& p : loads.bcs) {
      double x = p.velocity ? loads.value(p, T, total_time, m) : u[6 * p.node + p.dof];
      p.v0 = p.v1 = x; p.amp = -2;
    }
    freeze(loads.cloads); freeze(loads.temps); freeze(loads.fields);
    for (auto& f : loads.press) { double x = loads.value(f.v0, f.v1, f.amp, T, total_time, m); f.v0 = f.v1 = x; f.amp = -2; }
  }
  loads.ramp = ramp_default;
  loads.period = st.period;
  auto apply = [&](std::vector<Prescribed>& vec, const LoadBlock& lb, bool velocity, std::function<double(int, int)> cur) {
    std::map<long, int> idx;
    for (size_t i = 0; i < vec.size(); ++i) idx[(long)vec[i].node * 8 + vec[i].dof] = (int)i;
    for (auto& nv : lb.nodal) {
      long key = (long)nv.node * 8 + nv.dof;
      auto it = idx.find(key);
      if (it != idx.end()) {
        Prescribed& p = vec[it->second];  // v0: value at the end of the previous step (frozen)
        p.v1 = nv.value;
        p.amp = lb.amplitude;
        p.velocity = velocity;
      } else {
        Prescribed p{nv.node, nv.dof, cur(nv.node, nv.dof), nv.value, lb.amplitude, velocity};
        idx[key] = (int)vec.size();
        vec.push_back(p);
      }
    }
  };
  // OP=NEW removes the loads of the category first. Removed concentrated and distributed loads go to zero with the
  // default amplitude of the step: ramped down over a static step, removed at once in explicit dynamics.
  for (auto& lb : st.loads) {
    if (!lb.op_new) continue;
    const int removal = -1;
    if (lb.kind == LoadBlock::BOUNDARY || lb.kind == LoadBlock::VELOCITY_BC || lb.kind == LoadBlock::ACCELERATION_BC) loads.bcs.clear();
    else if (lb.kind == LoadBlock::CLOAD) for (auto& p : loads.cloads) { p.v1 = 0; p.amp = removal; }
    else if (lb.kind == LoadBlock::DSLOAD) for (auto& f : loads.press) { f.v1 = 0; f.amp = removal; }
    else if (lb.kind == LoadBlock::TEMPERATURE) loads.temps.clear();
    else if (lb.kind == LoadBlock::FIELD) loads.fields.clear();
  }
  for (auto& lb : st.loads) {
    switch (lb.kind) {
      case LoadBlock::BOUNDARY:
      case LoadBlock::VELOCITY_BC:
      case LoadBlock::ACCELERATION_BC:
        apply(loads.bcs, lb, lb.kind == LoadBlock::VELOCITY_BC, [&](int n, int d) { return u[6 * n + d]; });
        if (lb.kind == LoadBlock::ACCELERATION_BC)
          for (auto& nv : lb.nodal)
            for (auto& p : loads.bcs) if (p.node == nv.node && p.dof == nv.dof) { p.acceleration = true; p.velocity = false; }
        break;
      case LoadBlock::CLOAD:
        apply(loads.cloads, lb, false, [&](int, int) { return 0.0; });
        if (lb.follower)
          for (auto& nv : lb.nodal) {
            for (auto& p : loads.cloads) if (p.node == nv.node && p.dof == nv.dof) p.follower = true;
            if (rot_track.empty()) rot_track.assign(nn, -1);
            if (rot_track[nv.node] < 0) { rot_track[nv.node] = (int)node_R.size(); node_R.push_back(Mat3::identity()); }
          }
        break;
      case LoadBlock::TEMPERATURE: apply(loads.temps, lb, false, [&](int n, int) { return temp[n]; }); break;
      case LoadBlock::FIELD: apply(loads.fields, lb, false, [&](int n, int) { return fv1[n]; }); break;
      case LoadBlock::DSLOAD: {
        std::map<long, int> idx;
        for (size_t i = 0; i < loads.press.size(); ++i) idx[(long)loads.press[i].face.elem * 4096 + loads.press[i].face.face + 2048] = (int)i;
        for (auto& fv : lb.faces) {
          long key = (long)fv.face.elem * 4096 + fv.face.face + 2048;
          auto it = idx.find(key);
          if (it != idx.end()) {
            FaceLoad& f = loads.press[it->second];  // v0 frozen at the end of the previous step
            f.v1 = fv.value; f.amp = lb.amplitude;
          } else {
            idx[key] = (int)loads.press.size();
            loads.press.push_back({fv.face, 0.0, fv.value, lb.amplitude});
          }
        }
        break;
      }
    }
  }
  update_fixed(st);
  activate_links(st.proc == Step::STATIC);
}

void Analysis::update_fixed(const Step&) {
  std::fill(fixed.begin(), fixed.end(), 0);
  local_bcs.clear();
  report_rf.assign(6 * nn, 0);
  auto set = [&](int node, int dof, double value) {
    if (!m.transforms.count(node)) { fixed[6 * node + dof] = 1; return; }
    local_bcs.push_back({node, dof, value});
    for (int k = 0; k < 3; ++k) report_rf[6 * node + (dof < 3 ? 0 : 3) + k] = 1;
  };
  for (auto& b : m.model_bcs) set(b.node, b.dof, b.value);
  for (auto& b : loads.bcs) set(b.node, b.dof, b.v1);
}

Analysis::ElemSum Analysis::eval_elements(double dt, double bv1, double bv2, bool dynamic, bool commit, std::vector<double>& f) {
  ElemSum S;
  std::fill(f.begin(), f.end(), 0.0);
  if (spring_k.size() != groups.size()) {
    spring_k.resize(groups.size());
    for (size_t g = 0; g < groups.size(); ++g) spring_k[g].assign(groups[g]->size(), 0.0);
  }
  for (size_t g = 0; g < groups.size(); ++g) {
    if (!group_mask.empty() && !group_mask[g]) continue;
    ElementGroup& G = *groups[g];
    const int nen = G.nen;
    double* EF = efor[g].data();
    double dtmin = 1e30, ev = 0, ep = 0, ea = 0, ec = 0;
    int dte = -1, nfail = 0;
#pragma omp parallel num_threads(nthreads)
    {
      double Xe[3 * 8], xo[3 * 8], due[6 * 8], t0[8], t1[8], fve[8];
      double my_dt = 1e30;
      int my_e = -1;
#pragma omp for reduction(+ : ev, ep, ea, ec, nfail) schedule(static)
      for (int e = 0; e < G.size(); ++e) {
        const int* cn = &G.conn[e * nen];
        if (!group_rigid[g].empty() && group_rigid[g][e]) {  // element of a rigid body: mass only
          std::fill(EF + (size_t)e * nen * 6, EF + (size_t)(e + 1) * nen * 6, 0.0);
          continue;
        }
        for (int I = 0; I < nen; ++I) {
          int n = cn[I];
          for (int k = 0; k < 3; ++k) { Xe[3 * I + k] = X[3 * n + k]; xo[3 * I + k] = X[3 * n + k] + u[6 * n + k]; }
          for (int k = 0; k < 6; ++k) due[6 * I + k] = du[6 * n + k];
          t0[I] = temp_prev[n]; t1[I] = temp[n]; fve[I] = fv1[n];
        }
        ElemIn in;
        in.X = Xe; in.xo = xo; in.du = due;
        if (has_temp) { in.temp0 = t0; in.temp1 = t1; }
        in.fv1 = fve;
        in.dt = dt; in.bv1 = bv1; in.bv2 = bv2; in.dynamic = dynamic; in.nlgeom = nlgeom; in.initial = initial_pass;
        ElemOut out;
        out.f = EF + (size_t)e * nen * 6;
        G.compute(e, in, out, commit);
        ev += out.e_visc; ep += out.e_plastic; ea += out.e_art; ec += out.e_corr;
        if (out.failed) ++nfail;
        double dte_e = out.dt_stable;
        if (dte_e < 0) {  // spring-like element: stiffness handled per node below
          if (commit) G.dt_last[e] = dte_e;
          spring_k[g][e] = -dte_e;
          dte_e = 1e30;
        } else {
          if (commit) G.dt_last[e] = dte_e / std::sqrt(G.msf[e]);
          spring_k[g][e] = 0;
          if (dte_e < my_dt) { my_dt = dte_e; my_e = e; }
        }
      }
#pragma omp critical
      {
        if (my_dt < dtmin) { dtmin = my_dt; dte = my_e; }
      }
    }
    // scatter (serial, deterministic)
    for (int e = 0; e < G.size(); ++e) {
      const int* cn = &G.conn[e * nen];
      const double* fe = EF + (size_t)e * nen * 6;
      for (int I = 0; I < nen; ++I) {
        double* fn = &f[6 * cn[I]];
        for (int k = 0; k < 6; ++k) fn[k] += fe[6 * I + k];
      }
    }
    S.e_visc += ev; S.e_plastic += ep; S.e_art += ea; S.e_corr += ec; S.failed += nfail;
    if (dtmin < S.dt_min) { S.dt_min = dtmin; S.dt_elem = dte >= 0 ? G.ids[dte] : -1; }
    double& dcl = is_solid_type(G.type) ? S.dt_solid : S.dt_other;
    dcl = std::min(dcl, dtmin);
  }
  // springs / connectors: their stiffness adds to the frequency of the nodes they load
  bool any_spring = false;
  for (auto& v : spring_k) for (double k : v) if (k > 0) { any_spring = true; break; }
  if (any_spring) {
    node_om2.assign(nn, 0.0);
    node_ks.assign(nn, 0.0);
    std::vector<int> node_el(nn, -1);
    for (size_t g = 0; g < groups.size(); ++g) {
      ElementGroup& G = *groups[g];
      for (int e = 0; e < G.size(); ++e) {
        if (spring_k[g][e] > 0) {
          for (int I = 0; I < G.nen; ++I) {
            int n = G.conn[e * G.nen + I];
            node_ks[n] += spring_k[g][e];
            if (node_el[n] < 0) node_el[n] = G.ids[e];
          }
        } else if (G.dt_last[e] < 1e29) {
          double d = G.dt_last[e] * std::sqrt(G.msf[e]);
          double om2 = 4.0 / (d * d);
          for (int I = 0; I < G.nen; ++I) {
            int n = G.conn[e * G.nen + I];
            node_om2[n] = std::max(node_om2[n], om2);
          }
        }
      }
    }
    for (int n = 0; n < nn; ++n) {
      if (node_ks[n] <= 0 || mt[n] <= 0) continue;
      if (fixed[6 * n] && fixed[6 * n + 1] && fixed[6 * n + 2]) continue;
      double om2 = node_om2[n] + node_ks[n] / mt[n];
      double d = 2.0 / std::sqrt(om2);
      if (d < S.dt_min) { S.dt_min = d; S.dt_elem = node_el[n]; }
      S.dt_other = std::min(S.dt_other, d);
    }
  }
  return S;
}

namespace {
// outward-ordered face node lists (local indices)
const int HEXF[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
const int TETF[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
}  // namespace

void Analysis::commit_rotations() {
  for (int n = 0; n < (int)rot_track.size(); ++n)
    if (rot_track[n] >= 0) node_R[rot_track[n]] = rotation_from_vector(load3(&du[6 * n + 3])) * node_R[rot_track[n]];
}

void Analysis::external_forces(double st, double tt, std::vector<double>& f) {
  std::fill(f.begin(), f.end(), 0.0);
  for (auto& c : loads.cloads) {
    const double v = loads.value(c, st, tt, m);
    auto tr = m.transforms.find(c.node);
    if (tr != m.transforms.end() && !c.follower) {  // load given in the local axes of the node
      const int base = c.dof < 3 ? 0 : 3;
      add3(&f[6 * c.node + base], tr->second.row(c.dof - base) * v);
      continue;
    }
    if (!c.follower || rot_track.empty() || rot_track[c.node] < 0) { f[6 * c.node + c.dof] += v; continue; }
    // follower: the reference direction rotated with the node (rotation at the start of the increment, times the
    // rotation of the increment)
    const Mat3 R = rotation_from_vector(load3(&du[6 * c.node + 3])) * node_R[rot_track[c.node]];
    const int base = c.dof < 3 ? 0 : 3;
    add3(&f[6 * c.node + base], R.col(c.dof - base) * v);
  }
  if (loads.press.empty()) return;
  auto xcur = [&](int n) { return Vec3(X[3 * n] + u[6 * n] + du[6 * n], X[3 * n + 1] + u[6 * n + 1] + du[6 * n + 1], X[3 * n + 2] + u[6 * n + 2] + du[6 * n + 2]); };
  for (auto& pl : loads.press) {
    double p = loads.value(pl.v0, pl.v1, pl.amp, st, tt, m);
    if (p == 0) continue;
    const Element& el = m.elems[pl.face.elem];
    if (pl.face.face >= 1001 && pl.face.face <= 1003) {  // gravity: lumped nodal masses x acceleration component
      const int k = pl.face.face - 1001, off = elem_mass_off.empty() ? -1 : elem_mass_off[pl.face.elem];
      if (off < 0) continue;
      for (size_t I = 0; I < el.nodes.size(); ++I) f[6 * el.nodes[I] + k] += elem_mass[off + I] * p;
      continue;
    }
    int nodes[4];
    int nf = 0;
    double sgn = 1;
    switch (el.type) {
      case ElType::S4R: case ElType::S4: case ElType::S3R: case ElType::S3RS:
        nf = (int)el.nodes.size();
        for (int i = 0; i < nf; ++i) nodes[i] = el.nodes[i];
        sgn = pl.face.face < 0 ? -1.0 : 1.0;  // SPOS face: outward normal = +n
        break;
      case ElType::C3D8R: case ElType::C3D8: case ElType::C3D8I: {
        int fi = std::max(1, pl.face.face) - 1;
        nf = 4;
        for (int i = 0; i < 4; ++i) nodes[i] = el.nodes[HEXF[fi][i]];
        break;
      }
      case ElType::C3D4: {
        int fi = std::max(1, pl.face.face) - 1;
        nf = 3;
        for (int i = 0; i < 3; ++i) nodes[i] = el.nodes[TETF[fi][i]];
        break;
      }
      default: continue;
    }
    Vec3 An;  // area-weighted outward normal
    if (nf == 4) An = cross(xcur(nodes[2]) - xcur(nodes[0]), xcur(nodes[3]) - xcur(nodes[1])) * 0.5;
    else An = cross(xcur(nodes[1]) - xcur(nodes[0]), xcur(nodes[2]) - xcur(nodes[0])) * 0.5;
    Vec3 F = An * (-p * sgn / nf);  // positive pressure pushes into the face
    for (int i = 0; i < nf; ++i) add3(&f[6 * nodes[i]], F);
  }
}

void Analysis::save_state(BinOut& o) const {
  o.tag("analysis");
  o.pod(nn);
  o.vec(u); o.vec(v); o.vec(a); o.vec(du); o.vec(fint); o.vec(fext); o.vec(fcon); o.vec(rf);
  o.vec(mt); o.vec(mr); o.vec(cdamp); o.vec(temp); o.vec(temp_prev); o.vec(fv1); o.vec(fixed);
  o.vec(loads.bcs); o.vec(loads.cloads); o.vec(loads.temps); o.vec(loads.fields); o.vec(loads.press);
  o.pod(loads.ramp); o.pod(loads.period);
  o.pod(en); o.pod(total_time); o.pod(total_increments); o.pod(nlgeom);
  o.tag("groups");
  uint64_t ng = groups.size();
  o.pod(ng);
  for (auto& g : groups) { o.pod(g->type); o.pod(g->size()); g->save(o); }
}

void Analysis::load_state(BinIn& in) {
  in.tag("analysis");
  int n;
  in.pod(n);
  if (n != nn) throw std::runtime_error("restart file: node count differs from the model");
  in.vec(u); in.vec(v); in.vec(a); in.vec(du); in.vec(fint); in.vec(fext); in.vec(fcon); in.vec(rf);
  in.vec(mt); in.vec(mr); in.vec(cdamp); in.vec(temp); in.vec(temp_prev); in.vec(fv1); in.vec(fixed);
  in.vec(loads.bcs); in.vec(loads.cloads); in.vec(loads.temps); in.vec(loads.fields); in.vec(loads.press);
  in.pod(loads.ramp); in.pod(loads.period);
  in.pod(en); in.pod(total_time); in.pod(total_increments); in.pod(nlgeom);
  in.tag("groups");
  uint64_t ng;
  in.pod(ng);
  if (ng != groups.size()) throw std::runtime_error("restart file: element groups differ from the model");
  for (auto& g : groups) {
    ElType t; int sz;
    in.pod(t); in.pod(sz);
    if (t != g->type || sz != g->size()) throw std::runtime_error("restart file: element group mismatch");
    g->load(in);
  }
}

}  // namespace abq
