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

// Tied translations move together: common acceleration (r_a + r_b) / (m_a + m_b); a prescribed dof of one
// node holds the other one, whose force then adds to the reaction.
void Analysis::tie_accelerations() {
  for (auto [na, nb] : ties)
    for (int k = 0; k < 3; ++k) {
      int ia = 6 * na + k, ib = 6 * nb + k;
      double ra = fixed[ia] ? -rf[ia] : mt[na] * a[ia], rb = fixed[ib] ? -rf[ib] : mt[nb] * a[ib];
      if (fixed[ia] || fixed[ib]) {
        if (fixed[ia] && !fixed[ib]) { rf[ia] -= rb; a[ib] = a[ia]; }
        else if (fixed[ib] && !fixed[ia]) { rf[ib] -= ra; a[ia] = a[ib]; }
        continue;
      }
      double mm = mt[na] + mt[nb];
      double ac = mm > 0 ? (ra + rb) / mm : 0.0;
      a[ia] = a[ib] = ac;
    }
}
Analysis::~Analysis() = default;

void Analysis::setup() {
  nn = m.nn();
  X.resize(3 * nn);
  for (int i = 0; i < nn; ++i) store3(&X[3 * i], m.X[i]);
  u.assign(6 * nn, 0); v.assign(6 * nn, 0); a.assign(6 * nn, 0); du.assign(6 * nn, 0);
  fint.assign(6 * nn, 0); fext.assign(6 * nn, 0); fcon.assign(6 * nn, 0); rf.assign(6 * nn, 0);
  mt.assign(nn, 0); mr.assign(nn, 0); cdamp.assign(nn, 0);
  temp.assign(nn, 0); temp_prev.assign(nn, 0); fv1.assign(nn, 0);
  fixed.assign(6 * nn, 0);
  node_rot.assign(nn, 0);
  groups = make_element_groups(m);
  ties.clear();
  for (const auto& el : m.elems)
    if (el.type == ElType::CONN3D2 && el.nodes.size() == 2 && el.section >= 0 &&
        m.sections[el.section].ctype == Section::JOIN)
      ties.push_back({el.nodes[0], el.nodes[1]});
  efor.resize(groups.size());
  for (size_t g = 0; g < groups.size(); ++g) {
    auto& G = *groups[g];
    efor[g].assign((size_t)G.size() * G.nen * 6, 0.0);
    if (G.rot && G.type != ElType::CONN3D2)
      for (int n : G.conn) node_rot[n] = 1;
    int cnt = G.size();
    info("  %-8s %8d elements", eltype_name(G.type), cnt);
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
    if (lb.kind == LoadBlock::BOUNDARY || lb.kind == LoadBlock::VELOCITY_BC) loads.bcs.clear();
    else if (lb.kind == LoadBlock::CLOAD) for (auto& p : loads.cloads) { p.v1 = 0; p.amp = removal; }
    else if (lb.kind == LoadBlock::DSLOAD) for (auto& f : loads.press) { f.v1 = 0; f.amp = removal; }
    else if (lb.kind == LoadBlock::TEMPERATURE) loads.temps.clear();
    else if (lb.kind == LoadBlock::FIELD) loads.fields.clear();
  }
  for (auto& lb : st.loads) {
    switch (lb.kind) {
      case LoadBlock::BOUNDARY:
      case LoadBlock::VELOCITY_BC:
        apply(loads.bcs, lb, lb.kind == LoadBlock::VELOCITY_BC, [&](int n, int d) { return u[6 * n + d]; });
        break;
      case LoadBlock::CLOAD: apply(loads.cloads, lb, false, [&](int, int) { return 0.0; }); break;
      case LoadBlock::TEMPERATURE: apply(loads.temps, lb, false, [&](int n, int) { return temp[n]; }); break;
      case LoadBlock::FIELD: apply(loads.fields, lb, false, [&](int n, int) { return fv1[n]; }); break;
      case LoadBlock::DSLOAD: {
        std::map<long, int> idx;
        for (size_t i = 0; i < loads.press.size(); ++i) idx[(long)loads.press[i].face.elem * 16 + loads.press[i].face.face + 8] = (int)i;
        for (auto& fv : lb.faces) {
          long key = (long)fv.face.elem * 16 + fv.face.face + 8;
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
}

void Analysis::update_fixed(const Step&) {
  std::fill(fixed.begin(), fixed.end(), 0);
  for (auto& b : m.model_bcs) fixed[6 * b.node + b.dof] = 1;
  for (auto& b : loads.bcs) fixed[6 * b.node + b.dof] = 1;
}

Analysis::ElemSum Analysis::eval_elements(double dt, double bv1, double bv2, bool dynamic, bool commit, std::vector<double>& f) {
  ElemSum S;
  std::fill(f.begin(), f.end(), 0.0);
  if (spring_k.size() != groups.size()) {
    spring_k.resize(groups.size());
    for (size_t g = 0; g < groups.size(); ++g) spring_k[g].assign(groups[g]->size(), 0.0);
  }
  for (size_t g = 0; g < groups.size(); ++g) {
    ElementGroup& G = *groups[g];
    const int nen = G.nen;
    double* EF = efor[g].data();
    double dtmin = 1e30, ev = 0, ep = 0, ea = 0;
    int dte = -1, nfail = 0;
#pragma omp parallel num_threads(nthreads)
    {
      double Xe[3 * 8], xo[3 * 8], due[6 * 8], t0[8], t1[8], fve[8];
      double my_dt = 1e30;
      int my_e = -1;
#pragma omp for reduction(+ : ev, ep, ea, nfail) schedule(static)
      for (int e = 0; e < G.size(); ++e) {
        const int* cn = &G.conn[e * nen];
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
        ev += out.e_visc; ep += out.e_plastic; ea += out.e_art;
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
    S.e_visc += ev; S.e_plastic += ep; S.e_art += ea; S.failed += nfail;
    if (dtmin < S.dt_min) { S.dt_min = dtmin; S.dt_elem = dte >= 0 ? G.ids[dte] : -1; }
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
    }
  }
  return S;
}

namespace {
// outward-ordered face node lists (local indices)
const int HEXF[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
const int TETF[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
}  // namespace

void Analysis::external_forces(double st, double tt, std::vector<double>& f) {
  std::fill(f.begin(), f.end(), 0.0);
  for (auto& c : loads.cloads) f[6 * c.node + c.dof] += loads.value(c, st, tt, m);
  if (loads.press.empty()) return;
  auto xcur = [&](int n) { return Vec3(X[3 * n] + u[6 * n] + du[6 * n], X[3 * n + 1] + u[6 * n + 1] + du[6 * n + 1], X[3 * n + 2] + u[6 * n + 2] + du[6 * n + 2]); };
  for (auto& pl : loads.press) {
    double p = loads.value(pl.v0, pl.v1, pl.amp, st, tt, m);
    if (p == 0) continue;
    const Element& el = m.elems[pl.face.elem];
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
