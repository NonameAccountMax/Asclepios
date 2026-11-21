// General contact for explicit dynamics.
//  - contact entities: shell facets (triangles, two-sided, offset by a reduced half thickness),
//    beam segments (radius), slave nodes (all nodes of the included elements); beyond a facet
//    boundary the contact is with the rounded edge (closest point), so no gap between facets;
//  - free shell edges are rounded within the shell: the nodes of a free edge are moved back into the
//    shell plane by the half thickness (perpendicular to the edge, along the bisector at a corner;
//    at most half the facet size and 0.6 x facet size x sine of the angle to each edge), so that
//    the contact surface does not extend beyond the edge;
//  - node-to-facet and node-to-segment penalty from the stable increment of the elements of each node: contact
//    frequency W^2 = max(1/(3 dt^2), 1.1 max_i 1/dt_i^2) over the pair, nodal stiffnesses that keep each node
//    stable at 1/W, pair stiffness = softer of slave and master, 9/11 of its smallest value since the contact
//    started, shared between the contacts of a node; the increment is bounded by 1/W of the active contacts
//    (0.3 m/dt^2 between two critical nodes, stiffer for light nodes of soft elements);
//  - facet half thickness h = min(t/2, 0.3 x smallest facet height), segment radius
//    r = min(r, 0.3 x segment length) (thickness reduced on small facets);
//  - normal viscous damping of 3 % of the critical damping of the contact spring with the lighter side;
//  - Coulomb friction, tangential penalty stiffness equal to the normal one;
//  - incremental normal force (no jump when k changes, when a node switches facet or when a
//    contact is found deep); initial overclosures act as offsets released as the nodes separate;
//  - candidate pairs from a bucket grid, rebuilt when nodes moved more than half the margin.
#include "contact/contact.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <unordered_map>

#include "core/log.hpp"
#include "solver/analysis.hpp"

namespace abq {

namespace {

struct Hit {
  int master = INT_MIN;  // >=0 triangle, <0: -1-segment, INT_MIN: none
  double pen = 0;
  double k = 0;          // penalty stiffness used
  Vec3 n;                // direction pushing the slave out
  double w[3] = {0, 0, 0};
  int side = 0;
};

// closest point on a triangle (Ericson), barycentric weights; interior = projection inside
Vec3 closest_on_tri(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c, double w[3], bool& interior) {
  Vec3 ab = b - a, ac = c - a, ap = p - a;
  double d1 = dot(ab, ap), d2 = dot(ac, ap);
  interior = false;
  if (d1 <= 0 && d2 <= 0) { w[0] = 1; w[1] = w[2] = 0; return a; }
  Vec3 bp = p - b;
  double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) { w[1] = 1; w[0] = w[2] = 0; return b; }
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) { double v = d1 / (d1 - d3); w[0] = 1 - v; w[1] = v; w[2] = 0; return a + ab * v; }
  Vec3 cp = p - c;
  double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) { w[2] = 1; w[0] = w[1] = 0; return c; }
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) { double v = d2 / (d2 - d6); w[0] = 1 - v; w[2] = v; w[1] = 0; return a + ac * v; }
  double va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    w[1] = 1 - v; w[2] = v; w[0] = 0;
    return b + (c - b) * v;
  }
  double den = 1.0 / (va + vb + vc);
  double v = vb * den, ww = vc * den;
  w[0] = 1 - v - ww; w[1] = v; w[2] = ww;
  interior = true;
  return a + ab * v + ac * ww;
}

}  // namespace

struct GeneralContact::Impl {
  std::vector<int> tri;       // 3 per triangle
  std::vector<double> tri_h;
  std::vector<int> seg;       // 2 per segment
  std::vector<double> seg_r;
  std::vector<int> slaves;    // node ids
  std::vector<double> slave_h;
  // free shell edges: element nodes of each shell element touching a free edge (4, -1 padded) and its size,
  // free edges (two nodes and the element), nodes on free edges with their edges, current offset of every node
  std::vector<int> fe_elem;
  std::vector<double> fe_size;
  std::vector<int> fedge;
  std::vector<int> fnode, fn_off, fn_edge;
  std::vector<Vec3> shift;
  std::vector<int> slave_of;  // node -> slave index or -1
  std::vector<int> nb_off, nb;  // per slave: sorted nodes of its elements (facets touching them are excluded)
  std::vector<int> cand_off, cand;
  std::vector<double> xref;   // node positions at the last candidate search
  double margin = 0, cell = 0;
  std::vector<int> last_master;
  std::vector<signed char> last_side;
  std::vector<Vec3> ft;
  bool first_eval = true;
  // incremental normal force: force and penetration of the previous increment (pen < 0: no contact)
  std::vector<double> fn_t, pen_t, fn_s, pen_s;
  std::vector<Vec3> last_n;  // direction pushing the slave out at the previous increment
  // penalty: frequency bound of every node (1/dt_e^2 of its stiffest element), per slave and pass the contact
  // frequency, pair stiffness, damping mass and smallest pair stiffness since the contact started
  std::vector<double> node_w2;
  std::vector<double> cw2_t, cw2_s, md_t, md_s, kmin_t, kmin_s;
  double w2_max = 0;  // largest contact frequency of the increment (stable increment bound)
  double mu = 0;
  std::vector<Hit> hit_t, hit_s;
  std::vector<Vec3> force_t, force_s;
  int rebuilds = 0;
  // diagnostics (ABQ_CONTACT_DIAG)
  int diag = -1;
  long d_deep_n = 0;
};

GeneralContact::GeneralContact(Analysis& A) : A_(A) {}
GeneralContact::~GeneralContact() { delete p_; }

bool GeneralContact::define(const Step& st) {
  const Model& m = A_.m;
  delete p_;
  p_ = new Impl();
  Impl& P = *p_;
  std::vector<char> incl(m.ne(), 0);
  if (st.contact_all_exterior) std::fill(incl.begin(), incl.end(), 1);
  for (auto& sn : st.contact_surfaces) {
    const Surface& s = m.surfaces.at(sn);
    for (auto& f : s.faces) incl[f.elem] = 1;
  }
  int nsolid = 0;
  std::vector<double> node_h(m.nn(), 1e30);
  std::vector<std::array<int, 4>> shells;  // included shell elements (nodes, -1 padded)
  for (auto& gp : A_.groups) {
    ElementGroup& G = *gp;
    for (int e = 0; e < G.size(); ++e) {
      int ge = G.ids[e];
      if (!incl[ge]) continue;
      const int* cn = &G.conn[e * G.nen];
      if (G.type == ElType::S4R || G.type == ElType::S4 || G.type == ElType::S3R || G.type == ElType::S3RS) {
        Vec3 p[4];
        for (int I = 0; I < G.nen; ++I) p[I] = m.X[cn[I]];
        // smallest height of the triangular facets (a quadrilateral is split along its 0-2 diagonal)
        auto tri_height = [](const Vec3& a, const Vec3& b, const Vec3& c) {
          double l = std::max(norm(b - a), std::max(norm(c - b), norm(a - c)));
          return l > 0 ? norm(cross(b - a, c - a)) / l : 0.0;
        };
        double hmin = G.nen == 4 ? std::min(tri_height(p[0], p[1], p[2]), tri_height(p[0], p[2], p[3]))
                                 : tri_height(p[0], p[1], p[2]);
        double h = std::min(0.5 * G.contact_thickness(e), 0.3 * hmin);
        auto add = [&](int a, int b, int c) { P.tri.insert(P.tri.end(), {a, b, c}); P.tri_h.push_back(h); };
        if (G.nen == 4) { add(cn[0], cn[1], cn[2]); add(cn[0], cn[2], cn[3]); }
        else add(cn[0], cn[1], cn[2]);
        for (int I = 0; I < G.nen; ++I) node_h[cn[I]] = std::min(node_h[cn[I]], h);
        shells.push_back({cn[0], cn[1], cn[2], G.nen == 4 ? cn[3] : -1});
      } else if (G.type == ElType::B31) {
        double L = norm(m.X[cn[1]] - m.X[cn[0]]);
        double r = std::min(0.5 * G.contact_thickness(e), 0.3 * L);
        P.seg.insert(P.seg.end(), {cn[0], cn[1]});
        P.seg_r.push_back(r);
        for (int I = 0; I < 2; ++I) node_h[cn[I]] = std::min(node_h[cn[I]], r);
      } else if (G.type != ElType::CONN3D2) {
        ++nsolid;
      }
    }
  }
  {  // free shell edges: element edges used by a single included shell element
    std::vector<std::pair<uint64_t, int>> edges;  // (sorted node pair, element)
    for (int k = 0; k < (int)shells.size(); ++k) {
      const auto& q = shells[k];
      int nv = q[3] >= 0 ? 4 : 3;
      for (int i = 0; i < nv; ++i) {
        uint64_t a = (uint64_t)q[i], b = (uint64_t)q[(i + 1) % nv];
        edges.push_back({std::min(a, b) << 32 | std::max(a, b), k});
      }
    }
    std::sort(edges.begin(), edges.end());
    std::vector<int> elem_id(shells.size(), -1);
    std::vector<std::vector<int>> node_edges(m.nn());
    for (size_t i = 0; i < edges.size(); ++i) {
      bool dup = (i > 0 && edges[i - 1].first == edges[i].first) ||
                 (i + 1 < edges.size() && edges[i + 1].first == edges[i].first);
      if (dup) continue;
      int k = edges[i].second;
      if (elem_id[k] < 0) {
        elem_id[k] = (int)P.fe_size.size();
        const auto& q = shells[k];
        P.fe_elem.insert(P.fe_elem.end(), q.begin(), q.end());
        Vec3 a = m.X[q[0]], b = m.X[q[1]], c = m.X[q[2]];
        double area = q[3] >= 0 ? 0.5 * norm(cross(c - a, m.X[q[3]] - b)) : 0.5 * norm(cross(b - a, c - a));
        P.fe_size.push_back(std::sqrt(area));
      }
      int a = (int)(edges[i].first >> 32), b = (int)(uint32_t)edges[i].first;
      int id = (int)P.fedge.size() / 3;
      P.fedge.insert(P.fedge.end(), {a, b, elem_id[k]});
      node_edges[a].push_back(id);
      node_edges[b].push_back(id);
    }
    P.fn_off.push_back(0);
    for (int n = 0; n < m.nn(); ++n) {
      if (node_edges[n].empty()) continue;
      P.fnode.push_back(n);
      P.fn_edge.insert(P.fn_edge.end(), node_edges[n].begin(), node_edges[n].end());
      P.fn_off.push_back((int)P.fn_edge.size());
    }
    P.shift.assign(m.nn(), Vec3());
  }
  if (nsolid) warn("general contact: %d solid elements included; solid faces are not supported yet and are ignored", nsolid);
  P.slave_of.assign(m.nn(), -1);
  for (int n = 0; n < m.nn(); ++n)
    if (node_h[n] < 1e29) { P.slave_of[n] = (int)P.slaves.size(); P.slaves.push_back(n); P.slave_h.push_back(node_h[n]); }
  int ns = (int)P.slaves.size();
  {  // topological neighbourhood: nodes sharing a facet or segment with the slave
    std::vector<std::vector<int>> nbs(ns);
    auto link = [&](const int* v, int k) {
      for (int i = 0; i < k; ++i) {
        int si = P.slave_of[v[i]];
        if (si < 0) continue;
        for (int j = 0; j < k; ++j) nbs[si].push_back(v[j]);
      }
    };
    for (size_t t = 0; t < P.tri_h.size(); ++t) link(&P.tri[3 * t], 3);
    for (size_t sg = 0; sg < P.seg_r.size(); ++sg) link(&P.seg[2 * sg], 2);
    P.nb_off.assign(ns + 1, 0);
    for (int si = 0; si < ns; ++si) {
      auto& v = nbs[si];
      v.push_back(P.slaves[si]);
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      P.nb_off[si + 1] = P.nb_off[si] + (int)v.size();
    }
    P.nb.resize(P.nb_off[ns]);
    for (int si = 0; si < ns; ++si) std::copy(nbs[si].begin(), nbs[si].end(), P.nb.begin() + P.nb_off[si]);
  }
  P.last_master.assign(ns, INT_MIN);
  P.last_side.assign(ns, 0);
  P.ft.assign(ns, Vec3());
  P.fn_t.assign(ns, 0.0); P.fn_s.assign(ns, 0.0); P.pen_t.assign(ns, -1.0); P.pen_s.assign(ns, -1.0);
  P.last_n.assign(ns, Vec3());
  for (auto* v : {&P.cw2_t, &P.cw2_s, &P.md_t, &P.md_s, &P.kmin_t, &P.kmin_s}) v->assign(ns, 0.0);
  P.hit_t.resize(ns); P.hit_s.resize(ns); P.force_t.resize(ns); P.force_s.resize(ns);
  P.mu = st.contact_property >= 0 ? m.interactions[st.contact_property].friction : 0.0;
  double lsum = 0, hmax = 0;
  int cnt = 0;
  for (size_t t = 0; t < P.tri_h.size(); ++t) {
    lsum += norm(m.X[P.tri[3 * t + 1]] - m.X[P.tri[3 * t]]); ++cnt;
    hmax = std::max(hmax, P.tri_h[t]);
  }
  for (size_t s = 0; s < P.seg_r.size(); ++s) { lsum += norm(m.X[P.seg[2 * s + 1]] - m.X[P.seg[2 * s]]); ++cnt; hmax = std::max(hmax, P.seg_r[s]); }
  double lavg = cnt ? lsum / cnt : 1.0;
  P.margin = std::max(0.25 * lavg, hmax);
  P.cell = 2 * lavg + 2 * hmax + P.margin;
  active_ = !P.slaves.empty() && (!P.tri.empty() || !P.seg.empty());
  info("  general contact: %zu facets, %zu segments, %d nodes, friction %.3g", P.tri_h.size(), P.seg_r.size(), ns, P.mu);
  return active_;
}

void GeneralContact::compute(double dt, std::vector<double>& f, Energies& en, int increment) {
  std::fill(f.begin(), f.end(), 0.0);
  ncontacts_ = 0;
  if (!active_ || (dt <= 1e-20 && !p_->first_eval)) return;  // zero increment: only the initial adjustment
  Impl& P = *p_;
  const Analysis& A = A_;
  const int nn = A.nn;
  auto X = [&](int n) {
    return Vec3(A.X[3 * n] + A.u[6 * n] + A.du[6 * n], A.X[3 * n + 1] + A.u[6 * n + 1] + A.du[6 * n + 1],
                A.X[3 * n + 2] + A.u[6 * n + 2] + A.du[6 * n + 2]);
  };
  const int ns = (int)P.slaves.size();
  const int ntri = (int)P.tri_h.size(), nseg = (int)P.seg_r.size();
  // free shell edges moved back into the shell by the half thickness (current configuration)
  const int nfn = (int)P.fnode.size();
  for (int i = 0; i < nfn; ++i) {
    const int n = P.fnode[i];
    const int si = P.slave_of[n];
    const double h = si >= 0 ? P.slave_h[si] : 0.0;
    const Vec3 xn = X(n);
    Vec3 inw[2], along[2];
    double size = 1e300;
    const int ne = P.fn_off[i + 1] - P.fn_off[i];
    for (int j = 0; j < std::min(ne, 2); ++j) {
      const int* ed = &P.fedge[3 * P.fn_edge[P.fn_off[i] + j]];
      const int* q = &P.fe_elem[4 * ed[2]];
      const Vec3 a = X(q[0]), b = X(q[1]), c = X(q[2]);
      Vec3 nrm, ctr;
      if (q[3] >= 0) { const Vec3 d = X(q[3]); nrm = cross(c - a, d - b); ctr = (a + b + c + d) * 0.25; }
      else { nrm = cross(b - a, c - a); ctr = (a + b + c) * (1.0 / 3); }
      const Vec3 other = X(ed[0] == n ? ed[1] : ed[0]);
      along[j] = normalized(other - xn);
      Vec3 w = normalized(cross(nrm, along[j]));
      if (dot(w, ctr - xn) < 0) w = w * -1.0;
      inw[j] = w;
      size = std::min(size, P.fe_size[ed[2]]);
    }
    Vec3 d;
    if (ne == 1) d = inw[0] * h;
    else if (ne == 2) {
      const double c = dot(inw[0], inw[1]);
      if (c > -0.9) d = (inw[0] + inw[1]) * (h / (1 + c));  // both edges moved back by h
    }
    // at most half the facet size, and 0.6 x facet size x sine of the angle between the offset and each edge
    double dmax = 0.5 * size, dl = norm(d);
    if (dl > 0) {
      for (int j = 0; j < std::min(ne, 2); ++j) {
        const double ca = dot(d, along[j]) / dl;
        if (ca > 0) dmax = std::min(dmax, 0.6 * size * std::sqrt(std::max(0.0, 1 - ca * ca)));
      }
      if (dl > dmax) d = d * (dmax / dl);
    }
    P.shift[n] = d;
  }
  auto Xs = [&](int n) { return X(n) + P.shift[n]; };

  // ---- candidate search when nodes moved more than half the margin
  bool rebuild = P.xref.empty();
  if (!rebuild) {
    double dmax = 0;
    for (int si = 0; si < ns; ++si) {
      int n = P.slaves[si];
      dmax = std::max(dmax, norm2(Xs(n) - load3(&P.xref[3 * n])));
    }
    rebuild = std::sqrt(dmax) > 0.5 * P.margin;
  }
  if (rebuild) {
    ++P.rebuilds;
    P.xref.assign(3 * nn, 0.0);
    for (int si = 0; si < ns; ++si) store3(&P.xref[3 * P.slaves[si]], Xs(P.slaves[si]));
    const double cs = P.cell;
    auto key = [](int64_t i, int64_t j, int64_t k) { return ((i + (1 << 20)) << 42) | ((j + (1 << 20)) << 21) | (k + (1 << 20)); };
    std::vector<std::pair<int64_t, int>> cells;
    auto insert_box = [&](Vec3 lo, Vec3 hi, int id) {
      int64_t i0 = (int64_t)std::floor(lo.x / cs), i1 = (int64_t)std::floor(hi.x / cs);
      int64_t j0 = (int64_t)std::floor(lo.y / cs), j1 = (int64_t)std::floor(hi.y / cs);
      int64_t k0 = (int64_t)std::floor(lo.z / cs), k1 = (int64_t)std::floor(hi.z / cs);
      // a facet spanning a huge number of cells (or non-finite coordinates) means a diverging solution: stop
      // instead of filling the memory
      if (!(std::isfinite(lo.x + lo.y + lo.z + hi.x + hi.y + hi.z)) ||
          double(i1 - i0 + 1) * double(j1 - j0 + 1) * double(k1 - k0 + 1) > 1e6)
        throw std::runtime_error("contact: a surface facet spans more than 10^6 search cells (diverging solution)");
      for (int64_t i = i0; i <= i1; ++i)
        for (int64_t j = j0; j <= j1; ++j)
          for (int64_t k = k0; k <= k1; ++k) cells.push_back({key(i, j, k), id});
    };
    double hs_max = 0;
    for (double h : P.slave_h) hs_max = std::max(hs_max, h);
    for (int t = 0; t < ntri; ++t) {
      Vec3 a = Xs(P.tri[3 * t]), b = Xs(P.tri[3 * t + 1]), c = Xs(P.tri[3 * t + 2]);
      double inf = P.tri_h[t] + hs_max + P.margin;
      Vec3 lo(std::min({a.x, b.x, c.x}) - inf, std::min({a.y, b.y, c.y}) - inf, std::min({a.z, b.z, c.z}) - inf);
      Vec3 hi(std::max({a.x, b.x, c.x}) + inf, std::max({a.y, b.y, c.y}) + inf, std::max({a.z, b.z, c.z}) + inf);
      insert_box(lo, hi, t);
    }
    for (int s = 0; s < nseg; ++s) {
      Vec3 a = X(P.seg[2 * s]), b = X(P.seg[2 * s + 1]);
      double inf = P.seg_r[s] + hs_max + P.margin;
      Vec3 lo(std::min(a.x, b.x) - inf, std::min(a.y, b.y) - inf, std::min(a.z, b.z) - inf);
      Vec3 hi(std::max(a.x, b.x) + inf, std::max(a.y, b.y) + inf, std::max(a.z, b.z) + inf);
      insert_box(lo, hi, -1 - s);
    }
    std::sort(cells.begin(), cells.end());
    std::vector<std::vector<int>> lists(ns);
#pragma omp parallel for num_threads(A.nthreads) schedule(dynamic, 64)
    for (int si = 0; si < ns; ++si) {
      int n = P.slaves[si];
      Vec3 x = Xs(n);
      int64_t kk = key((int64_t)std::floor(x.x / cs), (int64_t)std::floor(x.y / cs), (int64_t)std::floor(x.z / cs));
      auto lo = std::lower_bound(cells.begin(), cells.end(), std::make_pair(kk, INT_MIN));
      for (auto it = lo; it != cells.end() && it->first == kk; ++it) {
        int id = it->second;
        auto near = [&](int node) {
          return std::binary_search(P.nb.begin() + P.nb_off[si], P.nb.begin() + P.nb_off[si + 1], node);
        };
        if (id >= 0) {
          const int* tn = &P.tri[3 * id];
          if (near(tn[0]) || near(tn[1]) || near(tn[2])) continue;  // facets touching the node's elements
          double w[3]; bool in;
          Vec3 q = closest_on_tri(x, Xs(tn[0]), Xs(tn[1]), Xs(tn[2]), w, in);
          if (norm(x - q) > P.slave_h[si] + P.tri_h[id] + P.margin) continue;
        } else {
          int s = -1 - id;
          if (near(P.seg[2 * s]) || near(P.seg[2 * s + 1])) continue;
        }
        lists[si].push_back(id);
      }
    }
    P.cand_off.assign(ns + 1, 0);
    for (int si = 0; si < ns; ++si) P.cand_off[si + 1] = P.cand_off[si] + (int)lists[si].size();
    P.cand.resize(P.cand_off[ns]);
    for (int si = 0; si < ns; ++si) std::copy(lists[si].begin(), lists[si].end(), P.cand.begin() + P.cand_off[si]);
  }

  // ---- penalty stiffness of a pair, from the frequency bound w_i^2 = 1/dt_e^2 of each node (stable increment
  // of its stiffest element) and the contact frequency W^2 = max(w_ref^2/3, 1.1 max_i w_i^2) of the pair: each
  // node gets k_i = 2 m_i W^2 (1 - r_i)/(1 - r_i/2), r_i = w_i^2/W^2, so that its own stiffness plus the
  // contact spring stays stable at the increment 1/W; the pair stiffness is the softer of the slave and of the
  // master facet (nodal flexibilities weighted by the contact point coordinates). A massive node stiffer
  // than the contact frequency is not limiting; light nodes of soft elements get a stiffer spring.
  const double wref2 = dt > 0 ? 1.0 / (dt * dt) : 0.0;
  P.node_w2.assign(nn, 0.0);
  for (auto& gp : A.groups) {
    const ElementGroup& G = *gp;
    for (int e = 0; e < G.size(); ++e) {
      const double d = G.dt_last[e] * std::sqrt(G.msf[e]);
      if (!(d > 0 && d < 1e29)) continue;
      for (int I = 0; I < G.nen; ++I) {
        const int nd = G.conn[e * G.nen + I];
        P.node_w2[nd] = std::max(P.node_w2[nd], 1.0 / (d * d));
      }
    }
  }
  auto node_k = [&](int nd, double w2) {
    const double mi = A.mt[nd];
    if (!(mi > 0)) return 1e300;
    const double r = std::min(P.node_w2[nd] / w2, 1.0 / 1.1);
    return 2 * mi * w2 * (1 - r) / (1 - 0.5 * r);
  };

  // ---- detection: deepest facet hit and deepest segment hit per slave node
  const bool first = P.first_eval;
  P.first_eval = false;
  int ninit = 0;
#pragma omp parallel for num_threads(A.nthreads) reduction(+ : ninit) schedule(dynamic, 64)
  for (int si = 0; si < ns; ++si) {
    int n = P.slaves[si];
    Vec3 x = Xs(n);
    double hs = P.slave_h[si];
    Hit bt, bs;
    for (int c = P.cand_off[si]; c < P.cand_off[si + 1]; ++c) {
      int id = P.cand[c];
      Hit h;
      h.master = id;
      if (id >= 0) {
        const int* tn = &P.tri[3 * id];
        Vec3 a = Xs(tn[0]), b = Xs(tn[1]), cc = Xs(tn[2]);
        bool in;
        Vec3 q = closest_on_tri(x, a, b, cc, h.w, in);
        Vec3 nrm = normalized(cross(b - a, cc - a));
        Vec3 d = x - q;
        double s = dot(d, nrm);
        double lat = in ? 0.0 : norm(d - nrm * s);  // distance beyond the facet boundary, in its plane
        const double hsum = hs + P.tri_h[id];
        if (lat <= 0.1 * hsum) {  // above the facet (or its edge): two-sided normal contact
          int side = s >= 0 ? 1 : -1;
          if (P.pen_t[si] >= 0) {  // in contact at the previous increment: keep the side of the surface
            double dn = dot(nrm, P.last_n[si]);
            if (std::fabs(dn) > 0.3) side = dn > 0 ? 1 : -1;
          }
          h.side = side;
          h.pen = hsum - side * s;
          h.n = nrm * side;
        } else {  // beyond an edge or a vertex: rounded edge (distance to the closest point)
          double dist = norm(d);
          h.side = 0;
          h.pen = hsum - dist;
          h.n = d / dist;
        }
      } else {
        int sg = -1 - id;
        Vec3 a = X(P.seg[2 * sg]), b = X(P.seg[2 * sg + 1]);
        Vec3 ab = b - a;
        double l = std::clamp(dot(x - a, ab) / std::max(norm2(ab), 1e-300), 0.0, 1.0);
        Vec3 d = x - (a + ab * l);
        double dist = norm(d);
        if (dist < 1e-14) continue;
        h.w[0] = 1 - l; h.w[1] = l;
        h.pen = hs + P.seg_r[sg] - dist;
        h.n = d / dist;
      }
      if (h.pen <= 0) continue;
      if (id >= 0) { if (h.pen > bt.pen) bt = h; }
      else if (h.pen > bs.pen) bs = h;
    }
    for (int pass = 0; pass < 2; ++pass) {
      Hit& h = pass ? bs : bt;
      if (h.master == INT_MIN) continue;
      int nm = h.master >= 0 ? 3 : 2;
      const int* mn = h.master >= 0 ? &P.tri[3 * h.master] : &P.seg[2 * (-1 - h.master)];
      double mm = 0;
      for (int i = 0; i < nm; ++i) mm += h.w[i] * A.mt[mn[i]];
      double msl = A.mt[n];
      double meff = std::min(msl > 0 ? msl : 1e300, mm > 0 ? mm : 1e300);
      if (meff >= 1e299) { h = Hit(); continue; }
      // contact frequency: the stiffest node of the pair, at least a third of the increment frequency
      double w2 = wref2 / 3;
      w2 = std::max(w2, 1.1 * P.node_w2[n]);
      for (int i = 0; i < nm; ++i) w2 = std::max(w2, 1.1 * P.node_w2[mn[i]]);
      const double ks = node_k(n, w2);
      double flex = 0;
      for (int i = 0; i < nm; ++i) if (h.w[i] > 0) flex += h.w[i] / node_k(mn[i], w2);
      const double kc = std::min(ks, flex > 0 ? 1.0 / flex : 1e300);
      if (!(kc > 0 && kc < 1e299)) { h = Hit(); continue; }
      h.k = kc;
      (pass ? P.cw2_s : P.cw2_t)[si] = w2;
      (pass ? P.md_s : P.md_t)[si] = meff;
    }
    P.hit_t[si] = bt;
    P.hit_s[si] = bs;
    if (first && (bt.master != INT_MIN || bs.master != INT_MIN)) ++ninit;
  }
  if (first && ninit) {
    // strain-free adjustment of the initial overclosures (as Abaqus does by default): the overclosed nodes are moved
    // out along the contact normal before any force is computed. Every node is both slave and master, so each side
    // of an overlap moves by half the penetration. The move is a displacement (reported in U, as in Abaqus) that
    // does not strain the elements (their strains come from the displacement increments).
    int nadj = 0;
    for (int si = 0; si < ns; ++si) {
      Hit& bt = P.hit_t[si];
      Hit& bs = P.hit_s[si];
      const Hit& h = (bt.master != INT_MIN && (bs.master == INT_MIN || bt.pen >= bs.pen)) ? bt : bs;
      if (h.master == INT_MIN || !(h.pen > 0)) continue;
      const int n = P.slaves[si];
      for (int k = 0; k < 3; ++k) A_.u[6 * n + k] += 0.5 * h.pen * h.n[k];
      bt = Hit(); bs = Hit();
      ++nadj;
    }
    warn("general contact: %d nodes initially overclosed, moved out without strain (half the penetration each side)",
         nadj);
  }
  if (dt <= 1e-20) {  // initial pass: no forces (all hits are cleared for the first real increment)
    for (int si = 0; si < ns; ++si) { P.hit_t[si] = Hit(); P.hit_s[si] = Hit(); }
    return;
  }

  // stiffness shared by several contacts on the same node: utilisation u of each node = sum over its contacts of
  // k_pair / k_node (weighted by the contact point coordinates on the master side); a node used beyond 1 has its
  // nodal stiffness divided by u in the pair stiffness
  std::vector<double> ksum(nn, 0.0);
  for (int si = 0; si < ns; ++si) {
    int n = P.slaves[si];
    for (int pass = 0; pass < 2; ++pass) {
      const Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
      if (h.master == INT_MIN) continue;
      int nm = h.master >= 0 ? 3 : 2;
      const int* mn = h.master >= 0 ? &P.tri[3 * h.master] : &P.seg[2 * (-1 - h.master)];
      const double w2 = (pass ? P.cw2_s : P.cw2_t)[si];
      ksum[n] += h.k / node_k(n, w2);
      for (int i = 0; i < nm; ++i) ksum[mn[i]] += h.w[i] * h.k / node_k(mn[i], w2);
    }
  }

  // ---- response: incremental normal force F = min(F_prev + k dpen, k pen) >= 0. Changes of k (with dt),
  // switches between facets and contacts found already deep do not create force jumps (no energy is
  // injected); a new contact only loads the approach of the current increment, so initial or lateral
  // overclosures act as offsets that are released as the nodes separate.
  double efric = 0, edamp = 0, w2mx = 0;
  int ncont = 0;
  long ndeep = 0;
#pragma omp parallel for num_threads(A.nthreads) reduction(+ : efric, edamp, ncont, ndeep) reduction(max : w2mx) \
    schedule(dynamic, 64)
  for (int si = 0; si < ns; ++si) {
    int n = P.slaves[si];
    const Vec3 dus = load3(&A.du[6 * n]);
    for (int pass = 0; pass < 2; ++pass) {
      Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
      double& fn = pass ? P.fn_s[si] : P.fn_t[si];
      double& pp = pass ? P.pen_s[si] : P.pen_t[si];
      Vec3& F = pass ? P.force_s[si] : P.force_t[si];
      const double fd_old = h.master == INT_MIN ? 0.0 : dot(F, h.n) - fn;  // damping force of the previous increment
      F = Vec3();
      if (h.master == INT_MIN) { pp = -1; fn = 0; continue; }
      int nm = h.master >= 0 ? 3 : 2;
      const int* mn = h.master >= 0 ? &P.tri[3 * h.master] : &P.seg[2 * (-1 - h.master)];
      // smallest pair stiffness since the contact started (the penalty never stiffens during a contact), 9/11 of
      // it, reduced by the sharing of the nodes
      const double w2 = (pass ? P.cw2_s : P.cw2_t)[si];
      double& kmin = (pass ? P.kmin_s : P.kmin_t)[si];
      kmin = pp >= 0 && kmin > 0 ? std::min(kmin, h.k) : h.k;
      double flex = 0;
      for (int i = 0; i < nm; ++i)
        if (h.w[i] > 0) flex += h.w[i] * std::max(1.0, ksum[mn[i]]) / node_k(mn[i], w2);
      const double kshared = std::min(node_k(n, w2) / std::max(1.0, ksum[n]), flex > 0 ? 1.0 / flex : 1e300);
      const double k = 9.0 / 11.0 * kmin * std::min(1.0, kshared / h.k);
      const double mdamp = (pass ? P.md_s : P.md_t)[si];
      w2mx = std::max(w2mx, w2);
      Vec3 dum;
      for (int i = 0; i < nm; ++i) dum += load3(&A.du[6 * mn[i]]) * h.w[i];
      double Fn;
      const double appr = std::max(0.0, -dot(dus - dum, h.n));  // relative approach in this increment
      if (pp >= 0) {
        // the force grows at most with the actual approach (a jump of the penetration after a change of
        // master facet or of contact geometry is not loaded)
        Fn = std::min(fn + k * std::min(h.pen - pp, appr), k * h.pen);
      } else {
        Fn = k * std::min(h.pen, appr);
        double hm = h.master >= 0 ? P.tri_h[h.master] : P.seg_r[-1 - h.master];
        if (!first && h.pen > 0.3 * (P.slave_h[si] + hm)) ++ndeep;
      }
      Fn = std::max(Fn, 0.0);
      fn = Fn;
      pp = h.pen;
      // normal viscous damping, 3 % of the critical damping of the contact spring with the lighter side of the
      // pair, c = 2 zeta sqrt(k m); the total normal force stays compressive
      double Ft = Fn;
      if (dt > 0) {
        const double c = 2 * 0.03 * std::sqrt(k * mdamp);
        Ft = std::max(0.0, Fn - c * dot(dus - dum, h.n) / dt);
        // dissipated by the damper, integrated as the work of the contact forces (trapezoidal rule)
        edamp += -0.5 * (fd_old + Ft - Fn) * dot(dus - dum, h.n);
      }
      F = h.n * Ft;
      ++ncont;
      const bool fric = pass == 0 || P.hit_t[si].master == INT_MIN;
      if (fric && P.mu > 0) {
        Vec3 dr = dus - dum;
        Vec3 dtg = dr - h.n * dot(dr, h.n);
        Vec3 f0 = P.last_master[si] == h.master ? P.ft[si] : Vec3();
        f0 -= h.n * dot(f0, h.n);
        Vec3 ftr = f0 - dtg * k;
        double fmax = P.mu * Fn, fm = norm(ftr);
        if (fm > fmax) {
          if (k > 0) efric += fmax * (fm - fmax) / k;
          ftr = fm > 0 ? ftr * (fmax / fm) : Vec3();
        }
        P.ft[si] = ftr;
        F += ftr;
      }
    }
    const Hit& bt = P.hit_t[si];
    const Hit& bs = P.hit_s[si];
    if (bt.master != INT_MIN) { P.last_master[si] = bt.master; P.last_side[si] = (signed char)bt.side; P.last_n[si] = bt.n; }
    else if (bs.master != INT_MIN) { P.last_master[si] = bs.master; P.last_side[si] = 0; P.last_n[si] = bs.n; }
    else { P.last_master[si] = INT_MIN; P.last_side[si] = 0; P.ft[si] = Vec3(); }
  }
  for (int si = 0; si < ns; ++si) {  // serial scatter
    int n = P.slaves[si];
    for (int pass = 0; pass < 2; ++pass) {
      const Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
      if (h.master == INT_MIN) continue;
      const Vec3& F = pass ? P.force_s[si] : P.force_t[si];
      add3(&f[6 * n], F);
      int nm = h.master >= 0 ? 3 : 2;
      const int* mn = h.master >= 0 ? &P.tri[3 * h.master] : &P.seg[2 * (-1 - h.master)];
      for (int i = 0; i < nm; ++i) add3(&f[6 * mn[i]], F * (-h.w[i]));
    }
  }
  en.friction += efric;
  P.w2_max = w2mx;
  // the damping part of the contact work is dissipation: reported with the viscous energy (as the contact
  // damping energy in ALLVD), and removed from the contact work so that the energy balance is unchanged
  en.viscous += edamp;
  en.contact -= edamp;
  ncontacts_ = ncont;
  P.d_deep_n += ndeep;
  if (P.diag < 0) P.diag = std::getenv("ABQ_CONTACT_DIAG") ? 1 : 0;
  if (P.diag > 0 && increment % 500 == 0)
    std::fprintf(stderr, "CDIAG inc %d contacts %d deep_new %ld  Wcontact %.4e fric %.4e\n", increment, ncont, P.d_deep_n, en.contact, en.friction);
  (void)increment;
}

double GeneralContact::dt_limit() const {
  return p_ && p_->w2_max > 0 ? 1.0 / std::sqrt(p_->w2_max) : 1e30;
}

void GeneralContact::save(BinOut& o) const {
  o.tag("contact.v3");
  bool has = p_ != nullptr;
  o.pod(has);
  if (!has) return;
  const Impl& P = *p_;
  o.vec(P.last_master); o.vec(P.last_side); o.vec(P.ft); o.vec(P.hit_t); o.vec(P.hit_s); o.vec(P.force_t); o.vec(P.force_s);
  o.vec(P.fn_t); o.vec(P.pen_t); o.vec(P.fn_s); o.vec(P.pen_s); o.vec(P.last_n);
  o.vec(P.kmin_t); o.vec(P.kmin_s);
  o.pod(P.first_eval);
}

void GeneralContact::load(BinIn& in) {
  std::string tg = in.read_tag();
  if (tg != "contact.v3" && tg != "contact.v2" && tg != "contact") throw std::runtime_error("restart file: contact section expected");
  bool has;
  in.pod(has);
  if (!has) return;
  if (!p_) throw std::runtime_error("restart: contact state without a contact definition");
  Impl& P = *p_;
  in.vec(P.last_master); in.vec(P.last_side); in.vec(P.ft); in.vec(P.hit_t); in.vec(P.hit_s); in.vec(P.force_t); in.vec(P.force_s);
  if (tg != "contact") {
    in.vec(P.fn_t); in.vec(P.pen_t); in.vec(P.fn_s); in.vec(P.pen_s); in.vec(P.last_n);
    if (tg == "contact.v3") { in.vec(P.kmin_t); in.vec(P.kmin_s); }
  } else {  // first format (total penalty with offsets per slave/master pair): convert
    struct KV { int64_t k; double v; };
    std::vector<KV> ip;
    in.vec(ip);
    std::unordered_map<int64_t, double> off;
    for (auto& kv : ip) off[kv.k] = kv.v;
    const int ns = (int)P.slaves.size();
    for (int si = 0; si < ns; ++si)
      for (int pass = 0; pass < 2; ++pass) {
        const Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
        double& fn = pass ? P.fn_s[si] : P.fn_t[si];
        double& pp = pass ? P.pen_s[si] : P.pen_t[si];
        if (h.master == INT_MIN) { fn = 0; pp = -1; continue; }
        auto it = off.find((int64_t)si * 4194304LL + (h.master + 2097152LL));
        pp = h.pen + (it != off.end() ? it->second : 0.0);
        fn = std::max(0.0, dot(pass ? P.force_s[si] : P.force_t[si], h.n));
        if (!pass) P.last_n[si] = h.n;
        else if (P.hit_t[si].master == INT_MIN) P.last_n[si] = h.n;
      }
  }
  in.pod(P.first_eval);
  P.xref.clear();  // candidate search redone at the first evaluation
}

}  // namespace abq
