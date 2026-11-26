// General contact for explicit dynamics.
//  - contact entities: shell facets (triangles, two-sided, offset by a reduced half thickness),
//    exterior faces of solid elements (one-sided, outward normal, no thickness), rigid elements (two-sided, no
//    thickness, master only), beam segments
//    (radius), slave nodes (all nodes of the included shells and beams, the nodes on the exterior
//    faces of the included solids); beyond a facet
//    boundary the contact is with the rounded edge (closest point), so no gap between facets;
//    beam nodes contact facets and segments, shell nodes only facets; the free shell edges contact
//    the beam segments (edge to segment, closest points of the two segments);
//  - beam segments contact each other edge to edge (closest points of the two segments, inside both; each pair once,
//    equal weighting), the ends of the segments by their nodes;
//  - free shell edges are rounded within the shell: the nodes of a free edge are moved back into the
//    shell plane by the half thickness (perpendicular to the edge, along the bisector at a corner;
//    at most half the facet size and 0.6 x facet size x sine of the angle to each edge), so that
//    the contact surface does not extend beyond the edge;
//  - node-to-facet and node-to-segment penalty from the stable increment of the elements of each node: contact
//    frequency W^2 = max(1/(3 (1.2 dt)^2), 1.1 max_i 1/dt_i^2) over the pair, nodal stiffnesses that keep each node
//    stable at 1/W, pair stiffness = softer of slave and master, 9/11 of its smallest value since the contact
//    started, shared between the contacts of a node; the increment is bounded by 1/W of the active contacts
//    (0.3 m/dt^2 between two critical nodes, stiffer for light nodes of soft elements);
//  - facet half thickness h = min(t/2, 0.3 x smallest facet height), segment radius
//    r = min(r, 0.3 x segment length) (thickness reduced on small facets);
//  - normal viscous damping of 3 % of the critical damping of the contact spring (slave mass, or reduced mass of the
//    slave and of the master facet when the master side sets the stiffness), not
//    bounded by the spring force (slightly tensile total force when overlapping nodes separate fast);
//  - softened contact (*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE = TABULAR / LINEAR / EXPONENTIAL): pressure of the
//    overclosure on the tributary area of the slave node, weight 0.7 per pass when the reverse pass is active too,
//    3 % damping on its stiffness, which also bounds the increment;
//  - Coulomb friction, tangential penalty stiffness equal to the normal one (stored force kept while the master stays
//    the same facet);
//  - incremental normal force (no jump when k changes, when a node switches facet or when a
//    contact is found deep); initial overclosures act as offsets released as the nodes separate;
//  - candidate pairs from a bucket grid, rebuilt when nodes moved more than half the margin;
//  - *CONTACT EXCLUSIONS: no contact between the nodes of one surface of an excluded pair and the facets or segments
//    of the other (both ways), nor between the elements of the two surfaces for edge contacts.
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

#include "contact/geom.hpp"
#include "core/log.hpp"
#include "solver/analysis.hpp"

namespace abq {

namespace {

struct Hit {
  int master = INT_MIN;  // >=0 triangle, <0: -1-segment, INT_MIN: none
  double pen = 0;
  double k = 0;          // penalty stiffness used
  Vec3 n;                // direction pushing the slave out
  double w[3] = {0, 0, 0};       // contact point on the master triangle or segment (stiffness, sharing)
  double wf[4] = {0, 0, 0, 0};   // weights of the nodes carrying the force (bilinear on quadrilaterals)
  int side = 0;
  double ws = 0;                 // edge contacts: contact point along the slave edge (0 at its first node)
};
struct HitV3 {  // restart files before contact.v4
  int master = INT_MIN;
  double pen = 0, k = 0;
  Vec3 n;
  double w[3] = {0, 0, 0};
  int side = 0;
};

// closest points of the segments p1-q1 and p2-q2 (Ericson): parameters s on the first, t on the second
void closest_seg_seg(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2, double& s, double& t) {
  const Vec3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
  const double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
  if (a <= 1e-300 && e <= 1e-300) { s = t = 0; return; }
  if (a <= 1e-300) { s = 0; t = std::clamp(f / e, 0.0, 1.0); return; }
  const double c = dot(d1, r);
  if (e <= 1e-300) { t = 0; s = std::clamp(-c / a, 0.0, 1.0); return; }
  const double b = dot(d1, d2), den = a * e - b * b;
  s = den > 1e-14 * a * e ? std::clamp((b * f - c * e) / den, 0.0, 1.0) : 0.0;
  t = (b * s + f) / e;
  if (t < 0) { t = 0; s = std::clamp(-c / a, 0.0, 1.0); }
  else if (t > 1) { t = 1; s = std::clamp((b - c) / a, 0.0, 1.0); }
}

}  // namespace

struct GeneralContact::Impl {
  std::vector<int> tri;       // 3 per triangle
  std::vector<double> tri_h;
  // per triangle, the nodes that carry its contact forces: the four nodes of its quadrilateral (bilinear weights,
  // half 0 = nodes 0-1-2, half 1 = nodes 0-2-3) or its own three nodes (4th = -1)
  std::vector<int> fac;
  std::vector<signed char> fac_half;
  std::vector<char> tri_solid;  // exterior face of a solid: one-sided (outward normal), no thickness
  std::vector<double> tri_depth;  // solid faces: deepest penetration still detected (beyond: other side of the solid)
  std::vector<int> st_off, st;    // per node: the exterior solid facets that contain it
  std::vector<Vec3> snorm;        // per node: outward normal of its solid surface (mean of its facets), zero if none
  // analytical rigid surfaces: facets (first triangle of each) around each of their nodes, and the current vertex
  // normals (mean of the facet normals: exact on an arc), interpolated over a facet as the contact normal
  std::vector<int> an_nodes, an_off, an;
  std::vector<Vec3> an_vn;
  // nodes of the master triangle or segment of a hit (weights w) and their number
  const int* stiff_nodes(const Hit& h, int& nm) const {
    nm = h.master >= 0 ? 3 : 2;
    return h.master >= 0 ? &tri[3 * h.master] : &seg[2 * (-1 - h.master)];
  }
  // nodes carrying the force of a hit (weights wf) and their number
  const int* master_nodes(const Hit& h, int& nm) const {
    if (h.master >= 0) { const int* q = &fac[4 * h.master]; nm = q[3] >= 0 ? 4 : 3; return q; }
    nm = 2;
    return &seg[2 * (-1 - h.master)];
  }
  std::vector<int> seg;       // 2 per segment
  std::vector<double> seg_r;
  std::vector<int> slaves;    // node ids
  std::vector<double> slave_h;
  // free shell edges: element nodes of each shell element touching a free edge (4, -1 padded) and its size,
  // free edges (two nodes and the element), nodes on free edges with their edges, current offset of every node
  std::vector<int> fe_elem, fe_gelem;  // fe_gelem: global element
  std::vector<double> fe_size;
  std::vector<int> fedge;
  std::vector<int> fnode, fn_off, fn_edge;
  std::vector<Vec3> shift;
  std::vector<int> slave_of;  // node -> slave index or -1
  std::vector<char> beam_node;  // node of a beam element (only these contact beam segments node-to-segment)
  std::vector<char> sfm_node;   // node of a surface element (SFM3D4R, SFM3D3)
  // softened pressure-overclosure law (null: hard contact, penalty), its clearance, contact area of each entity
  const SurfaceInteraction* soft = nullptr;
  double clear = 0;
  std::vector<double> area;
  std::vector<int> tri_elem, seg_elem;  // element of each triangle and segment
  // excluded surface pairs: elements and nodes of each side (self-contact exclusion: the same surface twice)
  struct Exclusion { std::vector<char> ea, eb, na, nb; };
  std::vector<Exclusion> excl;
  bool excluded(int node, int elem) const {
    for (const auto& x : excl) if ((x.na[node] && x.eb[elem]) || (x.nb[node] && x.ea[elem])) return true;
    return false;
  }
  bool excluded_elems(int e1, int e2) const {
    for (const auto& x : excl) if ((x.ea[e1] && x.eb[e2]) || (x.eb[e1] && x.ea[e2])) return true;
    return false;
  }
  // contact entities: the slave nodes (0 .. ns - 1), then, when there are beam segments, the free shell edges
  // (edge-to-segment contact: the shell perimeter against the beams, Abaqus/Explicit general contact), then the beam
  // segments (edge to edge between beams)
  int nent = 0, nfe = 0;
  std::vector<int> efedge;  // edges of the edge entities (two nodes and the element of fe_elem)
  std::vector<int> vowner;  // per node: the edge entity that owns its vertex contacts (lowest index)
  std::vector<int> bedge;   // segment of each beam edge entity
  // the two nodes of an edge entity
  const int* edge_nodes(int ei) const {
    const int k = ei - (int)slaves.size();
    return k < nfe ? &efedge[3 * k] : &seg[2 * bedge[k - nfe]];
  }
  int slave_nodes(int ei, const Hit& h, int sn[2], double sw[2]) const {
    if (ei < (int)slaves.size()) { sn[0] = slaves[ei]; sw[0] = 1; return 1; }
    const int* ed = edge_nodes(ei);
    sn[0] = ed[0]; sn[1] = ed[1]; sw[0] = 1 - h.ws; sw[1] = h.ws;
    return 2;
  }
  double ent_h(int ei) const {  // half thickness of a contact entity
    if (ei < (int)slaves.size()) return slave_h[ei];
    const int k = ei - (int)slaves.size();
    if (k >= nfe) return seg_r[bedge[k - nfe]];
    const int* ed = &efedge[3 * k];
    return std::min(slave_h[slave_of[ed[0]]], slave_h[slave_of[ed[1]]]);
  }
  void resize_state(int n) {
    last_master.resize(n, INT_MIN); last_side.resize(n, 0); ft.resize(n, Vec3());
    behind.resize(n, 0);
    fn_t.resize(n, 0.0); fn_s.resize(n, 0.0); pen_t.resize(n, -1.0); pen_s.resize(n, -1.0); last_n.resize(n, Vec3());
    for (auto* v : {&cw2_t, &cw2_s, &md_t, &md_s, &mm_t, &mm_s, &kmin_t, &kmin_s, &cdt_t, &cdt_s}) v->resize(n, 0.0);
    hit_t.resize(n); hit_s.resize(n); force_t.resize(n); force_s.resize(n);
  }
  std::vector<int> nb_off, nb;  // per slave: sorted nodes of its elements (facets touching them are excluded)
  std::vector<int> cand_off, cand;
  std::vector<double> xref;   // node positions at the last candidate search
  std::vector<int> tracked;   // nodes whose motion triggers a new search: slaves, facet and segment nodes
  double margin = 0, cell = 0, lavg = 1;
  std::vector<int> last_master;
  std::vector<signed char> last_side;
  std::vector<Vec3> ft;
  bool first_eval = true;
  // incremental normal force: force and penetration of the previous increment (pen < 0: no contact)
  std::vector<double> fn_t, pen_t, fn_s, pen_s;
  std::vector<char> behind;  // node behind the one-sided rigid surfaces while out of contact
  std::vector<Vec3> last_n;  // direction pushing the slave out at the previous increment
  // penalty: frequency bound of every node (1/dt_e^2 of its stiffest element), per slave and pass the contact
  // frequency, pair stiffness, damping mass and smallest pair stiffness since the contact started
  std::vector<double> node_w2;
  std::vector<double> node_g;  // per node: global-to-element ratio of the increment of the element setting node_w2
  std::vector<double> cdt_t, cdt_s;  // per slave and pass: stable increment bound of the contact
  std::vector<double> cw2_t, cw2_s, md_t, md_s, mm_t, mm_s, kmin_t, kmin_s;  // md: slave mass, mm: master mass
  double dt_min = 0;  // smallest stable increment bound of the contacts of the increment
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
  std::vector<double> node_h(m.nn(), 1e30);
  std::vector<double> node_area(m.nn(), 0.0);  // contact area of each node (share of its facets)
  std::vector<std::array<int, 4>> sfaces;  // faces of the included solids, oriented outward (4th node -1: triangle)
  std::vector<std::vector<int>> solid_elems;  // nodes of the included solid elements
  P.beam_node.assign(m.nn(), 0);
  std::vector<std::array<int, 4>> shells;  // included shell elements (nodes, -1 padded)
  std::vector<int> shell_elem, sface_elem;  // their elements, and those of the solid faces
  for (auto& gp : A_.groups) {
    ElementGroup& G = *gp;
    for (int e = 0; e < G.size(); ++e) {
      int ge = G.ids[e];
      if (!incl[ge]) continue;
      const int* cn = &G.conn[e * G.nen];
      if (is_shell_type(G.type)) {
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
        auto add = [&](int a, int b, int c, int half) {
          P.tri.insert(P.tri.end(), {a, b, c});
          P.tri_elem.push_back(ge);
          P.tri_h.push_back(h);
          P.fac.insert(P.fac.end(), {cn[0], cn[1], cn[2], G.nen == 4 ? cn[3] : -1});
          P.fac_half.push_back((signed char)half);
          P.tri_solid.push_back(0);
          P.tri_depth.push_back(0.0);
        };
        if (G.nen == 4) { add(cn[0], cn[1], cn[2], 0); add(cn[0], cn[2], cn[3], 1); }
        else add(cn[0], cn[1], cn[2], 0);
        for (int I = 0; I < G.nen; ++I) node_h[cn[I]] = std::min(node_h[cn[I]], h);
        {
          const double ar = G.nen == 4 ? 0.5 * norm(cross(p[2] - p[0], p[3] - p[1])) : 0.5 * norm(cross(p[1] - p[0], p[2] - p[0]));
          for (int I = 0; I < G.nen; ++I) node_area[cn[I]] += ar / G.nen;
        }
        shells.push_back({cn[0], cn[1], cn[2], G.nen == 4 ? cn[3] : -1});
        shell_elem.push_back(ge);
      } else if (G.type == ElType::B31) {
        double L = norm(m.X[cn[1]] - m.X[cn[0]]);
        double r = std::min(0.5 * G.contact_thickness(e), 0.3 * L);
        P.seg.insert(P.seg.end(), {cn[0], cn[1]});
        P.seg_elem.push_back(ge);
        P.seg_r.push_back(r);
        for (int I = 0; I < 2; ++I) {
          node_h[cn[I]] = std::min(node_h[cn[I]], r);
          P.beam_node[cn[I]] = 1;
          node_area[cn[I]] += L * r;  // half the length times the diameter
        }
      } else if (G.type == ElType::C3D8R || G.type == ElType::C3D8 || G.type == ElType::C3D8I ||
                 G.type == ElType::C3D4 || G.type == ElType::C3D6) {
        // faces in local node numbers, oriented by the element centroid below
        static const int hexf[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
        static const int tetf[4][4] = {{0, 1, 2, -1}, {0, 1, 3, -1}, {1, 2, 3, -1}, {2, 0, 3, -1}};
        static const int wedf[5][4] = {{0, 1, 2, -1}, {3, 4, 5, -1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
        // prisms computed as bricks with a doubled node (explicit): faces of the prism
        const bool degen = G.input_type == ElType::C3D6 && G.nen == 8;
        int c6[6];
        if (degen) { static const int k6[6] = {0, 2, 3, 4, 6, 7}; for (int I = 0; I < 6; ++I) c6[I] = cn[k6[I]]; cn = c6; }
        const int nen = degen ? 6 : G.nen;
        const bool hex = nen == 8, wedge = nen == 6;
        Vec3 ctr;
        for (int I = 0; I < nen; ++I) ctr += m.X[cn[I]] * (1.0 / nen);
        for (int f = 0; f < (hex ? 6 : wedge ? 5 : 4); ++f) {
          const int* lf = hex ? hexf[f] : wedge ? wedf[f] : tetf[f];
          std::array<int, 4> q = {cn[lf[0]], cn[lf[1]], cn[lf[2]], lf[3] >= 0 ? cn[lf[3]] : -1};
          const Vec3 a = m.X[q[0]], b = m.X[q[1]], c = m.X[q[2]];
          const Vec3 nrm = q[3] >= 0 ? cross(c - a, m.X[q[3]] - b) : cross(b - a, c - a);
          const Vec3 fc = q[3] >= 0 ? (a + b + c + m.X[q[3]]) * 0.25 : (a + b + c) * (1.0 / 3);
          if (dot(nrm, fc - ctr) < 0) std::swap(q[1], q[q[3] >= 0 ? 3 : 2]);
          sfaces.push_back(q);
          sface_elem.push_back(ge);
        }
        solid_elems.emplace_back(cn, cn + nen);
      }
    }
  }
  // rigid elements (R3D4, R3D3) and surface elements (SFM3D4R, SFM3D3): two-sided facets without thickness; facets of
  // analytical rigid surfaces: one-sided (normal to the left of the profile direction), a node coming from behind
  // never contacts them; the nodes of rigid facets carry no mass and are never slaves (always the master), those of
  // surface elements are slaves too
  P.sfm_node.assign(m.nn(), 0);
  for (int ge = 0; ge < m.ne(); ++ge)
    if (is_surface_type(m.elems[ge].type)) for (int n : m.elems[ge].nodes) P.sfm_node[n] = 1;
  std::vector<char> one_sided(m.ne(), 0);  // facets of analytical rigid surfaces: one-sided
  for (const auto& cp : m.couplings) if (!cp.analytic.empty()) for (int e : cp.rigid_elems) one_sided[e] = 1;
  for (int ge = 0; ge < m.ne(); ++ge) {
    const Element& el = m.elems[ge];
    if (!(is_rigid_type(el.type) || is_surface_type(el.type)) || !incl[ge]) continue;
    double lmax = 0;
    for (size_t I = 0; I < el.nodes.size(); ++I) lmax = std::max(lmax, norm(m.X[el.nodes[(I + 1) % el.nodes.size()]] - m.X[el.nodes[I]]));
    const int* cn = el.nodes.data();
    const int nv = (int)el.nodes.size();
    auto add = [&](int a, int b, int c, int half) {
      P.tri.insert(P.tri.end(), {a, b, c});
      P.tri_elem.push_back(ge);
      P.tri_h.push_back(0.0);
      P.fac.insert(P.fac.end(), {cn[0], cn[1], cn[2], nv == 4 ? cn[3] : -1});
      P.fac_half.push_back((signed char)half);
      P.tri_solid.push_back(one_sided[ge] ? 2 : 0);
      P.tri_depth.push_back(one_sided[ge] ? 0.5 * lmax : 0.0);
    };
    if (nv == 4) { add(cn[0], cn[1], cn[2], 0); add(cn[0], cn[2], cn[3], 1); }
    else if (nv == 3) add(cn[0], cn[1], cn[2], 0);
    if (is_surface_type(el.type)) {  // slaves too, no thickness
      const Vec3 a = m.X[cn[0]], b = m.X[cn[1]], c = m.X[cn[2]];
      const double ar = nv == 4 ? 0.5 * norm(cross(c - a, m.X[cn[3]] - b)) : 0.5 * norm(cross(b - a, c - a));
      for (int I = 0; I < nv; ++I) { node_h[cn[I]] = 0.0; node_area[cn[I]] += ar / nv; }
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
        P.fe_gelem.push_back(shell_elem[k]);
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
  {  // exterior faces of the solids: faces used by a single included solid element
    std::vector<std::pair<std::array<int, 4>, int>> keyed;
    for (int k = 0; k < (int)sfaces.size(); ++k) {
      std::array<int, 4> key = sfaces[k];
      std::sort(key.begin(), key.end());
      keyed.push_back({key, k});
    }
    std::sort(keyed.begin(), keyed.end());
    for (size_t i = 0; i < keyed.size(); ++i) {
      if ((i > 0 && keyed[i - 1].first == keyed[i].first) || (i + 1 < keyed.size() && keyed[i + 1].first == keyed[i].first))
        continue;
      const auto& q = sfaces[keyed[i].second];
      const int nv = q[3] >= 0 ? 4 : 3;
      // depth limit: half the smallest height of the face (a node further behind is on another side of the solid)
      double hmin = 1e300;
      for (int I = 0; I < nv; ++I) {
        const Vec3 a = m.X[q[I]], b = m.X[q[(I + 1) % nv]], c = m.X[q[(I + 2) % nv]];
        const double l = norm(b - a);
        if (l > 0) hmin = std::min(hmin, norm(cross(b - a, c - a)) / l);
      }
      for (int half = 0; half < nv - 2; ++half) {
        P.tri_depth.push_back(0.5 * hmin);
        P.tri.insert(P.tri.end(), {q[0], q[half + 1], q[half + 2]});
        P.tri_elem.push_back(sface_elem[keyed[i].second]);
        P.tri_h.push_back(0.0);
        P.fac.insert(P.fac.end(), {q[0], q[1], q[2], q[3]});
        P.fac_half.push_back((signed char)half);
        P.tri_solid.push_back(1);
      }
      for (int I = 0; I < nv; ++I) node_h[q[I]] = 0.0;
      const Vec3 a = m.X[q[0]], b = m.X[q[1]], c = m.X[q[2]];
      const double ar = nv == 4 ? 0.5 * norm(cross(c - a, m.X[q[3]] - b)) : 0.5 * norm(cross(b - a, c - a));
      for (int I = 0; I < nv; ++I) node_area[q[I]] += ar / nv;
    }
  }
  {  // facets of the solid surfaces around each node (folds of the surface: closest edge or vertex shared)
    std::vector<int> cnt(m.nn() + 1, 0);
    for (size_t t = 0; t < P.tri_h.size(); ++t)
      if (P.tri_solid[t] == 1) for (int k = 0; k < 3; ++k) ++cnt[P.tri[3 * t + k] + 1];
    for (int n = 0; n < m.nn(); ++n) cnt[n + 1] += cnt[n];
    P.st_off = cnt;
    P.st.assign(cnt[m.nn()], 0);
    for (size_t t = 0; t < P.tri_h.size(); ++t)
      if (P.tri_solid[t] == 1) for (int k = 0; k < 3; ++k) P.st[cnt[P.tri[3 * t + k]]++] = (int)t;
  }
  {
    std::vector<std::vector<int>> lists(m.nn());
    for (size_t t = 0; t < P.tri_h.size(); ++t) {
      if (P.tri_solid[t] != 2 || P.fac_half[t] != 0) continue;
      for (int k = 0; k < 4; ++k) if (P.fac[4 * t + k] >= 0) lists[P.fac[4 * t + k]].push_back((int)t);
    }
    P.an_off.assign(1, 0);
    for (int n = 0; n < m.nn(); ++n) {
      if (lists[n].empty()) continue;
      P.an_nodes.push_back(n);
      P.an.insert(P.an.end(), lists[n].begin(), lists[n].end());
      P.an_off.push_back((int)P.an.size());
    }
    P.an_vn.assign(m.nn(), Vec3());
  }
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
    for (const auto& el : solid_elems) link(el.data(), (int)el.size());
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
  int nedge = 0;  // free shell edges whose nodes are all slaves (edge-to-segment contacts)
  if (!P.seg.empty()) {
    std::vector<int> fe;
    for (size_t i = 0; i < P.fedge.size(); i += 3)
      if (P.slave_of[P.fedge[i]] >= 0 && P.slave_of[P.fedge[i + 1]] >= 0) fe.insert(fe.end(), {P.fedge[i], P.fedge[i + 1], P.fedge[i + 2]});
    nedge = (int)fe.size() / 3;
    P.efedge = fe;
    P.vowner.assign(m.nn(), INT_MAX);
    for (int k = 0; k < nedge; ++k)
      for (int j = 0; j < 2; ++j) P.vowner[fe[3 * k + j]] = std::min(P.vowner[fe[3 * k + j]], ns + k);
  }
  P.nfe = nedge;
  for (int sg = 0; sg < (int)P.seg_r.size(); ++sg)
    if (P.slave_of[P.seg[2 * sg]] >= 0 && P.slave_of[P.seg[2 * sg + 1]] >= 0) P.bedge.push_back(sg);
  P.nent = ns + nedge + (int)P.bedge.size();
  P.resize_state(P.nent);
  P.area.assign(P.nent, 0.0);
  for (int si = 0; si < ns; ++si) P.area[si] = node_area[P.slaves[si]];
  for (int ei = ns; ei < P.nent; ++ei) {  // edge contacts: length times the thickness (diameter)
    const int* ed = P.edge_nodes(ei);
    P.area[ei] = norm(m.X[ed[1]] - m.X[ed[0]]) * 2 * P.ent_h(ei);
  }
  {
    std::vector<char> tr(m.nn(), 0);
    for (int n : P.slaves) tr[n] = 1;
    for (int n : P.tri) tr[n] = 1;
    for (int n : P.seg) tr[n] = 1;
    for (int n = 0; n < m.nn(); ++n) if (tr[n]) P.tracked.push_back(n);
  }
  for (const auto& pr : st.contact_exclusions) {
    Impl::Exclusion x;
    x.ea.assign(m.ne(), 0); x.na.assign(m.nn(), 0);
    const Surface& sa = m.surfaces.at(pr.first);
    auto fill = [&](const Surface& S, std::vector<char>& em, std::vector<char>& nm) {
      for (int n : S.nodes) nm[n] = 1;
      for (const FaceRef& f : S.faces) {
        const Element& el = m.elems[f.elem];
        em[f.elem] = 1;
        int loc[4];
        const int k = face_local_nodes(el.type, f.face, (int)el.nodes.size(), loc);
        for (int i = 0; i < k; ++i) nm[el.nodes[loc[i]]] = 1;
      }
    };
    fill(sa, x.ea, x.na);
    if (pr.second.empty()) { x.eb = x.ea; x.nb = x.na; }
    else { x.eb.assign(m.ne(), 0); x.nb.assign(m.nn(), 0); fill(m.surfaces.at(pr.second), x.eb, x.nb); }
    P.excl.push_back(std::move(x));
  }
  P.mu = st.contact_property >= 0 ? m.interactions[st.contact_property].friction : 0.0;
  if (st.contact_property >= 0 && m.interactions[st.contact_property].law != SurfaceInteraction::HARD) {
    P.soft = &m.interactions[st.contact_property];
    P.clear = P.soft->clearance();
  }
  double lsum = 0, hmax = 0;
  int cnt = 0;
  for (size_t t = 0; t < P.tri_h.size(); ++t) {
    lsum += norm(m.X[P.tri[3 * t + 1]] - m.X[P.tri[3 * t]]); ++cnt;
    hmax = std::max(hmax, P.tri_h[t]);
  }
  for (size_t s = 0; s < P.seg_r.size(); ++s) { lsum += norm(m.X[P.seg[2 * s + 1]] - m.X[P.seg[2 * s]]); ++cnt; hmax = std::max(hmax, P.seg_r[s]); }
  double lavg = cnt ? lsum / cnt : 1.0;
  P.lavg = lavg;
  P.margin = std::max(0.25 * lavg, hmax) + P.clear;
  P.cell = 2 * lavg + 2 * hmax + P.margin;
  active_ = !P.slaves.empty() && (!P.tri.empty() || !P.seg.empty());
  info("  general contact: %zu facets, %zu segments, %d nodes, %d shell edges against segments, friction %.3g%s%s",
       P.tri_h.size(), P.seg_r.size(), ns, nedge, P.mu, P.excl.empty() ? "" : ", with exclusions",
       P.soft ? ", softened" : "");
  return active_;
}

// shared by compute() and implicit_points(): free shell edges moved back into the shell, vertex normals of the
// analytical surfaces, candidate search when the nodes moved more than half the margin (configuration X + u + du)
void GeneralContact::prepare() {
  Impl& P = *p_;
  const Analysis& A = A_;
  const int nn = A.nn;
  auto X = [&](int n) {
    return Vec3(A.X[3 * n] + A.u[6 * n] + A.du[6 * n], A.X[3 * n + 1] + A.u[6 * n + 1] + A.du[6 * n + 1],
                A.X[3 * n + 2] + A.u[6 * n + 2] + A.du[6 * n + 2]);
  };
  const int ns = (int)P.slaves.size(), nent = P.nent;
  const int ntri = (int)P.tri_h.size(), nseg = (int)P.seg_r.size();
  // outward normals of the slave nodes on solid surfaces (mean of the normals of their exterior faces)
  if (!P.st_off.empty()) {
    P.snorm.assign(nn, Vec3());
    for (int si = 0; si < ns; ++si) {
      const int n = P.slaves[si];
      Vec3 s;
      for (int j = P.st_off[n]; j < P.st_off[n + 1]; ++j) {
        const int* t = &P.tri[3 * P.st[j]];
        s += normalized(cross(X(t[1]) - X(t[0]), X(t[2]) - X(t[0])));
      }
      if (norm(s) > 1e-12) P.snorm[n] = normalized(s);
    }
  }
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
  for (size_t i = 0; i < P.an_nodes.size(); ++i) {  // vertex normals of the analytical rigid surfaces
    Vec3 s;
    for (int j = P.an_off[i]; j < P.an_off[i + 1]; ++j) {
      const int* tn = &P.tri[3 * P.an[j]];
      s += normalized(cross(X(tn[1]) - X(tn[0]), X(tn[2]) - X(tn[0])));
    }
    P.an_vn[P.an_nodes[i]] = normalized(s);
  }

  // ---- candidate search when nodes moved more than half the margin
  // (slave nodes and master-only nodes such as those of rigid facets)
  bool rebuild = P.xref.empty();
  if (!rebuild) {
    double dmax = 0;
    for (int n : P.tracked) dmax = std::max(dmax, norm2(Xs(n) - load3(&P.xref[3 * n])));
    rebuild = std::sqrt(dmax) > 0.5 * P.margin;
  }
  if (rebuild) {
    ++P.rebuilds;
    P.xref.assign(3 * nn, 0.0);
    for (int n : P.tracked) store3(&P.xref[3 * n], Xs(n));
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
          if (!P.excl.empty() && P.excluded(n, P.tri_elem[id])) continue;
          double w[3]; bool in;
          Vec3 q = closest_on_tri(x, Xs(tn[0]), Xs(tn[1]), Xs(tn[2]), w, in);
          if (norm(x - q) > P.slave_h[si] + P.tri_h[id] + P.clear + P.margin) continue;
        } else {
          int s = -1 - id;
          if (near(P.seg[2 * s]) || near(P.seg[2 * s + 1])) continue;
          if (!P.excl.empty() && P.excluded(n, P.seg_elem[s])) continue;
          // shell nodes do not contact beam segments (the shell perimeter does, edge to segment; free S4R sheet
          // dropped on a clamped B31 rod in Abaqus/Explicit: only the edges crossing the rod are stopped)
          if (!P.beam_node[n]) continue;
        }
        lists[si].push_back(id);
      }
    }
    lists.resize(nent);
#pragma omp parallel for num_threads(A.nthreads) schedule(dynamic, 64)
    for (int ei = ns; ei < nent; ++ei) {  // edges (free shell edges, beam segments): segments in the cells crossed
      const bool shell_edge = ei < ns + P.nfe;
      const int* ed = P.edge_nodes(ei);
      const int own = shell_edge ? -1 : P.bedge[ei - ns - P.nfe];  // segment of a beam edge
      const int* q = shell_edge ? &P.fe_elem[4 * ed[2]] : ed;
      const int nq = shell_edge ? 4 : 2;
      const int eown = shell_edge ? P.fe_gelem[ed[2]] : P.seg_elem[own];
      const Vec3 a = Xs(ed[0]), b = Xs(ed[1]);
      const int64_t i0 = (int64_t)std::floor(std::min(a.x, b.x) / cs), i1 = (int64_t)std::floor(std::max(a.x, b.x) / cs);
      const int64_t j0 = (int64_t)std::floor(std::min(a.y, b.y) / cs), j1 = (int64_t)std::floor(std::max(a.y, b.y) / cs);
      const int64_t k0 = (int64_t)std::floor(std::min(a.z, b.z) / cs), k1 = (int64_t)std::floor(std::max(a.z, b.z) / cs);
      if (double(i1 - i0 + 1) * double(j1 - j0 + 1) * double(k1 - k0 + 1) > 1e4) continue;
      auto& L = lists[ei];
      for (int64_t i = i0; i <= i1; ++i)
        for (int64_t j = j0; j <= j1; ++j)
          for (int64_t k = k0; k <= k1; ++k) {
            const int64_t kk = key(i, j, k);
            for (auto it = std::lower_bound(cells.begin(), cells.end(), std::make_pair(kk, INT_MIN));
                 it != cells.end() && it->first == kk; ++it) {
              const int id = it->second;
              if (id >= 0) continue;
              const int sg = -1 - id;
              if (sg <= own) continue;  // beam edges: each pair of segments once
              bool adj = false;  // segment attached to the element of the edge
              for (int v = 0; v < nq && q[v] >= 0; ++v) adj |= q[v] == P.seg[2 * sg] || q[v] == P.seg[2 * sg + 1];
              if (!adj && !(P.excl.size() && P.excluded_elems(eown, P.seg_elem[sg]))) L.push_back(id);
            }
          }
      std::sort(L.begin(), L.end());
      L.erase(std::unique(L.begin(), L.end()), L.end());
    }
    P.cand_off.assign(nent + 1, 0);
    for (int si = 0; si < nent; ++si) P.cand_off[si + 1] = P.cand_off[si] + (int)lists[si].size();
    P.cand.resize(P.cand_off[nent]);
    for (int si = 0; si < nent; ++si) std::copy(lists[si].begin(), lists[si].end(), P.cand.begin() + P.cand_off[si]);
  }
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
  const int ns = (int)P.slaves.size(), nent = P.nent;

  prepare();
  auto Xs = [&](int n) { return X(n) + P.shift[n]; };

  // ---- penalty stiffness of a pair, from the frequency bound w_i^2 = 1/dt_e^2 of each node (stable increment
  // of its stiffest element) and the contact frequency W^2 = max(w_ref^2/3, 1.1 max_i w_i^2) of the pair, w_ref =
  // 1/(1.2 dt) (floor identified on a soft sheet whose elements are far from critical): each
  // node gets k_i = 2 m_i W^2 (1 - r_i)/(1 - r_i/2), r_i = w_i^2/W^2, so that its own stiffness plus the
  // contact spring stays stable at the increment 1/W; the pair stiffness is the softer of the slave and of the
  // master facet (nodal flexibilities weighted by the contact point coordinates). A massive node stiffer
  // than the contact frequency is not limiting; light nodes of soft elements get a stiffer spring.
  // Solid elements run under the global estimate of the increment, longer than their own bound by a ratio g: a
  // slave node keeps the bound of its elements, the nodes of a master facet take the global one (w_i / g); when the
  // slave sets W the stable increment of the contact is g/W. Calibrated on C3D8R blocks dropped on C3D8R bases
  // and on an S4R plate in Abaqus/Explicit: same force histories increment by increment, increment shortened by
  // about 4 % during the contact.
  const double wref2 = dt > 0 ? 1.0 / (dt * dt) : 0.0;
  const double wfloor2 = wref2 / (3 * 1.2 * 1.2);
  P.node_w2.assign(nn, 0.0);
  P.node_g.assign(nn, 1.0);
  for (auto& gp : A.groups) {
    const ElementGroup& G = *gp;
    const double gs = is_solid_type(G.type) ? A.solid_dt_scale : 1.0;
    for (int e = 0; e < G.size(); ++e) {
      const double d = G.dt_last[e] * std::sqrt(G.msf[e]);
      if (!(d > 0 && d < 1e29)) continue;
      for (int I = 0; I < G.nen; ++I) {
        const int nd = G.conn[e * G.nen + I];
        if (1.0 / (d * d) > P.node_w2[nd]) { P.node_w2[nd] = 1.0 / (d * d); P.node_g[nd] = gs; }
      }
    }
  }
  // frequency bound of a node as slave (sl = true) or as a node of the master facet
  auto node_w = [&](int nd, bool sl) {
    const double g = P.node_g[nd];
    return sl ? P.node_w2[nd] : P.node_w2[nd] / (g * g);
  };
  // massless nodes: rigid facets are infinitely stiff, massless surface elements (SFM without DENSITY) carry no
  // contact stiffness (Abaqus/Explicit: a crimper of massless SFM3D4R elements passes through the tube)
  auto node_k = [&](int nd, double w2, bool sl) {
    const double mi = A.mt[nd];
    if (!(mi > 0)) return P.sfm_node[nd] ? 0.0 : 1e300;
    const double r = std::min(node_w(nd, sl) / w2, 1.0 / 1.1);
    return 2 * mi * w2 * (1 - r) / (1 - 0.5 * r);
  };

  // ---- detection: deepest facet hit and deepest segment hit per slave node
  const bool first = P.first_eval;
  P.first_eval = false;
  int ninit = 0;
#pragma omp parallel for num_threads(A.nthreads) reduction(+ : ninit) schedule(dynamic, 64)
  for (int si = 0; si < nent; ++si) {
    Hit bt, bs;
    if (si >= ns) {  // edge against beam segments: closest points of the two segments
      const bool shell_edge = si < ns + P.nfe;
      const int* ed = P.edge_nodes(si);
      const Vec3 p1 = Xs(ed[0]), q1 = Xs(ed[1]);
      const double he = P.ent_h(si);
      for (int c = P.cand_off[si]; c < P.cand_off[si + 1]; ++c) {
        const int id = P.cand[c], sg = -1 - id;
        const Vec3 a = X(P.seg[2 * sg]), b = X(P.seg[2 * sg + 1]);
        double se, te;
        closest_seg_seg(p1, q1, a, b, se, te);
        // shell edges, closest point at a vertex: one contact, held by one of the edges meeting there; beam edges:
        // inside both segments only (the nodes contact the segments)
        if (shell_edge && ((se <= 0 && P.vowner[ed[0]] != si) || (se >= 1 && P.vowner[ed[1]] != si))) continue;
        if (!shell_edge && !(se > 0 && se < 1 && te > 0 && te < 1)) continue;
        const Vec3 d = (p1 + (q1 - p1) * se) - (a + (b - a) * te);
        const double dist = norm(d);
        if (dist < 1e-14) continue;
        const double pen = he + P.seg_r[sg] + P.clear - dist;
        if (pen <= 0 || pen <= bs.pen) continue;
        Hit h;
        h.master = id;
        h.w[0] = h.wf[0] = 1 - te; h.w[1] = h.wf[1] = te;
        h.ws = se;
        h.pen = pen;
        h.n = d / dist;
        bs = h;
      }
    }
    const int n0 = si < ns ? P.slaves[si] : 0;
    const Vec3 x = Xs(n0);
    Hit bsol;
    int bsol_tier = 3;
    double bsol_appr = 0;
    double s_rigid = 1e300;  // signed distance to the nearest one-sided rigid facet (behind: < 0)
    const double hs = si < ns ? P.slave_h[si] : 0.0;
    for (int c = P.cand_off[si]; si < ns && c < P.cand_off[si + 1]; ++c) {
      int id = P.cand[c];
      Hit h;
      h.master = id;
      if (id >= 0) {
        const int* tn = &P.tri[3 * id];
        Vec3 a = Xs(tn[0]), b = Xs(tn[1]), cc = Xs(tn[2]);
        bool in;
        Vec3 q = closest_on_tri(x, a, b, cc, h.w, in);
        const double wc[3] = {h.w[0], h.w[1], h.w[2]};  // closest point: vertex (one weight) or edge (two)
        // weights of the master nodes from the coordinates of the contact point on the facet of the nodes
        // themselves (not the one drawn in by the rounding of free edges): the stiffness and its sharing use the
        // triangle, the force goes to the four nodes of a quadrilateral with the bilinear weights, (xi, eta) from
        // the corners (0,0) (1,0) (1,1) or (0,0) (1,1) (0,1) of the half (clamped S4R plate facet in
        // Abaqus/Explicit, small sheet dropped at 6 places: reactions of the plate nodes within 0.005, stiffness
        // of a light plate within 6 %)
        plane_coords(q, X(tn[0]), X(tn[1]), X(tn[2]), h.w);
        double sw = 0;
        for (int k = 0; k < 3; ++k) sw += (h.w[k] = std::max(h.w[k], 0.0));
        for (int k = 0; k < 3; ++k) h.w[k] /= sw;
        if (P.fac[4 * id + 3] >= 0) {
          const double xi = P.fac_half[id] ? h.w[1] : h.w[1] + h.w[2], eta = P.fac_half[id] ? h.w[1] + h.w[2] : h.w[2];
          h.wf[0] = (1 - xi) * (1 - eta); h.wf[1] = xi * (1 - eta); h.wf[2] = xi * eta; h.wf[3] = (1 - xi) * eta;
        } else {
          for (int k = 0; k < 3; ++k) h.wf[k] = h.w[k];
        }
        Vec3 nrm = normalized(cross(b - a, cc - a));
        Vec3 d = x - q;
        double s = dot(d, nrm);
        double lat = in ? 0.0 : norm(d - nrm * s);  // distance beyond the facet boundary, in its plane
        const double hsum = hs + P.tri_h[id] + P.clear;  // softened contact: from the clearance on
        // a node on the surface of a solid contacts only the solid faces turned towards it (its own outward
        // normal opposed to theirs, beyond 100 degrees): a node on the edge of a punch pressed into it is also on
        // the plane of the punch's side faces
        if (P.tri_solid[id] == 1 && !P.snorm.empty() && norm2(P.snorm[n0]) > 0 && dot(nrm, P.snorm[n0]) > -0.17) continue;
        if (P.tri_solid[id]) {  // solid face: one-sided, outward normal
          // no thickness: within the facet up to round-off (a node on the diagonal between the two triangles of a
          // quadrilateral face is otherwise outside both)
          double wp[3];
          plane_coords(x, a, b, cc, wp);
          if (lat <= 0.1 * hsum || std::min({wp[0], wp[1], wp[2]}) >= -1e-9) {
            if (P.tri_solid[id] == 2) {  // one-sided rigid facet: a node coming from behind never contacts it
              if (std::fabs(s) < std::fabs(s_rigid)) s_rigid = s;
              if (P.pen_t[si] < 0 && (s < 0 || P.behind[si])) continue;
            }
            if (s < -P.tri_depth[id]) continue;
            h.side = 1;
            h.pen = hsum - s;
            h.n = nrm;
            if (P.tri_solid[id] == 2) {  // analytical surface: normal of the smooth surface, from the vertex normals
              const Vec3 ni = P.an_vn[tn[0]] * h.w[0] + P.an_vn[tn[1]] * h.w[1] + P.an_vn[tn[2]] * h.w[2];
              if (norm(ni) > 0.5) h.n = normalized(ni);
            }
          } else if (s < 0) {
            // behind the face plane and beyond its boundary: inside the solid only under a concave fold of its
            // surface (a node below the fold line of a face bent into a valley projects outside both triangles),
            // i.e. behind the planes of all the facets that share the closest edge or vertex, with that edge or
            // vertex also their closest point; the penetration is the distance to that edge or vertex
            if (P.tri_solid[id] != 1) continue;
            const double dist = norm(d);
            if (!(dist > 1e-14) || dist > P.tri_depth[id]) continue;
            int fv[2], nf = 0;
            for (int k = 0; k < 3; ++k) if (wc[k] > 1e-12 && nf < 2) fv[nf++] = tn[k];
            if (nf == 0) continue;
            int nfold = 0;
            bool behind_all = true;
            for (int j = P.st_off[fv[0]]; j < P.st_off[fv[0] + 1] && behind_all; ++j) {
              const int t2 = P.st[j];
              if (t2 == id) continue;
              const int* t2n = &P.tri[3 * t2];
              if (nf == 2 && t2n[0] != fv[1] && t2n[1] != fv[1] && t2n[2] != fv[1]) continue;
              const Vec3 a2 = Xs(t2n[0]), b2 = Xs(t2n[1]), c2 = Xs(t2n[2]);
              const Vec3 n2 = normalized(cross(b2 - a2, c2 - a2));
              if (dot(x - a2, n2) >= 0) { behind_all = false; break; }
              double w2[3];
              bool in2;
              closest_on_tri(x, a2, b2, c2, w2, in2);
              for (int k = 0; k < 3; ++k)
                if (w2[k] > 1e-12 && t2n[k] != fv[0] && (nf < 2 || t2n[k] != fv[1])) behind_all = false;
              ++nfold;
            }
            if (!behind_all || nfold == 0) continue;
            h.side = 1;
            h.pen = hsum + dist;
            h.n = d * (-1.0 / dist);
          } else {  // beyond an edge, outside the face plane: rounded edge
            if (P.tri_solid[id] == 2 && P.pen_t[si] < 0 && P.behind[si]) continue;
            const double dist = norm(d);
            h.side = 0;
            h.pen = hsum - dist;
            h.n = d / dist;
          }
        } else if (lat <= 0.1 * hsum) {  // above the facet (or its edge): two-sided normal contact
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
        h.w[0] = h.wf[0] = 1 - l; h.w[1] = h.wf[1] = l;
        h.pen = hs + P.seg_r[sg] + P.clear - dist;
        h.n = d / dist;
      }
      if (h.pen <= 0) continue;
      if (id >= 0 && P.tri_solid[id]) {
        // inside a solid: the face of the contact of the previous increment while it is still hit, else a face
        // crossed during the increment (in front of it and within it at the start), else the shallowest (a node
        // just below a top face is also within the projection of the side faces, deeper; a node of another body
        // pushed into the bottom face right under a vertical edge is also slightly behind the side face)
        int tier = 2;
        double appr = 0;
        const int lm = P.last_master[si];
        if (P.pen_t[si] >= 0 && lm >= 0 && P.tri_solid[lm] &&
            std::equal(&P.fac[4 * lm], &P.fac[4 * lm] + 4, &P.fac[4 * id])) {
          tier = 0;
        } else {
          const int* tn = &P.tri[3 * id];
          const Vec3 a0 = Xs(tn[0]) - load3(&A.du[6 * tn[0]]), b0 = Xs(tn[1]) - load3(&A.du[6 * tn[1]]),
                     c0 = Xs(tn[2]) - load3(&A.du[6 * tn[2]]), x0 = x - load3(&A.du[6 * n0]);
          double wp0[3];
          plane_coords(x0, a0, b0, c0, wp0);
          const Vec3 n0 = normalized(cross(b0 - a0, c0 - a0));
          const double s0 = dot(x0 - a0, n0);
          if (s0 >= 0 && std::min({wp0[0], wp0[1], wp0[2]}) >= -1e-9) {
            tier = 1;
            appr = s0 + h.pen;  // motion into the face over the increment (signed distance s0, then -pen)
          }
        }
        // among faces crossed during the increment, the one the node moved into the most (a node on the edge of a
        // punch at the start, pushed into its bottom face, is also on the plane of its side face); otherwise the
        // shallowest
        bool better = bsol.master == INT_MIN || tier < bsol_tier;
        if (!better && tier == bsol_tier)
          better = tier == 1 && std::fabs(appr - bsol_appr) > 1e-9 * (std::fabs(appr) + std::fabs(bsol_appr)) ? appr > bsol_appr
                                                                                                         : h.pen < bsol.pen;
        if (better) { bsol = h; bsol_tier = tier; bsol_appr = appr; }
      }
      else if (id >= 0) { if (h.pen > bt.pen) bt = h; }
      else if (h.pen > bs.pen) bs = h;
    }
    if (si < ns && P.pen_t[si] < 0) {  // not in contact: remember a node behind the one-sided rigid surfaces
      if (s_rigid > 1e299 || s_rigid > hs) P.behind[si] = 0;
      else if (s_rigid < 0) P.behind[si] = 1;
    }
    if (bsol.master != INT_MIN && (bt.master == INT_MIN || bsol.pen > bt.pen)) bt = bsol;
    for (int pass = 0; pass < 2; ++pass) {
      Hit& h = pass ? bs : bt;
      if (h.master == INT_MIN) continue;
      int nm;
      const int* mn = P.stiff_nodes(h, nm);
      int sn[2];
      double sw[2];
      const int nsn = P.slave_nodes(si, h, sn, sw);
      double mm = 0, msl = 0;
      for (int i = 0; i < nm; ++i) mm += h.w[i] * A.mt[mn[i]];
      for (int i = 0; i < nsn; ++i) msl += sw[i] * A.mt[sn[i]];
      double meff = std::min(msl > 0 ? msl : 1e300, mm > 0 ? mm : 1e300);
      if (meff >= 1e299) { h = Hit(); continue; }
      // contact frequency: the stiffest node of the pair, at least a third of the squared frequency 1/(1.2 dt)
      double w2 = wfloor2, g = 1.0;
      for (int i = 0; i < nsn; ++i) if (1.1 * node_w(sn[i], true) > w2) { w2 = 1.1 * node_w(sn[i], true); g = P.node_g[sn[i]]; }
      for (int i = 0; i < nm; ++i) if (1.1 * node_w(mn[i], false) > w2) { w2 = 1.1 * node_w(mn[i], false); g = 1.0; }
      (pass ? P.cdt_s : P.cdt_t)[si] = std::max(g, 1.0) / std::sqrt(w2);
      double sflex = 0;  // slave side: a node, or the two nodes of an edge weighted by the contact point
      for (int i = 0; i < nsn; ++i) if (sw[i] > 0) sflex += sw[i] / node_k(sn[i], w2, true);
      const double ks = sflex > 0 ? 1.0 / sflex : 1e300;
      double flex = 0;
      for (int i = 0; i < nm; ++i) if (h.w[i] > 0) flex += h.w[i] / node_k(mn[i], w2, false);
      const double kc = std::min(ks, flex > 0 ? 1.0 / flex : 1e300);
      if (!(kc > 0 && kc < 1e299)) { h = Hit(); continue; }
      h.k = kc;
      (pass ? P.cw2_s : P.cw2_t)[si] = w2;
      (pass ? P.md_s : P.md_t)[si] = msl;
      (pass ? P.mm_s : P.mm_t)[si] = mm;  // zero: rigid master facet
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
      if (h.master == INT_MIN || !(h.pen - P.clear > 0)) continue;
      const int n = P.slaves[si];
      for (int k = 0; k < 3; ++k) A_.u[6 * n + k] += 0.5 * (h.pen - P.clear) * h.n[k];
      bt = Hit(); bs = Hit();
      ++nadj;
    }
    warn("general contact: %d nodes initially overclosed, moved out without strain (half the penetration each side)",
         nadj);
  }
  if (dt <= 1e-20) {  // initial pass: no forces (all hits are cleared for the first real increment)
    for (int si = 0; si < nent; ++si) { P.hit_t[si] = Hit(); P.hit_s[si] = Hit(); }
    return;
  }

  // stiffness shared by several contacts on the same node: utilisation u of each node = sum over its contacts of
  // k_pair / k_node (weighted by the contact point coordinates on the master side); a node used beyond 1 has its
  // nodal stiffness divided by u in the pair stiffness
  std::vector<double> ksum(nn, 0.0);
  // softened contact: nodes in contact as slaves this increment (the reverse pass of a two-way contact)
  std::vector<char> node_hit;
  if (P.soft) {
    node_hit.assign(nn, 0);
    for (int si = 0; si < ns; ++si)
      if (P.hit_t[si].master != INT_MIN || P.hit_s[si].master != INT_MIN) node_hit[P.slaves[si]] = 1;
  }
  for (int si = 0; si < nent; ++si) {
    for (int pass = 0; pass < 2; ++pass) {
      const Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
      if (h.master == INT_MIN) continue;
      int nm;
      const int* mn = P.stiff_nodes(h, nm);
      int sn[2];
      double sw[2];
      const int nsn = P.slave_nodes(si, h, sn, sw);
      const double w2 = (pass ? P.cw2_s : P.cw2_t)[si];
      for (int i = 0; i < nsn; ++i) ksum[sn[i]] += sw[i] * h.k / node_k(sn[i], w2, true);
      for (int i = 0; i < nm; ++i) ksum[mn[i]] += h.w[i] * h.k / node_k(mn[i], w2, false);
    }
  }

  // ---- response: incremental normal force F = min(F_prev + k dpen, k pen) >= 0. Changes of k (with dt),
  // switches between facets and contacts found already deep do not create force jumps (no energy is
  // injected); a new contact only loads the approach of the current increment, so initial or lateral
  // overclosures act as offsets that are released as the nodes separate.
  static const bool trace = std::getenv("ABQ_CONTACT_TRACE") != nullptr;
  double efric = 0, edamp = 0, dtmn = 1e30;
  int ncont = 0;
  long ndeep = 0;
#pragma omp parallel for num_threads(A.nthreads) reduction(+ : efric, edamp, ncont, ndeep) reduction(min : dtmn) \
    schedule(dynamic, 64)
  for (int si = 0; si < nent; ++si) {
    for (int pass = 0; pass < 2; ++pass) {
      Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
      double& fn = pass ? P.fn_s[si] : P.fn_t[si];
      double& pp = pass ? P.pen_s[si] : P.pen_t[si];
      Vec3& F = pass ? P.force_s[si] : P.force_t[si];
      const double fd_old = h.master == INT_MIN ? 0.0 : dot(F, h.n) - fn;  // damping force of the previous increment
      F = Vec3();
      if (h.master == INT_MIN) { pp = -1; fn = 0; continue; }
      int nm;
      const int* mn = P.stiff_nodes(h, nm);
      int sn[2];
      double sw[2];
      const int nsn = P.slave_nodes(si, h, sn, sw);
      Vec3 dus;
      for (int i = 0; i < nsn; ++i) dus += load3(&A.du[6 * sn[i]]) * sw[i];
      // smallest pair stiffness since the contact started (the penalty never stiffens during a contact), 9/11 of
      // it, reduced by the sharing of the nodes
      const double w2 = (pass ? P.cw2_s : P.cw2_t)[si];
      double& kmin = (pass ? P.kmin_s : P.kmin_t)[si];
      kmin = pp >= 0 && kmin > 0 ? std::min(kmin, h.k) : h.k;
      double flex = 0;
      for (int i = 0; i < nm; ++i)
        if (h.w[i] > 0) flex += h.w[i] * std::max(1.0, ksum[mn[i]]) / node_k(mn[i], w2, false);
      double sflex = 0;
      for (int i = 0; i < nsn; ++i)
        if (sw[i] > 0) sflex += sw[i] * std::max(1.0, ksum[sn[i]]) / node_k(sn[i], w2, true);
      const double kshared = std::min(sflex > 0 ? 1.0 / sflex : 1e300, flex > 0 ? 1.0 / flex : 1e300);
      const double k = 9.0 / 11.0 * kmin * std::min(1.0, kshared / h.k);
      // damping mass: the slave mass when the slave side sets the (shared) stiffness, the reduced mass of the
      // slave and of the master facet (weighted by the contact point) when the master side does
      const double msl = (pass ? P.md_s : P.md_t)[si], mms = (pass ? P.mm_s : P.mm_t)[si];
      const bool master_limits = flex > 0 && (sflex <= 0 || 1.0 / flex < 1.0 / sflex);
      const double mdamp = master_limits && mms > 0 ? msl * mms / (msl + mms) : msl;
      const bool rigid_master = !(mms > 0);
      if (!P.soft) dtmn = std::min(dtmn, (pass ? P.cdt_s : P.cdt_t)[si]);
      Vec3 dum;
      int nf;
      const int* mf = P.master_nodes(h, nf);
      for (int i = 0; i < nf; ++i) dum += load3(&A.du[6 * mf[i]]) * h.wf[i];
      double Fn;
      const double appr = std::max(0.0, -dot(dus - dum, h.n));  // relative approach in this increment
      double kd = k;  // stiffness of the damper
      if (P.soft) {
        // softened contact: pressure of the overclosure on the contact area of the slave entity (tributary area
        // of the node). When the reverse pass is active too (a node of the master facet in contact as a slave),
        // each pass carries 0.7 of it: total force 0.68 to 0.83 x (both passes at full area) in Abaqus/Explicit,
        // 1.0 with a single pass (sheets of one S4R element dropped on 2 x 2 plates of various sizes, per-node
        // contact forces and openings)
        double wp = 1.0;
        {
          int nf;
          const int* mf = P.master_nodes(h, nf);
          for (int i = 0; i < nf; ++i) if (node_hit[mf[i]]) wp = 0.7;
        }
        double dp;
        Fn = P.area[si] * wp * P.soft->pressure(h.pen - P.clear, &dp);
        kd = P.area[si] * wp * dp;
        double w2s = 0;
        for (int i = 0; i < nsn; ++i) w2s = std::max(w2s, node_w(sn[i], true));
        const double msl0 = (pass ? P.md_s : P.md_t)[si];
        if (kd > 0 && msl0 > 0) dtmn = std::min(dtmn, 1.0 / std::sqrt(w2s + kd / (4 * msl0)));
      } else if (pp >= 0) {
        // the force grows at most with the actual approach (a jump of the penetration after a change of
        // master facet or of contact geometry is not loaded)
        Fn = std::min(fn + k * std::min(h.pen - pp, appr), k * h.pen);
      } else {
        Fn = k * std::min(h.pen, appr);
        double hm = h.master >= 0 ? P.tri_h[h.master] : P.seg_r[-1 - h.master];
        if (!first && h.pen > 0.3 * (P.ent_h(si) + hm)) ++ndeep;
      }
      Fn = std::max(Fn, 0.0);
      fn = Fn;
      pp = h.pen;
      // normal viscous damping, 3 % of the critical damping of the contact spring with the damping mass,
      // c = 2 zeta sqrt(k m); between deformable bodies it is not bounded by the spring force: while the
      // spring is loaded and the nodes separate fast, the total normal force may become slightly tensile (as in
      // Abaqus/Explicit); no force within a released offset (spring force zero); against a rigid facet the total
      // force stays compressive
      double Ft = Fn;
      if (dt > 0) {
        const double c = 2 * 0.03 * std::sqrt(kd * mdamp);
        Ft = Fn > 0 ? Fn - c * dot(dus - dum, h.n) / dt : 0.0;
        if (rigid_master) Ft = std::max(Ft, 0.0);
        // dissipated by the damper, integrated as the work of the contact forces (trapezoidal rule)
        edamp += -0.5 * (fd_old + Ft - Fn) * dot(dus - dum, h.n);
      }
      F = h.n * Ft;
      if (trace)  // ABQ_CONTACT_TRACE: every contact force (stderr), for comparisons with nodal contact output
        std::fprintf(stderr, "CTRACE inc %d node %d pass %d k %.6e hk %.6e pen %.6e appr %.6e fn %.6e ft %.6e vn %.6e m %.6e master %d side %d n %.4f %.4f %.4f\n",
                     increment, si < ns ? P.slaves[si] : -1, pass, k, h.k, h.pen, appr, Fn, Ft,
                     dt > 0 ? dot(dus - dum, h.n) / dt : 0.0, mdamp, h.master, h.side, h.n.x, h.n.y, h.n.z);
      ++ncont;
      const bool fric = pass == 0 || P.hit_t[si].master == INT_MIN;
      if (fric && P.mu > 0) {
        Vec3 dr = dus - dum;
        Vec3 dtg = dr - h.n * dot(dr, h.n);
        // stored tangential force kept while the master stays the same facet (either triangle of a quadrilateral:
        // a node on the diagonal switches triangle without slipping) or segment
        const int lm = P.last_master[si];
        const bool same = lm == h.master || (lm >= 0 && h.master >= 0 && std::equal(&P.fac[4 * lm], &P.fac[4 * lm] + 4, &P.fac[4 * h.master]));
        Vec3 f0 = same ? P.ft[si] : Vec3();
        f0 -= h.n * dot(f0, h.n);
        Vec3 ftr = f0 - dtg * k;
        double fmax = P.mu * Fn, fm = norm(ftr);
        if (fm > fmax) {
          if (k > 0) efric += fmax * (fm - fmax) / k;
          ftr = fm > 0 ? ftr * (fmax / fm) : Vec3();
        }
        P.ft[si] = ftr;
        F += ftr;
        if (trace)
          std::fprintf(stderr, "FTRACE inc %d node %d master %d fn %.6e ft %.6e %.6e %.6e slip %.6e %.6e %.6e\n", increment,
                       si < ns ? P.slaves[si] : -1, h.master, Fn, ftr.x, ftr.y, ftr.z, dtg.x, dtg.y, dtg.z);
      }
    }
    const Hit& bt = P.hit_t[si];
    const Hit& bs = P.hit_s[si];
    if (bt.master != INT_MIN) { P.last_master[si] = bt.master; P.last_side[si] = (signed char)bt.side; P.last_n[si] = bt.n; }
    else if (bs.master != INT_MIN) { P.last_master[si] = bs.master; P.last_side[si] = 0; P.last_n[si] = bs.n; }
    else { P.last_master[si] = INT_MIN; P.last_side[si] = 0; P.ft[si] = Vec3(); }
  }
  for (int si = 0; si < nent; ++si) {  // serial scatter
    for (int pass = 0; pass < 2; ++pass) {
      const Hit& h = pass ? P.hit_s[si] : P.hit_t[si];
      if (h.master == INT_MIN) continue;
      const Vec3& F = pass ? P.force_s[si] : P.force_t[si];
      int sn[2];
      double sw[2];
      const int nsn = P.slave_nodes(si, h, sn, sw);
      for (int i = 0; i < nsn; ++i) add3(&f[6 * sn[i]], F * sw[i]);
      int nm;
      const int* mn = P.master_nodes(h, nm);
      for (int i = 0; i < nm; ++i) add3(&f[6 * mn[i]], F * (-h.wf[i]));
    }
  }
  en.friction += efric;
  P.dt_min = dtmn;
  // the damping part of the contact work is dissipation: reported with the viscous energy (as the contact
  // damping energy in ALLVD), and removed from the contact work so that the energy balance is unchanged
  en.viscous += edamp;
  en.contact -= edamp;
  ncontacts_ = ncont;
  P.d_deep_n += ndeep;
  if (P.diag < 0) P.diag = std::getenv("ABQ_CONTACT_DIAG") ? std::max(1, std::atoi(std::getenv("ABQ_CONTACT_DIAG"))) : 0;
  if (P.diag > 0 && increment % P.diag == 0)  // ABQ_CONTACT_DIAG = interval in increments
    std::fprintf(stderr, "CDIAG inc %d contacts %d deep_new %ld  Wcontact %.4e fric %.4e\n", increment, ncont, P.d_deep_n, en.contact, en.friction);
  (void)increment;
}

double GeneralContact::dt_limit() const {
  return p_ && p_->dt_min > 0 ? p_->dt_min : 1e30;
}

void GeneralContact::save(BinOut& o) const {
  o.tag("contact.v5");
  bool has = p_ != nullptr;
  o.pod(has);
  if (!has) return;
  const Impl& P = *p_;
  o.vec(P.last_master); o.vec(P.last_side); o.vec(P.ft); o.vec(P.hit_t); o.vec(P.hit_s); o.vec(P.force_t); o.vec(P.force_s);
  o.vec(P.fn_t); o.vec(P.pen_t); o.vec(P.fn_s); o.vec(P.pen_s); o.vec(P.last_n);
  o.vec(P.kmin_t); o.vec(P.kmin_s);
  o.vec(P.behind);
  o.pod(P.first_eval);
}

void GeneralContact::load(BinIn& in) {
  std::string tg = in.read_tag();
  if (tg != "contact.v5" && tg != "contact.v4" && tg != "contact.v3" && tg != "contact.v2" && tg != "contact")
    throw std::runtime_error("restart file: contact section expected");
  bool has;
  in.pod(has);
  if (!has) return;
  if (!p_) throw std::runtime_error("restart: contact state without a contact definition");
  Impl& P = *p_;
  in.vec(P.last_master); in.vec(P.last_side); in.vec(P.ft);
  if (tg == "contact.v4" || tg == "contact.v5") {
    in.vec(P.hit_t); in.vec(P.hit_s);
  } else {  // three master weights (triangle facets): bilinear quadrilateral weights from the same point
    for (auto* hv : {&P.hit_t, &P.hit_s}) {
      std::vector<HitV3> old;
      in.vec(old);
      hv->assign(old.size(), Hit());
      for (size_t i = 0; i < old.size(); ++i) {
        Hit& h = (*hv)[i];
        h.master = old[i].master; h.pen = old[i].pen; h.k = old[i].k; h.n = old[i].n; h.side = old[i].side;
        for (int k = 0; k < 3; ++k) h.w[k] = h.wf[k] = old[i].w[k];
        if (h.master >= 0 && P.fac[4 * h.master + 3] >= 0) {
          const bool half = P.fac_half[h.master];
          const double xi = half ? h.w[1] : h.w[1] + h.w[2], eta = half ? h.w[1] + h.w[2] : h.w[2];
          h.wf[0] = (1 - xi) * (1 - eta); h.wf[1] = xi * (1 - eta); h.wf[2] = xi * eta; h.wf[3] = (1 - xi) * eta;
        }
      }
    }
  }
  in.vec(P.force_t); in.vec(P.force_s);
  if (tg != "contact") {
    in.vec(P.fn_t); in.vec(P.pen_t); in.vec(P.fn_s); in.vec(P.pen_s); in.vec(P.last_n);
    if (tg != "contact.v2") { in.vec(P.kmin_t); in.vec(P.kmin_s); }
    if (tg == "contact.v5") in.vec(P.behind);
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
  P.resize_state(P.nent);  // files without the edge entities
  P.xref.clear();  // candidate search redone at the first evaluation
}

}  // namespace abq

namespace abq {

double GeneralContact::friction() const { return p_ ? p_->mu : 0.0; }
const SurfaceInteraction* GeneralContact::softened() const { return p_ ? p_->soft : nullptr; }
double GeneralContact::clearance() const { return p_ ? p_->clear : 0.0; }
double GeneralContact::mean_facet_length() const { return p_ ? p_->lavg : 1.0; }

void GeneralContact::implicit_points(std::vector<Point>& pts, double near) {
  pts.clear();
  if (!active_) return;
  prepare();
  Impl& P = *p_;
  const Analysis& A = A_;
  auto X = [&](int n) {
    return Vec3(A.X[3 * n] + A.u[6 * n] + A.du[6 * n], A.X[3 * n + 1] + A.u[6 * n + 1] + A.du[6 * n + 1],
                A.X[3 * n + 2] + A.u[6 * n + 2] + A.du[6 * n + 2]);
  };
  auto Xs = [&](int n) { return X(n) + P.shift[n]; };
  // facet normal; quadrilaterals: one normal for their two triangles (cross product of the diagonals)
  auto facet_normal = [&](int id) {
    const int* tn = &P.tri[3 * id];
    Vec3 nrm = normalized(cross(Xs(tn[1]) - Xs(tn[0]), Xs(tn[2]) - Xs(tn[0])));
    const int* qf = &P.fac[4 * id];
    if (qf[3] >= 0) {
      const Vec3 nq = normalized(cross(Xs(qf[2]) - Xs(qf[0]), Xs(qf[3]) - Xs(qf[1])));
      nrm = dot(nq, nrm) >= 0 ? nq : nq * -1.0;
    }
    return nrm;
  };
  auto Xo = [&](int n) {  // start of the increment
    return Vec3(A.X[3 * n] + A.u[6 * n], A.X[3 * n + 1] + A.u[6 * n + 1], A.X[3 * n + 2] + A.u[6 * n + 2]) + P.shift[n];
  };
  const int ns = (int)P.slaves.size();
  std::vector<Point> found(ns);
  std::vector<char> ok(ns, 0);
#pragma omp parallel for num_threads(A.nthreads) schedule(dynamic, 64)
  for (int si = 0; si < ns; ++si) {
    const int n0 = P.slaves[si];
    const Vec3 x = Xs(n0);
    const double hs = P.slave_h[si];
    Point bshell{}, bsolid{};
    bool have_shell = false, have_solid = false;
    for (int c = P.cand_off[si]; c < P.cand_off[si + 1]; ++c) {
      const int id = P.cand[c];
      if (id < 0) continue;  // beam segments: no implicit contact
      const int* tn = &P.tri[3 * id];
      const Vec3 a = Xs(tn[0]), b = Xs(tn[1]), cc = Xs(tn[2]);
      double w[3];
      bool in;
      const Vec3 q = closest_on_tri(x, a, b, cc, w, in);
      Vec3 nrm = facet_normal(id);
      const bool solid = P.tri_solid[id] != 0;
      if (solid && !P.st_off.empty()) {
        // solid faces: normal interpolated from the vertex normals (mean of the adjacent exterior faces within 30
        // degrees), continuous from one facet to the next: a slave node at a shared vertex does not switch between
        // the tilted normals of the faces around it from one iteration to the next
        double wq[3];
        plane_coords(q, a, b, cc, wq);
        Vec3 ns;
        for (int k = 0; k < 3; ++k) {
          Vec3 vn;
          for (int j = P.st_off[tn[k]]; j < P.st_off[tn[k] + 1]; ++j) {
            const Vec3 nf = facet_normal(P.st[j]);
            if (dot(nf, nrm) > 0.866) vn = vn + nf;
          }
          ns = ns + (norm(vn) > 0 ? normalized(vn) : nrm) * std::max(wq[k], 0.0);
        }
        if (norm(ns) > 0.5) nrm = normalized(ns);
      }
      // solid surfaces: only the faces turned towards the node (see compute())
      if (P.tri_solid[id] == 1 && !P.snorm.empty() && norm2(P.snorm[n0]) > 0 && dot(nrm, P.snorm[n0]) > -0.17) continue;
      const Vec3 d = x - q;
      const double s = dot(d, nrm);
      const double lat = in ? 0.0 : norm(d - nrm * s);
      const double hsum = hs + P.tri_h[id] + P.clear;
      double pen;
      Vec3 nv;
      if (lat <= 0.1 * hsum) {
        if (solid) {
          if (s < -P.tri_depth[id]) continue;
          pen = hsum - s;
          nv = nrm;
        } else {  // two-sided shell facet: the side where the node was at the start of the increment
          const Vec3 a0 = Xo(tn[0]), b0 = Xo(tn[1]), c0 = Xo(tn[2]);
          const double s0 = dot(Xo(n0) - a0, normalized(cross(b0 - a0, c0 - a0)));
          const int side = s0 >= 0 ? 1 : -1;
          pen = hsum - side * s;
          nv = nrm * side;
        }
      } else {  // rounded edge or vertex (shell facets only: on solid faces the normal would switch between the face
                // and its edge from one iteration to the next)
        if (solid) continue;
        const double dist = norm(d);
        if (dist < 1e-14) continue;
        pen = hsum - dist;
        nv = d / dist;
      }
      if (pen <= -near || (near == 0 && pen <= 0)) continue;
      Point pt{};
      pt.slave = n0;
      pt.master = id;
      pt.pen = pen;
      pt.n = nv;
      pt.area = P.area.empty() ? 0.0 : P.area[si];
      double wp[3];
      plane_coords(q, X(tn[0]), X(tn[1]), X(tn[2]), wp);
      double sw = 0;
      for (int k = 0; k < 3; ++k) sw += (wp[k] = std::max(wp[k], 0.0));
      for (int k = 0; k < 3; ++k) wp[k] /= sw;
      const int* qf = &P.fac[4 * id];
      if (qf[3] >= 0) {  // force on the four nodes of the quadrilateral, bilinear weights
        const double xi = P.fac_half[id] ? wp[1] : wp[1] + wp[2], eta = P.fac_half[id] ? wp[1] + wp[2] : wp[2];
        pt.nm = 4;
        const double wf[4] = {(1 - xi) * (1 - eta), xi * (1 - eta), xi * eta, (1 - xi) * eta};
        for (int k = 0; k < 4; ++k) { pt.mn[k] = qf[k]; pt.w[k] = wf[k]; }
      } else {
        pt.nm = 3;
        for (int k = 0; k < 3; ++k) { pt.mn[k] = tn[k]; pt.w[k] = wp[k]; }
      }
      // shell facets: the deepest; inside a solid: the face nearest to the node (shallowest penetration); a face
      // actually penetrated wins over faces only near contact (a node at the corner of a punch is also on the plane
      // of its side faces)
      if (solid) {
        const bool better = !have_solid || (pen > 0 && (bsolid.pen <= 0 || pen < bsolid.pen)) || (pen <= 0 && bsolid.pen <= 0 && pen > bsolid.pen);
        if (better) { bsolid = pt; have_solid = true; }
      } else if (!have_shell || pen > bshell.pen) { bshell = pt; have_shell = true; }
    }
    if (have_solid && (!have_shell || bsolid.pen > bshell.pen)) { found[si] = bsolid; ok[si] = 1; }
    else if (have_shell) { found[si] = bshell; ok[si] = 1; }
  }
  for (int si = 0; si < ns; ++si) if (ok[si]) pts.push_back(found[si]);
}

}  // namespace abq
