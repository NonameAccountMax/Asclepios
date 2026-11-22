// Linear dependencies between degrees of freedom: u[d] = sum_j c_j u[j].
//  - surface ties (*TIE): every slave node is paired once, at the start, with its closest point on the master
//    surface (triangles: exact closest point; quadrilaterals: closest point of the bilinear patch); its translations
//    (and its rotations when both sides carry rotational degrees of freedom) are the shape-function interpolation
//    of those of the facet nodes (explicit steps); implicit steps use a surface formulation, the master motion
//    averaged over the share of the slave faces of each slave node. ADJUST=YES moves the slave node onto that point before the analysis (no strain).
//    Slave nodes farther than the position tolerance stay free (default tolerance: 5 % of the mean diagonal of
//    the master facets, plus the half thicknesses of shell sides);
//  - JOIN connectors: the translations of the second node are those of the first one (reversed when only the
//    second one is prescribed);
//  - *COUPLING KINEMATIC and *RIGID BODY: the coupling nodes follow the rigid motion of the reference node (the
//    elements of a rigid body only carry mass); *COUPLING DISTRIBUTING: the reference node follows the weighted
//    mean motion of the coupling nodes (centroid translation and least-squares rotation; surface weights =
//    tributary areas, times (1 - r/r_max)^p for the LINEAR, QUADRATIC and CUBIC methods); *EQUATION: the first
//    term is eliminated; *MPC BEAM (rigid link), TIE, PIN: the first node follows the second one.
// A prescribed dependent dof keeps its boundary condition (the link is dropped for the step). Chains (a term that
// is itself dependent) are replaced by its terms. Implicit: the dependent dofs are eliminated (u = T u_r). Explicit: the
// forces and masses of a dependent dof are moved onto its terms with the coefficients (row-sum lumping of T^T M T),
// its velocity and acceleration interpolate those of its terms. An offset slave node with tied rotations is
// carried rigidly by the interpolated rotation (u_s = u_p + theta_p x d, linear in the rotations); its mass then
// adds the inertia m d^2 to the rotations of the master nodes.
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "core/log.hpp"
#include "solver/analysis.hpp"

namespace abq {

namespace {

// nodes of one face of an element (local indices), in a closed loop
int face_local_nodes(ElType t, int face, int nen, int* out) {
  static const int hex[6][4] = {{0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};
  static const int tet[4][3] = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {2, 3, 0}};
  if (t == ElType::C3D8R || t == ElType::C3D8 || t == ElType::C3D8I) {
    if (face < 1 || face > 6) return 0;
    for (int i = 0; i < 4; ++i) out[i] = hex[face - 1][i];
    return 4;
  }
  if (t == ElType::C3D4) {
    if (face < 1 || face > 4) return 0;
    for (int i = 0; i < 3; ++i) out[i] = tet[face - 1][i];
    return 3;
  }
  if (t == ElType::S4R || t == ElType::S4 || t == ElType::S3R || t == ElType::S3RS) {
    const int n = std::min(nen, 4);
    for (int i = 0; i < n; ++i) out[i] = i;
    return n;
  }
  return 0;
}

bool is_shell(ElType t) { return t == ElType::S4R || t == ElType::S4 || t == ElType::S3R || t == ElType::S3RS; }

struct Facet {
  int n = 0;
  int nodes[4];
  double half_t = 0;  // half thickness of a shell facet
};

// closest point of a triangle; barycentric weights
double closest_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c, double* w) {
  const Vec3 ab = b - a, ac = c - a, ap = p - a;
  const double d1 = dot(ab, ap), d2 = dot(ac, ap);
  auto set = [&](double u, double v, double x) { w[0] = u; w[1] = v; w[2] = x; };
  if (d1 <= 0 && d2 <= 0) set(1, 0, 0);
  else {
    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    const double vc = d1 * d4 - d3 * d2, vb = d5 * d2 - d1 * d6, va = d3 * d6 - d5 * d4;
    if (d3 >= 0 && d4 <= d3) set(0, 1, 0);
    else if (d6 >= 0 && d5 <= d6) set(0, 0, 1);
    else if (vc <= 0 && d1 >= 0 && d3 <= 0) { const double s = d1 / (d1 - d3); set(1 - s, s, 0); }
    else if (vb <= 0 && d2 >= 0 && d6 <= 0) { const double s = d2 / (d2 - d6); set(1 - s, 0, s); }
    else if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
      const double s = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      set(0, 1 - s, s);
    } else {
      const double den = 1.0 / (va + vb + vc);
      set(va * den, vb * den, vc * den);
    }
  }
  const Vec3 q = a * w[0] + b * w[1] + c * w[2];
  return norm(p - q);
}

// closest point of a bilinear quadrilateral (Gauss-Newton on the natural coordinates, kept in [-1, 1]^2)
double closest_on_quad(const Vec3& p, const Vec3* x, double* w) {
  double r = 0, s = 0;
  auto shape = [](double r, double s, double* N) {
    N[0] = 0.25 * (1 - r) * (1 - s); N[1] = 0.25 * (1 + r) * (1 - s);
    N[2] = 0.25 * (1 + r) * (1 + s); N[3] = 0.25 * (1 - r) * (1 + s);
  };
  for (int it = 0; it < 30; ++it) {
    double N[4];
    shape(r, s, N);
    Vec3 q = x[0] * N[0] + x[1] * N[1] + x[2] * N[2] + x[3] * N[3];
    const Vec3 xr = (x[1] - x[0]) * (0.25 * (1 - s)) + (x[2] - x[3]) * (0.25 * (1 + s));
    const Vec3 xs = (x[3] - x[0]) * (0.25 * (1 - r)) + (x[2] - x[1]) * (0.25 * (1 + r));
    const Vec3 d = q - p;
    const double g0 = dot(xr, d), g1 = dot(xs, d);
    const double h00 = dot(xr, xr), h01 = dot(xr, xs), h11 = dot(xs, xs);
    const double det = h00 * h11 - h01 * h01;
    if (!(det > 0)) break;
    const double dr = -(h11 * g0 - h01 * g1) / det, ds = -(h00 * g1 - h01 * g0) / det;
    const double r1 = std::clamp(r + dr, -1.0, 1.0), s1 = std::clamp(s + ds, -1.0, 1.0);
    const double step = std::fabs(r1 - r) + std::fabs(s1 - s);
    r = r1; s = s1;
    if (step < 1e-13) break;
  }
  shape(r, s, w);
  const Vec3 q = x[0] * w[0] + x[1] * w[1] + x[2] * w[2] + x[3] * w[3];
  return norm(p - q);
}

}  // namespace

namespace {

// closest point of the facets of a master surface: facet index (-1: none within tol), weights, distance and the
// distance along the surface from the facet (0 when the point projects inside the facet)
struct FacetHit {
  int facet = -1;
  double w[4] = {0, 0, 0, 0};
  double dist = 1e300, lateral = 0;
};
FacetHit closest_facet(const Vec3& p, const std::vector<Facet>& facets, const std::vector<Vec3>& lo,
                       const std::vector<Vec3>& hi, const std::vector<Vec3>& X, double tol, double eps, int skip_node) {
  FacetHit h;
  for (size_t k = 0; k < facets.size(); ++k) {
    if (p.x < lo[k].x - tol || p.y < lo[k].y - tol || p.z < lo[k].z - tol || p.x > hi[k].x + tol ||
        p.y > hi[k].y + tol || p.z > hi[k].z + tol)
      continue;
    const Facet& F = facets[k];
    if (skip_node >= 0 && std::find(F.nodes, F.nodes + F.n, skip_node) != F.nodes + F.n) { h.facet = -2; return h; }
    double w[4] = {0, 0, 0, 0};
    double d;
    Vec3 nrm;
    if (F.n == 3) {
      d = closest_on_triangle(p, X[F.nodes[0]], X[F.nodes[1]], X[F.nodes[2]], w);
      nrm = cross(X[F.nodes[1]] - X[F.nodes[0]], X[F.nodes[2]] - X[F.nodes[0]]);
    } else {
      const Vec3 x[4] = {X[F.nodes[0]], X[F.nodes[1]], X[F.nodes[2]], X[F.nodes[3]]};
      d = closest_on_quad(p, x, w);
      nrm = cross(x[2] - x[0], x[3] - x[1]);
    }
    if (d < h.dist - eps) {
      h.dist = d;
      h.facet = (int)k;
      std::copy(w, w + 4, h.w);
      Vec3 q;
      for (int i = 0; i < F.n; ++i) q += X[F.nodes[i]] * w[i];
      const Vec3 n = normalized(nrm), e = p - q;
      h.lateral = norm(e - n * dot(n, e));
    }
  }
  if (h.dist > tol && h.facet >= 0) h.facet = -1;
  return h;
}

// nodes of the faces of an element-based surface (shell faces: half thickness in *half_t)
std::vector<Facet> surface_facets(const Model& m, const Surface& S) {
  std::vector<Facet> out;
  for (const FaceRef& f : S.faces) {
    const Element& el = m.elems[f.elem];
    int loc[4];
    const int n = face_local_nodes(el.type, f.face, (int)el.nodes.size(), loc);
    if (n < 3) continue;
    Facet F;
    F.n = n;
    for (int i = 0; i < n; ++i) F.nodes[i] = el.nodes[loc[i]];
    if (n == 4 && F.nodes[3] == F.nodes[2]) F.n = 3;  // degenerate quadrilateral
    if (is_shell(el.type) && el.section >= 0) F.half_t = 0.5 * m.sections[el.section].thickness;
    out.push_back(F);
  }
  return out;
}

}  // namespace

void Analysis::tie_surfaces() {
  node_ties.clear();
  for (const SurfaceTie& t : m.surface_ties) {
    const Surface& S = m.surfaces.at(t.slave);
    const Surface& M = m.surfaces.at(t.master);
    if (M.node_based) { warn("*TIE %s: node-based master surface %s not supported: ignored", t.name.c_str(), t.master.c_str()); continue; }
    const std::vector<Facet> facets = surface_facets(m, M), sfaces = surface_facets(m, S);
    if (facets.empty()) { warn("*TIE %s: master surface %s has no supported facet: ignored", t.name.c_str(), t.master.c_str()); continue; }
    double diag = 0, tm = 0, ts = 0;
    std::vector<Vec3> lo(facets.size()), hi(facets.size());
    for (size_t k = 0; k < facets.size(); ++k) {
      const Facet& F = facets[k];
      tm = std::max(tm, F.half_t);
      double dk = 0;
      lo[k] = hi[k] = m.X[F.nodes[0]];
      for (int i = 0; i < F.n; ++i) {
        const Vec3& x = m.X[F.nodes[i]];
        for (int c = 0; c < 3; ++c) { lo[k][c] = std::min(lo[k][c], x[c]); hi[k][c] = std::max(hi[k][c], x[c]); }
        for (int j = i + 1; j < F.n; ++j) dk = std::max(dk, norm(m.X[F.nodes[j]] - x));
      }
      diag += dk;
    }
    diag /= facets.size();
    // slave nodes (and the half thickness of shell slave faces)
    std::set<int> snodes;
    if (S.node_based) snodes.insert(S.nodes.begin(), S.nodes.end());
    for (const Facet& F : sfaces) {
      snodes.insert(F.nodes, F.nodes + F.n);
      ts = std::max(ts, F.half_t);
    }
    const double tol = t.tolerance >= 0 ? t.tolerance : 0.05 * diag + ts + tm;
    const bool shells = ts > 0 || tm > 0;
    const double eps = 1e-12 * diag;
    int ntied = 0, nfar = 0, nmoved = 0;
    std::map<int, int> tie_of;  // slave node -> index in node_ties
    // point formulation: closest point of each slave node
    for (int sn : snodes) {
      const FacetHit h = closest_facet(m.X[sn], facets, lo, hi, m.X, tol, eps, sn);
      if (h.facet == -2) continue;  // node of the master surface itself
      if (h.facet < 0) { ++nfar; continue; }
      const Facet& F = facets[h.facet];
      NodeTie nt;
      nt.slave = sn;
      nt.rotations = t.rotations;
      Vec3 q;
      for (int i = 0; i < F.n; ++i) { nt.masters.push_back(F.nodes[i]); nt.w_pt.push_back(h.w[i]); q += m.X[F.nodes[i]] * h.w[i]; }
      if (t.adjust && !shells && h.dist > eps) { m.X[sn] = q; ++nmoved; }
      nt.off_pt = m.X[sn] - q;
      tie_of[sn] = (int)node_ties.size();
      node_ties.push_back(nt);
      ++ntied;
    }
    // surface formulation (implicit steps): each slave node follows the master surface averaged over its share of
    // the slave faces with dual weights psi_s, c_m = int(psi_s N_m) / int(psi_s) over the part of the slave faces
    // that projects inside the master facets (4 x 4 cells of 2 x 2 Gauss points per slave face)
    std::map<int, double> den;
    std::map<int, std::map<int, double>> num;
    const double g = 0.5 / std::sqrt(3.0);
    for (const Facet& F : sfaces) {
      Vec3 xs[4];
      for (int i = 0; i < 4; ++i) xs[i] = m.X[F.nodes[F.n == 3 && i == 3 ? 2 : i]];
      // dual shape functions of the slave face, psi_a = A_ab N_b with int(psi_a N_b) = delta_ab int(N_b): the
      // averaged constraint stays exact for a linear field (the slave nodes of matching meshes sit on the master)
      double Me[4][4] = {}, De[4] = {};
      for (int gi = 0; gi < 3; ++gi)
        for (int gj = 0; gj < 3; ++gj) {
          static const double gp[3] = {-0.7745966692414834, 0.0, 0.7745966692414834}, gw[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
          const double r = gp[gi], s = gp[gj];
          const double N[4] = {0.25 * (1 - r) * (1 - s), 0.25 * (1 + r) * (1 - s), 0.25 * (1 + r) * (1 + s), 0.25 * (1 - r) * (1 + s)};
          const Vec3 xr = (xs[1] - xs[0]) * (0.25 * (1 - s)) + (xs[2] - xs[3]) * (0.25 * (1 + s));
          const Vec3 xq = (xs[3] - xs[0]) * (0.25 * (1 - r)) + (xs[2] - xs[1]) * (0.25 * (1 + r));
          const double wgt = norm(cross(xr, xq)) * gw[gi] * gw[gj];
          for (int a = 0; a < 4; ++a) { De[a] += N[a] * wgt; for (int b = 0; b < 4; ++b) Me[a][b] += N[a] * N[b] * wgt; }
        }
      double A[4][4] = {};
      const int nd = F.n == 3 ? 3 : 4;  // triangles: repeated fourth node merged into the third
      if (nd == 3) {
        for (int a = 0; a < 4; ++a) { Me[a][2] += Me[a][3]; Me[a][3] = 0; }
        for (int b = 0; b < 4; ++b) { Me[2][b] += Me[3][b]; Me[3][b] = 0; }
        De[2] += De[3]; De[3] = 0;
      }
      {  // A = diag(De) Me^-1 (Gauss-Jordan on the nd x nd block)
        double Mi[4][8] = {};
        for (int a = 0; a < nd; ++a) { for (int b = 0; b < nd; ++b) Mi[a][b] = Me[a][b]; Mi[a][nd + a] = 1; }
        for (int c = 0; c < nd; ++c) {
          int piv = c;
          for (int a = c + 1; a < nd; ++a) if (std::fabs(Mi[a][c]) > std::fabs(Mi[piv][c])) piv = a;
          for (int b = 0; b < 2 * nd; ++b) std::swap(Mi[c][b], Mi[piv][b]);
          const double dd = Mi[c][c];
          if (!(std::fabs(dd) > 0)) break;
          for (int b = 0; b < 2 * nd; ++b) Mi[c][b] /= dd;
          for (int a = 0; a < nd; ++a)
            if (a != c) { const double f = Mi[a][c]; for (int b = 0; b < 2 * nd; ++b) Mi[a][b] -= f * Mi[c][b]; }
        }
        for (int a = 0; a < nd; ++a) for (int b = 0; b < nd; ++b) A[a][b] = De[a] * Mi[a][nd + b];
      }
      const int nc = 4;
      for (int ci = 0; ci < nc; ++ci)
        for (int cj = 0; cj < nc; ++cj)
          for (int gi = 0; gi < 2; ++gi)
            for (int gj = 0; gj < 2; ++gj) {
              const double r = -1 + (2.0 * ci + 1 + (gi ? 2 : -2) * g) / nc, s = -1 + (2.0 * cj + 1 + (gj ? 2 : -2) * g) / nc;
              const double N[4] = {0.25 * (1 - r) * (1 - s), 0.25 * (1 + r) * (1 - s), 0.25 * (1 + r) * (1 + s), 0.25 * (1 - r) * (1 + s)};
              const Vec3 xr = (xs[1] - xs[0]) * (0.25 * (1 - s)) + (xs[2] - xs[3]) * (0.25 * (1 + s));
              const Vec3 xq = (xs[3] - xs[0]) * (0.25 * (1 - r)) + (xs[2] - xs[1]) * (0.25 * (1 + r));
              const double wgt = norm(cross(xr, xq)) / (nc * nc);
              if (!(wgt > 0)) continue;
              const Vec3 p = xs[0] * N[0] + xs[1] * N[1] + xs[2] * N[2] + xs[3] * N[3];
              const FacetHit h = closest_facet(p, facets, lo, hi, m.X, tol, eps, -1);
              if (h.facet < 0 || h.lateral > 1e-6 * diag) continue;
              const Facet& Mf = facets[h.facet];
              double Ns[4] = {N[0], N[1], N[2], N[3]};
              if (nd == 3) { Ns[2] += Ns[3]; Ns[3] = 0; }
              for (int a = 0; a < nd; ++a) {
                const int sa = F.nodes[a];
                if (!tie_of.count(sa)) continue;
                double psi = 0;
                for (int b = 0; b < nd; ++b) psi += A[a][b] * Ns[b];
                den[sa] += psi * (Ns[0] + Ns[1] + Ns[2] + Ns[3]) * wgt;
                for (int b = 0; b < Mf.n; ++b) if (h.w[b] != 0) num[sa][Mf.nodes[b]] += psi * h.w[b] * wgt;
              }
            }
    }
    for (auto& [sn, d] : den) {
      if (!(d > 0)) continue;
      NodeTie& nt = node_ties[tie_of[sn]];
      nt.w_avg.assign(nt.masters.size(), 0.0);
      for (auto [mn, v] : num[sn]) {
        auto it = std::find(nt.masters.begin(), nt.masters.end(), mn);
        if (it == nt.masters.end()) { nt.masters.push_back(mn); nt.w_pt.push_back(0.0); nt.w_avg.push_back(v / d); }
        else nt.w_avg[it - nt.masters.begin()] = v / d;
      }
      Vec3 q;
      for (size_t i = 0; i < nt.masters.size(); ++i) q += m.X[nt.masters[i]] * nt.w_avg[i];
      nt.off_avg = m.X[sn] - q;
    }
    if (t.adjust && shells) warn_once_flag("TIEADJSHELL", "*TIE ADJUST=YES with shell surfaces: slave nodes not moved");
    info("  tie %s: %d slave nodes tied (%d moved onto the master surface), %d beyond the position tolerance %.4g",
         t.name.c_str(), ntied, nmoved, nfar, tol);
  }
}

// *COUPLING, *RIGID BODY, *EQUATION, *MPC: dependencies that do not change between steps (setup, after the
// element groups)
void Analysis::build_constraints() {
  constraint_links.clear();
  auto push = [&](int dof, std::vector<LinkTerm> terms) { constraint_links.push_back({dof, std::move(terms)}); };
  // rigid offset: translation k of a node at offset d from node r: u_r,k + (theta_r x d)_k
  auto rigid_terms = [&](int r, const Vec3& d, int k) {
    const double C[3][3] = {{0, d.z, -d.y}, {-d.z, 0, d.x}, {d.y, -d.x, 0}};
    std::vector<LinkTerm> t{{6 * r + k, 1.0, 1.0}};
    for (int j = 0; j < 3; ++j) if (C[k][j] != 0) t.push_back({6 * r + 3 + j, C[k][j], C[k][j] * C[k][j]});
    return t;
  };
  for (const NodeCoupling& cp : m.couplings) {
    std::vector<int> nodes = cp.nodes;
    std::map<int, double> area;  // tributary area of the surface nodes
    if (!cp.surface.empty()) {
      const Surface& S = m.surfaces.at(cp.surface);
      if (S.node_based) nodes.insert(nodes.end(), S.nodes.begin(), S.nodes.end());
      for (const Facet& F : surface_facets(m, S)) {
        Vec3 xs[4];
        for (int i = 0; i < 4; ++i) xs[i] = m.X[F.nodes[F.n == 3 && i == 3 ? 2 : i]];
        for (int gi = 0; gi < 2; ++gi)
          for (int gj = 0; gj < 2; ++gj) {
            const double r = gi ? 0.5773502691896258 : -0.5773502691896258, q = gj ? 0.5773502691896258 : -0.5773502691896258;
            const double N[4] = {0.25 * (1 - r) * (1 - q), 0.25 * (1 + r) * (1 - q), 0.25 * (1 + r) * (1 + q), 0.25 * (1 - r) * (1 + q)};
            const Vec3 xr = (xs[1] - xs[0]) * (0.25 * (1 - q)) + (xs[2] - xs[3]) * (0.25 * (1 + q));
            const Vec3 xq = (xs[3] - xs[0]) * (0.25 * (1 - r)) + (xs[2] - xs[1]) * (0.25 * (1 + r));
            const double w = norm(cross(xr, xq));
            for (int a = 0; a < 4; ++a) area[F.nodes[F.n == 3 && a == 3 ? 2 : a]] += N[a] * w;
          }
        nodes.insert(nodes.end(), F.nodes, F.nodes + F.n);
      }
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    nodes.erase(std::remove(nodes.begin(), nodes.end(), cp.ref), nodes.end());
    if (nodes.empty()) { warn("coupling %s: no coupling node", cp.name.c_str()); continue; }
    const Vec3 xr = m.X[cp.ref];
    if (cp.kind == NodeCoupling::KINEMATIC) {
      // every coupling node follows the rigid motion of the reference node
      for (int n : nodes) {
        const Vec3 d = m.X[n] - xr;
        for (int k = 0; k < 3; ++k) if (cp.dof[k]) push(6 * n + k, rigid_terms(cp.ref, d, k));
        // rotations too, also on nodes without rotational stiffness (solids): reported as the reference rotation
        for (int k = 3; k < 6; ++k) if (cp.dof[k]) push(6 * n + k, {{6 * cp.ref + k, 1.0, 1.0}});
      }
    } else {
      // the reference node follows the weighted mean motion of the coupling nodes: translation of the weighted
      // centroid plus the least-squares rotation, theta = J^-1 sum w_i r_i x u_i, J = sum w_i (r_i^2 I - r_i r_i^T)
      std::vector<double> w(nodes.size(), 1.0);
      double rmax = 0;
      Vec3 xc;
      double W = 0;
      for (size_t i = 0; i < nodes.size(); ++i) {
        if (!area.empty()) w[i] = area.count(nodes[i]) ? area[nodes[i]] : 0.0;
        W += w[i];
        xc += m.X[nodes[i]] * w[i];
      }
      if (!(W > 0)) { warn("coupling %s: zero weights", cp.name.c_str()); continue; }
      xc = xc / W;
      for (int n : nodes) rmax = std::max(rmax, norm(m.X[n] - xc));
      if (cp.weighting > 0)
        for (size_t i = 0; i < nodes.size(); ++i) {
          const double rho = rmax > 0 ? norm(m.X[nodes[i]] - xc) / rmax : 0.0;
          w[i] *= std::pow(1 - rho, cp.weighting) + 1e-12;
        }
      W = 0;
      xc = Vec3();
      for (size_t i = 0; i < nodes.size(); ++i) { W += w[i]; xc += m.X[nodes[i]] * w[i]; }
      xc = xc / W;
      double J[3][3] = {};
      for (size_t i = 0; i < nodes.size(); ++i) {
        const Vec3 r = m.X[nodes[i]] - xc;
        const double r2 = dot(r, r);
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) J[a][b] += w[i] * ((a == b ? r2 : 0.0) - r[a] * r[b]);
      }
      double tr = J[0][0] + J[1][1] + J[2][2];
      for (int a = 0; a < 3; ++a) J[a][a] += 1e-12 * tr;  // collinear coupling nodes: no rotation about the line
      Mat3 Jm;
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) Jm.m[a][b] = J[a][b];
      const Mat3 Ji = Jm.inverse();
      const Vec3 e = xr - xc;
      const double E[3][3] = {{0, -e.z, e.y}, {e.z, 0, -e.x}, {-e.y, e.x, 0}};  // E t = e x t
      std::vector<LinkTerm> terms[6];
      for (size_t i = 0; i < nodes.size(); ++i) {
        const Vec3 r = m.X[nodes[i]] - xc;
        const double R[3][3] = {{0, -r.z, r.y}, {r.z, 0, -r.x}, {-r.y, r.x, 0}};  // R u = r x u
        double A[3][3];  // theta = sum A_i u_i
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) {
            A[a][b] = 0;
            for (int c = 0; c < 3; ++c) A[a][b] += Ji.m[a][c] * R[c][b] * w[i];
          }
        for (int l = 0; l < 3; ++l) {
          for (int k = 0; k < 3; ++k) {  // u_ref = mean u + theta x e = mean u - e x theta
            double c = (k == l ? w[i] / W : 0.0);
            for (int j = 0; j < 3; ++j) c -= E[k][j] * A[j][l];
            if (c != 0) terms[k].push_back({6 * nodes[i] + l, c, std::fabs(c)});
          }
          for (int k = 0; k < 3; ++k)
            if (A[k][l] != 0) terms[3 + k].push_back({6 * nodes[i] + l, A[k][l], std::fabs(A[k][l])});
        }
      }
      for (int k = 0; k < 6; ++k) if (cp.dof[k] && !terms[k].empty()) push(6 * cp.ref + k, terms[k]);
    }
    for (int e : cp.rigid_elems) elem_rigid[e] = 1;
  }
  for (const LinearConstraint& eq : m.equations) {
    if (eq.node.size() < 2 || eq.coef[0] == 0) continue;
    std::vector<LinkTerm> t;
    for (size_t i = 1; i < eq.node.size(); ++i) {
      const double c = -eq.coef[i] / eq.coef[0];
      t.push_back({6 * eq.node[i] + eq.dof[i], c, std::fabs(c)});
    }
    push(6 * eq.node[0] + eq.dof[0], t);
  }
  for (const NodePairMpc& mp : m.mpcs) {
    const Vec3 d = m.X[mp.a] - m.X[mp.b];
    for (int k = 0; k < 3; ++k)
      push(6 * mp.a + k, mp.kind == NodePairMpc::BEAM ? rigid_terms(mp.b, d, k) : std::vector<LinkTerm>{{6 * mp.b + k, 1.0, 1.0}});
    if (mp.kind != NodePairMpc::PIN && node_rot[mp.a] && node_rot[mp.b])
      for (int k = 3; k < 6; ++k) push(6 * mp.a + k, {{6 * mp.b + k, 1.0, 1.0}});
  }
}

void Analysis::warn_once_flag(const std::string& key, const std::string& msg) {
  if (warned_keys.insert(key).second) warn("%s", msg.c_str());
}

void Analysis::activate_links(bool implicit) {
  links.clear();
  link_terms.clear();
  link_of.assign(6 * nn, -1);
  link_term_dofs.clear();
  std::map<int, std::vector<LinkTerm>> dep;  // dependent dof -> terms
  auto add = [&](int d, std::vector<LinkTerm> terms) {
    if (dep.count(d)) { warn_once_flag("LINKTWICE", "a degree of freedom is constrained twice (ties/JOIN): first constraint kept"); return; }
    dep[d] = std::move(terms);
  };
  for (const NodeTie& nt : node_ties) {
    bool rot = nt.rotations && node_rot[nt.slave] && !node_aniso[nt.slave];
    for (int mn : nt.masters) rot = rot && node_rot[mn] && !node_aniso[mn];
    // implicit steps use the surface formulation (when the slave faces gave one), explicit steps the point one
    const bool avg = implicit && !nt.w_avg.empty();
    const std::vector<double>& wv = avg ? nt.w_avg : nt.w_pt;
    // with tied rotations, an offset slave node is carried rigidly: u_s = u_p + theta_p x d (d = X_s - X_p),
    // C[k][j] theta_j = (theta x d)_k
    const Vec3& d = avg ? nt.off_avg : nt.off_pt;
    const double C[3][3] = {{0, d.z, -d.y}, {-d.z, 0, d.x}, {d.y, -d.x, 0}};
    for (int k = 0; k < (rot ? 6 : 3); ++k) {
      const int dof = 6 * nt.slave + k;
      if (fixed[dof]) continue;
      std::vector<LinkTerm> terms;
      for (size_t i = 0; i < nt.masters.size(); ++i) {
        const double w = wv[i];
        if (w == 0) continue;
        terms.push_back({6 * nt.masters[i] + k, w, w});
        if (rot && k < 3)
          for (int j = 0; j < 3; ++j)
            if (C[k][j] != 0) terms.push_back({6 * nt.masters[i] + 3 + j, w * C[k][j], w * C[k][j] * C[k][j]});
      }
      add(dof, std::move(terms));
    }
  }
  for (const auto& [dof, terms] : constraint_links)
    if (!fixed[dof]) add(dof, terms);
  // boundary conditions in the local axes of transformed nodes (zero values): the motion is kept in the plane
  // (one condition), along the free axis (two) or fixed (three)
  {
    std::map<std::pair<int, int>, std::vector<int>> by;  // (node, translations 0 / rotations 3) -> local axes held
    for (const LocalBC& b : local_bcs) {
      if (b.value != 0) warn_once_flag("LOCALBCVAL", "non-zero boundary condition in transformed axes: held at zero");
      by[{b.node, b.dof < 3 ? 0 : 3}].push_back(b.dof % 3);
    }
    for (auto& [key, ax] : by) {
      std::sort(ax.begin(), ax.end());
      ax.erase(std::unique(ax.begin(), ax.end()), ax.end());
      const Mat3& R = m.transforms.at(key.first);
      const int d0 = 6 * key.first + key.second;
      auto argmax = [](const Vec3& v) { int k = 0; for (int i = 1; i < 3; ++i) if (std::fabs(v[i]) > std::fabs(v[k])) k = i; return k; };
      if (ax.size() == 1) {  // e.u = 0: the largest component depends on the others
        const Vec3 e = R.row(ax[0]);
        const int k = argmax(e);
        std::vector<LinkTerm> t;
        for (int l = 0; l < 3; ++l) if (l != k && e[l] != 0) t.push_back({d0 + l, -e[l] / e[k], std::fabs(e[l] / e[k])});
        if (t.empty()) fixed[d0 + k] = 1; else add(d0 + k, t);
      } else if (ax.size() == 2) {  // u = s f along the free axis f
        const Vec3 f = R.row(3 - ax[0] - ax[1]);
        const int k = argmax(f);
        for (int l = 0; l < 3; ++l)
          if (l != k) { if (f[l] == 0) fixed[d0 + l] = 1; else add(d0 + l, {{d0 + k, f[l] / f[k], std::fabs(f[l] / f[k])}}); }
      } else {
        for (int l = 0; l < 3; ++l) fixed[d0 + l] = 1;
      }
    }
  }
  for (auto [na, nb] : joins)
    for (int k = 0; k < 3; ++k) {
      const int ia = 6 * na + k, ib = 6 * nb + k;
      if (fixed[ib] && fixed[ia]) continue;
      if (fixed[ib]) add(ia, {{ib, 1.0, 1.0}});
      else add(ib, {{ia, 1.0, 1.0}});
    }
  // replace the terms that are themselves dependent (chains), drop self references
  for (int pass = 0; pass < 8; ++pass) {
    bool changed = false;
    for (auto& [d, terms] : dep) {
      std::vector<LinkTerm> out;
      for (const LinkTerm& t : terms) {
        auto it = dep.find(t.dof);
        if (it != dep.end() && t.dof != d) {
          for (const LinkTerm& u : it->second) out.push_back({u.dof, t.c * u.c, t.mc * u.mc});
          changed = true;
        } else out.push_back(t);
      }
      terms = std::move(out);
    }
    if (!changed) break;
  }
  std::set<int> tset;
  for (auto& [d, terms] : dep) {
    std::map<int, std::pair<double, double>> merged;
    for (const LinkTerm& t : terms)
      if (t.dof != d && !dep.count(t.dof)) { merged[t.dof].first += t.c; merged[t.dof].second += t.mc; }
    if (merged.empty()) continue;
    DofLink L{d, (int)link_terms.size(), 0};
    for (auto& [j, cm] : merged) {
      link_terms.push_back({j, cm.first, cm.second});
      ++L.count;
      tset.insert(j);
    }
    link_of[d] = (int)links.size();
    links.push_back(L);
  }
  link_term_dofs.assign(tset.begin(), tset.end());
  for (const DofLink& L : links)
    if (L.count > 64) warn_once_flag("LINKTERMS", "a tied degree of freedom depends on more than 64 others: truncated in implicit steps");
  if (!links.empty()) res.assign(6 * nn, 0.0);
}

void Analysis::gather_links(std::vector<double>& f) const {
  for (const DofLink& L : links) {
    const double fs = f[L.dof];
    if (fs == 0) continue;
    for (int t = L.first; t < L.first + L.count; ++t) f[link_terms[t].dof] += link_terms[t].c * fs;
    f[L.dof] = 0;
  }
}

void Analysis::spread_links(std::vector<double>& x) const {
  for (const DofLink& L : links) {
    double s = 0;
    for (int t = L.first; t < L.first + L.count; ++t) s += link_terms[t].c * x[link_terms[t].dof];
    x[L.dof] = s;
  }
}

void Analysis::link_accelerations(const std::vector<double>& r) {
  // residual forces and masses of the term dofs, with the shares of their dependent dofs
  link_r.assign(r.begin(), r.end());
  link_m.resize(6 * nn);
  for (int n = 0; n < nn; ++n)
    for (int k = 0; k < 6; ++k) link_m[6 * n + k] = k < 3 ? mt[n] : mr[n];
  for (const DofLink& L : links) {
    const double rs = link_r[L.dof], ms = link_m[L.dof];
    for (int t = L.first; t < L.first + L.count; ++t) {
      const LinkTerm& lt = link_terms[t];
      link_r[lt.dof] += lt.c * rs;
      link_m[lt.dof] += lt.mc * ms;  // interpolation weight (translations), offset inertia m d^2 (rotations)
    }
  }
  for (int j : link_term_dofs) {
    if (fixed[j]) rf[j] = -link_r[j];
    else a[j] = link_m[j] > 0 ? link_r[j] / link_m[j] : 0.0;
  }
  for (const DofLink& L : links) {
    double s = 0;
    for (int t = L.first; t < L.first + L.count; ++t) s += link_terms[t].c * a[link_terms[t].dof];
    a[L.dof] = s;
    rf[L.dof] = 0;
  }
}

}  // namespace abq
