// Contact pairs for explicit dynamics (*CONTACT PAIR), kinematic hard contact:
//  - at the end of each increment the next one is predicted without contact (current velocities and accelerations);
//    the slave nodes that penetrate the master surface in the predicted configuration are moved back onto it by
//    corrections of the current accelerations that conserve momentum (as Abaqus/Explicit reports them: the
//    correction belongs to the accelerations, hence to the kinetic energy, of the current time): the
//    resisting force of a slave node is the one that would have stopped it exactly on the master surface in the
//    increment (its mass times the penetration over dt^2), distributed to the master nodes with the bilinear weights
//    of the contact point, together with the slave mass; the master nodes get the correction of these forces over
//    their mass plus the distributed slave masses, and each slave node the correction that puts it on the corrected
//    master surface (normal direction only: no friction);
//  - balanced master-slave weighting between two deformable surfaces (default 0.5, WEIGHT = weight of the first
//    surface as master): the corrections of both passes (each surface as slave) are averaged, then a second
//    correction of both passes, not weighted, removes the remaining penetration; a rigid or massless surface is
//    always the master;
//  - contact surfaces offset by the half thickness of the shells on the side of the surface (SPOS / SNEG), the
//    exterior faces of solids without offset;
//  - the corrections are returned as nodal contact forces (mass times correction; on prescribed degrees of freedom
//    the force taken by the reaction), whose work is the contact work (plastic impact: energy is absorbed).
#include "contact/pairs.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>

#include "contact/geom.hpp"
#include "core/log.hpp"
#include "solver/analysis.hpp"

namespace abq {

namespace {

struct Side {
  std::vector<std::array<int, 4>> fac;  // facet nodes (4th = -1: triangle), oriented so that s n is the surface side
  std::vector<double> fh;               // offset of the contact surface from the facet (shell half thickness)
  std::vector<signed char> fs;          // side of the surface: +1 positive normal, -1 negative
  std::vector<int> nodes;               // nodes of the surface (slaves when this side is the slave)
  std::vector<double> nh;               // offset at each node (smallest of its facets)
  std::vector<int> nf_off, nf;          // facets around each node of the surface
  bool rigid = false;                   // rigid or massless surface: always the master
};

struct Pair {
  Side a, b;  // first and second surface
  double f = 0.5;  // weight of the first surface as master
  double weight = -1;  // WEIGHT given (< 0: default)
  bool resolved = false;  // masters chosen (needs the nodal masses)
  std::string name;
};

// contact found in a pass: slave node, master facet nodes and bilinear weights, normal pushing the slave out,
// penetration
struct Hit {
  int s, nm;
  int mn[4];
  double w[4];
  Vec3 n;
  double d;
};

}  // namespace

struct ContactPairs::Impl {
  std::vector<Pair> pairs;
};

ContactPairs::ContactPairs(Analysis& A) : A_(A) {}
ContactPairs::~ContactPairs() { delete p_; }

namespace {

Side build_side(const Model& m, const Surface& S) {
  Side sd;
  std::unordered_map<int, int> idx;
  int nrigid = 0;
  for (const FaceRef& f : S.faces) {
    const Element& el = m.elems[f.elem];
    int loc[4];
    const int k = face_local_nodes(el.type, f.face, (int)el.nodes.size(), loc);
    if (k < 3) continue;
    std::array<int, 4> q = {el.nodes[loc[0]], el.nodes[loc[1]], el.nodes[loc[2]], k == 4 ? el.nodes[loc[3]] : -1};
    if (q[3] >= 0 && q[3] == q[2]) q[3] = -1;
    double h = 0;
    signed char s = 1;
    const bool shell = is_shell_type(el.type);
    if (shell || is_rigid_type(el.type) || is_surface_type(el.type)) {
      if (f.face == 0) {
        static bool once = false;
        if (!once) warn("contact pair on a shell surface without side (SPOS / SNEG): positive side used");
        once = true;
      }
      s = f.face < 0 ? -1 : 1;
      if (shell && el.section >= 0) h = 0.5 * m.sections[el.section].thickness;
    } else {  // solid face: oriented outward (normal away from the element centroid)
      Vec3 c;
      for (int n : el.nodes) c += m.X[n] * (1.0 / el.nodes.size());
      const Vec3 a = m.X[q[0]], b = m.X[q[1]], cc = m.X[q[2]];
      const Vec3 nr = q[3] >= 0 ? cross(cc - a, m.X[q[3]] - b) : cross(b - a, cc - a);
      const Vec3 fc = q[3] >= 0 ? (a + b + cc + m.X[q[3]]) * 0.25 : (a + b + cc) * (1.0 / 3);
      if (dot(nr, fc - c) < 0) std::swap(q[1], q[q[3] >= 0 ? 3 : 2]);
    }
    if (is_rigid_type(el.type)) ++nrigid;
    sd.fac.push_back(q);
    sd.fh.push_back(h);
    sd.fs.push_back(s);
  }
  std::vector<std::vector<int>> around;
  for (size_t t = 0; t < sd.fac.size(); ++t)
    for (int j = 0; j < 4 && sd.fac[t][j] >= 0; ++j) {
      const int n = sd.fac[t][j];
      auto it = idx.find(n);
      if (it == idx.end()) {
        it = idx.emplace(n, (int)sd.nodes.size()).first;
        sd.nodes.push_back(n);
        sd.nh.push_back(sd.fh[t]);
        around.emplace_back();
      }
      sd.nh[it->second] = std::min(sd.nh[it->second], sd.fh[t]);
      around[it->second].push_back((int)t);
    }
  sd.nf_off.assign(1, 0);
  for (auto& v : around) { sd.nf.insert(sd.nf.end(), v.begin(), v.end()); sd.nf_off.push_back((int)sd.nf.size()); }
  sd.rigid = nrigid > 0 && nrigid == (int)sd.fac.size();
  return sd;
}

}  // namespace

bool ContactPairs::define(const Step& st) {
  delete p_;
  p_ = new Impl();
  const Model& m = A_.m;
  for (const auto& cp : st.contact_pairs) {
    Pair pr;
    pr.a = build_side(m, m.surfaces.at(cp.first));
    pr.b = build_side(m, m.surfaces.at(cp.second));
    if (pr.a.fac.empty() || pr.b.fac.empty()) {
      warn("contact pair %s - %s: a surface without element faces, ignored", cp.first.c_str(), cp.second.c_str());
      continue;
    }
    pr.weight = cp.weight;
    pr.name = cp.first + " - " + cp.second;
    if (cp.interaction >= 0) {
      const SurfaceInteraction& si = m.interactions[cp.interaction];
      if (si.friction > 0) warn("contact pair %s - %s: friction not supported in kinematic contact pairs (frictionless)", cp.first.c_str(), cp.second.c_str());
      if (si.law != SurfaceInteraction::HARD) warn("contact pair %s - %s: softened kinematic contact not supported (hard contact)", cp.first.c_str(), cp.second.c_str());
    }
    info("  contact pair %s - %s: %zu / %zu facets", cp.first.c_str(), cp.second.c_str(), pr.a.fac.size(), pr.b.fac.size());
    p_->pairs.push_back(std::move(pr));
  }
  active_ = !p_->pairs.empty();
  return active_;
}

namespace {

// contacts of the slave nodes of side sl against the facets of side ms in the configuration X + u + du
void find_hits(const Analysis& A, const Side& sl, const Side& ms, const std::vector<double>& du, std::vector<Hit>& out) {
  out.clear();
  auto X = [&](int n) {
    return Vec3(A.X[3 * n] + A.u[6 * n] + du[6 * n], A.X[3 * n + 1] + A.u[6 * n + 1] + du[6 * n + 1],
                A.X[3 * n + 2] + A.u[6 * n + 2] + du[6 * n + 2]);
  };
  const int nfac = (int)ms.fac.size();
  std::vector<Vec3> fn(nfac);
  double lsum = 0, hmax = 0;
  for (int t = 0; t < nfac; ++t) {
    const auto& q = ms.fac[t];
    const Vec3 a = X(q[0]), b = X(q[1]), c = X(q[2]);
    fn[t] = normalized(q[3] >= 0 ? cross(c - a, X(q[3]) - b) : cross(b - a, c - a));
    lsum += norm(b - a);
    hmax = std::max(hmax, ms.fh[t]);
  }
  double hsl = 0;
  for (double h : sl.nh) hsl = std::max(hsl, h);
  const double lavg = nfac ? lsum / nfac : 1.0;
  const double infl = hmax + hsl + 0.25 * lavg;
  const double cs = 2 * lavg + 2 * infl;
  // bucket grid of the master facets (bounding boxes inflated by the offsets and a margin)
  auto key = [](int64_t i, int64_t j, int64_t k) { return ((i + (1 << 20)) << 42) | ((j + (1 << 20)) << 21) | (k + (1 << 20)); };
  std::vector<std::pair<int64_t, int>> cells;
  for (int t = 0; t < nfac; ++t) {
    Vec3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
    for (int j = 0; j < 4 && ms.fac[t][j] >= 0; ++j) {
      const Vec3 x = X(ms.fac[t][j]);
      lo = Vec3(std::min(lo.x, x.x), std::min(lo.y, x.y), std::min(lo.z, x.z));
      hi = Vec3(std::max(hi.x, x.x), std::max(hi.y, x.y), std::max(hi.z, x.z));
    }
    const int64_t i0 = (int64_t)std::floor((lo.x - infl) / cs), i1 = (int64_t)std::floor((hi.x + infl) / cs);
    const int64_t j0 = (int64_t)std::floor((lo.y - infl) / cs), j1 = (int64_t)std::floor((hi.y + infl) / cs);
    const int64_t k0 = (int64_t)std::floor((lo.z - infl) / cs), k1 = (int64_t)std::floor((hi.z + infl) / cs);
    if (double(i1 - i0 + 1) * double(j1 - j0 + 1) * double(k1 - k0 + 1) > 1e6) continue;  // diverging facet
    for (int64_t i = i0; i <= i1; ++i)
      for (int64_t j = j0; j <= j1; ++j)
        for (int64_t k = k0; k <= k1; ++k) cells.push_back({key(i, j, k), t});
  }
  std::sort(cells.begin(), cells.end());
  for (size_t si = 0; si < sl.nodes.size(); ++si) {
    const int s = sl.nodes[si];
    if (!(A.mt[s] > 0)) continue;
    // outward normal of the slave surface at the node (mean of its facets) and its contact point
    Vec3 ns;
    for (int j = sl.nf_off[si]; j < sl.nf_off[si + 1]; ++j) {
      const auto& q = sl.fac[sl.nf[j]];
      const Vec3 a = X(q[0]), b = X(q[1]), c = X(q[2]);
      ns += normalized(q[3] >= 0 ? cross(c - a, X(q[3]) - b) : cross(b - a, c - a)) * (double)sl.fs[sl.nf[j]];
    }
    ns = normalized(ns);
    const Vec3 p = X(s) + ns * sl.nh[si];
    const int64_t kk = key((int64_t)std::floor(p.x / cs), (int64_t)std::floor(p.y / cs), (int64_t)std::floor(p.z / cs));
    Hit best;
    best.s = -1;
    best.d = 1e300;
    for (auto it = std::lower_bound(cells.begin(), cells.end(), std::make_pair(kk, INT_MIN)); it != cells.end() && it->first == kk; ++it) {
      const int t = it->second;
      const auto& q = ms.fac[t];
      bool own = false;
      for (int j = 0; j < 4 && q[j] >= 0; ++j) own |= q[j] == s;
      if (own) continue;
      const Vec3 nu = fn[t] * (double)ms.fs[t];  // outward normal of the master surface
      if (dot(ns, nu) > -0.1) continue;          // the surfaces must face each other
      const Vec3 off = nu * ms.fh[t];
      const Vec3 a = X(q[0]) + off, b = X(q[1]) + off, c = X(q[2]) + off;
      double w[4] = {0, 0, 0, 0};
      int nm = 3;
      double wt[3];
      bool in;
      Vec3 qp = closest_on_tri(p, a, b, c, wt, in);
      int half = 0;
      if (q[3] >= 0) {
        nm = 4;
        const Vec3 d4 = X(q[3]) + off;
        double wt2[3];
        bool in2;
        const Vec3 qp2 = closest_on_tri(p, a, c, d4, wt2, in2);
        if (norm2(p - qp2) < norm2(p - qp)) { qp = qp2; in = in2; half = 1; std::copy(wt2, wt2 + 3, wt); }
      }
      // inside the facet (projection), up to round-off
      double wp[3];
      if (half == 0) plane_coords(p, a, b, c, wp);
      else plane_coords(p, a, c, X(q[3]) + off, wp);
      if (!in && std::min({wp[0], wp[1], wp[2]}) < -1e-6) continue;
      const double g = dot(p - qp, nu);  // signed gap
      const double d = -g;
      const double lf = norm(b - a);
      if (d <= 0 || d > sl.nh[si] + ms.fh[t] + 0.25 * lf || d >= best.d) continue;
      Hit h;
      h.s = s;
      h.nm = nm;
      for (int j = 0; j < nm; ++j) h.mn[j] = q[j];
      if (nm == 4) {  // bilinear weights from the corners (0,0) (1,0) (1,1) of half 0 or (0,0) (1,1) (0,1) of half 1
        const double xi = half ? wt[1] : wt[1] + wt[2], eta = half ? wt[1] + wt[2] : wt[2];
        w[0] = (1 - xi) * (1 - eta); w[1] = xi * (1 - eta); w[2] = xi * eta; w[3] = (1 - xi) * eta;
      } else {
        for (int j = 0; j < 3; ++j) w[j] = wt[j];
      }
      std::copy(w, w + 4, h.w);
      h.n = nu;
      h.d = d;
      best = h;
    }
    if (best.s >= 0) out.push_back(best);
  }
}

// acceleration corrections (3 per node, added to da) of one pure master-slave pass, scaled by wgt (a correction x of
// the acceleration moves the predicted position by c x); forces taken by prescribed dofs added to fkin
void pass_corrections(const Analysis& A, const std::vector<Hit>& hits, double c, double wgt, std::vector<double>& da,
                      std::vector<double>& fkin) {
  if (hits.empty() || wgt == 0) return;
  const double idt2 = 1.0 / c;
  std::unordered_map<int, std::array<double, 4>> master;  // force (3) and added mass of the master nodes
  for (const Hit& h : hits) {
    const double f = A.mt[h.s] * h.d * idt2;
    for (int j = 0; j < h.nm; ++j) {
      auto& M = master[h.mn[j]];
      for (int k = 0; k < 3; ++k) M[k] -= h.w[j] * f * h.n[k];
      M[3] += h.w[j] * A.mt[h.s];
    }
  }
  std::unordered_map<int, Vec3> dam;  // corrections of the master nodes
  for (const auto& kv : master) {
    const int n = kv.first;
    const auto& M = kv.second;
    const double mm = A.mt[n];
    Vec3 a;
    for (int k = 0; k < 3; ++k) {
      if (A.fixed[6 * n + k] || !(mm > 0)) { fkin[6 * n + k] += wgt * M[k]; continue; }  // reaction / rigid
      a[k] = M[k] / (mm + M[3]);
    }
    dam[n] = a;
  }
  for (const Hit& h : hits) {
    Vec3 am;
    for (int j = 0; j < h.nm; ++j) am += dam[h.mn[j]] * h.w[j];
    const Vec3 as = h.n * (h.d * idt2 + dot(h.n, am));
    for (int k = 0; k < 3; ++k)
      if (!A.fixed[6 * h.s + k]) da[3 * h.s + k] += wgt * as[k];
  }
  for (const auto& kv : dam)
    for (int k = 0; k < 3; ++k) da[3 * kv.first + k] += wgt * kv.second[k];
}

}  // namespace

void ContactPairs::correct(double c, std::vector<double>& du, std::vector<double>& a, std::vector<double>& fk) {
  ncontacts_ = 0;
  if (!active_ || !(c > 0)) return;
  const Analysis& A = A_;
  const int nn = A.nn;
  std::vector<double> da(3 * nn, 0.0), dat(3 * nn, 0.0);
  std::vector<Hit> h1, h2;
  auto apply = [&]() {  // corrections moved into the predicted increment and accumulated
    for (int n = 0; n < nn; ++n)
      for (int k = 0; k < 3; ++k) {
        const double x = da[3 * n + k];
        if (x == 0) continue;
        du[6 * n + k] += c * x;
        dat[3 * n + k] += x;
      }
    std::fill(da.begin(), da.end(), 0.0);
  };
  for (Pair& pr : p_->pairs) {
    if (!pr.resolved) {  // massless surfaces (rigid facets, surface elements without density) are masters
      auto massless = [&](const Side& s) {
        for (int n : s.nodes) if (A.mt[n] > 0) return false;
        return true;
      };
      pr.a.rigid = pr.a.rigid || massless(pr.a);
      pr.b.rigid = pr.b.rigid || massless(pr.b);
      pr.f = pr.a.rigid && pr.b.rigid ? -1.0 : pr.a.rigid ? 1.0 : pr.b.rigid ? 0.0
             : (pr.weight >= 0 ? std::clamp(pr.weight, 0.0, 1.0) : 0.5);
      if (pr.f < 0) warn("contact pair %s: both surfaces rigid or massless, ignored", pr.name.c_str());
      pr.resolved = true;
    }
    if (pr.f < 0) continue;
    // first surface as slave (weight 1 - f), second surface as slave (weight f)
    if (pr.f < 1) find_hits(A, pr.a, pr.b, du, h1); else h1.clear();
    if (pr.f > 0) find_hits(A, pr.b, pr.a, du, h2); else h2.clear();
    ncontacts_ += (int)(h1.size() + h2.size());
    pass_corrections(A, h1, c, 1 - pr.f, da, fk);
    pass_corrections(A, h2, c, pr.f, da, fk);
    apply();
    if (pr.f > 0 && pr.f < 1) {  // balanced: second correction of the remaining penetration, not weighted
      find_hits(A, pr.a, pr.b, du, h1);
      find_hits(A, pr.b, pr.a, du, h2);
      pass_corrections(A, h1, c, 1.0, da, fk);
      pass_corrections(A, h2, c, 1.0, da, fk);
      apply();
    }
  }
  for (int n = 0; n < nn; ++n)  // corrected accelerations, and the corresponding contact forces
    for (int k = 0; k < 3; ++k) {
      const double x = dat[3 * n + k];
      if (x == 0) continue;
      a[6 * n + k] += x;
      fk[6 * n + k] += A.mt[n] * x;
    }
}

}  // namespace abq
