#include "model/model.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <unordered_map>

#include "core/log.hpp"
#include "inp/parser.hpp"

namespace abq {

const char* eltype_name(ElType t) {
  switch (t) {
    case ElType::C3D8R: return "C3D8R";
    case ElType::C3D8: return "C3D8";
    case ElType::C3D8I: return "C3D8I";
    case ElType::C3D4: return "C3D4";
    case ElType::S4R: return "S4R";
    case ElType::S4: return "S4";
    case ElType::S3R: return "S3R";
    case ElType::S3RS: return "S3RS";
    case ElType::B31: return "B31";
    case ElType::T3D2: return "T3D2";
    case ElType::CONN3D2: return "CONN3D2";
    default: return "UNKNOWN";
  }
}

ElType eltype_from(const std::string& n0, int& nn) {
  std::string n = upper(n0);
  if (n == "C3D8R") { nn = 8; return ElType::C3D8R; }
  if (n == "C3D8") { nn = 8; return ElType::C3D8; }
  if (n == "C3D8I") { nn = 8; return ElType::C3D8I; }
  if (n == "C3D4") { nn = 4; return ElType::C3D4; }
  if (n == "S4R" || n == "S4RS" || n == "S4R5") { nn = 4; return ElType::S4R; }
  if (n == "S4") { nn = 4; return ElType::S4; }
  if (n == "S3R" || n == "S3") { nn = 3; return ElType::S3R; }
  if (n == "S3RS") { nn = 3; return ElType::S3RS; }
  if (n == "B31" || n == "B31H") { nn = 2; return ElType::B31; }
  if (n == "T3D2") { nn = 2; return ElType::T3D2; }
  if (n == "CONN3D2") { nn = 2; return ElType::CONN3D2; }
  nn = 0;
  return ElType::UNKNOWN;
}

double Amplitude::eval(double st, double tt) const {
  double tm = total_time ? tt : st;
  if (t.empty()) return 1.0;
  const auto& T = t;
  if (tm <= T.front()) return v.front();
  if (tm >= T.back()) return v.back();
  size_t i = std::upper_bound(T.begin(), T.end(), tm) - T.begin();
  if (i == 0) i = 1;
  double t0 = T[i - 1], t1 = T[i], dt = t1 - t0;
  double xi = dt > 0 ? (tm - t0) / dt : 1.0;
  if (kind == SMOOTH_STEP) xi = xi * xi * xi * (10.0 - 15.0 * xi + 6.0 * xi * xi);
  return v[i - 1] + (v[i] - v[i - 1]) * xi;
}

double Material::yield(double ep, double* H) const {
  if (plastic_sy.empty()) { if (H) *H = 0; return 1e300; }
  const auto& s = plastic_sy;
  const auto& e = plastic_ep;
  if (s.size() == 1 || ep <= e.front()) { if (H) *H = s.size() > 1 ? (s[1] - s[0]) / std::max(e[1] - e[0], 1e-30) : 0; return s.front(); }
  if (ep >= e.back()) { if (H) *H = 0; return s.back(); }
  size_t i = std::upper_bound(e.begin(), e.end(), ep) - e.begin();
  double h = (s[i] - s[i - 1]) / std::max(e[i] - e[i - 1], 1e-30);
  if (H) *H = h;
  return s[i - 1] + h * (ep - e[i - 1]);
}

double Material::shear_modulus() const {
  if (elastic) return E / (2 * (1 + nu));
  if (fung) {
    double g = 0;
    for (int i = 3; i < 6; ++i) g = std::max(g, fung_b[i][i]);
    return fung_c * g;
  }
  if (neohooke) return 2 * C10;
  return 0;
}
double Material::bulk_modulus() const {
  if (elastic) return E / (3 * (1 - 2 * nu));
  if (fung) return fung_D > 0 ? 2.0 / fung_D : 20 * shear_modulus();
  if (neohooke) return D1 > 0 ? 2.0 / D1 : 20 * shear_modulus();
  return 0;
}

double ConnectorComponent::force(double urel, double fv1, double* dFdu) const {
  if (!nonlinear) {
    double ke = k;
    if (kf.size() > 1) {  // stiffness interpolated in the field variable (constant beyond the table)
      std::vector<std::pair<double, double>> t;
      for (size_t i = 0; i < kf.size(); ++i) t.push_back({fv[i], kf[i]});
      std::sort(t.begin(), t.end());
      if (fv1 <= t.front().first) ke = t.front().second;
      else if (fv1 >= t.back().first) ke = t.back().second;
      else {
        size_t i = 1;
        while (t[i].first < fv1) ++i;
        ke = t[i - 1].second + (t[i].second - t[i - 1].second) * (fv1 - t[i - 1].first) / (t[i].first - t[i - 1].first);
      }
    }
    if (dFdu) *dFdu = ke;
    return ke * urel;
  }
  // table: groups of equal field value, each sorted by u
  std::vector<double> fvals;
  for (double f : fv) if (std::find(fvals.begin(), fvals.end(), f) == fvals.end()) fvals.push_back(f);
  std::sort(fvals.begin(), fvals.end());
  auto interp_group = [&](double fval, double* d) {
    std::vector<std::pair<double, double>> pts;
    for (size_t i = 0; i < u.size(); ++i) if (fv[i] == fval) pts.push_back({u[i], F[i]});
    std::sort(pts.begin(), pts.end());
    if (pts.empty()) { *d = 0; return 0.0; }
    if (pts.size() == 1) { *d = 0; return pts[0].second; }
    if (urel <= pts.front().first) { *d = 0; return pts.front().second; }
    if (urel >= pts.back().first) { *d = 0; return pts.back().second; }
    size_t i = 1;
    while (i < pts.size() && pts[i].first < urel) ++i;
    double s = (pts[i].second - pts[i - 1].second) / (pts[i].first - pts[i - 1].first);
    *d = s;
    return pts[i - 1].second + s * (urel - pts[i - 1].first);
  };
  double d = 0;
  if (fvals.size() <= 1 || ndep == 0) {
    double f = interp_group(fvals.empty() ? 0 : fvals[0], &d);
    if (dFdu) *dFdu = d;
    return f;
  }
  if (fv1 <= fvals.front()) { double f = interp_group(fvals.front(), &d); if (dFdu) *dFdu = d; return f; }
  if (fv1 >= fvals.back()) { double f = interp_group(fvals.back(), &d); if (dFdu) *dFdu = d; return f; }
  size_t j = 1;
  while (j < fvals.size() && fvals[j] < fv1) ++j;
  double d0, d1;
  double f0 = interp_group(fvals[j - 1], &d0), f1 = interp_group(fvals[j], &d1);
  double w = (fv1 - fvals[j - 1]) / (fvals[j] - fvals[j - 1]);
  if (dFdu) *dFdu = (1 - w) * d0 + w * d1;
  return (1 - w) * f0 + w * f1;
}

std::string Model::node_name(int i) const {
  int in = node_instance[i];
  return (in >= 0 && !instances[in].empty() ? instances[in] + "." : std::string()) + std::to_string(node_label[i]);
}
std::string Model::elem_name(int e) const {
  int in = elems[e].instance;
  return (in >= 0 && !instances[in].empty() ? instances[in] + "." : std::string()) + std::to_string(elems[e].label);
}

namespace {

struct Scope {
  std::string prefix;  // "MV." or ""
  int inst = -1;
  Vec3 trans;
  bool rotate = false;
  Vec3 ra, rb;
  double angle = 0;
};

struct PendingSection {
  Section s;
  std::string prefix;
  std::string orient_name;
  std::string mat_name;
  std::string behavior_name;
  std::string controls_name;
  std::string where;
};

class Builder {
 public:
  Model m;
  std::vector<Card> cards;
  std::map<std::string, std::pair<size_t, size_t>> parts;  // name -> [first card, end card)
  std::unordered_map<std::string, int> node_key, elem_key;
  std::vector<PendingSection> pend_sections;
  struct ControlsCard { double hg_scale = 1.0; bool enhanced = false, face_dt = true, relax = true; };
  std::map<std::string, ControlsCard> section_controls;  // *SECTION CONTROLS by name
  std::set<std::string> warned;
  int cur_material = -1;
  int cur_behavior = -1;
  int cur_interaction = -1;
  enum { CTX_NONE, CTX_MATERIAL, CTX_BEHAVIOR, CTX_INTERACTION } ctx = CTX_NONE;

  void warn_once(const std::string& key, const std::string& msg) {
    if (warned.insert(key).second) warn("%s", msg.c_str());
  }

  std::string resolve_name(const std::string& raw, const std::string& prefix) {
    std::string n = upper(trim(raw));
    if (n.rfind("ASSEMBLY.", 0) == 0) n = n.substr(9);
    return n;
  }
  const std::vector<int>* find_set(std::map<std::string, std::vector<int>>& sets, const std::string& raw,
                                   const std::string& prefix) {
    std::string n = resolve_name(raw, prefix);
    auto it = sets.find(prefix + n);
    if (it != sets.end()) return &it->second;
    it = sets.find(n);
    if (it != sets.end()) return &it->second;
    return nullptr;
  }
  int find_node(const std::string& tok, const std::string& prefix) {
    std::string t = upper(trim(tok));
    if (t.rfind("ASSEMBLY.", 0) == 0) t = t.substr(9);
    auto it = node_key.find(is_number(t) ? prefix + std::to_string(to_int(t)) : t);
    if (it == node_key.end()) {
      // "INST.label" with integer label formatting
      size_t dot = t.rfind('.');
      if (dot != std::string::npos && is_number(t.substr(dot + 1))) {
        it = node_key.find(t.substr(0, dot + 1) + std::to_string(to_int(t.substr(dot + 1))));
        if (it != node_key.end()) return it->second;
      }
      return -1;
    }
    return it->second;
  }
  int find_elem(const std::string& tok, const std::string& prefix) {
    std::string t = upper(trim(tok));
    if (t.rfind("ASSEMBLY.", 0) == 0) t = t.substr(9);
    auto it = elem_key.find(is_number(t) ? prefix + std::to_string(to_int(t)) : t);
    return it == elem_key.end() ? -1 : it->second;
  }
  // node list from a token which is either a node label or a node set name
  std::vector<int> nodes_of(const std::string& tok, const std::string& prefix, const Card& c) {
    std::string t = trim(tok);
    if (is_number(t) || (t.find('.') != std::string::npos && is_number(t.substr(t.rfind('.') + 1)))) {
      int n = find_node(t, prefix);
      if (n >= 0) return {n};
    }
    auto* s = find_set(m.nsets, t, prefix);
    if (s) return *s;
    int n = find_node(t, prefix);
    if (n >= 0) return {n};
    throw std::runtime_error("unknown node or node set '" + t + "' at " + c.where());
  }
  std::vector<int> elems_of(const std::string& tok, const std::string& prefix, const Card& c) {
    std::string t = trim(tok);
    if (is_number(t)) {
      int e = find_elem(t, prefix);
      if (e >= 0) return {e};
    }
    auto* s = find_set(m.elsets, t, prefix);
    if (s) return *s;
    int e = find_elem(t, prefix);
    if (e >= 0) return {e};
    throw std::runtime_error("unknown element or element set '" + t + "' at " + c.where());
  }

  // -------------------------------------------------------------- geometry cards
  void card_node(const Card& c, const Scope& sc) {
    std::string nset = c.has("NSET") ? sc.prefix + upper(c.get("NSET")) : "";
    for (auto& d : c.data) {
      if (d.empty() || d[0].empty()) continue;
      int label = to_int(d[0], &c);
      Vec3 x(d.size() > 1 ? to_double(d[1], &c) : 0, d.size() > 2 ? to_double(d[2], &c) : 0,
             d.size() > 3 ? to_double(d[3], &c) : 0);
      if (sc.rotate) x = rotate_point(x, sc);
      x += sc.trans;
      std::string key = sc.prefix + std::to_string(label);
      int idx;
      auto it = node_key.find(key);
      if (it != node_key.end()) { idx = it->second; m.X[idx] = x; }
      else {
        idx = m.nn();
        node_key[key] = idx;
        m.X.push_back(x);
        m.node_label.push_back(label);
        m.node_instance.push_back(sc.inst);
      }
      if (!nset.empty()) m.nsets[nset].push_back(idx);
    }
  }
  Vec3 rotate_point(const Vec3& x, const Scope& sc) {
    Vec3 axis = normalized(sc.rb - sc.ra);
    Mat3 R = rotation_from_vector(axis * (sc.angle * M_PI / 180.0));
    return sc.ra + R * (x - sc.ra);
  }

  void card_element(const Card& c, const Scope& sc) {
    int nn = 0;
    std::string tname = c.get("TYPE");
    ElType t = eltype_from(tname, nn);
    if (t == ElType::UNKNOWN) {
      warn_once("ELTYPE" + upper(tname), "element type " + upper(tname) + " not supported, elements ignored (" + c.where() + ")");
      return;
    }
    std::string elset = c.has("ELSET") ? sc.prefix + upper(c.get("ELSET")) : "";
    std::vector<std::string> buf;
    auto flush = [&](std::vector<std::string>& toks) {
      Element e;
      e.type = t;
      e.label = to_int(toks[0], &c);
      e.instance = sc.inst;
      for (size_t i = 1; i < toks.size() && (int)e.nodes.size() < nn; ++i) {
        if (toks[i].empty()) continue;
        int n = find_node(toks[i], sc.prefix);
        if (n < 0) throw std::runtime_error("element " + toks[0] + ": unknown node " + toks[i] + " at " + c.where());
        e.nodes.push_back(n);
      }
      if ((int)e.nodes.size() < nn) {
        if (t == ElType::CONN3D2 && e.nodes.size() == 1) {
          warn_once("CONNGROUND", "grounded CONN3D2 (single node) not supported, ignored");
          return;
        }
        throw std::runtime_error("element " + toks[0] + ": not enough nodes at " + c.where());
      }
      if ((t == ElType::S4R || t == ElType::S4) && e.nodes[2] == e.nodes[3]) { e.type = ElType::S3R; e.nodes.pop_back(); }
      std::string key = sc.prefix + std::to_string(e.label);
      int idx = m.ne();
      elem_key[key] = idx;
      m.elems.push_back(std::move(e));
      if (!elset.empty()) m.elsets[elset].push_back(idx);
    };
    for (auto& d : c.data) {
      for (auto& tok : d) buf.push_back(tok);
      // count node tokens (non-empty after the label)
      int cnt = 0;
      for (size_t i = 1; i < buf.size(); ++i) if (!buf[i].empty()) ++cnt;
      if (cnt >= nn || (t == ElType::CONN3D2 && !buf.empty() && cnt >= 1 && d.back() != "")) {
        flush(buf);
        buf.clear();
      } else if (!buf.empty() && buf.back().empty()) {
        buf.pop_back();  // trailing comma continuation
      }
    }
    if (!buf.empty()) flush(buf);
  }

  // a set holds each member once (a member listed again, or a set defined again, adds nothing)
  static void dedupe(std::vector<int>& v) {
    if (v.empty()) return;
    std::vector<char> seen((size_t)*std::max_element(v.begin(), v.end()) + 1, 0);
    size_t k = 0;
    for (int x : v)
      if (x < 0 || !seen[x]) { if (x >= 0) seen[x] = 1; v[k++] = x; }
    v.resize(k);
  }

  void card_set(const Card& c, const Scope& sc, bool is_node) {
    std::string name = upper(c.get(is_node ? "NSET" : "ELSET"));
    std::string target_prefix = sc.prefix;
    if (c.has("INSTANCE")) target_prefix = upper(c.get("INSTANCE")) + ".";
    std::string full = sc.prefix + name;
    auto& sets = is_node ? m.nsets : m.elsets;
    auto& out = sets[full];
    if (is_node && c.has("ELSET")) {  // *NSET, ELSET=<element set>: the nodes of its elements
      auto it = m.elsets.find(target_prefix + upper(c.get("ELSET")));
      if (it == m.elsets.end()) it = m.elsets.find(upper(c.get("ELSET")));
      if (it == m.elsets.end()) throw std::runtime_error("*NSET: unknown element set " + c.get("ELSET") + " (" + c.where() + ")");
      for (int e : it->second) out.insert(out.end(), m.elems[e].nodes.begin(), m.elems[e].nodes.end());
      dedupe(out);
      return;
    }
    bool gen = c.has("GENERATE");
    for (auto& d : c.data) {
      if (gen) {
        if (d.size() < 2) continue;
        int a = to_int(d[0], &c), b = to_int(d[1], &c), s = d.size() > 2 && !d[2].empty() ? to_int(d[2], &c) : 1;
        for (int k = a; k <= b; k += s) {
          int idx = is_node ? find_node(std::to_string(k), target_prefix) : find_elem(std::to_string(k), target_prefix);
          if (idx >= 0) out.push_back(idx);
        }
        continue;
      }
      for (auto& tok : d) {
        if (tok.empty()) continue;
        std::vector<int> v = is_node ? nodes_of(tok, target_prefix, c) : elems_of(tok, target_prefix, c);
        out.insert(out.end(), v.begin(), v.end());
      }
    }
    dedupe(out);
  }

  void card_surface(const Card& c, const Scope& sc) {
    Surface s;
    s.name = sc.prefix + upper(c.get("NAME"));
    std::string type = c.getu("TYPE", "ELEMENT");
    if (type == "NODE") {
      s.node_based = true;
      for (auto& d : c.data) {
        if (d.empty() || d[0].empty()) continue;
        auto v = nodes_of(d[0], sc.prefix, c);
        s.nodes.insert(s.nodes.end(), v.begin(), v.end());
      }
    } else {
      for (auto& d : c.data) {
        if (d.empty() || d[0].empty()) continue;
        int face = 0;
        std::string f = d.size() > 1 ? upper(d[1]) : "";
        if (f == "SPOS") face = 1;
        else if (f == "SNEG") face = -1;
        else if (f.size() == 2 && f[0] == 'S' && std::isdigit((unsigned char)f[1])) face = f[1] - '0';
        else if (!f.empty()) warn_once("FACE" + f, "surface face label " + f + " not supported (" + c.where() + ")");
        // could also be the name of another surface
        std::vector<int> es;
        auto* sub = find_set(m.elsets, d[0], sc.prefix);
        if (sub || is_number(trim(d[0]))) es = elems_of(d[0], sc.prefix, c);
        else {
          std::string n = resolve_name(d[0], sc.prefix);
          auto it = m.surfaces.find(sc.prefix + n);
          if (it == m.surfaces.end()) it = m.surfaces.find(n);
          if (it == m.surfaces.end()) throw std::runtime_error("unknown set/surface '" + d[0] + "' at " + c.where());
          s.faces.insert(s.faces.end(), it->second.faces.begin(), it->second.faces.end());
          continue;
        }
        for (int e : es) s.faces.push_back({e, face});
      }
    }
    m.surfaces[s.name] = std::move(s);
  }

  void card_orientation(const Card& c, const Scope& sc) {
    Orientation o;
    o.name = sc.prefix + upper(c.get("NAME"));
    std::string sys = c.getu("SYSTEM", "RECTANGULAR");
    if (sys != "RECTANGULAR") warn_once("ORISYS" + sys, "orientation system " + sys + " treated as RECTANGULAR");
    if (!c.data.empty() && c.data[0].size() >= 6) {
      auto& d = c.data[0];
      Vec3 a(to_double(d[0]), to_double(d[1]), to_double(d[2]));
      Vec3 b(to_double(d[3]), to_double(d[4]), to_double(d[5]));
      Vec3 orig;
      if (d.size() >= 9 && !d[6].empty()) orig = Vec3(to_double(d[6]), to_double(d[7]), to_double(d[8]));
      Vec3 e1 = normalized(a - orig);
      Vec3 e3 = normalized(cross(a - orig, b - orig));
      Vec3 e2 = cross(e3, e1);
      o.R = Mat3::rows(e1, e2, e3);
      if (c.data.size() > 1 && c.data[1].size() >= 2 && !c.data[1][1].empty()) {
        int ax = to_int(c.data[1][0]);
        double ang = to_double(c.data[1][1]) * M_PI / 180.0;
        if (ang != 0 && ax >= 1 && ax <= 3) {
          Vec3 axis = o.R.row(ax - 1);
          Mat3 Rr = rotation_from_vector(axis * ang);
          Vec3 r0 = Rr * o.R.row(0), r1 = Rr * o.R.row(1), r2 = Rr * o.R.row(2);
          o.R = Mat3::rows(r0, r1, r2);
        }
      }
    }
    m.orient_index[o.name] = (int)m.orients.size();
    m.orients.push_back(o);
  }

  void card_section(const Card& c, const Scope& sc) {
    PendingSection p;
    p.prefix = sc.prefix;
    p.where = c.where();
    Section& s = p.s;
    s.elset = upper(c.get("ELSET"));
    p.mat_name = upper(c.get("MATERIAL"));
    p.orient_name = upper(c.get("ORIENTATION"));
    p.controls_name = upper(c.get("CONTROLS"));
    if (c.keyword == "SHELL SECTION") {
      s.kind = Section::SHELL;
      if (c.has("COMPOSITE")) warn_once("COMPOSITE", "composite shell sections not supported");
      if (!c.data.empty()) {
        s.thickness = to_double(c.data[0][0], &c);
        if (c.data[0].size() > 1 && !c.data[0][1].empty()) s.nip = to_int(c.data[0][1], &c);
      }
    } else if (c.keyword == "SOLID SECTION") {
      s.kind = Section::SOLID;
    } else if (c.keyword == "BEAM SECTION") {
      s.kind = Section::BEAM;
      s.shape = c.getu("SECTION");
      s.poisson = c.getd("POISSON", 0.0);
      s.rotary_isotropic = c.getu("ROTARY INERTIA") == "ISOTROPIC";
      if (!c.data.empty()) for (auto& t : c.data[0]) if (!t.empty()) s.dims.push_back(to_double(t, &c));
      if (c.data.size() > 1 && c.data[1].size() >= 3)
        s.n1 = Vec3(to_double(c.data[1][0]), to_double(c.data[1][1]), to_double(c.data[1][2]));
      if (s.shape != "CIRC" && s.shape != "RECT" && s.shape != "PIPE")
        throw std::runtime_error("beam section " + s.shape + " not supported at " + c.where());
    } else if (c.keyword == "CONNECTOR SECTION") {
      s.kind = Section::CONNECTOR;
      p.behavior_name = upper(c.get("BEHAVIOR"));
      std::string t = c.data.empty() ? "" : upper(c.data[0][0]);
      if (t == "AXIAL") s.ctype = Section::AXIAL;
      else if (t == "CARTESIAN") s.ctype = Section::CARTESIAN;
      else if (t == "JOIN") s.ctype = Section::JOIN;
      else { s.ctype = Section::OTHER; warn_once("CONN" + t, "connector type " + t + " not supported (" + c.where() + ")"); }
      if (c.data.size() > 1 && !c.data[1].empty() && !c.data[1][0].empty()) p.orient_name = upper(c.data[1][0]);
    }
    pend_sections.push_back(p);
  }

  // -------------------------------------------------------------- materials
  void card_material_sub(const Card& c) {
    if (cur_material < 0) return;
    Material& mat = m.materials[cur_material];
    const std::string& k = c.keyword;
    if (k == "DENSITY") mat.density = to_double(c.data.at(0).at(0), &c);
    else if (k == "ELASTIC") {
      std::string t = c.getu("TYPE", "ISOTROPIC");
      if (t != "ISOTROPIC") throw std::runtime_error("*ELASTIC TYPE=" + t + " not supported at " + c.where());
      mat.elastic = true;
      mat.E = to_double(c.data.at(0).at(0), &c);
      mat.nu = c.data[0].size() > 1 ? to_double(c.data[0][1], &c) : 0.0;
      if (c.data.size() > 1) warn_once("ELTEMP" + mat.name, "temperature-dependent elasticity: first row used (" + mat.name + ")");
    } else if (k == "PLASTIC") {
      std::string h = c.getu("HARDENING", "ISOTROPIC");
      if (h != "ISOTROPIC") warn_once("HARD" + h, "*PLASTIC HARDENING=" + h + " treated as ISOTROPIC");
      for (auto& d : c.data) {
        mat.plastic_sy.push_back(to_double(d.at(0), &c));
        mat.plastic_ep.push_back(d.size() > 1 ? to_double(d[1], &c) : 0.0);
      }
    } else if (k == "EXPANSION") {
      mat.alpha = to_double(c.data.at(0).at(0), &c);
      if (c.has("ZERO") && c.getd("ZERO", 0) != 0) warn_once("EXPZERO", "*EXPANSION ZERO ignored (reference 0)");
    } else if (k == "DAMPING") {
      mat.damp_alpha = c.getd("ALPHA", 0);
      mat.damp_beta = c.getd("BETA", 0);
    } else if (k == "ANISOTROPIC HYPERELASTIC") {
      if (!c.has("FUNG-ANISOTROPIC") && !c.has("FUNG-ORTHOTROPIC"))
        throw std::runtime_error("only FUNG-ANISOTROPIC / FUNG-ORTHOTROPIC supported at " + c.where());
      std::vector<double> v;
      for (auto& d : c.data) for (auto& t : d) if (!t.empty()) v.push_back(to_double(t, &c));
      mat.fung = true;
      // Voigt index of tensor pairs: 11->0, 22->1, 33->2, 12->3, 13->4, 23->5
      auto set = [&](int I, int J, double x) { mat.fung_b[I][J] = mat.fung_b[J][I] = x; };
      if (c.has("FUNG-ANISOTROPIC")) {
        if (v.size() < 23) throw std::runtime_error("FUNG-ANISOTROPIC needs 23 constants at " + c.where());
        // order: b1111 b1122 b2222 b1133 b2233 b3333 b1112 b2212 b3312 b1212 b1113 b2213 b3313 b1213 b1313
        //        b1123 b2223 b3323 b1223 b1323 b2323 c D
        const int I[21] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5};
        const int J[21] = {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5};
        for (int i = 0; i < 21; ++i) set(I[i], J[i], v[i]);
        mat.fung_c = v[21];
        mat.fung_D = v[22];
      } else {
        if (v.size() < 11) throw std::runtime_error("FUNG-ORTHOTROPIC needs 11 constants at " + c.where());
        // b1111 b1122 b2222 b1133 b2233 b3333 b1212 b1313 b2323 c D
        set(0, 0, v[0]); set(0, 1, v[1]); set(1, 1, v[2]); set(0, 2, v[3]); set(1, 2, v[4]); set(2, 2, v[5]);
        set(3, 3, v[6]); set(4, 4, v[7]); set(5, 5, v[8]);
        mat.fung_c = v[9];
        mat.fung_D = v[10];
      }
    } else if (k == "HYPERELASTIC") {
      if (!c.has("NEO HOOKE")) throw std::runtime_error("*HYPERELASTIC: only NEO HOOKE supported at " + c.where());
      mat.neohooke = true;
      mat.C10 = to_double(c.data.at(0).at(0), &c);
      mat.D1 = c.data[0].size() > 1 ? to_double(c.data[0][1], &c) : 0.0;
    } else {
      warn_once("MATKW" + k, "material option *" + k + " ignored (" + c.where() + ")");
    }
  }

  void card_behavior_sub(const Card& c) {
    if (cur_behavior < 0) return;
    ConnectorBehavior& b = m.conn_behaviors[cur_behavior];
    if (c.keyword == "CONNECTOR ELASTICITY") {
      if (c.has("RIGID")) { warn_once("CONNRIGID", "rigid connector elasticity not supported"); return; }
      int comp = c.geti("COMPONENT", 1);
      if (comp < 1 || comp > 6) return;
      ConnectorComponent& cc = b.comp[comp - 1];
      cc.elastic = true;
      cc.nonlinear = c.has("NONLINEAR");
      cc.ndep = c.geti("DEPENDENCIES", 0);
      if (!cc.nonlinear) {
        cc.k = to_double(c.data.at(0).at(0), &c);
        if (cc.ndep > 0)  // rows: stiffness, frequency, temperature, field 1
          for (auto& d : c.data) {
            cc.kf.push_back(to_double(d.at(0), &c));
            cc.fv.push_back(d.size() > 3 && !d[3].empty() ? to_double(d[3], &c) : 0.0);
          }
      }
      else {
        for (auto& d : c.data) {
          cc.F.push_back(to_double(d.at(0), &c));
          cc.u.push_back(to_double(d.at(1), &c));
          cc.fv.push_back(cc.ndep > 0 && d.size() > 3 ? to_double(d[3], &c) : 0.0);
        }
      }
    } else if (c.keyword == "CONNECTOR CONSTITUTIVE REFERENCE") {
      if (!c.data.empty())
        for (size_t i = 0; i < c.data[0].size() && i < 6; ++i)
          if (!c.data[0][i].empty()) { b.comp[i].has_ref = true; b.comp[i].ref = to_double(c.data[0][i], &c); }
    } else {
      warn_once("CONNKW" + c.keyword, "connector option *" + c.keyword + " ignored");
    }
  }

  void card_amplitude(const Card& c) {
    Amplitude a;
    a.name = upper(c.get("NAME"));
    std::string def = c.getu("DEFINITION", "TABULAR");
    if (def == "SMOOTH STEP") a.kind = Amplitude::SMOOTH_STEP;
    else if (def != "TABULAR") throw std::runtime_error("amplitude definition " + def + " not supported at " + c.where());
    a.total_time = c.getu("TIME", "STEP TIME") == "TOTAL TIME";
    std::vector<double> v;
    for (auto& d : c.data) for (auto& t : d) if (!t.empty()) v.push_back(to_double(t, &c));
    for (size_t i = 0; i + 1 < v.size(); i += 2) { a.t.push_back(v[i]); a.v.push_back(v[i + 1]); }
    m.amplitude_index[a.name] = (int)m.amplitudes.size();
    m.amplitudes.push_back(a);
  }

  int amplitude_of(const Card& c) {
    if (!c.has("AMPLITUDE")) return -1;
    auto it = m.amplitude_index.find(upper(c.get("AMPLITUDE")));
    if (it == m.amplitude_index.end()) throw std::runtime_error("unknown amplitude at " + c.where());
    return it->second;
  }

  // BC data line -> list of (node, dof, value)
  void bc_line(const Card& c, const std::vector<std::string>& d, std::vector<NodalValue>& out) {
    auto nodes = nodes_of(d.at(0), "", c);
    std::string t1 = d.size() > 1 ? upper(d[1]) : "";
    std::vector<int> dofs;
    std::vector<double> vals;
    if (!t1.empty() && !is_number(t1)) {
      if (t1 == "ENCASTRE") dofs = {0, 1, 2, 3, 4, 5};
      else if (t1 == "PINNED") dofs = {0, 1, 2};
      else if (t1 == "XSYMM") dofs = {0, 4, 5};
      else if (t1 == "YSYMM") dofs = {1, 3, 5};
      else if (t1 == "ZSYMM") dofs = {2, 3, 4};
      else if (t1 == "XASYMM") dofs = {1, 2, 3};
      else if (t1 == "YASYMM") dofs = {0, 2, 4};
      else if (t1 == "ZASYMM") dofs = {0, 1, 5};
      else throw std::runtime_error("unknown boundary type " + t1 + " at " + c.where());
      vals.assign(dofs.size(), 0.0);
    } else {
      int a = to_int(d.at(1), &c);
      int b = d.size() > 2 && !d[2].empty() ? to_int(d[2], &c) : a;
      double v = d.size() > 3 && !d[3].empty() ? to_double(d[3], &c) : 0.0;
      for (int k = a; k <= b; ++k) { if (k >= 1 && k <= 6) { dofs.push_back(k - 1); vals.push_back(v); } }
    }
    for (int n : nodes)
      for (size_t i = 0; i < dofs.size(); ++i) out.push_back({n, dofs[i], vals[i]});
  }

  // -------------------------------------------------------------- general contact
  // *CONTACT, *CONTACT INCLUSIONS, *CONTACT PROPERTY ASSIGNMENT, in the model data (active in every step) or in a step
  Step model_contact;
  bool contact_keyword(const std::string& k, const Card& c, Step& st) {
    if (k == "CONTACT") {
      st.contact_defined = true;
    } else if (k == "CONTACT INCLUSIONS") {
      st.contact_defined = true;
      if (c.has("ALL EXTERIOR")) st.contact_all_exterior = true;
      for (auto& d : c.data) {
        for (auto& t : d) if (!t.empty()) {
          std::string n = resolve_name(t, "");
          if (!m.surfaces.count(n)) throw std::runtime_error("unknown contact surface '" + t + "' at " + c.where());
          st.contact_surfaces.push_back(n);
        }
      }
    } else if (k == "CONTACT PROPERTY ASSIGNMENT") {
      for (auto& d : c.data) {
        if (d.size() >= 3 && d[0].empty() && d[1].empty()) {
          auto it = m.interaction_index.find(upper(d[2]));
          if (it == m.interaction_index.end()) throw std::runtime_error("unknown surface interaction at " + c.where());
          st.contact_property = it->second;
        } else warn_once("CPA", "surface-specific contact property assignment not supported (global only)");
      }
    } else {
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------- steps
  void parse_step(size_t& i) {
    const Card& sc = cards[i];
    Step st;
    st.contact_defined = model_contact.contact_defined;  // model-level general contact applies to every step
    st.contact_all_exterior = model_contact.contact_all_exterior;
    st.contact_surfaces = model_contact.contact_surfaces;
    st.contact_property = model_contact.contact_property;
    st.name = sc.get("NAME", "Step-" + std::to_string(m.steps.size() + 1));
    const bool nlg_given = sc.has("NLGEOM");
    st.nlgeom = sc.getu("NLGEOM", "NO") == "YES" || (!nlg_given && !m.steps.empty() && m.steps.back().nlgeom);
    st.amp_ramp = sc.has("AMPLITUDE") ? (sc.getu("AMPLITUDE") == "RAMP" ? 1 : 0) : -1;
    st.maxinc = sc.geti("INC", 100);
    bool have_proc = false;
    int last_output = 0;  // 1 field, 2 history
    int hist_freq = 0;    // FREQUENCY of the current history request
    for (++i; i < cards.size(); ++i) {
      const Card& c = cards[i];
      const std::string& k = c.keyword;
      if (k == "END STEP") break;
      if (k == "STATIC") {
        st.proc = Step::STATIC;
        have_proc = true;
        if (!c.data.empty()) {
          auto& d = c.data[0];
          st.dt0 = d.size() > 0 && !d[0].empty() ? to_double(d[0], &c) : 1.0;
          st.period = d.size() > 1 && !d[1].empty() ? to_double(d[1], &c) : 1.0;
          st.dtmin = d.size() > 2 && !d[2].empty() ? to_double(d[2], &c) : std::min(st.dt0, 1e-5 * st.period);
          st.dtmax = d.size() > 3 && !d[3].empty() ? to_double(d[3], &c) : st.period;
        } else { st.dt0 = st.period = st.dtmax = 1.0; st.dtmin = 1e-5; }
      } else if (k == "DYNAMIC") {
        if (!c.has("EXPLICIT")) throw std::runtime_error("implicit *DYNAMIC not supported (use EXPLICIT) at " + c.where());
        st.proc = Step::EXPLICIT;
        have_proc = true;
        if (c.has("DIRECT USER CONTROL") || c.has("FIXED TIME INCREMENTATION")) st.fixed_dt = true;
        st.dt_scale = c.getd("SCALE FACTOR", 1.0);
        st.element_by_element = c.has("ELEMENT BY ELEMENT");
        if (!c.data.empty()) {
          auto& d = c.data[0];
          if (d.size() > 1 && !d[1].empty()) st.period = to_double(d[1], &c);
          if (d.size() > 0 && !d[0].empty()) st.user_dt = to_double(d[0], &c);
        }
      } else if (k == "BULK VISCOSITY") {
        if (!c.data.empty()) {
          if (!c.data[0].empty() && !c.data[0][0].empty()) st.bv1 = to_double(c.data[0][0], &c);
          if (c.data[0].size() > 1 && !c.data[0][1].empty()) st.bv2 = to_double(c.data[0][1], &c);
        }
      } else if (k == "VARIABLE MASS SCALING" || k == "FIXED MASS SCALING") {
        st.mscale.on = true;
        st.mscale.dt = c.getd("DT", 0);
        st.mscale.factor = c.getd("FACTOR", 0);
        std::string t = c.getu("TYPE", "BELOW MIN");
        st.mscale.below_min = (t == "BELOW MIN");
        if (t != "BELOW MIN" && t != "UNIFORM") warn_once("MSTYPE" + t, "mass scaling TYPE=" + t + " treated as BELOW MIN");
        st.mscale.frequency = k == "VARIABLE MASS SCALING" ? c.geti("FREQUENCY", 0) : 0;
        st.mscale.number_interval = c.geti("NUMBER INTERVAL", 0);
        if (c.has("ELSET")) warn_once("MSELSET", "mass scaling ELSET ignored: applied to the whole model");
      } else if (k == "BOUNDARY" || k == "CLOAD" || k == "TEMPERATURE" || k == "FIELD") {
        LoadBlock lb;
        lb.kind = k == "BOUNDARY" ? (c.getu("TYPE", "DISPLACEMENT") == "VELOCITY" ? LoadBlock::VELOCITY_BC : LoadBlock::BOUNDARY)
                  : k == "CLOAD" ? LoadBlock::CLOAD : k == "TEMPERATURE" ? LoadBlock::TEMPERATURE : LoadBlock::FIELD;
        lb.op_new = c.getu("OP", "MOD") == "NEW";
        lb.amplitude = amplitude_of(c);
        lb.field_var = c.geti("VARIABLE", 1);
        for (auto& d : c.data) {
          if (d.empty() || d[0].empty()) continue;
          if (lb.kind == LoadBlock::BOUNDARY || lb.kind == LoadBlock::VELOCITY_BC) bc_line(c, d, lb.nodal);
          else if (lb.kind == LoadBlock::CLOAD) {
            auto nodes = nodes_of(d[0], "", c);
            int dof = to_int(d.at(1), &c) - 1;
            double v = d.size() > 2 ? to_double(d[2], &c) : 0.0;
            for (int n : nodes) lb.nodal.push_back({n, dof, v});
          } else {
            auto nodes = nodes_of(d[0], "", c);
            double v = d.size() > 1 ? to_double(d[1], &c) : 0.0;
            for (int n : nodes) lb.nodal.push_back({n, 0, v});
          }
        }
        st.loads.push_back(std::move(lb));
      } else if (k == "DSLOAD" || k == "DLOAD") {
        LoadBlock lb;
        lb.kind = LoadBlock::DSLOAD;
        lb.op_new = c.getu("OP", "MOD") == "NEW";
        lb.amplitude = amplitude_of(c);
        for (auto& d : c.data) {
          if (d.empty() || d[0].empty()) continue;
          std::string lt = d.size() > 1 ? upper(d[1]) : "P";
          double v = d.size() > 2 ? to_double(d[2], &c) : 0.0;
          if (k == "DSLOAD") {
            if (lt != "P") { warn_once("DSL" + lt, "*DSLOAD type " + lt + " not supported"); continue; }
            std::string n = resolve_name(d[0], "");
            auto it = m.surfaces.find(n);
            if (it == m.surfaces.end()) throw std::runtime_error("unknown surface '" + d[0] + "' at " + c.where());
            for (auto& f : it->second.faces) lb.faces.push_back({f, v});
          } else {
            int face;
            if (lt == "P") face = 1;
            else if (lt == "PNEG") face = -1;
            else if (lt.size() == 2 && lt[0] == 'P' && std::isdigit((unsigned char)lt[1])) face = lt[1] - '0';
            else { warn_once("DL" + lt, "*DLOAD type " + lt + " not supported"); continue; }
            for (int e : elems_of(d[0], "", c)) lb.faces.push_back({{e, face}, v});
          }
        }
        st.loads.push_back(std::move(lb));
      } else if (contact_keyword(k, c, st)) {
      } else if (k == "OUTPUT") {
        if (c.has("FIELD")) {
          last_output = 1;
          if (!st.output_defined) { st.output = OutputRequest(); st.output_defined = true; }
          st.output.field_intervals = c.geti("NUMBER INTERVAL", st.output.field_intervals);
          st.output.field_time_interval = c.getd("TIME INTERVAL", st.output.field_time_interval);
          st.output.field_frequency = c.geti("FREQUENCY", st.output.field_frequency);
          if (c.has("TIME MARKS")) st.output.time_marks = c.getu("TIME MARKS") == "YES";
          if (st.output.field_intervals == 0 && st.output.field_time_interval == 0 && st.output.field_frequency == 0)
            st.output.field_intervals = 20;
        } else {
          last_output = 2;
          hist_freq = c.geti("FREQUENCY", 0);
          if (c.has("NUMBER INTERVAL")) st.output.history_intervals = c.geti("NUMBER INTERVAL", 200);
          else if (c.has("TIME INTERVAL") && st.period > 0)
            st.output.history_intervals = std::max(1, (int)std::lround(st.period / c.getd("TIME INTERVAL", st.period)));
        }
      } else if (k == "NODE OUTPUT" || k == "ELEMENT OUTPUT") {
        if (last_output == 2 && k == "NODE OUTPUT" && c.has("NSET") && hist_freq > 0) {
          bool rf = false;
          for (auto& d : c.data) for (auto& t : d) { std::string u = upper(t); rf |= u == "RF" || u == "RF1" || u == "RF2" || u == "RF3"; }
          if (rf) { st.output.rf_hist_set = upper(c.get("NSET")); st.output.rf_hist_freq = hist_freq; }
          bool mot = false;
          for (auto& d : c.data) for (auto& t : d) { std::string u = upper(t); mot |= u == "U" || u == "U1" || u == "UR" || u == "UR1"; }
          if (mot && !rf) { st.output.mon_hist_set = upper(c.get("NSET")); st.output.mon_hist_freq = hist_freq; }
        }
        if (last_output == 1)
          for (auto& d : c.data) for (auto& t : d) if (!t.empty()) (k == "NODE OUTPUT" ? st.output.node_vars : st.output.elem_vars).push_back(upper(t));
      } else if (k == "RESTART" || k == "DIAGNOSTICS" || k == "INCREMENTATION OUTPUT" || k == "ENERGY OUTPUT" ||
                 k == "MONITOR" || k == "CONTACT OUTPUT" || k == "PRINT" || k == "NODE PRINT" || k == "EL PRINT" ||
                 k == "NODE FILE" || k == "EL FILE" || k == "CONTROLS" || k == "CONTACT CONTROLS ASSIGNMENT") {
        // not relevant
      } else {
        warn_once("STEPKW" + k, "step option *" + k + " ignored (" + c.where() + ")");
      }
    }
    if (!have_proc) throw std::runtime_error("step " + st.name + " has no procedure");
    if (!nlg_given && st.proc == Step::EXPLICIT) st.nlgeom = true;  // Abaqus/Explicit default
    if (st.proc == Step::STATIC && st.mscale.on) st.mscale.on = false;
    m.steps.push_back(std::move(st));
  }

  // -------------------------------------------------------------- dispatcher
  bool geometry_card(const Card& c, const Scope& sc) {
    const std::string& k = c.keyword;
    if (k == "NODE") card_node(c, sc);
    else if (k == "ELEMENT") card_element(c, sc);
    else if (k == "NSET") card_set(c, sc, true);
    else if (k == "ELSET") card_set(c, sc, false);
    else if (k == "SURFACE") card_surface(c, sc);
    else if (k == "ORIENTATION") card_orientation(c, sc);
    else if (k == "SHELL SECTION" || k == "SOLID SECTION" || k == "BEAM SECTION" || k == "CONNECTOR SECTION") card_section(c, sc);
    else if (k == "TRANSVERSE SHEAR" || k == "TRANSVERSE SHEAR STIFFNESS") {
      if (!pend_sections.empty() && !c.data.empty()) {
        Section& s = pend_sections.back().s;
        s.has_ts = true;
        s.ts_k11 = to_double(c.data[0].at(0), &c);
        s.ts_k22 = c.data[0].size() > 1 && !c.data[0][1].empty() ? to_double(c.data[0][1], &c) : s.ts_k11;
      }
    } else return false;
    return true;
  }

  void process_range(size_t a, size_t b, const Scope& sc) {
    for (size_t i = a; i < b; ++i) {
      const Card& c = cards[i];
      if (!geometry_card(c, sc)) {
        if (c.keyword == "END PART" || c.keyword == "END INSTANCE") continue;
        warn_once("PARTKW" + c.keyword, "part option *" + c.keyword + " ignored (" + c.where() + ")");
      }
    }
  }

  void run(const std::string& path) {
    cards = parse_inp(path);
    Scope root;
    bool in_assembly = false;
    for (size_t i = 0; i < cards.size(); ++i) {
      const Card& c = cards[i];
      const std::string& k = c.keyword;
      // material / behavior sub-options
      if (ctx == CTX_MATERIAL && (k == "DENSITY" || k == "ELASTIC" || k == "PLASTIC" || k == "EXPANSION" || k == "DAMPING" ||
                                  k == "ANISOTROPIC HYPERELASTIC" || k == "HYPERELASTIC" || k == "USER MATERIAL" ||
                                  k == "DEPVAR" || k == "VISCOELASTIC" || k == "CONDUCTIVITY" || k == "SPECIFIC HEAT")) {
        card_material_sub(c);
        continue;
      }
      if (ctx == CTX_BEHAVIOR && k.rfind("CONNECTOR", 0) == 0 && k != "CONNECTOR SECTION" && k != "CONNECTOR BEHAVIOR") {
        card_behavior_sub(c);
        continue;
      }
      if (ctx == CTX_INTERACTION && (k == "FRICTION" || k == "SURFACE BEHAVIOR" || k == "CONTACT DAMPING")) {
        if (k == "FRICTION" && !c.data.empty() && !c.data[0].empty() && !c.data[0][0].empty())
          m.interactions[cur_interaction].friction = to_double(c.data[0][0], &c);
        else if (k != "FRICTION") warn_once("INTKW" + k, "interaction option *" + k + " ignored");
        continue;
      }
      ctx = CTX_NONE;
      if (k == "HEADING") {
        for (auto& d : c.data) for (auto& t : d) m.heading += t + " ";
      } else if (k == "PART") {
        std::string name = upper(c.get("NAME"));
        size_t a = i + 1;
        while (i < cards.size() && cards[i].keyword != "END PART") ++i;
        parts[name] = {a, i};
      } else if (k == "ASSEMBLY") {
        in_assembly = true;
      } else if (k == "END ASSEMBLY") {
        in_assembly = false;
      } else if (k == "INSTANCE") {
        Scope sc;
        std::string iname = upper(c.get("NAME"));
        std::string pname = upper(c.get("PART"));
        sc.prefix = iname + ".";
        sc.inst = (int)m.instances.size();
        m.instances.push_back(c.get("NAME"));
        if (!c.data.empty() && c.data[0].size() >= 3) sc.trans = Vec3(to_double(c.data[0][0]), to_double(c.data[0][1]), to_double(c.data[0][2]));
        if (c.data.size() > 1 && c.data[1].size() >= 7) {
          auto& d = c.data[1];
          sc.rotate = true;
          sc.ra = Vec3(to_double(d[0]), to_double(d[1]), to_double(d[2]));
          sc.rb = Vec3(to_double(d[3]), to_double(d[4]), to_double(d[5]));
          sc.angle = to_double(d[6]);
        }
        auto it = parts.find(pname);
        if (it == parts.end()) throw std::runtime_error("instance of unknown part " + pname + " at " + c.where());
        process_range(it->second.first, it->second.second, sc);
        // cards up to *END INSTANCE belong to the instance
        size_t a = i + 1;
        while (i < cards.size() && cards[i].keyword != "END INSTANCE") ++i;
        process_range(a, i, sc);
      } else if (geometry_card(c, root)) {
        (void)in_assembly;
      } else if (k == "MATERIAL") {
        Material mat;
        mat.name = upper(c.get("NAME"));
        cur_material = (int)m.materials.size();
        m.material_index[mat.name] = cur_material;
        m.materials.push_back(mat);
        ctx = CTX_MATERIAL;
      } else if (k == "CONNECTOR BEHAVIOR") {
        ConnectorBehavior b;
        b.name = upper(c.get("NAME"));
        cur_behavior = (int)m.conn_behaviors.size();
        m.conn_behavior_index[b.name] = cur_behavior;
        m.conn_behaviors.push_back(b);
        ctx = CTX_BEHAVIOR;
      } else if (k == "SURFACE INTERACTION") {
        SurfaceInteraction si;
        si.name = upper(c.get("NAME"));
        cur_interaction = (int)m.interactions.size();
        m.interaction_index[si.name] = cur_interaction;
        m.interactions.push_back(si);
        ctx = CTX_INTERACTION;
      } else if (k == "AMPLITUDE") {
        card_amplitude(c);
      } else if (contact_keyword(k, c, model_contact)) {
      } else if (k == "SECTION CONTROLS") {
        double s = 1.0;
        if (!c.data.empty() && !c.data[0].empty() && !c.data[0][0].empty()) s = to_double(c.data[0][0], &c);
        ControlsCard cc;
        cc.hg_scale = s;
        const std::string hg = c.getu("HOURGLASS");
        cc.enhanced = hg == "ENHANCED";
        cc.relax = hg.empty() || hg == "RELAX STIFFNESS";
        if (hg == "VISCOUS" || hg == "COMBINED") warn_once("HGCTRL" + hg, "*SECTION CONTROLS: HOURGLASS=" + hg + " not supported, relaxed stiffness used");
        cc.face_dt = c.getu("IMPROVED DT METHOD") != "NO";
        section_controls[upper(c.get("NAME"))] = cc;
      } else if (k == "INITIAL CONDITIONS") {
        std::string t = c.getu("TYPE");
        for (auto& d : c.data) {
          if (d.empty() || d[0].empty()) continue;
          auto nodes = nodes_of(d[0], "", c);
          if (t == "TEMPERATURE") for (int n : nodes) m.init_temperature.push_back({n, 0, to_double(d.at(1), &c)});
          else if (t == "VELOCITY") for (int n : nodes) m.init_velocity.push_back({n, to_int(d.at(1), &c) - 1, to_double(d.at(2), &c)});
          else { warn_once("IC" + t, "*INITIAL CONDITIONS TYPE=" + t + " ignored"); break; }
        }
      } else if (k == "BOUNDARY") {
        for (auto& d : c.data) if (!d.empty() && !d[0].empty()) bc_line(c, d, m.model_bcs);
      } else if (k == "STEP") {
        parse_step(i);
      } else if (k == "PHYSICAL CONSTANTS" || k == "PREPRINT" || k == "END PART" || k == "END INSTANCE" || k == "RESTART") {
      } else {
        warn_once("KW" + k, "keyword *" + k + " ignored (" + c.where() + ")");
      }
    }
    assign_sections();
  }

  void assign_sections() {
    for (auto& p : pend_sections) {
      Section s = p.s;
      if (!p.mat_name.empty()) {
        auto it = m.material_index.find(p.mat_name);
        if (it == m.material_index.end()) throw std::runtime_error("unknown material " + p.mat_name + " at " + p.where);
        s.material = it->second;
      }
      if (!p.orient_name.empty()) {
        auto it = m.orient_index.find(p.prefix + p.orient_name);
        if (it == m.orient_index.end()) it = m.orient_index.find(p.orient_name);
        if (it == m.orient_index.end()) throw std::runtime_error("unknown orientation " + p.orient_name + " at " + p.where);
        s.orientation = it->second;
      }
      if (!p.behavior_name.empty()) {
        auto it = m.conn_behavior_index.find(p.behavior_name);
        if (it == m.conn_behavior_index.end()) throw std::runtime_error("unknown connector behavior " + p.behavior_name + " at " + p.where);
        s.behavior = it->second;
      }
      if (!p.controls_name.empty() && section_controls.count(p.controls_name)) {
        const ControlsCard& cc = section_controls[p.controls_name];
        s.hg_scale = cc.hg_scale;
        s.hg_enhanced = cc.enhanced;
        s.hg_relax = cc.relax;
        s.face_dt = cc.face_dt;
      }
      const std::vector<int>* es = find_set(m.elsets, s.elset, p.prefix);
      if (!es) throw std::runtime_error("unknown element set " + s.elset + " for section at " + p.where);
      int idx = (int)m.sections.size();
      m.sections.push_back(s);
      for (int e : *es) m.elems[e].section = idx;
    }
    int missing = 0;
    for (auto& e : m.elems) if (e.section < 0) ++missing;
    if (missing) warn("%d elements have no section and will be ignored", missing);
  }
};

}  // namespace

Model build_model(const std::string& inp_path) {
  Builder b;
  b.run(inp_path);
  return std::move(b.m);
}

}  // namespace abq
