#include <map>

#include "core/log.hpp"
#include "element/groups.hpp"

namespace abq {

const char* out_channel_name(int c) {
  static const char* n[N_OUT] = {"S_Mises", "S11", "S22", "S33", "S12", "S13", "S23", "PEEQ", "LE_max", "STH", "SF", "SDV"};
  return n[c];
}

std::vector<std::unique_ptr<ElementGroup>> make_element_groups(const Model& m) {
  std::map<ElType, std::unique_ptr<ElementGroup>> g;
  for (int e = 0; e < m.ne(); ++e) {
    const Element& el = m.elems[e];
    if (el.section < 0) continue;
    const Section& s = m.sections[el.section];
    ElType t = el.type;
    auto it = g.find(t);
    if (it == g.end()) {
      std::unique_ptr<ElementGroup> grp;
      switch (t) {
        case ElType::C3D8R: case ElType::C3D8: case ElType::C3D8I: case ElType::C3D4: grp = make_solid_group(); break;
        case ElType::S4R: case ElType::S4: case ElType::S3R: case ElType::S3RS: grp = make_shell_group(); break;
        case ElType::B31: grp = make_beam_group(); break;
        case ElType::CONN3D2: grp = make_connector_group(); break;
        case ElType::T3D2: grp = make_truss_group(); break;
        default: continue;
      }
      grp->type = t;
      it = g.emplace(t, std::move(grp)).first;
    }
    bool kind_ok = (s.kind == Section::SOLID && (t == ElType::C3D8R || t == ElType::C3D8 || t == ElType::C3D8I || t == ElType::C3D4 || t == ElType::T3D2)) ||
                   (s.kind == Section::SHELL && (t == ElType::S4R || t == ElType::S4 || t == ElType::S3R || t == ElType::S3RS)) ||
                   (s.kind == Section::BEAM && t == ElType::B31) || (s.kind == Section::CONNECTOR && t == ElType::CONN3D2);
    if (!kind_ok) throw std::runtime_error("element " + m.elem_name(e) + ": section type does not match element type");
    if (s.kind != Section::CONNECTOR && s.material < 0) throw std::runtime_error("element " + m.elem_name(e) + ": section without material");
    it->second->ids.push_back(e);
    for (int n : el.nodes) it->second->conn.push_back(n);
  }
  std::vector<std::unique_ptr<ElementGroup>> out;
  for (auto& kv : g) {
    auto& grp = kv.second;
    grp->nen = grp->conn.size() / grp->ids.size();
    grp->msf.assign(grp->size(), 1.0);
    grp->dt_last.assign(grp->size(), 1e30);
    grp->out.assign(grp->size() * N_OUT, 0.f);
    grp->setup(m);
    out.push_back(std::move(grp));
  }
  return out;
}

}  // namespace abq
