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
  // Abaqus/Explicit computes the 6-node prism as a reduced-integration brick whose first node is doubled (bricks
  // 1 1 2 3 4 4 5 6: uniform strain, hourglass forces and stable increment of that brick, equal nodal masses; single
  // C3D6 elements with prescribed hourglass motions, four triangle shapes)
  bool xpl = !m.steps.empty();
  for (const auto& st : m.steps) xpl = xpl && st.proc == Step::EXPLICIT;
  for (int e = 0; e < m.ne(); ++e) {
    const Element& el = m.elems[e];
    if (el.section < 0) continue;
    const Section& s = m.sections[el.section];
    ElType t = el.type;
    auto it = g.find(t);
    if (it == g.end()) {
      std::unique_ptr<ElementGroup> grp;
      switch (t) {
        case ElType::C3D8R: case ElType::C3D8: case ElType::C3D8I: case ElType::C3D6: case ElType::C3D4:
        case ElType::C3D10: case ElType::C3D10M: grp = make_solid_group(); break;
        case ElType::S4R: case ElType::S4: case ElType::S3R: case ElType::S3RS:
        case ElType::M3D4R: case ElType::M3D4: case ElType::M3D3: grp = make_shell_group(); break;
        case ElType::B31: grp = make_beam_group(); break;
        case ElType::B32: grp = make_beam3_group(); break;
        case ElType::CONN3D2: grp = make_connector_group(); break;
        case ElType::T3D2: grp = make_truss_group(); break;
        case ElType::SPRINGA: case ElType::DASHPOTA: case ElType::MASS: case ElType::SFM3D4: case ElType::SFM3D3:
          grp = make_discrete_group(); break;
        default: continue;
      }
      grp->type = t == ElType::C3D6 && xpl ? ElType::C3D8R : t;
      grp->input_type = t;
      it = g.emplace(t, std::move(grp)).first;
    }
    bool kind_ok = (s.kind == Section::SOLID && (is_solid_type(t) || t == ElType::T3D2)) ||
                   (s.kind == Section::SHELL && is_shell_type(t) && s.membrane == is_membrane_type(t)) ||
                   (s.kind == Section::BEAM && (t == ElType::B31 || t == ElType::B32)) || (s.kind == Section::CONNECTOR && t == ElType::CONN3D2) ||
                   (s.kind == Section::DISCRETE && (t == ElType::SPRINGA || t == ElType::DASHPOTA || t == ElType::MASS)) ||
                   (s.kind == Section::SURFACE && is_surface_type(t));
    if (!kind_ok) throw std::runtime_error("element " + m.elem_name(e) + ": section type does not match element type");
    if (s.kind != Section::CONNECTOR && s.kind != Section::DISCRETE && s.kind != Section::SURFACE && s.material < 0) throw std::runtime_error("element " + m.elem_name(e) + ": section without material");
    it->second->ids.push_back(e);
    if (t == ElType::C3D6 && it->second->type == ElType::C3D8R) {
      static const int brick[8] = {0, 0, 1, 2, 3, 3, 4, 5};
      for (int k : brick) it->second->conn.push_back(el.nodes[k]);
    } else {
      for (int n : el.nodes) it->second->conn.push_back(n);
    }
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
