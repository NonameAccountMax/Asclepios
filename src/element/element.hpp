// Element group interface. One group per element family; elements are processed independently
// (incremental, corotational formulation shared by the explicit and the implicit drivers).
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "core/binio.hpp"
#include "model/model.hpp"

namespace abq {

// Per-element output channels written to VTK (cell data).
enum OutChannel { OUT_MISES = 0, OUT_S11, OUT_S22, OUT_S33, OUT_S12, OUT_S13, OUT_S23, OUT_PEEQ, OUT_LEMAX, OUT_STH, OUT_SF, OUT_SDV, N_OUT };
const char* out_channel_name(int c);

struct ElemIn {
  const double* X = nullptr;      // nen*3 reference coordinates
  const double* xo = nullptr;     // nen*3 coordinates at the start of the increment
  const double* du = nullptr;     // nen*6 increments (u, rotation vector)
  const double* temp0 = nullptr;  // nen temperatures at start (nullptr: none)
  const double* temp1 = nullptr;  // nen temperatures at end
  const double* fv1 = nullptr;    // nen field variable 1 at end
  double dt = 0;                  // time increment
  double dt_mat = 0;              // time increment seen by time-dependent materials (0: instantaneous response)
  double bv1 = 0, bv2 = 0;        // bulk viscosity (explicit only)
  bool dynamic = true;            // explicit: rate-dependent artificial terms on
  bool nlgeom = true;
  bool initial = false;           // increment-0 pass on the initial velocities (linear bulk viscosity only)
};

struct ElemOut {
  double* f = nullptr;     // nen*6 internal force, global axes
  double dt_stable = 1e30; // stable increment for unit mass scale
  double e_visc = 0;       // viscous (bulk viscosity) dissipation increment
  double e_plastic = 0;    // plastic dissipation increment
  double e_art = 0;        // artificial (hourglass/drill) energy increment
  // internal energy at the material points minus the work of the element forces over the increment: the thermal
  // strain work is removed, and for linear elasticity with J2 plasticity the change of the elastic strain energy
  // plus the plastic dissipation replaces the trapezoidal sigma : d(eps) (as ALLIE = ALLSE + ALLPD)
  double e_corr = 0;
  bool failed = false;     // material point did not converge
};

class ElementGroup {
 public:
  ElType type = ElType::UNKNOWN;        // formulation computed by the group
  ElType input_type = ElType::UNKNOWN;  // element type of the input (C3D6 in explicit procedures: computed as C3D8R)
  int nen = 0;
  bool rot = false;                // uses rotational dofs
  std::vector<int> ids;            // global element indices
  std::vector<int> conn;           // size()*nen global node indices
  std::vector<double> msf;         // mass scaling factor (density multiplier) per element
  std::vector<double> dt_last;     // last stable increment estimate (unscaled)
  std::vector<float> out;          // size()*N_OUT output channels
  const Model* model = nullptr;
  bool physical_rotary = false;    // implicit dynamics: physical rotary inertia (no explicit stability scaling)

  virtual ~ElementGroup() = default;
  int size() const { return (int)ids.size(); }
  virtual void setup(const Model& m) = 0;
  // lumped translational mass and isotropic rotary inertia per element node (unscaled density)
  virtual void mass(int e, const double* X, double* mt, double* mr) const = 0;
  // rotary inertia per element node that mass proportional damping (alpha) multiplies; default: the lumped one, mr
  virtual void damping_rotary(int e, const double* mr, double* mrd) const {
    for (int I = 0; I < nen; ++I) mrd[I] = mr[I];
  }
  // anisotropic rotary inertia, added to the isotropic part of mass(): symmetric tensor (xx, yy, zz, xy, xz, yz)
  // carried by every node of element e in the current configuration (x: current nodal coordinates), unscaled
  // density. Returns false when the element has none.
  virtual bool rotary_tensor(int e, const double* x, double* J6) const { return false; }
  // internal forces for the increment; when commit is true the new state is stored
  virtual void compute(int e, const ElemIn& in, ElemOut& out, bool commit) = 0;
  // material Rayleigh damping coefficient (mass proportional)
  virtual double alpha_damping(int e) const { return 0; }
  // element with nodal pressure unknowns in the slot of the first rotation (C3D4H, implicit procedures)
  virtual bool hybrid(int e) const { return false; }
  // contact thickness (shells/beams) for general contact
  virtual double contact_thickness(int e) const { return 0; }
  // strain energy currently stored (for energy balance; optional)
  virtual double strain_energy(int e) const { return 0; }
  // restart: element state (derived groups append their own state after the base vectors)
  virtual void save(BinOut& o) const { o.vec(msf); o.vec(dt_last); o.vec(out); }
  virtual void load(BinIn& in) { in.vec(msf); in.vec(dt_last); in.vec(out); }
};

std::vector<std::unique_ptr<ElementGroup>> make_element_groups(const Model& m);

}  // namespace abq
