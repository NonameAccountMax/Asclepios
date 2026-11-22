// Shared analysis state for the explicit and implicit drivers.
#pragma once
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "element/element.hpp"
#include "model/model.hpp"

namespace abq {

struct Prescribed {
  int node, dof;
  double v0, v1;  // value at step start / given value
  int amp;        // amplitude index, -1 default, -2 constant (propagated)
  bool velocity = false;
};
struct FaceLoad {
  FaceRef face;
  double v0, v1;
  int amp;
};

struct LoadState {
  std::vector<Prescribed> bcs, cloads, temps, fields;
  std::vector<FaceLoad> press;
  bool ramp = false;  // default amplitude of the current step
  double period = 1;
  double value(double v0, double v1, int amp, double st, double tt, const Model& m) const;
  double value(const Prescribed& p, double st, double tt, const Model& m) const { return value(p.v0, p.v1, p.amp, st, tt, m); }
};

class GeneralContact;
class OutputWriter;

struct Energies {
  double kinetic = 0, internal = 0, external = 0, viscous = 0, plastic = 0, artificial = 0, damping = 0,
         contact = 0, friction = 0, mass_scaling_work = 0;
  double total(double ke) const { return ke + internal + viscous + damping + contact - external; }
};

struct Analysis {
  Model& m;
  std::string job;
  std::vector<std::unique_ptr<ElementGroup>> groups;
  int nn = 0;
  std::vector<double> X;                          // 3nn reference coordinates
  std::vector<double> u, v, a, du;                // 6nn
  std::vector<double> fint, fext, fcon, rf;       // 6nn
  std::vector<double> mt, mr, cdamp;              // nn masses, rotary inertia, alpha*m
  std::vector<double> temp, temp_prev, fv1;       // nn
  std::vector<char> node_rot;                     // node carries rotational stiffness/inertia
  // anisotropic rotary inertia (beams with the default exact inertia): per node symmetric tensor (xx, yy, zz, xy,
  // xz, yz) added to mr, in the current configuration; node_aniso flags the nodes that carry one
  std::vector<double> jrot;
  std::vector<char> node_aniso;
  bool any_aniso = false;
  void rotary_tensors(bool with_increment);       // jrot at X + u (+ du)
  // rotational accelerations of an anisotropic node: J a = r - w x (J w) on the free rotational dofs; r and the
  // result on the fixed dofs are the reactions (Abaqus sign: rf = -(r - J a))
  void rotary_accelerations(int n, const double* r, const double* w, double* acc, double* react) const;
  double rotary_energy(int n, const double* w) const;  // 0.5 w.(mr + J) w
  std::vector<char> fixed;                        // 6nn prescribed dof (current step)
  // linear dependencies between dofs (links.cpp): u[dof] = sum c u[term], from surface ties and JOIN connectors
  struct LinkTerm { int dof; double c, mc; };    // coefficient, mass share (c^2 / interpolation weight)
  struct NodeTie {                                // slave node following master facet nodes
    int slave;
    std::vector<int> masters;
    std::vector<double> w_pt, w_avg;              // closest-point interpolation; surface average (implicit)
    bool rotations;
    Vec3 off_pt, off_avg;                         // X_slave - interpolated master point (after ADJUST)
  };
  std::vector<NodeTie> node_ties;
  std::vector<std::pair<int, int>> joins;         // CONN3D2 JOIN: the translations of node b are those of node a
  struct DofLink { int dof, first, count; };
  std::vector<DofLink> links;                     // active in the current step (dependent dof not prescribed)
  std::vector<LinkTerm> link_terms;
  std::vector<int> link_of;                       // 6nn: index in links, -1 for an independent dof
  std::vector<int> link_term_dofs;                // dofs appearing as terms
  std::vector<double> link_r, link_m, res;        // explicit work arrays (res: residual forces of all dofs)
  void tie_surfaces();                            // setup, before the element groups: pairing, ADJUST
  void build_constraints();                       // setup: couplings, rigid bodies, equations, MPCs
  std::vector<std::pair<int, std::vector<LinkTerm>>> constraint_links;
  std::vector<char> elem_rigid;                   // per model element: part of a rigid body (not evaluated)
  std::vector<std::vector<char>> group_rigid;     // the same per group element
  void activate_links(bool implicit);             // start of a step, after the boundary conditions
  void gather_links(std::vector<double>& f) const;  // forces of the dependent dofs moved onto their terms
  void spread_links(std::vector<double>& x) const;  // dependent values from their terms
  void link_accelerations(const std::vector<double>& r);  // explicit: accelerations and reactions with the links
  std::set<std::string> warned_keys;
  void warn_once_flag(const std::string& key, const std::string& msg);
  std::vector<std::vector<double>> efor;          // per group element force buffers
  std::vector<std::vector<double>> spring_k;      // per group: spring stiffness of connector-like elements
  std::vector<double> node_om2, node_ks;          // nodal frequency bound work arrays
  bool has_temp = false;
  LoadState loads;
  Energies en;
  double total_time = 0;
  int total_increments = 0;
  std::unique_ptr<GeneralContact> contact;
  std::unique_ptr<OutputWriter> writer;
  int nthreads = 1;
  // restart: periodic save (total time interval, 0 = off) and resume from a restart file
  double save_interval = 0, next_save = 0;
  std::string restart_in;  // file to resume from ("" = none)
  int resume_step = -1;    // step index stored in restart_in
  void save_state(BinOut& o) const;
  void load_state(BinIn& in);
  bool nlgeom = true;

  explicit Analysis(Model& model);
  ~Analysis();
  void setup();
  void lump_masses();
  void begin_step(const Step& st, int step_index, bool ramp_default);
  void update_fixed(const Step& st);
  // element loop: computes fint from u (start) + du; returns min stable dt and accumulates energies
  struct ElemSum {
    double dt_min = 1e30;
    int dt_elem = -1;
    double dt_solid = 1e30, dt_other = 1e30;  // smallest increment of the solid elements and of all the others
    double e_visc = 0, e_plastic = 0, e_art = 0;
    int failed = 0;
  };
  ElemSum eval_elements(double dt, double bv1, double bv2, bool dynamic, bool commit, std::vector<double>& f);
  bool initial_pass = false;  // eval_elements: increment-0 viscous forces of the initial velocities
  std::vector<char> group_mask;  // eval_elements: groups evaluated (empty: all)
  double solid_dt_scale = 1.0;   // global estimate of the increment: ratio applied to the solid elements
  void external_forces(double st, double tt, std::vector<double>& f);  // at configuration X+u+du
  void temperatures(double st, double tt);
};

}  // namespace abq
