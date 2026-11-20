// In-memory finite element model built from an Abaqus-style input deck (flattened assembly).
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "core/math.hpp"

namespace abq {

enum class ElType { C3D8R, C3D8, C3D8I, C3D4, S4R, S4, S3R, S3RS, B31, T3D2, CONN3D2, UNKNOWN };
const char* eltype_name(ElType t);
ElType eltype_from(const std::string& abaqus_name, int& nnodes);

struct Orientation {
  std::string name;
  Mat3 R = Mat3::identity();  // rows = local axes 1, 2, 3 in global coordinates
};

struct Amplitude {
  std::string name;
  enum Kind { TABULAR, SMOOTH_STEP } kind = TABULAR;
  bool total_time = false;
  std::vector<double> t, v;
  double eval(double step_time, double total_time) const;
};

struct Material {
  std::string name;
  double density = 0;
  bool elastic = false;
  double E = 0, nu = 0;
  std::vector<double> plastic_sy, plastic_ep;  // isotropic hardening table
  double alpha = 0;                            // isotropic thermal expansion
  double damp_alpha = 0, damp_beta = 0;        // Rayleigh
  bool fung = false;                           // generalized Fung, fully anisotropic
  double fung_b[6][6] = {};                    // Voigt (11,22,33,12,13,23), tensor components
  double fung_c = 0, fung_D = 0;
  bool neohooke = false;
  double C10 = 0, D1 = 0;
  bool plastic() const { return !plastic_sy.empty(); }
  double yield(double ep, double* H = nullptr) const;  // flow stress and hardening slope
  // small-strain moduli used for stable time / hourglass / transverse shear defaults
  double shear_modulus() const;
  double bulk_modulus() const;
};

struct ConnectorComponent {
  bool elastic = false;
  bool nonlinear = false;
  double k = 0;
  int ndep = 0;                    // number of field-variable dependencies
  std::vector<double> F, u, fv;    // table rows (force, relative motion, field var 1)
  std::vector<double> kf;          // linear elasticity depending on field variable 1: stiffness per row (with fv)
  bool has_ref = false;
  double ref = 0;                  // constitutive reference (e.g. free length)
  double force(double urel, double fv1, double* dFdu) const;
};
struct ConnectorBehavior {
  std::string name;
  ConnectorComponent comp[6];
};

struct Section {
  enum Kind { SOLID, SHELL, BEAM, CONNECTOR } kind = SOLID;
  std::string elset;
  int material = -1;
  int orientation = -1;
  // shell
  double thickness = 0;
  int nip = 5;
  bool has_ts = false;
  double ts_k11 = 0, ts_k22 = 0;
  // beam
  std::string shape;  // CIRC, RECT, PIPE
  std::vector<double> dims;
  Vec3 n1{0, 0, -1};
  double poisson = -1;
  bool rotary_isotropic = false;  // *BEAM SECTION, ROTARY INERTIA=ISOTROPIC (default EXACT)
  // connector
  enum ConnType { AXIAL, CARTESIAN, JOIN, OTHER } ctype = AXIAL;
  int behavior = -1;
  // hourglass control (from *SECTION CONTROLS): scale of the default stiffness, 1 = default; enhanced = stiffness
  // of the hourglass modes from the physical (enhanced strain) element instead of the default 0.005 G
  double hg_scale = 1.0;
  bool hg_enhanced = false;
  // explicit procedures: the hourglass stiffness relaxes (HOURGLASS=RELAX STIFFNESS, the Abaqus/Explicit default);
  // HOURGLASS=STIFFNESS keeps it constant
  bool hg_relax = true;
  // element length of the explicit stable increment of C3D8R: larger of the uniform-gradient bound and a length
  // built on the largest face (*SECTION CONTROLS, IMPROVED DT METHOD=YES, the default); NO keeps the bound alone
  bool face_dt = true;
};

struct Element {
  ElType type = ElType::UNKNOWN;
  int label = 0;
  int instance = -1;
  std::vector<int> nodes;  // global node indices
  int section = -1;
};

struct FaceRef {
  int elem;  // global element index
  int face;  // shells: +1 SPOS, -1 SNEG, 0 both/unspecified; solids: 1..6 (S1..S6); beams: 0
};
struct Surface {
  std::string name;
  bool node_based = false;
  std::vector<FaceRef> faces;
  std::vector<int> nodes;
};

struct SurfaceInteraction {
  std::string name;
  double friction = 0;
};

// One prescribed quantity on one node degree of freedom (BC, concentrated load, temperature, field).
struct NodalValue {
  int node;
  int dof;  // 0..5 for BC/CLOAD; 0 for temperature/field
  double value;
};
struct FaceValue {
  FaceRef face;
  double value;  // pressure
};

struct LoadBlock {  // one keyword occurrence inside a step
  enum Kind { BOUNDARY, CLOAD, DSLOAD, TEMPERATURE, FIELD, VELOCITY_BC } kind;
  bool op_new = false;
  int amplitude = -1;
  int field_var = 1;
  std::vector<NodalValue> nodal;
  std::vector<FaceValue> faces;
};

struct OutputRequest {
  int field_intervals = 0;         // NUMBER INTERVAL
  double field_time_interval = 0;  // TIME INTERVAL
  int field_frequency = 0;         // FREQUENCY (Abaqus/Standard): a frame every n increments and at the step end
  // TIME MARKS=YES: explicit increments are shortened to reach the field output times exactly; by default (NO, as
  // Abaqus/Explicit) a frame is written at the end of the first increment that reaches or passes the output time
  bool time_marks = false;
  std::vector<std::string> node_vars, elem_vars;
  int history_intervals = 200;
  // *OUTPUT, HISTORY, FREQUENCY=n + *NODE OUTPUT, NSET=set with RF: total reaction of the set every n increments
  std::string rf_hist_set;
  int rf_hist_freq = 0;
  // same with U / UR: displacements and rotations of the nodes of the set every n increments
  std::string mon_hist_set;
  int mon_hist_freq = 0;
};

struct Step {
  std::string name;
  enum Proc { EXPLICIT, STATIC } proc = EXPLICIT;
  bool nlgeom = false;
  int amp_ramp = -1;  // *STEP AMPLITUDE: -1 default (RAMP static, STEP explicit), 0 STEP, 1 RAMP
  double period = 1.0;
  // static incrementation
  double dt0 = 1.0, dtmin = 1e-5, dtmax = 1.0;
  int maxinc = 100;
  // explicit
  double bv1 = 0.06, bv2 = 1.2;
  double dt_scale = 1.0;
  bool element_by_element = false;  // *DYNAMIC, EXPLICIT, ELEMENT BY ELEMENT: no global estimate of the increment
  bool fixed_dt = false;
  double user_dt = 0;
  struct MassScaling {
    bool on = false;
    bool below_min = true;  // TYPE=BELOW MIN (else UNIFORM)
    double dt = 0;
    double factor = 0;
    int frequency = 0;        // every N increments (0 = only at step start)
    int number_interval = 0;
  } mscale;
  // loads
  std::vector<LoadBlock> loads;
  // general contact (explicit)
  bool contact_defined = false;
  bool contact_all_exterior = false;
  std::vector<std::string> contact_surfaces;
  int contact_property = -1;
  // output
  OutputRequest output;
  bool output_defined = false;
};

struct Model {
  std::string heading;
  std::string job;
  // nodes
  std::vector<Vec3> X;
  std::vector<int> node_label;
  std::vector<int> node_instance;
  std::vector<std::string> instances;  // instance names ("" for flat models)
  // elements
  std::vector<Element> elems;
  // sets and surfaces (names upper case, instance-qualified "INST.NAME" for part sets)
  std::map<std::string, std::vector<int>> nsets, elsets;
  std::map<std::string, Surface> surfaces;
  std::vector<Orientation> orients;
  std::map<std::string, int> orient_index;
  std::vector<Material> materials;
  std::map<std::string, int> material_index;
  std::vector<Section> sections;
  std::vector<ConnectorBehavior> conn_behaviors;
  std::map<std::string, int> conn_behavior_index;
  std::vector<Amplitude> amplitudes;
  std::map<std::string, int> amplitude_index;
  std::vector<SurfaceInteraction> interactions;
  std::map<std::string, int> interaction_index;
  // initial / model-level data
  std::vector<NodalValue> init_temperature;
  std::vector<NodalValue> init_velocity;
  std::vector<NodalValue> model_bcs;  // model-level *BOUNDARY (always active)
  std::vector<Step> steps;

  int nn() const { return (int)X.size(); }
  int ne() const { return (int)elems.size(); }
  std::string node_name(int i) const;
  std::string elem_name(int e) const;
};

Model build_model(const std::string& inp_path);

}  // namespace abq
