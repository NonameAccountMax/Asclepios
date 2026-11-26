// In-memory finite element model built from an Abaqus-style input deck (flattened assembly).
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "core/math.hpp"

namespace abq {

enum class ElType { C3D8R, C3D8, C3D8I, C3D6, C3D4, C3D10, C3D10M, DC3D8, S4R, S4, S3R, S3RS, B31, B32, T3D2, CONN3D2, R3D4, R3D3, SFM3D4, SFM3D3, M3D4R, M3D4, M3D3, SPRINGA, DASHPOTA, MASS, UNKNOWN };
const char* eltype_name(ElType t);
inline bool is_rigid_type(ElType t) { return t == ElType::R3D4 || t == ElType::R3D3; }
// surface elements (SFM3D4R, SFM3D3): facets without stiffness or mass, moved by their boundary conditions
inline bool is_surface_type(ElType t) { return t == ElType::SFM3D4 || t == ElType::SFM3D3; }
// membranes (M3D4R, M3D4, M3D3): in-plane forces only, no rotational degrees of freedom
inline bool is_membrane_type(ElType t) { return t == ElType::M3D4R || t == ElType::M3D4 || t == ElType::M3D3; }
// shells and membranes (elements of the shell group)
inline bool is_shell_type(ElType t) {
  return t == ElType::S4R || t == ElType::S4 || t == ElType::S3R || t == ElType::S3RS || is_membrane_type(t);
}
inline bool is_solid_type(ElType t) {
  return t == ElType::C3D8R || t == ElType::C3D8 || t == ElType::C3D8I || t == ElType::C3D6 || t == ElType::C3D4 ||
         t == ElType::C3D10 || t == ElType::C3D10M;
}
ElType eltype_from(const std::string& abaqus_name, int& nnodes);

struct Orientation {
  std::string name;
  Mat3 R = Mat3::identity();  // rows = local axes 1, 2, 3 in global coordinates
  std::vector<Vec3> local_dirs;  // LOCAL DIRECTIONS (e.g. fibre directions), in the local axes
  // axis of the additional rotation (second data line, default the local 1-axis): on shells and membranes the
  // following axis in the cyclic order 1 -> 2 -> 3 -> 1 is projected onto the surface as the material 1-direction
  int rot_axis = 1;
};

struct Amplitude {
  std::string name;
  enum Kind { TABULAR, SMOOTH_STEP, PERIODIC, DECAY } kind = TABULAR;
  bool total_time = false;
  std::vector<double> t, v;
  // PERIODIC: a = A0 + sum_n A_n cos(n w (t - t0)) + B_n sin(n w (t - t0)) for t >= t0 (A0 before);
  // DECAY: a = A0 + A exp(-(t - t0) / td) for t >= t0 (A0 before)
  double w = 0, t0 = 0, a0 = 0, ad = 0, td = 1;
  std::vector<double> an, bn;
  double eval(double step_time, double total_time) const;
};

struct Material {
  std::string name;
  double density = 0;
  bool elastic = false;
  double E = 0, nu = 0;
  std::vector<double> plastic_sy, plastic_ep;  // isotropic hardening table
  double alpha = 0;                            // isotropic thermal expansion
  std::vector<double> el_T, el_E, el_nu;       // temperature-dependent elasticity (rows E, nu, T), empty otherwise
  void elastic_at(double T, double& E, double& nu) const;
  double damp_alpha = 0, damp_beta = 0;        // Rayleigh
  bool fung = false;                           // generalized Fung, fully anisotropic
  double fung_b[6][6] = {};                    // Voigt (11,22,33,12,13,23), tensor components
  double fung_c = 0, fung_D = 0;
  bool neohooke = false;
  double C10 = 0, D1 = 0;
  // other isotropic hyperelastic forms: Mooney-Rivlin (C10, C01, D1), Yeoh (C10, C20, C30, D1, D2, D3), Ogden
  // (mu_i, alpha_i, D_i)
  bool mooney = false, yeoh = false, ogden = false;
  double C01 = 0, C20 = 0, C30 = 0, D2 = 0, D3 = 0;
  std::vector<double> og_mu, og_alpha, og_D;
  // orthotropic linear elasticity (*ELASTIC, TYPE=ENGINEERING CONSTANTS): stiffness in the material axes, tensor
  // strains (11 22 33 12 13 23; shear terms already doubled)
  bool ortho = false;
  double Cortho[6][6] = {};
  // linear kinematic hardening (*PLASTIC, HARDENING=KINEMATIC): yield stress and slope of the backstress
  bool kinematic = false;
  double kin_C = 0;
  // Holzapfel-Gasser-Ogden (*ANISOTROPIC HYPERELASTIC, HOLZAPFEL): C10, D, k1, k2, kappa, fibre families
  bool hgo = false;
  double hgo_C10 = 0, hgo_D = 0, hgo_k1 = 0, hgo_k2 = 0, hgo_kappa = 0;
  int hgo_nfib = 0;
  bool hyper() const { return fung || neohooke || hgo || mooney || yeoh || ogden; }
  // superelastic shape memory alloy (material named ABQ_SUPER_ELASTIC..., *USER MATERIAL constants): austenite and
  // martensite moduli, transformation strain, stress-temperature slopes, transformation stresses at T0
  bool sma = false;
  double sma_EA = 0, sma_nuA = 0, sma_EM = 0, sma_nuM = 0, sma_epsL = 0, sma_CL = 0, sma_sLs = 0, sma_sLe = 0,
         sma_T0 = 0, sma_CU = 0, sma_sUs = 0, sma_sUe = 0, sma_sCLs = 0, sma_epsV = 0;
  // linear viscoelasticity (*VISCOELASTIC, TIME=PRONY) on isotropic *ELASTIC: Prony terms (shear ratio g, bulk ratio
  // k, relaxation time tau). Once the input is read, E and nu are the instantaneous moduli (*ELASTIC gives the
  // long-term ones unless MODULI=INSTANTANEOUS).
  std::vector<double> pr_g, pr_k, pr_tau;
  bool el_instantaneous = false;
  bool visco() const { return !pr_tau.empty(); }
  // rate-dependent yield (*RATE DEPENDENT, TYPE=POWER LAW): flow stress times 1 + (plastic strain rate / D)^(1/n)
  double rate_D = 0, rate_n = 1;
  bool rate_dependent() const { return rate_D > 0; }
  // heat transfer: isotropic conductivity and specific heat (temperature-independent)
  double conductivity = 0, specific_heat = 0;
  // ductile damage (*DAMAGE INITIATION, CRITERION=DUCTILE): equivalent plastic strain at the onset of damage as a
  // function of the stress triaxiality; *DAMAGE EVOLUTION with linear softening: equivalent plastic displacement at
  // failure (TYPE=DISPLACEMENT) or fracture energy (TYPE=ENERGY, u_f = 2 G_f / yield stress at the onset)
  std::vector<double> dmg_eps, dmg_eta;
  double dmg_uf = 0, dmg_gf = 0;
  bool damage() const { return !dmg_eps.empty(); }
  bool damage_evolves() const { return damage() && (dmg_uf > 0 || dmg_gf > 0); }
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
  double damping = 0;              // *CONNECTOR DAMPING: linear coefficient on the relative velocity
  bool stop = false;               // *CONNECTOR STOP: admissible range of the relative position (AXIAL: length)
  double stop_lo = -1e300, stop_hi = 1e300;
  double force(double urel, double fv1, double* dFdu) const;
};
struct ConnectorBehavior {
  std::string name;
  ConnectorComponent comp[6];
};

struct Section {
  enum Kind { SOLID, SHELL, BEAM, CONNECTOR, DISCRETE, SURFACE } kind = SOLID;
  double area_density = 0;  // *SURFACE SECTION, DENSITY: mass per unit area of surface elements
  // discrete elements (*SPRING, *DASHPOT, *MASS): stiffness or table (force, relative displacement), damping
  // coefficient, point mass
  double spring_k = 0, dashpot_c = 0, point_mass = 0;
  std::vector<double> spring_F, spring_u;
  std::string elset;
  int material = -1;
  int orientation = -1;
  // solid section of trusses: cross-sectional area (first data value, default 1)
  double area = 1;
  // shell
  double thickness = 0;
  int nip = 5;
  // OFFSET: position of the reference surface (the nodes) from the mid-surface, in thicknesses, positive towards
  // the positive face (SPOS = 0.5, SNEG = -0.5)
  double offset = 0;
  // COMPOSITE: layers from the negative face up (thickness, section points, material; orientation angle unused)
  struct Layer { double t = 0; int nip = 3; int material = -1; };
  std::vector<Layer> layers;
  bool membrane = false;  // *MEMBRANE SECTION
  bool has_ts = false;
  double ts_k11 = 0, ts_k22 = 0;
  // beam
  std::string shape;  // CIRC, RECT, PIPE, BOX, I
  std::vector<double> dims;
  Vec3 n1{0, 0, -1};
  double poisson = -1;
  bool rotary_isotropic = false;  // *BEAM SECTION, ROTARY INERTIA=ISOTROPIC (default EXACT)
  // connector
  enum ConnType { AXIAL, CARTESIAN, JOIN, HINGE, BEAM_LINK, BUSHING, OTHER } ctype = AXIAL;
  int behavior = -1;
  // hourglass control (from *SECTION CONTROLS): scale of the default stiffness, 1 = default; enhanced = stiffness
  // of the hourglass modes from the physical (enhanced strain) element instead of the default 0.005 G
  double hg_scale = 1.0;
  bool hg_enhanced = false;
  // no HOURGLASS= given: enhanced control for hyperelastic materials (the Abaqus/Explicit and Abaqus/Standard default)
  bool hg_default = true;
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
  bool heat = false;       // temperature degrees of freedom (DC3D8: heat transfer only; C3D8T: coupled, mechanics of C3D8)
};

struct FaceRef {
  int elem;  // global element index
  int face;  // shells: +1 SPOS, -1 SNEG, 0 both/unspecified; solids: 1..6 (S1..S6); beams: 0
};
// nodes of one face of an element (local indices, closed loop; shells, rigid and surface elements: all their nodes,
// beams and trusses: both nodes); returns their number, 0 for an invalid face
int face_local_nodes(ElType t, int face, int nen, int* out);

struct Surface {
  std::string name;
  bool node_based = false;
  std::vector<FaceRef> faces;
  std::vector<int> nodes;
};

// *TIE: the nodes of the slave surface follow the master surface (interpolated at their closest point)
struct SurfaceTie {
  std::string name, slave, master;  // surface names
  bool adjust = true;               // ADJUST=YES: slave nodes moved onto the master surface (strain free)
  double tolerance = -1;            // POSITION TOLERANCE (< 0: default)
  bool rotations = true;            // rotations tied too when both sides carry them (NO ROTATION: off)
};

// *COUPLING (KINEMATIC / DISTRIBUTING) and *RIGID BODY: coupling nodes tied to a reference node
struct NodeCoupling {
  std::string name;
  int ref = -1;                  // reference node
  std::string surface;           // coupling nodes: nodes of this surface...
  std::vector<int> nodes;        // ...or these nodes
  enum Kind { KINEMATIC, DISTRIBUTING } kind = KINEMATIC;
  bool dof[6] = {true, true, true, true, true, true};  // constrained degrees of freedom
  int weighting = 0;             // DISTRIBUTING: UNIFORM 0, LINEAR 1, QUADRATIC 2, CUBIC 3
  std::vector<int> rigid_elems;  // *RIGID BODY: elements moving rigidly (mass only, not evaluated)
  std::string analytic;          // *RIGID BODY, ANALYTICAL SURFACE=: surface discretised into rigid facets
};
// analytical rigid surface (*SURFACE, TYPE=CYLINDER): profile in the local (x, y) plane swept along the generator
struct AnalyticSurface {
  std::string name;
  Vec3 a, b, c;  // origin, point on local x, point on the generator (local y = x cross generator)
  struct Seg { int kind; double x, y, cx, cy; };  // kind 0 START, 1 LINE, 2 CIRCL (end point, centre)
  std::vector<Seg> segs;
};
// *EQUATION: sum coef_i u(node_i, dof_i) = 0, the first term is eliminated
struct LinearConstraint {
  std::vector<int> node, dof;
  std::vector<double> coef;
};
// *MPC BEAM / TIE / PIN: node a follows node b
struct NodePairMpc {
  enum Kind { BEAM, TIE, PIN } kind = BEAM;
  int a = -1, b = -1;
};

struct SurfaceInteraction {
  std::string name;
  double friction = 0;
  // *SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE: HARD (penalty), or a softened law p(h) of the overclosure h (negative:
  // clearance): TABULAR (points (p, h), linear beyond the last one), EXPONENTIAL (zero at the clearance c0, p0 at
  // contact), LINEAR (slope k)
  enum Law { HARD, TABULAR, EXPONENTIAL, LINEAR } law = HARD;
  std::vector<double> tab_p, tab_h;
  double c0 = 0, p0 = 0, slope = 0;
  // clearance at which the pressure starts (>= 0)
  double clearance() const;
  // contact pressure at the overclosure h and its slope dp/dh
  double pressure(double h, double* dpdh) const;
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
  enum Kind { BOUNDARY, CLOAD, DSLOAD, TEMPERATURE, FIELD, VELOCITY_BC, ACCELERATION_BC } kind;
  bool op_new = false;
  int amplitude = -1;
  int field_var = 1;
  bool follower = false;  // *CLOAD, FOLLOWER: the force rotates with the node
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
  enum Proc { EXPLICIT, STATIC, FREQUENCY, BUCKLE, MODAL_DYNAMIC, STEADY_STATE, HEAT } proc = EXPLICIT;
  // linear perturbation step (*STEP, PERTURBATION; *FREQUENCY, *BUCKLE, *STEADY STATE DYNAMICS always): response about
  // the base state, loads and boundary conditions of the step only, not carried over
  bool perturbation = false;
  int neig = 0;                      // *FREQUENCY / *BUCKLE: number of eigenvalues requested
  double eig_max = 0;                // *BUCKLE: largest eigenvalue of interest (0: none); *FREQUENCY: highest frequency
  double eig_min = 0;                // *FREQUENCY: lowest frequency of interest (cycles per time)
  struct ModalDamping { int first = 1, last = 1; double ratio = 0; };
  std::vector<ModalDamping> modal_damping;  // *MODAL DAMPING (fraction of critical damping per mode range)
  double ssd_f0 = 0, ssd_f1 = 0, ssd_bias = 3; int ssd_points = 20;  // *STEADY STATE DYNAMICS frequency range
  // arc length (*STATIC, RIKS): initial, minimum, maximum arc length increments and total arc length (period);
  // end at a load proportionality factor or a displacement of a node
  bool riks = false;
  double riks_lpf_max = 0;           // 0: none
  int riks_node = -1, riks_dof = 0;
  double riks_u = 0;
  bool nlgeom = false;
  int amp_ramp = -1;  // *STEP AMPLITUDE: -1 default (RAMP static, STEP explicit), 0 STEP, 1 RAMP
  double period = 1.0;
  // static incrementation (also implicit dynamics)
  double dt0 = 1.0, dtmin = 1e-5, dtmax = 1.0;
  // implicit dynamics (*DYNAMIC without EXPLICIT): HHT-alpha integration in the implicit driver
  bool implicit_dynamic = false;
  double hht_alpha = -0.05;
  // time-dependent material response in a static step: *VISCO (creep and relaxation over the step time); *STATIC
  // responds instantaneously, or with the long-term moduli under LONG TERM
  bool visco = false, long_term = false;
  // heat transfer (*HEAT TRANSFER: proc HEAT) and fully coupled temperature-displacement (*COUPLED TEMPERATURE-
  // DISPLACEMENT: proc STATIC with coupled): steady state or transient (backward Euler), largest temperature change
  // of an increment (DELTMX, 0: none)
  bool coupled = false, heat_steady = false;
  double deltmx = 0;
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
  // *CONTACT EXCLUSIONS: surface pairs (second name empty: self-contact of the first)
  std::vector<std::pair<std::string, std::string>> contact_exclusions;
  // *CONTACT PAIR (explicit, kinematic hard contact): first and second surface, interaction, weighting factor of the
  // first surface as master (WEIGHT; < 0: default, balanced 0.5 between deformable surfaces)
  struct ContactPair { std::string first, second; int interaction = -1; double weight = -1; };
  std::vector<ContactPair> contact_pairs;
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
  std::vector<SurfaceTie> surface_ties;
  std::vector<NodeCoupling> couplings;
  std::map<std::string, AnalyticSurface> analytic_surfaces;
  std::vector<LinearConstraint> equations;
  std::vector<NodePairMpc> mpcs;
  // initial / model-level data
  std::vector<NodalValue> init_temperature;
  std::vector<NodalValue> init_velocity;
  struct ElemStress { int elem; double s[6]; };
  std::vector<ElemStress> init_stress;  // *INITIAL CONDITIONS, TYPE=STRESS (global components 11 22 33 12 13 23)
  std::vector<NodalValue> model_bcs;  // model-level *BOUNDARY (always active)
  std::map<int, Mat3> transforms;      // *TRANSFORM (rectangular): node -> local axes (rows) of its dofs
  std::vector<Step> steps;

  int nn() const { return (int)X.size(); }
  int ne() const { return (int)elems.size(); }
  std::string node_name(int i) const;
  std::string elem_name(int e) const;
};

Model build_model(const std::string& inp_path);

}  // namespace abq
