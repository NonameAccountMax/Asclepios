#pragma once
#include "solver/analysis.hpp"

#include <Eigen/Sparse>

namespace abq {
void run_explicit_step(Analysis& A, const Step& st, int si);
void run_static_step(Analysis& A, const Step& st, int si);
// arc length (*STATIC, RIKS), riks.cpp
void run_riks_step(Analysis& A, const Step& st, int si);
// linear perturbation procedures (linear.cpp): *FREQUENCY, *BUCKLE, *MODAL DYNAMIC, *STEADY STATE DYNAMICS, *STATIC
void run_linear_step(Analysis& A, const Step& st, int si);
// heat transfer (thermal.cpp): *HEAT TRANSFER steps; temperatures of an increment (dt = 0: steady state) into A.temp
// from A.temp_prev, prescribed temperatures at step time st / total time tt; returns the largest change of a free node
void run_heat_step(Analysis& A, const Step& st, int si);
double heat_increment(Analysis& A, double dt, double st, double tt);

// ---- shared by the implicit procedures (static.cpp)
struct DofMap {
  std::vector<int> eq;  // 6nn -> equation or -1 (fixed, inactive or dependent dofs)
  int neq = 0;
};
// active dofs (translations of nodes with elements, rotations of shell and beam nodes) not fixed; dependent dofs of
// linear constraints eliminated
DofMap make_dofmap(const Analysis& A);
// finite-difference step of the element tangents (from the shortest element edge)
double fd_step(const Analysis& A);
// tangent stiffness of all elements at the current state and increment A.du, in the equations of dm (constraint terms
// distributed); nlgeom of A
void assemble_tangent(Analysis& A, const DofMap& dm, double dt, double hfd, std::vector<Eigen::Triplet<double>>& trip);
// lumped masses per dof (physical rotary inertia), dependent dofs lumped onto their terms
void dof_masses(Analysis& A, std::vector<double>& mdof);
// connector stops (stops.cpp): presence of supported stops; kinematic correction of the accelerations a so that the
// predicted configuration X + u + dpred respects them (c: displacement per unit acceleration correction), forces in fk
bool has_connector_stops(const Model& m);
void connector_stops(const Analysis& A, double c, std::vector<double>& dpred, std::vector<double>& a, std::vector<double>& fk);
}  // namespace abq
