// Constitutive updates. Stresses are Cauchy, components in a corotational frame.
#pragma once
#include "core/math.hpp"
#include "model/model.hpp"

namespace abq {

// ---- incremental (hypoelastic) isotropic elasticity + J2 plasticity, isotropic hardening
// deps: strain increment (tensor components, Abaqus order 11,22,33,12,13,23), thermal part already removed.
// Returns plastic dissipation increment per unit volume.
double j2_update_3d(const Material& m, const Sym6& deps, Sym6& sig, double& ep);

// Plane stress (sigma33 = sigma13 = sigma23 = 0), in-plane increments d11, d22, d12 (tensor shear).
// Returns plastic dissipation; deps33 receives the thickness strain increment.
double j2_update_ps(const Material& m, double d11, double d22, double d12, double sig[3], double& ep, double& deps33);

// Beam stress state (sigma11, tau12, tau13), engineering shear increments.
double j2_update_beam(const Material& m, double G, double d11, double dg12, double dg13, double& s11, double& t12,
                      double& t13, double& ep);

// ---- hyperelasticity, total formulation
// Fung (generalized, anisotropic) and neo-Hooke. C in the material (orientation) basis.
// 3D: Cauchy stress (material basis) from the deformation gradient F (material basis).
void hyper_cauchy_3d(const Material& m, const Mat3& F, Sym6& sigma, double* energy = nullptr);
// Plane stress membrane: in-plane C (C11, C22, C12) in material basis; solves for C33 with S33 = 0.
// Outputs 2nd Piola-Kirchhoff in-plane components S (S11, S22, S12) and the thickness stretch; for Fung,
// *qexp receives the exponent Q (tracks the stiffening between tangent updates), else -1.
bool hyper_ps(const Material& m, const double C[3], double& lam3, double S[3], double* energy = nullptr, double* qexp = nullptr);

// Effective dilatational modulus (lambda + 2 mu) for stable time increment estimates (initial state).
double dilatational_modulus(const Material& m);
double plane_stress_modulus(const Material& m);

}  // namespace abq
