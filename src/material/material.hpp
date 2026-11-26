// Constitutive updates. Stresses are Cauchy, components in a corotational frame.
#pragma once
#include "core/math.hpp"
#include "model/model.hpp"

namespace abq {

// ---- incremental (hypoelastic) isotropic elasticity + J2 plasticity, isotropic hardening
// deps: strain increment (tensor components, Abaqus order 11,22,33,12,13,23), thermal part already removed.
// Returns plastic dissipation increment per unit volume.
// dt: time increment of rate-dependent yield (<= 0: rate-independent).
double j2_update_3d(const Material& m, const Sym6& deps, Sym6& sig, double& ep, double dt = 0);

// Linear kinematic hardening (Prager), backstress alpha (deviatoric tensor components).
double j2_kinematic_3d(const Material& m, const Sym6& deps, Sym6& sig, Sym6& alpha, double& ep);

// Plane stress (sigma33 = sigma13 = sigma23 = 0), in-plane increments d11, d22, d12 (tensor shear).
// Returns plastic dissipation; deps33 receives the thickness strain increment.
double j2_update_ps(const Material& m, double d11, double d22, double d12, double sig[3], double& ep, double& deps33);

// Beam stress state (sigma11, tau12, tau13), engineering shear increments.
double j2_update_beam(const Material& m, double G, double d11, double dg12, double dg13, double& s11, double& t12,
                      double& t13, double& ep);

// ---- superelastic shape memory alloy (Auricchio-Taylor, single martensite fraction xi, linear kinetics)
// eps: total strain (corotational, tensor components), etr: transformation strain, xi: martensite fraction (both
// updated), T: temperature, sig_old: stress at the start of the increment, sig: new stress.
void sma_update(const Material& m, const Sym6& eps, Sym6& etr, double& xi, double T, const Sym6& sig_old, Sym6& sig);
// transformation energy density of a change dxi of the martensite fraction between the stresses s0 and s1
double sma_transformation_energy(const Material& m, const Sym6& s0, const Sym6& s1, double dxi);

// ---- hyperelasticity, total formulation
// Fung (generalized, anisotropic) and neo-Hooke. C in the material (orientation) basis.
// 3D: Cauchy stress (material basis) from the deformation gradient F (material basis).
// fib: Holzapfel-Gasser-Ogden fibre directions (unit vectors, reference configuration, basis of F), 3 per family.
void hyper_cauchy_3d(const Material& m, const Mat3& F, Sym6& sigma, double* energy = nullptr, const double* fib = nullptr);
// Fully integrated elements with a mean volume change (volume ratio jbar of the element, J of the point): isotropic
// models (neo-Hooke, Mooney-Rivlin, Yeoh, Ogden) evaluated with the invariants of the volume-averaged gradient
// Fbar = (jbar / J)^(1/3) F normalised by the point's J, i.e. b* = jbar^(2/3) J^(-4/3) F F^T; with
// tau^ = 2 b* dW/db*: sigma = dev(tau^) / jbar + U'(jbar) I, and pproj = tr(tau^) / (3 jbar) + U'(jbar) is the
// part of the stress that works on the element's mean volume change. Reproduces the linear stiffness of
// Abaqus/Standard C3D8 and C3D8H with these materials (a stiffness 5/3 mu on the variation of the volume change
// between integration points, beyond the B-bar element). Returns false for other models (Fung, HGO).
bool hyper_cauchy_mean_dilatation(const Material& m, const Mat3& F, double jbar, Sym6& sigma, double& pproj);
// Holzapfel-Gasser-Ogden: largest tangent modulus of the fibres in their direction (stable increment estimate).
double hgo_fibre_modulus(const Material& m, const Mat3& F, const double* fib);
// Plane stress membrane: in-plane C (C11, C22, C12) in material basis; solves for C33 with S33 = 0.
// Outputs 2nd Piola-Kirchhoff in-plane components S (S11, S22, S12) and the thickness stretch; for Fung,
// *qexp receives the exponent Q (tracks the stiffening between tangent updates), else -1.
bool hyper_ps(const Material& m, const double C[3], double& lam3, double S[3], double* energy = nullptr, double* qexp = nullptr);

// Effective dilatational modulus (lambda + 2 mu) for stable time increment estimates (initial state).
double dilatational_modulus(const Material& m);
double plane_stress_modulus(const Material& m);

}  // namespace abq
