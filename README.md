# Asclepios — an open-source finite element solver

Asclepios (`asclepios`) is a lightweight finite element solver written in C++17. It reads a subset of the Abaqus keyword
input format (`.inp`), so existing decks can run unchanged, and implements only the modules needed for
nonlinear structural dynamics of soft tissues and slender devices: explicit dynamics, static implicit
analysis, shells, beams, connectors, hyperelastic and elasto-plastic materials, and general contact.

It is a research code: small (≈ 5,300 lines), readable, and easy to modify, rather than a general-purpose
replacement for a commercial package.

## Features

| Area | Supported |
|---|---|
| Input | `*INCLUDE`, `*PARAMETER` (expressions), `*PART` / `*ASSEMBLY` / `*INSTANCE` (translation, rotation), `Assembly.Instance.set` names, `*NODE`, `*ELEMENT`, `*NSET` / `*ELSET` (GENERATE, INSTANCE, `*NSET, ELSET=`), `*SURFACE` (ELEMENT SPOS/SNEG/S1…S6, NODE), `*ORIENTATION` (rectangular) |
| Elements | C3D8R (Flanagan-Belytschko with hourglass control), C3D8 (B-bar), C3D8I (incompatible modes), C3D4, S4R (corotational Belytschko-Lin-Tsay), S4 (S4R kinematics with the stiffness of a physical element and output at four integration points), S3R (C0 triangle, stabilised shear), S3RS (C0 triangle, shear from the edge strains), B31 (corotational Timoshenko beam, CIRC / RECT / PIPE sections), CONN3D2 (AXIAL, CARTESIAN) |
| Materials | isotropic `*ELASTIC`, tabulated isotropic `*PLASTIC` (J2 in 3D, plane stress and beams), `*EXPANSION`, `*DENSITY`, `*DAMPING` ALPHA, `*ANISOTROPIC HYPERELASTIC` FUNG-ANISOTROPIC / FUNG-ORTHOTROPIC, `*HYPERELASTIC` NEO HOOKE, linear and field-dependent tabulated `*CONNECTOR ELASTICITY`, `*CONNECTOR CONSTITUTIVE REFERENCE` |
| Sections | `*SHELL SECTION` (Simpson integration, odd number of points), `*TRANSVERSE SHEAR`, `*BEAM SECTION`, `*SOLID SECTION`, `*CONNECTOR SECTION`, `*SECTION CONTROLS` (hourglass scaling, `HOURGLASS=ENHANCED`) |
| Steps | `*DYNAMIC, EXPLICIT` (NLGEOM by default; `DIRECT USER CONTROL`, `ELEMENT BY ELEMENT`, `SCALE FACTOR`), `*STATIC` (Newton-Raphson with automatic incrementation), `*BULK VISCOSITY`, `*VARIABLE / FIXED MASS SCALING` (BELOW MIN, UNIFORM) |
| Loads | `*BOUNDARY` (displacement, velocity, ENCASTRE…, OP=NEW), `*CLOAD`, `*DSLOAD` / `*DLOAD` P (follower pressure), `*TEMPERATURE`, `*FIELD`, `*AMPLITUDE` (TABULAR, SMOOTH STEP, step or total time); loads propagate between steps |
| Contact | explicit general contact: `*CONTACT INCLUSIONS` (surfaces or ALL EXTERIOR), global `*CONTACT PROPERTY ASSIGNMENT`, Coulomb `*FRICTION`; node-to-facet with rounded edges beyond the facet boundaries (free shell edges rounded within the shell), node-to-beam segment |
| Output | `*OUTPUT, FIELD` (NUMBER INTERVAL / TIME INTERVAL; FREQUENCY in static steps, ignored in explicit ones as Abaqus/Explicit does), `*OUTPUT, HISTORY` (number of intervals) |

Unknown keywords produce a single warning and are ignored.

## Build

Requirements: a C++17 compiler with OpenMP, and [Eigen](https://eigen.tuxfamily.org).
The Makefile targets macOS with Apple clang and Homebrew (`brew install libomp eigen`):

```sh
make            # builds bin/asclepios
make test       # verification suite (needs Python 3 and mpmath)
```

Paths can be overridden, for example `make OMP_DIR=/usr/local/opt/libomp EIGEN=/usr/include/eigen3`.
Other platforms need the OpenMP flags in the Makefile adjusted.

## Usage

```sh
bin/asclepios job=model cpus=8                                # reads model.inp
bin/asclepios model.inp -j 8 -o run1                          # same, job name run1
bin/asclepios job=model cpus=8 save=0.0005                    # restart file model.rst every 0.5 ms
bin/asclepios job=run2 input=model.inp restart=model.rst      # resume (same model, same binary)
```

A restart file holds the complete state (nodal fields, element states, contact, position within the explicit
step); a resumed run reproduces the uninterrupted run bit for bit. `ABQ_KEEP_RESTARTS=1` keeps one copy per
save time (`<job>_t<time>.rst`).

Output files:

- `<job>.pvd` and `<job>_vtk/*.vtu` for ParaView: U, V, RF, UR, NT11 at nodes; S_Mises, S11…S23, PEEQ,
  LE_max, STH, SF per element (shell stresses in the local shell axes: the section orientation if any, otherwise the
  projection of the global x axis on the surface, rotating with the element);
- `<job>_history.csv`: energies ALLKE, ALLIE, ALLWK, ALLVD (bulk viscosity and contact damping), ALLPD, ALLAE,
  ALLDMD, ALLFD, ETOTAL;
- `<job>.sta` (progress) and `<job>.dat` (messages).

Diagnostics on stderr: `ABQ_TRACE_NODE=Inst.123` (node state at each increment), `ABQ_TRACE_ELEM=<id>`
(shell state), `ABQ_TRACE_DT=<id>` (components of a shell's stable increment), `ABQ_CONTACT_DIAG=<n>`
(contact pairs and contact work every n increments). `ABQ_HISTORY_EVERY_INCREMENT=1` writes the energy history at every increment of
an explicit step (the `.sta` keeps the requested interval). `make elemtest` builds `bin/elemtest`; `bin/elemtest deck.inp [nl]` checks the first element of a deck: finite
difference stiffness, symmetry, eigenvalues, exact critical increment against the estimated one.
`tools/cmp_vtk.py` compares a result frame with a reference surface (legacy VTK).

## Numerical methods

- **Explicit dynamics**: central differences with half-step velocities and a lumped mass matrix. The stable
  increment is computed element by element, damping included: Δt = (2/ω)(√(1+ξ²) − ξ). Every stiffness
  (hourglass, drilling, connectors, contact) and every damping term enters this estimate. For models without
  rotational degrees of freedom, a global estimate (power iteration on M⁻¹K, K from directional derivatives of the
  internal forces) replaces it once converged: 0.917 times the converged estimate, kept only when at least 5 %
  above the element bound; `ELEMENT BY ELEMENT` turns it off.
  C3D8R element length: the larger of the uniform-gradient bound and a length built on the largest face
  (`*SECTION CONTROLS, IMPROVED DT METHOD`, YES by default); the initial estimate includes the linear bulk
  viscosity damping, and with initial velocities the quadratic damping at the initial compression rate; the
  reactions written at increment 0 include the linear bulk viscosity of the initial velocities (on the
  stable-increment length). Field frames are written at the end of the first increment past each output time; with
  `TIME MARKS=YES` the output times are reached exactly, the increments before each of them being made equal (the
  remainder spread over up to 10 increments).
- **C3D8R hourglass control**: Flanagan-Belytschko vectors. Explicit procedures: relaxed stiffness (`RELAX STIFFNESS`,
  the default), a standard linear solid on each mode: instantaneous stiffness 8/15·G·V·Σ|∇N|² on the current
  shape of the brick, long-term 0.1 of it, relaxation of the generalized hourglass strains with a time of 1000
  undamped stable increments of the initial element (force increments k·(Δq − 0.9·Δq_v) at the current stiffness);
  during plastic flow the force increments are scaled by max(0.1, H/(3G + H)); `HOURGLASS=STIFFNESS` keeps the
  instantaneous value. Implicit
  procedures: 0.005 G-based with a stiffness per spatial direction (sum of the squared mean gradients in that
  direction). Enhanced control (`HOURGLASS=ENHANCED`, and by default for hyperelastic materials): coupled stiffness of the 12 hourglass modes taken from the
  incompatible-mode brick (9 modes condensed, initial moduli) evaluated on the parallelepiped of the centre
  Jacobian.
- **Shells**: corotational formulation, through-thickness Simpson integration, linear bulk viscosity on membrane
  dilatation (b1·ρ·c·L_m) and curvature rates (b1·ρ·c·L_b·h²/12), with the plate wave speed c (Poisson ratio
  limited to 0.49 at fixed bulk modulus for nearly incompressible materials), L_m the element length of the stable
  increment and L_b the area over the longest midline (quadrilaterals) or the shortest edge (triangles). S4R: one-point membrane and bending; transverse shear from the covariant shear
  strains of the four edges (no shear for any Kirchhoff state), uniform part at the centre with the stiffness
  κGt/(1 + 0.25·10⁻⁴ A/t²) plus one stabilised mode (half the sum of the edge-strain variations, projected on the
  bilinear hourglass vector); hourglass stiffness of the membrane and rotation modes in Flanagan-Belytschko form
  (static: membrane 0.005 G, the documented default, rotations 0.18; explicit: instantaneous stiffnesses 0.2·G·t·A·ΣB²
  and 0.2·G·t³·A·ΣB², relaxed as for C3D8R with 500 undamped element increments), or the stiffness of an
  incompatible-mode element with `HOURGLASS=ENHANCED`; drilling stiffness G t A t²/(24 (t² + 0.1 A)) on the rotation
  of each corner against the spin of the bilinear field there (static) or against the mean spin of the element
  (explicit), with a warping correction. The stable increment uses the element's membrane and bending frequencies, with the tangent modulus of the outer
  section points for hyperelastic materials: membrane, the smaller of the element eigenvalue bound and
  L_m/c·(√(1 + b1²) − b1), L_m² = max(1/ΣB², 0.92·L_e²) (S3R: max(4/3/ΣB², 0.97·L_e²)), L_e from the exact one-point
  membrane eigenvalue; bending side of quadrilaterals, the larger of the uniform transverse shear (exact eigenvalue
  of its 2 × 2 mass-normalised matrix) and of the curvature, rotation hourglass and shear variation modes, which are
  orthogonal to it (thin shells: the increment does not depend on the thickness). S3R: constant membrane strain and curvature; transverse shear of the
  triangle seen as a degenerate quadrilateral whose two merged corners sit at its most isosceles vertex (homogeneous part from the covariant edge
  strains plus a shape-dependent share of their circulation, weak stiffness 0.001 on the circulation, same thin-shell
  factor as S4R), drilling on the corner rotations (explicit: against the element spin, the first node counted twice).
  S3RS: membrane, curvature and drilling as S3R; uniform transverse shear from the three edge strains (the shear field
  whose tangential components at the edge midpoints reproduce them, one point at the centroid, no stiffness on their
  circulation); full shear stiffness k G t on thick elements (t > L/5, L the shortest edge), a penalty
  k G t · t / (2.4 L) times the thin-shell factor on thin ones.
- **Shell rotary inertia** (explicit): the same about the three axes, nodal mass times A/12 + ((L²/16)⁴ + (t²/6)⁴)^¼,
  with L the shortest edge of a triangle or the area over the longest median of a quadrilateral (the distance
  between opposite edges of a parallelogram); S3RS: nodal mass times max(A/4, t²/6); rotations do not control the
  stable increment.
- **B31 rotary inertia**: by default the lumped inertia of the half element, anisotropic and carried by the beam
  frame, ρ·L/2·(I11 + I22, I11, I22) about (axis, n1, n2) per node; the rotational equations of these nodes include
  the gyroscopic moment ω × (J·ω). `ROTARY INERTIA=ISOTROPIC`: ρ·L/2·I_max·max((2 + 3u)/(2 + u), 6u/(1 + u))
  about every axis, u = A·L²/(24·I_min), independent of the elastic constants.
- **B31 bulk viscosity** (explicit): b1·ρ·c·L on the axial strain rate and half of it on the curvature rates
  (moments 0.5·b1·ρ·c·L·I·κ̇); none on transverse shear and torsion; not included in the section force output.
- **B31 kinematics**: section frame rotated by the mean nodal rotation, axial strain from the chord, transverse
  shear strains and forces along the rotated section directions; transverse shear stiffness k·G·A with the
  slenderness compensation 1/(1 + 0.25·L²·A/(12·I)) of each bending plane (I22 for shear along n1, I11 along
  n2); k = 0.85 for RECT, 0.89 for CIRC, 0.53 for PIPE.
- **B31 sections**: PIPE is thin-walled (A = 2πrt, I = πr³t, J = 2I on the mean radius r, section points on the
  mean radius); default section points, 8 circumferential points starting on the section axes; radially R/2 and R
  for CIRC (weights making area and second moment exact); 5 x 5 Simpson for RECT.
- **Plane-stress hyperelasticity**: Newton iteration on the thickness stretch (S33 = 0), with an exact
  derivative for the Fung model.
- **Contact**: penalty method. Each node has a frequency bound 1/Δt_e (stable increment of its stiffest element);
  a contact pair gets the frequency Ω² = max(1/(3Δt²), 1.1·max Ω_i²) of its nodes, nodal springs that keep each
  node stable at 1/Ω, and the softer of the slave and master springs (0.3·m/Δt² between two critical nodes,
  stiffer for light nodes of soft elements), never stiffened during a contact and shared between the contacts of
  a node; the increment is bounded by 1/Ω of the active contacts. Normal viscous damping of 3 % of critical
  (its dissipation is reported in ALLVD; the penalty energy only enters ETOTAL),
  shell half thickness limited to 0.3 times the smallest triangular facet height; free shell edges rounded within
  the shell (their nodes moved back into the shell plane by the half thickness, along the bisector at corners, at
  most half the facet size), so that the contact surface does not extend beyond the edge. The normal force is updated
  incrementally, so it does not jump when the stiffness changes, when a node moves to another facet, or when a contact
  is detected already penetrated. Initial overclosures are removed by a strain-free adjustment (half the
  penetration on each side, reported in the displacements). Coulomb friction.
- **Static implicit**: Newton-Raphson; element tangents by central finite differences of the internal force
  routine, so the same element code serves both procedures. Double precision throughout.

## Validation (`make validate`, `make parallel`)

181 small decks, run after every change of the solver:

- every supported element type (C3D8R, C3D8, C3D4, S4R, S3R, S3RS, B31, CONN3D2) under tension, bending, torsion, shear,
  large rotations, pressure, hyperelasticity, plasticity with load reversal, thermal loads, waves, vibrations and
  rigid rotation;
- classic benchmarks of the finite element literature, compared with their published values (below);
- solver features: contact between every pair of element types, constraints, connectors, other Abaqus element types,
  loads and boundary conditions, multi-step analyses, materials, options specific to each procedure;
- every case in both procedures (implicit and explicit), and parallel runs (results independent of the number of
  threads, speed-up).

Displacements, rotations, reaction forces and energy histories are compared over the whole load history: with the
closed-form solution wherever the nodal solution is known by hand (46 cases, within 0.1 %; most are exact), with
published values for the literature benchmarks (17 cases), and otherwise with reference results that each user
produces on a machine where a licence of the reference solver is available. Reference results of commercial solvers,
and comparisons with them, are not distributed with this repository.

Cases, compared quantities, tolerances and how to produce the references: [validation/README.md](validation/README.md).

### Benchmarks with published values

Refined meshes, linear static unless stated; error relative to the published converged value.

| Benchmark | Element, mesh | Quantity | Asclepios | Published | Error | Source |
|---|---|---|---|---|---|---|
| Scordelis-Lo roof | S4R, 16 × 16 | vertical displacement, middle of the free edge | 0.3043 | 0.3024 | 0.63 % | MacNeal and Harder (1985) |
| Scordelis-Lo roof | S3R, 48 × 48 × 2 | same | 0.2997 | 0.3024 | 0.89 % | MacNeal and Harder (1985) |
| Pinched cylinder with diaphragms | S4R, 32 × 32 | displacement under the load | 1.8205e-5 | 1.8248e-5 | 0.23 % | Belytschko et al. (1985) |
| Pinched hemisphere, 18° hole | S4R, 16 × 16 | displacement under the loads | 0.09383 | 0.094 | 0.18 % | MacNeal and Harder (1985) |
| Twisted beam, load along the clamped width | S4R, 24 × 4 | tip displacement | 1.7737e-3 | 1.754e-3 | 1.12 % | MacNeal and Harder (1985) |
| Twisted beam, load normal to the clamped end | S4R, 24 × 4 | tip displacement | 5.4263e-3 | 5.424e-3 | 0.04 % | MacNeal and Harder (1985) |
| Cook's membrane | S4R, 48 × 48 | vertical displacement, middle of the loaded edge | 24.30 | 23.96 | 1.42 % | Cook (1974) |
| Cook's membrane | S3R, 48 × 48 × 2 | same | 24.09 | 23.96 | 0.52 % | Cook (1974) |
| NAFEMS LE1 elliptic membrane | S4R, 48 × 24 | σyy at point D (extrapolated) | 92.19 MPa | 92.7 MPa | 0.55 % | NAFEMS (1990) |
| NAFEMS LE1 elliptic membrane | S3R, 48 × 24 × 2 | same | 93.97 MPa | 92.7 MPa | 1.36 % | NAFEMS (1990) |
| Clamped square plate, a / t = 100 | S4R, 16 × 16 | centre deflection | 0.013198 | 0.013159 | 0.29 % | Timoshenko and Woinowsky-Krieger (1959) |
| Clamped square plate, a / t = 100 | S3R, 32 × 32 × 2 | same | 0.013056 | 0.013159 | 0.78 % | Timoshenko and Woinowsky-Krieger (1959) |
| Clamped square plate, a / t = 100 | S3RS, 16 × 16 × 2 | same | 0.013238 | 0.013159 | 0.60 % | Timoshenko and Woinowsky-Krieger (1959) |
| Simply supported thick plate, a / t = 10 (Mindlin) | S3RS, 16 × 16 × 2 | centre deflection | 2.1963e-4 | 2.2219e-4 | 1.15 % | Wang, Reddy and Lee (2000) |
| Thick cylinder under internal pressure (Lamé) | C3D8R, 8 × 16 | inner / outer radial displacement | 9.0874e-3 / 5.7829e-3 | 9.0794e-3 / 5.7778e-3 | 0.09 % | Timoshenko and Goodier (1970) |
| Thick cylinder under internal pressure (Lamé) | C3D4, 16 × 32 × 6 | same | 9.0122e-3 / 5.7761e-3 | 9.0794e-3 / 5.7778e-3 | 0.74 % | Timoshenko and Goodier (1970) |
| Elastica, P L² / EI = 5 (large rotation) | B31, 40 | tip displacement across / along the beam | 71.368 / 38.700 | 71.379 / 38.763 | 0.16 % | Mattiasson (1981) |

Mesh-convergence curves and deformed shapes of every benchmark: [validation/LITERATURE.md](validation/LITERATURE.md) (`make literature`). Full references, tolerances and decks: `validation/gen_literature.py` (cases `lit_*`).

![Scordelis-Lo roof](validation/literature/scordelis_lo.png)

## Verification against closed-form solutions (`make test`)

| Case | Result | Reference |
|---|---|---|
| C3D8R / C3D8 bar, quasi-static explicit tension | ALLIE 0.99927 N·mm, ETOTAL/ALLWK < 1e-9 | 1.0 |
| C3D8 bar, linear static | reaction 200.000 N | 200 |
| C3D8R bar, static NLGEOM (λ = 1.1) | S11 = 19,061.6 MPa | E·ln 1.1 = 19,062 |
| S4R cantilever plate | 0.1983 mm | 0.2 (beam theory) |
| S3R cantilever plate | 0.1962 mm | 0.2 |
| B31 beam, CIRC / RECT | 0.021222 / 0.099994 mm | 0.021221 / 0.1 |
| B31 thermal expansion | −0.1000 mm | −0.1 |
| Plate impact on a fixed plate (general contact) | ETOTAL conserved to 1e-4 | energy balance |
| Fung anisotropic shell, biaxial stretch 1.15 × 1.10, S4R and S3R, material axes at 0° and 30° | σ_xx, σ_yy within 1e-6 | Fung energy differentiated with 40-digit arithmetic (mpmath) |

## Roadmap

Features exercised by the validation suite and not supported yet:

- Contact: solid element faces and nodes in general contact, beam-to-beam edge contact, contact pairs
  (`*CONTACT PAIR`), rigid surfaces (R3D4, `*RIGID BODY`), contact in implicit procedures, one-sided SPOS/SNEG
  contact.
- Constraints: `*TIE`, `*COUPLING` (kinematic and distributing), `*MPC`, `*EQUATION`, `*RIGID BODY`, connector types
  other than AXIAL, CARTESIAN and JOIN (HINGE, ...).
- Elements: C3D10M, C3D6, S4, M3D4R, T3D2, B32, SPRINGA, DASHPOTA, MASS.
- Loads and conditions: gravity (`GRAV`), follower forces, nodal transformations (`*TRANSFORM`), periodic amplitudes,
  initial stresses, temperature-dependent elastic properties.
- Materials: Mooney-Rivlin, kinematic hardening, viscoelasticity, superelastic Nitinol, user materials (UMAT/VUMAT).
- Procedures: implicit dynamics (`*DYNAMIC`), arc length (`*STATIC, RIKS`), natural frequencies (`*FREQUENCY`), linear
  buckling (`*BUCKLE`).

Accuracy and performance:

- Shells: drilling rotations and element stresses of static in-plane problems; S3 is read but computed with the S3R
  formulation.
- Dynamics: beam-to-shell impact (energy and reactions within 20 %, beam end rotations 16 %: the B31 stable
  increment, which sets the contact stiffness of the beam nodes, is 15 % below the reference one for this section)
  and self-contact of a coiled strip (first contact within 8 %; the later chatter of the coil on the clamped edge, which
  goes straight into the reactions, is not compared); viscous dissipation of
  shell bulk viscosity on strongly tapered quadrilaterals (up to 15 % high); stability of thin Fung shells closing with contact (valve model).
- Static: predictor from the previous increment and line search (plastic shells need many small increments).
- Energies: ALLIE from material-point quantities (strain energy plus plastic dissipation) instead of the work of nodal
  forces.
- Performance: vectorisation of the constitutive updates; parallel linear solver for the implicit procedure (speed-up
  1.4 on 8 threads); fewer threads on small models (OpenMP overhead).

## Relationship with Abaqus

Asclepios is not affiliated with, endorsed by or connected to Dassault Systèmes. Abaqus is a trademark of
Dassault Systèmes or its subsidiaries. Asclepios is an independent implementation and contains no Abaqus code.
It reads a subset of the Abaqus input format so that existing models can be run.

## Support and professional services

Questions and bug reports: GitHub issues. Product-specific simulations and custom validation studies can be
arranged: see [SUPPORT.md](SUPPORT.md).

## Citing

If you use Asclepios in published work, please cite it: [CITATION.cff](CITATION.cff) (GitHub button "Cite this
repository").

## License

Copyright 2026 NonameAccountMax. Licensed under the Apache License, Version 2.0: see [LICENSE](LICENSE) and
[NOTICE](NOTICE).
