# Validation

The validation suite is run after every change of the solver (`make validate`, `make parallel`). It contains 270
small decks (216 cases and 54 twins), each running in seconds:

- **elements and deformation modes**: every element type supported by asclepios (C3D8R, C3D8, C3D4, S4R, S3R, B31,
  CONN3D2) under tension, bending, torsion, shear, large rotations, pressure, hyperelasticity, plasticity with load
  reversal, thermal loads, waves, vibrations, rigid rotation;
- **patch tests and robustness**: membrane, plate bending and solid patch tests on distorted elements, very thin shells
  (locking), warped and trapezoidal elements, rotations beyond 360 degrees, near incompressibility, mesh convergence;
- **standard benchmarks**: Scordelis-Lo roof, pinched cylinder, pinched hemisphere, twisted beam, Cook's membrane,
  NAFEMS LE1 and LE10, Taylor-type impact;
- **input-file features**: parts and instances with transformations, `*INCLUDE`, `*PARAMETER`, sets with GENERATE,
  contact on selected surfaces, section controls, transverse shear stiffness, number of section points;
- **solver features**: contact between every pair of element types (shell, solid, beam, rigid and analytical
  surfaces, self-contact, friction, contact pairs, softened contact, exclusions, initial overclosures, crimping
  tools), constraints (TIE, couplings, MPC, equations, rigid bodies), connectors (JOIN, HINGE, BEAM, BUSHING,
  nonlinear, field-dependent, damping, stops), other Abaqus element types (C3D10, C3D10M, C3D8I, C3D8H, C3D4H, C3D6,
  S4, M3D4R, T3D2, B32, box and I beam sections, shell offsets and composite sections, springs, dashpots, point
  masses), loads and conditions (gravity, follower forces, local systems, amplitudes, prescribed accelerations,
  multi-step analyses with OP=NEW, symmetry conditions, initial stresses, temperature and field dependence),
  materials (Mooney-Rivlin, Ogden, Yeoh, orthotropic, Holzapfel-Gasser-Ogden, superelastic nitinol, kinematic
  hardening, rate dependence, damage and element deletion, viscoelasticity);
- **procedures**: every case compared with Abaqus exists in its native procedure and as a twin in the other one
  (static as quasi-static explicit, `_xpl`; explicit as implicit dynamic, `_std`), options specific to the explicit
  procedure (mass scaling, bulk viscosity, fixed increment, damping, adaptive meshing), implicit procedures (arc
  length, frequencies, buckling, stabilization, linear perturbation, `*VISCO`, modal and steady-state dynamics, heat
  transfer, coupled temperature-displacement), transfer between procedures (`*IMPORT`);
- **application models**: valve closure (Fung leaflets, flexible annulus, contact) and stent crimping (elastic-plastic
  beam ring, shrinking tube, beam-to-shell contact);
- **solver checks** (`make parallel`, `PARALLEL.md`): results independent of the number of threads, reproducibility,
  restart equivalence, run time against the model size (with Abaqus on 1 and 4 cpus, double and single precision).

Cases that asclepios does not support yet run in Abaqus only: they measure the coverage and list what remains to be
implemented (status NOT SUPPORTED, reason given in `RESULTS.md`).

## References

- **Published values** (15 cases, `lit_*`, `gen_literature.py`): classic benchmarks of the finite element literature
  (Scordelis-Lo roof, pinched cylinder, pinched hemisphere, twisted beam, Cook's membrane, NAFEMS LE1, clamped plate,
  Lame cylinder, elastica) on refined meshes, compared at the last frame with the converged value of the cited
  publication: nodal displacements, or a stress extrapolated to the published point from the element values next to
  it (quadratic through the first three elements). Tolerance per check (1 to 3 %): remaining discretisation error and
  spread of the published values. The values obtained are listed in the main README; `literature.py` (`make literature`) runs each benchmark on a series of meshes and writes the convergence curves and deformed shapes to `LITERATURE.md` and `literature/`.
- **Closed form** (56 cases): wherever the nodal solution of the discrete model is known by hand, it is the reference
  (ground truth): homogeneous states (tension, simple shear, hyperelastic membranes and blocks, plasticity with
  isotropic or kinematic hardening and load reversal, thermal stress, orthotropic elasticity, initial stress,
  temperature dependence), patch tests, constant moment (exact circular arcs, thin shells), uniform twisting, springs
  and connectors, mass-spring oscillator, rigid motions, slender beams with constraints (MPC, equations, pins, local
  systems, gravity) refined until the discretisation error is below 0.05 %, prescribed accelerations, amplitudes,
  steady conduction. Tolerance: 0.1 % on displacements, rotations, reaction forces, element stresses and energies.
  See `analytic.py`.
- **Reference solver** (165 cases): the same decks run in a commercial solver that reads the same input format
  (implicit or explicit according to the procedure), for the cases without a hand solution. These reference results
  are produced locally by each user who holds a licence of that solver; they and the comparisons with them are not
  distributed with the repository (licence terms of the reference solver).

Output frames are written at the exact requested times in both codes, so that no interpolation error enters the
comparison.

## Compared quantities

For every field output frame: displacements and rotations of the monitored node set `MON`; total reaction force and
moment of the node set `FIX`; for the monitored elements `MONE` (first, middle and last element of the main element
block): largest von Mises stress, von Mises stress of the mean (membrane for shells) stress, equivalent plastic
strain, axial force of beams; temperatures for heat transfer. Over time: the energies ALLIE, ALLKE and ALLWK. Errors
are relative to the peak reference value of each quantity over the whole history (numerical zeros are not compared):

- U, UR, NT: max over frames and monitored nodes of |x − x_ref| / max |x_ref|; rotations below 1e-3 rad are numerical
  zeros (for the case or at a monitored node), and large rotations are compared in the representation of Abaqus (half
  angle in (−π, π]). Explicit cases record U and UR of `MON` at every increment (history output, both sides) and are
  compared on these histories after the moving average described below (rotations below 1 rad);
- RF, RM: same on the total reaction force and moment. Explicit cases also record the total reaction of `FIX` at
  every increment (history output, both sides); their reaction is compared on these histories after a moving average
  over 30 reference increments, i.e. below about a tenth of the highest frequency of the mesh: in impacts and sudden
  loads the reactions and rotations ring at the highest plate modes (200-600 kHz for a 2.5 mm steel mesh), and frames
  sampled every few microseconds catch that ringing with an arbitrary phase (the spectra agree up to those modes);
- S: same on each element quantity, the largest of them is reported (the membrane stress of a shell is compared when
  it reaches 5 % of the peak stress, and not once the section has yielded);
- E: max over time of |E − E_ref| for ALLIE, ALLKE and ALLWK, divided by the largest of these energies, evaluated at
  the history times of the coarser of the two series (the finer one is interpolated).

Also reported, without pass/fail: the ratio of the numbers of increments (asclepios / reference). The Abaqus references
also store the stable time increment (DT), the mass change (DMASS), the hourglass energy (ALLAE) and the number of
iterations of the implicit procedure.

## Coverage

`coverage.py` (run by `make validate`, writes `COVERAGE.md`) lists the keywords read by asclepios that no case run by
asclepios exercises, the keywords used by such cases but not recognised by the reader, and the element types accepted by
asclepios, with those computed with the formulation of another type. Every new feature or fix of the solver comes with a
case, so that these lists stay empty or justified.

## Producing the Abaqus references

On a machine with an Abaqus licence, from the repository root (the cases run one after the other):

```sh
bash validation/abaqus/run_references.sh                 # all cases (ABAQUS=<command> to choose the Abaqus command)
bash validation/abaqus/run_references.sh s4r_pressure    # selected cases
DATACHECK=1 bash validation/abaqus/run_references.sh     # input check only (Abaqus pre-processor), a few minutes
FORCE=1 bash validation/abaqus/run_references.sh         # rerun cases that already have a reference
```

Each case writes `validation/reference/<case>.json` (frames, histories, numbers of increments and iterations, solver
version and precision) and `<case>.time` (wall time); cases named `par_*` also run on 4 cpus, cases marked
`** validation: single` also in single precision, and a case marked `** validation: oldjob=<case>` is run with the
results of that case (`*IMPORT`). The `reference/` folder, `RESULTS.md` and `PARALLEL.md` stay local (ignored by git):
check the licence terms of the reference solver before sharing any of them. The extraction script
(`abaqus/extract_odb.py`) runs with `abaqus python`.

## Running the comparison

```sh
make validate                          # all cases supported by asclepios; writes RESULTS.md (local) and COVERAGE.md
make parallel                          # threads, reproducibility, restart, run time; writes PARALLEL.md
python3 validation/compare.py s4r_rollup b31_rollup --jobs 4
python3 validation/gen_cases.py        # regenerate the decks and cases/manifest.json
```

Status per case: PASS, FAIL, RUN FAILED (asclepios did not finish), NO REFERENCE, NOT SUPPORTED, KNOWN FAIL (documented
limitation, listed at the end of `RESULTS.md`; it does not fail the run). A documented difference on a single quantity
(for example the energy) excludes only that quantity; the others are still checked. The command returns a non-zero
status when a case fails (`--strict`: also when a reference is missing).
