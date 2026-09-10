# QuakeCore — Phase 9L

**Research release, September 9, 2026.** Start with [the job interface](docs/JOB_FORMAT.md), the [soil and pile guide](docs/SOIL_FOUNDATIONS.md), the [Phase 9L verification report](docs/PHASE9L_SOIL_PILES_REVIEW.md), the [steel component guide](docs/STEEL_COMPONENTS.md), and the [wall element guide](docs/WALL_ELEMENTS.md). The JSON runner executes plane-frame, wall, steel-component, and vertical-pile NRHA and IDA. Independent OpenSees comparisons cover frame, wall, steel, damper, and the defined bilinear pile fixture. Berkeley runs remain synthetic-proxy studies; physical component validation and ASCE 41-23/AISC 342-22 compliance are not established.

Research prototype for a high-throughput nonlinear response-history analysis engine for buildings.

## Current objective

Test the central solver and compiler hypotheses before writing GPU code:

1. **Compile structural topology once** so constraints, sparse locations, nonlinear deformation operators, and output maps are frozen before time stepping.
2. **Avoid repeated global sparse factorizations** when building tangent changes are localized.
3. Reuse model/dt-dependent factorization work across an entire ground-motion suite.

The repository now contains a working CPU reference path from the original shear-building proof through compiled 2D frames and true 6-DOF/node 3D space-frame models with general MPC condensation.

**Phase 9D.5 architecture:** RC code-model iteration can now reuse one compiled 2D/3D frame while replacing only a state-layout-compatible nonlinear material field. Native RC demand recording also supports the objective finite-rotation corotational 3D path. See `PHASE9D5_IMMUTABLE_TOPOLOGY_COROTATIONAL_DEMAND.md`.

## Implemented

### Numerical core

- C++20 sparse matrix infrastructure (CSC)
- SuperLU sparse direct factorization/solve backend
- **same-pattern SuperLU refactorization** reusing column ordering and the elimination tree while permitting new row pivots
- exact eager Woodbury solver for `A + B diag(Δk) Bᵀ`
- exact **lazy Woodbury solver** that computes `A⁻¹b_j` only when hinge `j` first becomes nonlinear
- **sparse nonlinear update basis**: candidate hinge operators store only their actual nonzero DOF coefficients rather than a dense `n x m` matrix
- active-rank compression
- prepared Newmark/Woodbury workspace reusable across many records
- shared-memory multi-record suite runner with one prepared solver/cache per CPU worker
- Newmark average-acceleration NRHA
- full Newton, same-pattern direct refactorization, modified Newton, and exact low-rank solution strategies
- robust Newmark driver with rollback, backtracking line search, recursive time-step subdivision, and direct fallback
- material-neutral flat nonlinear state boundary plus generalized constitutive-bank API between component laws and the transient solver
- independent research **peak-oriented IMK-family** material bank with monotonic/cyclic deterioration state and event diagnostics
- parameter-driven **ASCE 41-style degrading hinge** with A-B-C-D-E-F envelope, IO/LS/CP tracking, explicit cyclic branch state, E lateral-loss and F effective/gravity-loss states
- **virgin-elastic active-front fast path** that bypasses full IMK reversal/deterioration logic until a component leaves its initial elastic domain
- model-specific Woodbury/direct **crossover calibration** on the actual compiled sparse matrix
- prepared **adaptive Newton linear solver**: rank-zero direct triangular solve, calibrated Woodbury range, same-pattern direct fallback beyond crossover
- exact story-block Schur factorization and substructured lazy Woodbury research backend (serial CPU path benchmarked but not default)
- Craig-Bampton/Guyan research reduction wrapper for static/modal/NRHA trade studies
- nonlinear active-front diagnostics (total vs non-initial-tangent component evaluations)
- termination taxonomy separating completed response, configured collapse/limit criteria, numerical failure, and initial instability; a configured criterion is not an independently verified physical-collapse diagnosis
- scalable sparse initial-tangent positive-definiteness certification plus near-zero mechanism tracking
- global robust-driver energy ledger (input/internal/damping/kinetic/balance)
- finite-rotation, energy-based **3D corotational reference element** with objective rigid-body kinematics and Euler-stability mesh-convergence tests
- exact generalized localized tangent solver for `A + U C(u) U^T`, including state-dependent geometric-update representation and cost-based direct fallback
- **IDA/collapse-suite runner** with record-level CPU parallelism, adaptive collapse-scale refinement, censoring status, and collapse-mechanism classification

### Compiled 2D/3D frame models

- 2D arbitrary node IDs and 3 DOF/node (`UX`, `UY`, `RZ`)
- fixed and equal-DOF constraint elimination with mass accumulation
- rigid-floor horizontal master/slave constraints
- 2D Euler-Bernoulli frame coordinate transformation
- in-plane **MVLEM wall element** with discrete concrete/steel fibers and independent nonlinear shear spring
- in-plane **SFI-MVLEM wall element** with plane-stress macro-panels, local zero-transverse-stress equilibrium, and exact condensed coupled tangent
- explicit elastic, Concrete01-compatible, Steel01-compatible, layered, and fixed-angle RC research material laws; full FSAM is not implemented
- wall self-mass, wall/frame/MPC assembly, wall response histories, NRHA and IDA integration
- constant-preload beam-column geometric stiffness for baseline P-Delta benchmarking
- state-dependent updated P-Delta reference path with consistent axial-lateral Newton Jacobian
- lumped nodal translational/rotational masses
- zero-length/generalized rotational bilinear hinges
- two-node 2D steel beams and columns with two statically condensed nonlinear end hinges
- bilinear, ASCE 41-style, or IMK-family end-hinge material selection with explicit parameter provenance
- zero-length rotational panel-zone surrogate with its own response classification and recorder
- small-displacement axial bilinear BRB component
- axial power-law viscous damper with a consistent velocity tangent and explicit regularization for sublinear exponents
- directional translational and relative-rotational zero-length soil springs and dashpots between coincident nodes
- rollback-safe asymmetric compression/uplift spring with permanent gap state
- explicit bilinear p-y, t-z, and q-z tributary-force adapters with mandatory provenance
- vertical elastic pile-line compilation with fixed soil nodes and depth-wise force/deformation recorders
- automatically generated hinge deformation basis `B`
- **compiled CSC scatter locations** for nonlinear tangent changes
- 3D 6 DOF/node (`UX`, `UY`, `UZ`, `RX`, `RY`, `RZ`) and 12-DOF space-frame elements
- 3D axial, torsional, and biaxial Euler-Bernoulli member stiffness
- 3D user reference vectors for local-axis orientation
- **general sparse linear MPC transformation** `u_full = T q`
- true horizontal rigid diaphragms with master `RZ` translation coupling
- non-diagonal condensed generalized mass and eccentric torsional inertia
- arbitrary generalized 3D spring deformation operators and rotational-vector hinges
- compiled low-rank basis
- horizontal base excitation
- Rayleigh damping
- roof and interstory response reducers

### Verification/tooling

- analytical 2D and biaxial 3D cantilever stiffness tests
- analytical MVLEM axial/flexural/shear compliance, rigid-body, symmetry, rollback, unit-scaling, and nonlinear tangent checks
- 740-step MVLEM/SFI-MVLEM cyclic static comparisons and 1,200-step three-story dynamic comparisons against OpenSeesPy 3.8.0
- wall JSON contract tests and synthetic wall NRHA/IDA integration fixtures
- 740-step steel-member, panel-zone, BRB, and power-law damper comparisons against OpenSeesPy 3.8.0
- linear viscous-damper SDOF transient-history comparison against OpenSeesPy 3.8.0
- integrated steel-frame NRHA full/direct-fallback equivalence, time-step sensitivity, nonlinear activation, and synthetic IDA checks
- p-y/t-z/q-z conversion, asymmetric contact/gap, and 200-element beam-on-Winkler closed-form checks
- full-history single-pile comparison against an independently assembled OpenSeesPy 3.8.0 model
- geometric-stiffness softening test
- rigid constraint/mass accumulation test
- true rigid-diaphragm kinematics, generalized mass coupling, and torsional base-load tests
- arbitrary hinge-basis/vector-operator compilation tests
- same-pattern direct-refactorization equivalence tests
- robust failed-step/subdivision regression test
- full-Newton vs Woodbury NRHA history equivalence tests
- prepared-suite equivalence test
- release and ASan/UBSan test configurations
- nonlinear shear, sparse-grid, frame, record-suite, rank-crossover, and component-mode-reduction benchmarks

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Sanitizers:

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DQUAKE_ENABLE_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

Benchmarks:

```bash
./build/low_rank_bench 45 45 100 5
./build/nrha_bench 120 1200 6
./build/frame_nrha_bench 20 2 600
./build/frame_suite_bench 20 2 12 600
./build/frame3d_nrha_bench 10 2 2 300 4.0
./build/frame3d_diaphragm_bench
./build/rank_crossover_bench 3000 300 20
./build/reduction_bench 8 12 300 8
./build/frame3d_suite_bench 10 8 300 4
./build/imk_active_front_bench 10 2 2 300 3.0
./build/calibration_bench 20 2 2 4
./build/substructure_bench 40 4 100 200
./build/stability_bench
./build/collapse_bench 6 20 1 2
./build/corotational_bench 1000
./build/corotational_buckling_bench
./build/geometric_update_bench
./build/ida_bench 4
```

## Current Phase 7 result

Phase 7 adds a finite-rotation 3D corotational reference path, exact generalized `U C U^T` localized tangent updates, and IDA/collapse-suite infrastructure. The corotational element is objective under arbitrary 3D rigid-body motion, recovers the linear frame at small response, and converges to Euler cantilever stability with mesh refinement (8 elements ~0.33% high; 16 elements ~0.08% high).

The generalized geometric-update formulation reproduces the Phase 6 updated-P-Delta tangent exactly, but the benchmark also establishes an important rejection condition: when all member geometric states evolve, the active generalized update dimension can rival or exceed the structural DOF count, making global Woodbury slower than same-pattern sparse refactorization. The model-aware fallback therefore remains authoritative; future geometric acceleration should be block-local/domain-parallel.

The new IDA engine runs records/scales in parallel, refines first-collapse brackets logarithmically, preserves physical-collapse vs numerical-failure vs censoring status, and records the collapse mechanism. A representative 8-record synthetic suite completed 89 nonlinear runs with zero numerical failures; four CPU workers reduced wall time from ~0.466 s to ~0.165 s while producing identical collapse brackets and causes. See [`docs/PHASE7_COROT_IDA.md`](docs/PHASE7_COROT_IDA.md).

## Current Phase 6 result

Phase 6 adds collapse-fidelity infrastructure: a parameter-driven ASCE 41-style degrading hinge, explicit IO/LS/CP/E/F event tracking, state-dependent P-Delta with a consistent Newton Jacobian, physical-collapse vs numerical-failure termination, energy-balance diagnostics, and a scalable initial-stability certification step.

A six-story degrading-hinge benchmark now progresses through acceptance limits and post-capping response without equating CP exceedance with collapse. Under one updated-P-Delta record, all three linear strategies reach IO/LS/CP/E/F at steps 140/152/168/195/219 and terminate at step 219 by the configured F effective/gravity-loss criterion. Under a stronger record, all three identify tangent-instability collapse at step 174 before any E/F loss. The full, same-pattern, and Woodbury reference paths have identical Newton iteration counts and peak responses in these checks.

The robust driver now certifies the starting tangent as positive definite using sparse symmetric elimination before integrating. Representative one-time certification costs are ~0.013 s at 780 DOF, ~0.026 s at 1,560 DOF, and ~0.059 s at 3,120 DOF in the current environment. Small validation matrices can additionally be cross-checked with dense LAPACK eigenvalues.

With state-dependent P-Delta enabled, geometric tangent changes currently force exact direct fallbacks; Phase 6 deliberately does not claim a Woodbury acceleration for that mode. With a fixed geometric baseline, the same strong collapse benchmark retains exact low-rank behavior (one baseline factorization versus 125 same-pattern and 256 fresh factorizations in the representative run).

The degrading-hinge cyclic implementation required an explicit committed branch state to eliminate a real near-collapse line-search pathology. After that rewrite the formerly failing record completes with no time-step subdivision, demonstrating that the failure was constitutive discontinuity rather than a Woodbury artifact. See [`docs/PHASE6_COLLAPSE.md`](docs/PHASE6_COLLAPSE.md).

## Current Phase 5 result

Phase 5 adds heterogeneous IMK-family hinges, a real active-front fast path, per-model solver calibration/adaptive switching, and exact story-block Schur substructuring. Representative internal CPU benchmarks (same physics and convergence tolerances in all compared paths) are:

### 10-story, 2x2-bay 3D IMK frame

```text
DOF:                         780
IMK hinges:                  240
steps:                       300

fresh full SuperLU:          ~1.61 s
same-pattern SuperLU:        ~0.70 s
lazy Woodbury:               ~0.068 s
Woodbury vs same-pattern:    ~10.4x

Newton iterations:           721 for all three
peak roof response:          0.075914 for full/Woodbury
max active tangent rank:     66
IMK virgin-elastic fast path ~85.2% of component evaluations
```

A prepared adaptive solver uses the direct factor at rank zero and Woodbury through the model-calibrated profitable range; this record ran in ~0.064 s after preparation. Its one-time calibration cost (~0.15 s in this small benchmark) is intended to be amortized across a ground-motion suite.

### 20-story, 2x2-bay 3D IMK frame

```text
DOF:                         1,560
IMK hinges:                  480
steps:                       250

fresh full SuperLU:          ~2.64 s
same-pattern SuperLU:        ~0.95 s
lazy Woodbury:               ~0.161 s
Woodbury vs same-pattern:    ~5.9x

Newton iterations:           572 for all three
peak roof response:          0.373069 for full/Woodbury
max active tangent rank:     204
IMK virgin-elastic fast path ~91.0% of component evaluations
```

Calibration on this actual 1,560-DOF matrix found Woodbury strongly faster for tangent-refresh ranks through 320 sampled directions and slower at rank 480. The strong-motion record only reached rank 204, so the adaptive policy correctly stayed in the low-rank path during nonlinear response rather than using a universal rank-fraction rule.

Exact story-block Schur substructuring is also implemented and verified against monolithic sparse solves and monolithic Woodbury. It reduces a 40-story research matrix to independent local interior blocks plus a sparse interface but is generally slower in the current **serial CPU** implementation. It is therefore retained as a domain-parallel/GPU backend and is not promoted to the default local CPU solver.

These remain **internal research benchmarks, not claims against OpenSees, xara, Perform-3D, LS-DYNA, or GPU hardware**. The IMK-family implementation is independently written and has not yet passed external OpenSees cyclic-history parity tests. See [`docs/PHASE5_IMK_ADAPTIVE.md`](docs/PHASE5_IMK_ADAPTIVE.md).

## Governing decomposition

For a fixed baseline effective tangent `A` and localized stiffness changes,

`K = A + B D Bᵀ`

where `D = diag(Δk)`. The exact Woodbury solve is

`K⁻¹r = y - W (I + D G)⁻¹ D Bᵀ y`

with

- `y = A⁻¹r`
- `W = A⁻¹B`
- `G = BᵀA⁻¹B`

The lazy implementation does not construct all of `W` and `G` up front. It grows the cached influence basis only for hinges whose tangent actually changes.

## Scope and limitations

This is research software, not design software. Current limitations include:

- state-dependent updated P-Delta exists as a consistent small-rotation reference path, but a true 3D corotational beam-column formulation is still pending;
- compiled 3D frames can mix bilinear and independent research IMK-family hinges, but IMK cyclic parity against an external reference implementation is still pending;
- panel zones and BRBs currently use engineer-supplied idealized backbones; geometry-derived AISC panel-zone properties, BRB fatigue/fracture, cumulative strain, and connection models are not implemented;
- viscous dampers are memoryless power-law devices; Maxwell series compliance, relief valves, gaps, stroke/force limits, and manufacturer qualification are not implemented;
- steel-member axial yielding, P-M interaction, local/lateral-torsional buckling, fracture, distributed plasticity, and evolving gravity/P-Delta are not implemented;
- no soil spring, pile, or foundation-failure component is exposed by the JSON job interface;
- modal analysis exists as a dense validation utility; production sparse eigensolvers are not yet implemented;
- external verification currently uses OpenSeesPy 3.8.0; physical specimen and full-building blind validation remain pending;
- same-pattern SuperLU is implemented, but a second symmetric-indefinite/LDL CPU backend would strengthen the baseline;
- no CUDA/cuDSS execution in this environment;
- no fiber/distributed-plasticity path.

## Immediate roadmap

1. Add nonlinear translational/rotational soil springs with gapping, uplift, compression-only, p-y/t-z/q-z backbone families, dashpots, and explicit depth/provenance metadata.
2. Compile pile lines as beam-column meshes with distributed soil springs and tributary-length conversion; verify free-field, single-pile, and group-constraint fixtures before any lateral pile claim.
3. Add external cyclic/history parity for the independent IMK-family bank and physical component/specimen benchmarks for steel, wall, BRB, panel-zone, and damper models.
4. Add steel axial yielding/P-M interaction, gravity-state transfer, evolving P-Delta, cyclic deterioration/fracture, and auditable ASCE 41/AISC parameter providers.
5. Extend low-rank updates to coupled member/device blocks where the measured cost model predicts a benefit; same-pattern direct refactorization remains authoritative otherwise.
6. Mirror representative full buildings in OpenSees/xara for modal, static, cyclic, NRHA, and collapse/IDA comparison.
7. Only after the CPU validation gate, begin CUDA/cuDSS; preserve local CPU record parallelism as the normal workstation path.

## Phase 8 RC1 — localized FSC column validation

Phase 8 adds a mechanics-based flexure-shear-critical column path for the UC Berkeley three-story nonductile RC benchmark. The accepted RC1 model uses 39-in clear flexible column geometry with rigid offsets, equilibrium-consistent initial-stress P-Delta, and separate degrading series shear springs with cyclic deterioration and E-state tracking. The experimental moving-P flexural wrapper remains disabled pending a true moving-surface return-mapping implementation.

Under the current deterministic 1.52g development proxy, RC1 completes with `T1=0.48 s`, peak base shear `41.00 kip`, and peak story drifts `[3.271, 4.588, 2.365]%`. The local B1 shear-degradation initiation state (`3.175%` drift, `9.792 kip` shear, `26.15 kip` compression) is close to the published Test-1 state, while global base shear remains too high. FullFactorization and SamePattern histories are identical in the release run. The proxy 5%-damped pseudo-spectral acceleration is `Sa(0.48s)=1.900g` and `Sa(0.34s)=2.423g`.

This is still **PRE-VALIDATION / PROXY INPUT**; final external DT1 parity requires the actual recorded shake-table acceleration. See [`docs/PHASE8_RELEASE_NOTES.md`](docs/PHASE8_RELEASE_NOTES.md).

## Phase 8.1 mechanism-validation note

Phase 8.1 separates dynamic inertial reaction from component-level first-story restoring shear, decouples ASCE41 B-C physical hardening from the 100x zero-length hinge penalty stiffness, and adds a displacement-controlled Berkeley-frame pushover. The pushover peaks at **23.56 kip**, compared with the independent NIST/FEMA P-2018 `Vy = 23.3 kip`; the proxy NRHA component-recovered first-story shear peaks at **28.80 kip**, versus 27.9 kip Perform3D and 29.8 kip measured, while the legacy inertial reaction is 38.53 kip and is no longer treated as mechanism capacity. See `PHASE8_1_RC1_RELEASE_NOTES.md` and `validation/uc_berkeley_3story/quakecore_phase8_1_release_result.json`.

## Phase 9A excitation/parity path

Phase 9 freezes the v0.8.1-rc1 structural model and moves validation effort to ground-motion provenance. See `PHASE9A_ALPHA1_RELEASE_NOTES.md` and `validation/uc_berkeley_3story/phase9/README.md`. The exact Llolleo Component 100 source is metadata-gated before any source-motion reconstruction is allowed.

## Phase 9D.6 prepared model-comparison harness

Phase 9D.6 adds a prepared robust NRHA driver and an order-invariant material-field comparison harness on top of the immutable compiled topology introduced in Phase 9D.5. Named variants restore the captured baseline before applying their nonlinear-material delta, allowing research, NIST/ASCE 41-17, and resolved production ASCE 41-23 / ACI 369 model fields to be compared with identical topology, excitation, solver controls, and recorders. The harness records aligned story response, hinge force-deformation-tangent histories, native column P/V histories, IO/LS/CP/E/F event timing, instantaneous tangent modes at story-peak states, and condensed story tangent stiffness. Prepared solver contexts are reused only while the exact effective-initial matrices remain unchanged; changes to hinge initial stiffness automatically invalidate the preparation. Production ASCE 41-23 / ACI 369 coefficients remain externally resolved and are not inferred by QuakeCore. See `PHASE9D6_PREPARED_COMPARISON_HARNESS.md`.
