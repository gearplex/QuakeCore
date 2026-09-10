# Next Development Milestones

## Milestone 2A: production-quality CPU baseline

Implemented: same-pattern SuperLU refactorization with retained ordering/permutations/factor structures, fixed CSC topology, and a direct comparison against the exact matrices used by the low-rank path.

Remaining:

- timing breakdown: model update, tangent scatter, factor/refactor, solve, residual, constitutive update;
- optional symmetric-indefinite/LDL backend where appropriate;
- benchmark against an external OpenSees/xara-class implementation.

## Milestone 2B: robustness hierarchy

Implemented: residual-based backtracking line search, time-step subdivision with rollback, solver-context caching by subdivision level, low-rank/direct fallback by active-rank threshold, and a deliberately failing coarse-step regression case.

Remaining:

- energy-based line-search option and residual-contraction policy;
- near-singular/negative-stiffness diagnostics and collapse classification;
- integrate calibrated direct/Woodbury switching into the robust subdivision driver; the normal prepared adaptive path already uses a measured model-specific crossover.

## Milestone 3: complete the 3D building compiler

Already implemented: 6 DOF/node, 12-DOF space frame, torsion, biaxial bending, local-axis reference vectors, generalized DOF springs, constant-preload geometric stiffness, and compiled sparse/nonlinear scatter.

Implemented beyond Phase 2: general sparse linear MPC transformation, true horizontal rigid diaphragms with master `RZ` coupling, non-diagonal condensed generalized mass, and arbitrary generalized/rotational-vector spring operators.

Remaining:

- optimized analytic corotational tangent validated against the Phase 7 finite-rotation energy reference (reference path implemented);
- member-local helper APIs that generate generalized hinge axes from element geometry;
- floor-level response reduction and torsional drift metrics.

## Milestone 4: constitutive models

Implemented as research paths: independently written peak-oriented IMK-family law, heterogeneous flat committed/trial state banks, event diagnostics, a virgin-elastic active-front fast path, concentrated-plasticity steel beams/columns with condensed end hinges, an axial bilinear BRB, a zero-length rotational panel-zone surrogate, and a power-law viscous damper with consistent velocity tangent. The defined bilinear steel mechanisms and damper law have deterministic OpenSeesPy 3.8.0 history parity; the degrading laws still require external cyclic parity before design use.

Remaining:

- external IMK cyclic-history/energy parity against OpenSees/xara or another validated implementation;
- BRB fatigue/fracture, cumulative strain, tension/compression asymmetry, casing and connection behavior;
- geometry-derived panel-zone properties and a finite joint-region formulation;
- Maxwell/relief-valve/gap/limit-state damper models;
- coupling-beam spring models;
- wall hinge/macro-element abstraction;
- soil/foundation springs and pile-line compilation;
- richer deterministic cyclic energy/deterioration reference tests.

## Milestone 5: external validation

- Mirror benchmark models in OpenSees/xara.
- Compare modal properties, static response, cyclic loops, and NRHA histories.
- Explicit tolerances for drift, residual drift, hinge deformation, base shear, acceleration, and dissipated energy.
- Add Perform/DRAIN comparison models where outputs are legally/operationally available.

## Milestone 6: GPU prototype

Requires a CUDA-capable external runner.

- CUDA SoA/AoSoA state kernels.
- Device-resident residual/tangent scatter.
- cuDSS backend behind the existing solver abstraction.
- Batched ground-motion state dimension.
- Active-record compaction.
- Nsight profiling.
- FP64 reference and guarded mixed-precision refinement.

## Decision gate

Do not expand into shells/solids before a real 3D concentrated-plasticity building demonstrates a material performance advantage against an optimized OpenSees/xara-class CPU baseline.

## Phase 9L: soil springs and piles

Implemented in the September 9, 2026 research release in three bounded layers:

1. A general zero-length translational/rotational spring and dashpot component with unilateral contact, gap/uplift, force caps, residual branches, and rollback-safe histories.
2. Named p-y, t-z, and q-z backbone adapters that convert distributed soil resistance to nodal springs using explicit tributary length and preserve the engineer's source/provenance. Do not embed copyrighted code tables or silently infer geotechnical parameters.
3. A pile-line compiler that meshes existing beam-column elements, attaches local-axis soil springs by depth, supports pile-head constraints, and exposes depth-wise force/deformation recorders. Validate spring laws first, then elastic beam-on-Winkler closed forms, then OpenSees single-pile histories, and only then pile groups/kinematic interaction.

The first release does not claim liquefaction, cyclic pore-pressure generation, group shadowing, free-field soil-column interaction, or SSI-compatible input motion. It provides bilinear tributary adapters rather than Simple1 internal laws. See `docs/PHASE9L_SOIL_PILES_REVIEW.md`.


## Phase 5 efficiency decisions

- Use the actual compiled sparse matrix to calibrate the direct-refactor/Woodbury crossover; rank zero uses the current direct factor, positive ranks use Woodbury only inside the measured profitable range.
- Preserve IMK active-front screening: representative strong records route roughly 85-91% of component evaluations through the virgin-elastic fast path.
- Exact story-block Schur substructuring is correct but not a serial CPU win at current block sizes; preserve it as a domain-parallel/GPU backend rather than making it the workstation default.
- Nonlinear-support DOFs need not be retained on the story-block Schur interface because low-rank updates act through exact baseline inverse operations.

See `docs/PHASE5_IMK_ADAPTIVE.md`.

## Phase 4 efficiency decisions

- Keep circuit-style localized tangent stamping and exact Woodbury as the concentrated-plasticity fast path.
- Calibrate direct-vs-low-rank crossover per model/topology; do not use a universal active-rank fraction.
- Keep global Craig-Bampton reduction as a validation/research utility; pursue local/story-block reduction that preserves sparse block structure.
- Use constitutive active-front diagnostics to justify screening only for expensive laws such as IMK/pinching/wall models.
- Make shared-memory ground-motion parallelism a first-class local workstation mode before requiring GPUs.

See `docs/PHASE4_EFFICIENCY.md`.


## Phase 7 decisions

- Keep the finite-rotation corotational element as the correctness reference; derive an analytic production tangent rather than optimizing numerical differentiation.
- Keep exact generalized `U C U^T` geometric updates, but let the cost model reject global Woodbury when active update dimension is large.
- Pursue geometric acceleration inside local/story blocks where update rank stays small and parallelism is available.
- Promote IDA/collapse-suite orchestration and record-level CPU parallelism.
- Next statistical gate: interval/right-censored fragility estimation rather than the current descriptive uncensored lognormal summary.
- Next credibility gate remains independent OpenSees/xara parity before GPU performance claims.
