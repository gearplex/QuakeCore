# Phase 7: corotational reference, localized geometric updates, and IDA

Phase 7 addresses three Phase 6 gaps: finite-rotation frame kinematics, exact localized representation of changing geometric tangents, and record/scale collapse-suite infrastructure. The design rule remains evidence-driven: an exact formulation is not promoted to the fast path unless it also wins the measured cost model.

## 1. Finite-rotation 3D corotational reference

`CorotationalFrame3D` is a correctness-first, energy-based 3D reference element. It updates the member chord and co-rotated frame from the current configuration, expresses end rotations relative to that frame, and obtains internal force and a symmetric consistent numerical tangent from the strain-energy potential.

The element is deliberately not advertised as algebraically identical to any production coordinate-transformation implementation. Its present role is a finite-rotation reference and a validation target for a later optimized analytic tangent.

Validation gates now include:

- arbitrary-axis rigid-body rotation plus translation: negligible basic deformation/internal force;
- small-response tangent recovery of the existing linear 12-DOF space-frame element;
- small-motion NRHA convergence to the linear solution;
- Euler cantilever stability under mesh refinement.

For the current cantilever stability benchmark, with exact Euler load `Pcr = 2313.1885`:

| Elements | Predicted Pcr / Euler Pcr | Error |
| ---: | ---: | ---: |
| 1 | 1.21615 | +21.62% |
| 2 | 1.05264 | +5.26% |
| 4 | 1.01316 | +1.32% |
| 8 | 1.00332 | +0.33% |
| 16 | 1.00081 | +0.08% |

The one-element result demonstrates why objectivity alone is not a sufficient validation criterion. Mesh convergence to the classical stability solution is the more meaningful gate.

## 2. Exact generalized localized tangent updates

Phase 7 generalizes the hinge-only Woodbury form

`A + B diag(dk) B^T`

to

`A + U C(u) U^T`,

where `U` is a fixed sparse basis and `C(u)` is a small state-dependent dense block matrix. Ordinary hinge tangent changes are the diagonal special case.

For the Phase 6 updated-P-Delta formulation, each element's fixed unit geometric matrix is decomposed into a compact local eigenbasis and combined with the axial-deformation direction. The resulting state-dependent tangent increment is represented exactly in `U C U^T`; direct and generalized-Woodbury solutions agree to numerical tolerance.

### Performance result: exact does not mean faster

When every column geometric state evolves, the active generalized space can become too large. In the six-story collapse model the generalized basis reaches 268 active directions for only 168 structural DOF. Forcing a global generalized Woodbury solve is therefore much slower than same-pattern sparse refactorization.

A 20-story column-chain microbenchmark similarly gives:

- structural DOF: 120;
- generalized basis/active dimension: 138;
- generalized Woodbury: ~0.00235 s;
- direct tangent solve: ~0.000054 s;
- generalized-Woodbury/direct speed ratio: ~0.023.

The production rule is therefore unchanged: the model-aware cost/fallback policy must be allowed to choose direct refactorization. The exact `U C U^T` representation is retained because it is useful for block-local/domain-parallel solvers and for cases where only a small subset of geometric states changes.

## 3. IDA and collapse-suite infrastructure

`run_ida_suite` runs many records over user-specified scale factors using the robust Phase 6 collapse driver. It supports:

- shared-memory parallelism by record;
- physical-collapse vs numerical-failure vs right-censored outcomes;
- collapse-mechanism classification;
- logarithmic bracket refinement after the first collapse/noncollapse pair;
- per-record last-noncollapse and first-collapse scales;
- collapse PGA summaries;
- descriptive lognormal median and log standard deviation for uncensored collapse PGAs.

The current lognormal summary is intentionally descriptive, not a censored-data maximum-likelihood fragility estimator.

Representative synthetic suite:

- 8 records;
- 89 total nonlinear analyses after adaptive scale refinement;
- 8 physical collapses;
- 0 numerical failures;
- median collapse PGA: 7.8845;
- descriptive beta_ln: 0.05136;
- collapse causes include tangent instability and F-loss.

Representative wall time:

- 1 worker: ~0.466 s;
- 4 workers: ~0.165 s;
- speedup: ~2.8x.

The collapse brackets and mechanisms are identical across worker counts.

## 4. Collapse mechanism taxonomy

`AnalysisResult` now carries an explicit `CollapseMechanism` independent of the human-readable termination reason:

- initial instability;
- tangent instability;
- drift limit;
- lateral-resistance loss (E);
- gravity/effective-resistance loss (F);
- configured beyond-CP component criterion;
- unrecoverable lateral-loss equilibrium;
- numerical failure.

This enables IDA/fragility studies to distinguish *how* a model terminated rather than treating every failed run as one undifferentiated collapse observation.

## 5. Phase 7 decisions

Promote:

- finite-rotation corotational reference path for fidelity/validation;
- explicit Euler-stability convergence tests;
- generalized exact `U C U^T` update algebra;
- IDA record/scale orchestration and cause-of-collapse classification;
- record-level CPU parallelism for local workstations.

Do not promote as default fast path:

- one giant global generalized-Woodbury solve for all changing geometric states;
- the current numerically differentiated corotational tangent as a production-speed element.

## 6. Next gates

1. derive/implement an optimized analytic corotational basic-to-global tangent and validate against the finite-difference energy reference;
2. exploit geometric updates inside story/local blocks rather than one global reduced update space;
3. add a statistically rigorous fragility fit that handles right-censored and interval-censored IDA observations;
4. external validation against OpenSees/xara: modal, static, buckling, cyclic hinge, single-record NRHA, and IDA collapse scale/cause;
5. add BRB, panel-zone, coupling-beam, wall, and foundation component families;
6. only after these CPU/fidelity gates, start the CUDA/cuDSS backend.
