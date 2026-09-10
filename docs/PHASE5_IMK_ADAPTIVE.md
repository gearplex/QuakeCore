# Phase 5: IMK active front, adaptive solver calibration, and exact story-block substructuring

Phase 5 connects four ideas that matter for production concentrated-plasticity NRHA: a deterioration-capable hinge law, event/frontier screening, model-specific linear-solver selection, and local exact substructuring that preserves the full-order equations.

## 1. Independent peak-oriented IMK-family bank

`IMKPeakOrientedMaterial` is an independently written research implementation using the familiar peak-oriented IMK parameter semantics: initial stiffness, positive/negative yield force, pre-capping and post-capping deformation capacity, residual strength, ultimate deformation, and energy/cyclic deterioration controls.

The implementation is deliberately not copied from OpenSees. Current verification is internal and invariant-based:

- monotonic elastic/post-yield/post-capping/residual/ultimate branches;
- reversal detection;
- unloading-stiffness deterioration;
- strength/capping deterioration state;
- mixed bilinear + IMK flat state banks;
- full/same-pattern/Woodbury NRHA equivalence because all linear strategies consume the same constitutive evaluation.

**Important limitation:** external cyclic-history parity against OpenSees/xara or another validated IMK implementation is still required before this law can be treated as a design-grade material model.

## 2. Active-front fast path and instrumentation

Each nonlinear trial reports event flags for elastic fast-path use, yielding/capping, reversal, deterioration, reloading, and failure. A virgin IMK component that remains inside its initial elastic domain bypasses the peak/reversal/deterioration machinery and evaluates only the elastic force/tangent.

The transient driver records:

- total component evaluations;
- fast-path vs full-state evaluations;
- tangent-active evaluations;
- transition, reversal, deterioration, and failure events;
- maximum active tangent rank.

Representative results:

| Model / motion | Component evaluations | Virgin-elastic fast path | Full IMK path | Max active rank |
| --- | ---: | ---: | ---: | ---: |
| 10-story, 2x2-bay, amp 3 | 173,040 | 85.17% | 14.83% | 66 |
| 20-story, 2x2-bay, amp 6 | 274,560 | 91.01% | 8.99% | 204 |

This supports a data-oriented material-bank design: most potential hinges can remain on a cheap screening path even during a nonlinear record. The current bilinear law remains fused/allocation-free because a more elaborate event scheduler would cost more than it saves for such a cheap constitutive update.

## 3. Per-model Woodbury/direct crossover calibration

`calibrate_solver_crossover()` benchmarks the actual compiled effective matrix rather than applying a fixed `rank / DOF` rule. It compares:

- a warmed lazy Woodbury solve for representative active update directions; and
- a same-sparsity-pattern numeric SuperLU refactorization on the same matrix topology.

For the 20-story, 1,560-DOF / 480-hinge frame, a representative four-repeat calibration produced:

| Active rank | Direct-refactor / Woodbury time ratio |
| ---: | ---: |
| 1 | 37.7x |
| 5 | 42.3x |
| 20 | 33.2x |
| 80 | 15.0x |
| 120 | 12.6x |
| 200 | 5.39x |
| 320 | 1.71x |
| 480 | 0.55x |

The recommended sampled Woodbury rank limit is therefore 320 for this matrix. This is a tangent-refresh calibration, not a claim that every direct triangular solve is slow: rank-zero response should use the already-current direct factorization.

### Prepared adaptive policy

`PreparedAdaptiveNewmark` owns both a lazy Woodbury solver and a same-pattern direct solver:

1. **rank 0:** direct triangular solve using the current baseline factor;
2. **0 < rank <= calibrated limit:** exact Woodbury solve;
3. **rank above calibrated limit:** same-pattern direct numeric refactorization/solve.

The preparation/calibration cost is meant to be amortized across a ground-motion suite. On short synthetic records it can exceed the runtime of one analysis and should not be hidden in single-record performance claims.

## 4. Exact story-block Schur substructuring

`BlockSchurFactor` decomposes the baseline sparse matrix into independent local interior blocks plus a sparse retained interface. For each local block,

`S = A_BB - A_BI A_II^-1 A_IB`

is assembled exactly, then the sparse global Schur complement is factorized. There is no modal truncation and no dynamic approximation.

A key Phase 5 correction is that nonlinear hinge-support DOFs do **not** need to remain on the Schur interface. Substructuring happens inside the baseline inverse operator `A^-1`; Woodbury can still compute `A^-1 b_j` exactly even when the nonzeros of `b_j` lie inside eliminated local blocks. This keeps the interface much smaller than the earlier conservative partition.

`SubstructuredLazyWoodburySolver` uses this exact block factor for both the residual baseline solve and lazy influence-column construction. Regression tests match monolithic SuperLU and monolithic Woodbury to tight floating-point tolerance.

### Performance decision

The current serial CPU implementation does not beat the monolithic sparse factor on the tested story blocks. A 40-story research matrix with 1,280 DOF and four-story blocks has about 320 interface and 960 interior DOF, but multiple local triangular solves plus the Schur solve add more overhead than they remove in serial. Repeated high-rank Woodbury comparisons hover around parity and are noisy rather than a reliable win.

**Decision:** retain exact story-block substructuring as a domain-parallel/GPU building block, but do not make it the default workstation CPU path. Future parallelization should be enabled only when a compiler cost model predicts blocks large enough to amortize scheduling overhead.

## 5. Representative Phase 5 NRHA benchmarks

### 10-story 3D IMK frame

- 780 DOF, 240 IMK hinges, 300 steps;
- full SuperLU: ~1.61 s;
- same-pattern SuperLU: ~0.70 s;
- lazy Woodbury: ~0.068 s;
- Woodbury / same-pattern: ~10.4x;
- 721 Newton iterations for all paths;
- same printed peak roof response: 0.075914;
- adaptive prepared solve after calibration: ~0.064 s.

### 20-story 3D IMK frame

- 1,560 DOF, 480 IMK hinges, 250 steps;
- full SuperLU: ~2.64 s;
- same-pattern SuperLU: ~0.95 s;
- lazy Woodbury: ~0.161 s;
- Woodbury / same-pattern: ~5.9x;
- 572 Newton iterations for all paths;
- same printed peak roof response: 0.373069;
- adaptive prepared solve after calibration: ~0.165 s.

All figures are internal research timings in this container and are not claims against OpenSees, xara, Perform-3D, or commercial/HPC hardware.

## 6. Production implications

Phase 5 reinforces a three-level CPU policy:

- exploit exact low-rank structure while nonlinear rank is modest;
- retain a high-quality same-pattern direct factorization as a deterministic safety path;
- calibrate the crossover on the actual building rather than guessing from model size.

It also reinforces the local-first deployment thesis. IMK event screening and adaptive low-rank/direct switching reduce work before GPU acceleration; story-block substructuring is being preserved for future shared-memory/GPU/domain execution rather than imposed on the serial workstation path.
