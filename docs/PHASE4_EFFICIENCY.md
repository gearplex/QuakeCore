# Phase 4: cross-domain efficiency experiments

Phase 4 tested four ideas imported from circuit simulation, component-mode synthesis, event-driven simulation, and shared-memory HPC. The goal was to keep only optimizations supported by measurement.

## 1. Circuit-style low-rank/direct crossover

A 3,000-DOF tridiagonal network with 300 two-DOF device-stamp directions was solved repeatedly while increasing the active tangent-update rank. The same-pattern SuperLU baseline and exact Woodbury solver agreed to printed roundoff.

Representative 20-repeat run:

| Active rank | Woodbury speedup vs same-pattern |
| ---: | ---: |
| 1 | 7.15x |
| 5 | 6.29x |
| 10 | 5.36x |
| 20 | 3.40x |
| 40 | 2.15x |
| 80 | 1.15x |
| 120 | 0.75x |
| 200 | 0.22x |

**Decision:** active rank alone is not a universal solver-selection criterion. The crossover depends strongly on sparse fill/topology. A real 3D frame remains strongly favorable to Woodbury at rank fractions for which the tridiagonal network has already crossed over. The production solver should therefore use a measured/calibrated policy rather than a fixed fraction-of-DOF threshold.

The benchmark is `rank_crossover_bench`.

## 2. Craig-Bampton / Guyan model reduction

A generic reduction wrapper now supports

`u_full = T q`

with static constraint modes plus optional fixed-interface modes. Nonlinear constitutive state remains unchanged and the hinge basis is transformed as

`B_r = T^T B`.

Verification includes exact retained-DOF static condensation and modal enrichment.

For an 8-story 3D frame (624 DOF, 192 nonlinear hinges):

- nonlinear-interface retention reduced the Guyan system to 336 DOF;
- adding 12 fixed-interface modes produced 348 reduced DOF;
- Guyan-only first-mode period error was about 7.5%;
- 12 fixed-interface modes reduced the first three period errors to roughly 1-3 ppm;
- the nonlinear roof-history max-norm error was about 0.096% in the tested strong motion.

However, the globally condensed reduced matrices became much denser. The Craig-Bampton model was slower than the unreduced sparse model in this prototype despite its excellent accuracy.

**Decision:** do not use global dense component-mode reduction as the production acceleration path. Preserve the research implementation, but pursue *local/story-block substructuring* that retains global sparse/block structure.

The benchmark is `reduction_bench`.

## 3. Event-driven nonlinear frontier diagnostics

The transient driver now records total nonlinear-component evaluations and the subset whose tangent differs from its initial elastic tangent. This measures the potential benefit of an active-front/material-screening strategy before implementing a complex event scheduler.

For the 10-story, 2x2-bay 3D frame (240 candidate hinges, 300 steps):

- moderate yielding, amplitude 3: about 5.5% of component evaluations active;
- stronger yielding, amplitude 6: about 13.0% active;
- very strong excitation, amplitude 12: about 20.6% active.

For the 20-story case at amplitude 6, only about 4.3% of all component evaluations were in a non-initial tangent state even though the instantaneous active rank reached 102.

**Decision:** do not complicate the bilinear law with event scheduling; its full update is already cheap. Preserve a constitutive-bank boundary and add active-front screening when IMK/pinching/BRB/wall laws make inactive-state evaluation materially more expensive.

The new material-neutral API `evaluate_nonlinear_deformations()` supports that direction, while normal frame residual evaluation remains fused/allocation-free for performance.

## 4. Shared-memory record parallelism

Independent earthquake records are embarrassingly parallel. `run_record_suite_parallel()` now creates one prepared Woodbury solver per CPU worker and dynamically schedules records from a shared queue. Lazy influence caches are worker-local, so there is no shared mutable solver state.

In the current 5-core container:

### 20-story 2D, 12 records x 600 steps

- fresh full-factorization suite: ~1.43 s;
- sequential prepared Woodbury: ~0.125 s (11.4x vs fresh full);
- 4-worker prepared Woodbury wall time: ~0.051 s;
- additional parallel speedup: ~2.46x;
- wall-time speedup vs fresh full: ~28x.

A favorable repeat produced ~3.1x parallel scaling and ~35.7x vs the fresh baseline, so CPU scheduling/noise is visible at these very short runtimes.

### 10-story 3D, 8 records x 300 steps

- sequential prepared Woodbury: ~0.42-0.44 s;
- 4-worker wall time: ~0.14-0.16 s;
- parallel speedup: ~2.6-3.1x.

**Decision:** local CPU suite parallelism is a first-class deployment feature, not a fallback. GPU/cloud acceleration should be an additional throughput tier rather than a requirement.

## Current representative single-record result

After restoring a fused allocation-free constitutive/scatter loop, a 10-story, 2x2-bay 3D frame at amplitude 6 gives approximately:

- fresh SuperLU: ~1.5 s;
- optimized same-pattern SuperLU: ~0.46 s;
- sparse lazy Woodbury: ~0.062 s;
- Woodbury vs same-pattern: ~7.4x;
- identical 677 Newton iterations and identical printed peak roof response.

A 20-story, 1,560-DOF case remains roughly 6x-class faster than same-pattern SuperLU in representative runs.

These are internal algorithm benchmarks, not OpenSees/xara/Perform claims.

## Negative optimization retained as evidence

Replacing the tiny reduced solve with the container's generic LAPACK `DGESV` did not improve the rank-crossover benchmark; transpose/call overhead and the available BLAS implementation offset the benefit. The hand-written tiny dense solve remains for now. This should be revisited when benchmarking against MKL/OpenBLAS/vendor BLAS on the target Dell and GPU dense kernels later.

## Next efficiency experiments

1. Compile local/story-block static condensation or component-mode reduction without destroying global sparsity.
2. Add model-specific direct-vs-Woodbury calibration using actual factor fill and observed active-rank distribution.
3. Add timing decomposition for constitutive update, sparse residual, baseline solve, influence-cache growth, reduced solve, and output reduction.
4. Add optimized CPU sparse backends when available (KLU/LDL/PARDISO) and compare them on identical matrices.
5. Apply event/frontier screening only after an independently implemented IMK-family material bank exists.
6. Keep multi-record CPU parallelism as the local default; later map the same record dimension to GPU batching.
