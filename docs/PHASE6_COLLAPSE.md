# Phase 6: degrading hinges, near-collapse robustness, and stability certification

Phase 6 shifts the project from performance-first benchmarking toward collapse-fidelity infrastructure. The objective is not to declare an ASCE 41-compliant design engine; it is to provide a research implementation that can represent degrading concentrated-plasticity components, continue through post-capping response, and distinguish physical instability from numerical failure.

## 1. ASCE 41-style degrading hinge envelope

`ASCE41HingeMaterial` is a parameter-driven generalized component model. It does not embed tabulated acceptance/modeling values. The user provides the component-specific yield strength, plastic-deformation capacities, residual-strength ratio, optional IO/LS/CP limits, and cyclic-degradation controls.

The current envelope is:

- A-B: elastic response to effective yield;
- B-C: optional hardening to capping;
- C-D: finite negative post-capping branch to residual strength;
- D-E: residual branch followed by a short configurable residual-to-zero ramp;
- E-F: zero lateral/seismic resistance while effective/gravity capacity may remain;
- beyond F: configured component effective/gravity failure.

The finite C-D and pre-E ramps avoid artificial force cliffs. Cyclic response uses explicit committed branch state (envelope, unload-to-zero, reload-to-target, E-loss, F-loss) so nearby Newton trial points cannot infer different branches from the same committed history.

**Validation status:** the law has internal backbone, reversal, continuity, deterioration, E/F persistence, and NRHA solver-equivalence tests. External cyclic-history parity remains required before design use.

## 2. Performance levels are not collapse states

The robust driver records the first IO, LS, CP, E, and F events independently. CP exceedance does not terminate the analysis unless the user explicitly configures a beyond-CP component-count criterion.

`AnalysisTermination` distinguishes:

- `Completed`
- `PhysicalCollapse`
- `NumericalFailure`
- `InitialInstability`

Physical-collapse criteria can be combined independently:

- maximum story drift ratio;
- maximum number of components beyond CP;
- maximum number reaching E lateral-resistance loss;
- maximum number reaching F effective/gravity loss;
- near-zero/negative tangent-stability criterion;
- unrecoverable loss of equilibrium coincident with E-loss during repeated time-step subdivision.

This avoids the common but unsafe shortcut `nonconvergence == collapse`.

## 3. Explicit cyclic branch state fixed a real convergence pathology

An early Phase 6 degrading-hinge implementation could fail at one strong-motion step under all three linear solver strategies (fresh factorization, same-pattern factorization, and Woodbury). Increasing Newton iterations, line-search attempts, and time-step subdivision did not fix it.

Tracing the hinge path showed that the constitutive law inferred unloading/reloading branch from trial-point direction and peak history. Nearby Newton trials could therefore switch branch discontinuously.

The material state was rewritten to store the active branch, zero-force intercept, reload target, and last converged tangent explicitly. After the rewrite the same formerly pathological record completed with zero subdivisions and zero line-search failures under all three linear strategies. This is an important negative/positive result: the problem was constitutive continuity, not the low-rank solver.

## 4. Updated P-Delta reference path

The 3D frame can now update geometric stiffness from state-dependent axial compression. The internal-force formulation and Newton tangent include the derivative of axial force with respect to displacement; a finite-difference Jacobian regression test verifies consistency.

This remains a small-rotation/state-dependent geometric formulation, not a full corotational beam-column element. Because geometric tangent changes are not currently represented in the hinge-only Woodbury update basis, the robust Woodbury strategy deliberately falls back to direct tangent refactorization when updated P-Delta is active. Correctness takes priority over reporting an artificial low-rank speedup.

## 5. Initial stability and evolving mechanism tracking are separate

Phase 6 now uses two different tools for two different mathematical questions.

### Initial positive-definiteness certification

Before robust NRHA, the starting structural tangent is certified using a sparse symmetric elimination/Cholesky test after reverse-Cuthill-McKee ordering. A symmetric matrix is positive definite if and only if elimination proceeds with strictly positive pivots. The production path does not form an `n x n` dense matrix.

For small validation matrices, an optional LAPACK dense eigenvalue cross-check reports negative/near-zero eigenvalue counts and the minimum eigenvalue.

Representative one-time startup timings on stable 2x2-bay 3D frames:

| Stories | DOF | Initial K nnz | Sparse factor nnz | SPD check |
| ---: | ---: | ---: | ---: | ---: |
| 10 | 780 | 6,480 | 41,608 | ~0.013 s |
| 20 | 1,560 | 13,140 | 87,428 | ~0.026 s |
| 40 | 3,120 | 26,460 | 179,068 | ~0.059 s |

The robust driver performs this check by default and terminates at time zero as `InitialInstability` if it fails.

### Near-zero tangent tracking during response

After the initial tangent is certified SPD, inverse iteration is used only to track the eigenvalue nearest zero. The accepted-state tangent vector is carried explicitly so the stability check observes the actual converged branch tangent, not an unloading tangent obtained by re-evaluating the material at zero increment.

This is a mechanism-tracking metric, not a complete inertia calculation.

## 6. Energy ledger

For every accepted robust step (including subdivided substeps), the driver accumulates:

- earthquake input work;
- internal work;
- damping dissipation;
- kinetic energy;
- global balance residual.

Representative degrading-hinge sweeps in the six-story test model retain relative balance errors on the order of a few tenths of one percent near strong response. Energy balance is diagnostic rather than a collapse criterion, but it provides an independent numerical-health signal.

## 7. Representative collapse results

### Updated P-Delta, amplitude 15 synthetic record

All three linear strategies produce the same result:

- termination: component F effective/gravity-loss criterion;
- first IO: step 140;
- first LS: step 152;
- first CP: step 168;
- first E: step 195;
- first F: step 219;
- termination step: 219;
- max roof response: 0.750825;
- max story drift ratio: 0.059836;
- Newton iterations: 678.

### Updated P-Delta, amplitude 20 synthetic record

All three linear strategies again agree:

- termination: tangent-instability criterion;
- first IO: step 83;
- first LS: step 150;
- first CP: step 160;
- E/F not yet reached at termination;
- termination step: 174;
- max roof response: 0.533968;
- max story drift ratio: 0.033396;
- converged tangent near-zero estimate becomes negative;
- Newton iterations: 516.

This is an important case: physical instability can precede local E/F component loss.

### Fixed geometric baseline, amplitude 20

With state-dependent geometric updates disabled, the hinge-only tangent remains compatible with the exact low-rank solver:

- full / same-pattern / Woodbury termination: step 195 by the same 6% drift criterion;
- identical response and 452 Newton iterations;
- fresh full factorizations: 256;
- same-pattern factorizations: 125;
- Woodbury baseline factorizations: 1;
- representative elapsed times: ~0.033 s / 0.014 s / 0.009 s.

These are internal research timings, not external software claims.

## 8. Current limitations and next gate

Phase 6 does **not** establish design-grade collapse prediction. Major remaining items include:

1. external cyclic and NRHA parity against OpenSees/xara/other validated references;
2. a true 3D corotational beam-column formulation;
3. richer ASCE 41 component-specific hinge families and calibrated cyclic rules;
4. BRB, panel-zone, coupling-beam, wall, and foundation failure abstractions;
5. recovery of low-rank acceleration for state-dependent geometric stiffness, likely through localized element/block updates rather than hinge-only `B`;
6. probabilistic collapse/IDA validation over many records and scale factors;
7. output-level collapse taxonomy and audit trail suitable for engineering review.

The Phase 6 design principle is: **preserve the distinction between constitutive degradation, acceptance-limit exceedance, loss of component resistance, global instability, and numerical failure.**
