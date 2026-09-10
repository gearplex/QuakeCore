# Phase 3: General MPCs, production CPU baseline, and robust transient stepping

## Purpose

Phase 3 removes three prototype shortcuts that could otherwise inflate the apparent advantage of the low-rank solver or prevent realistic building modeling:

1. equality-only 3D constraints;
2. fresh sparse symbolic/order work on every conventional tangent factorization;
3. an all-or-nothing Newton step with no rollback/subdivision path.

It also decouples constitutive state from the transient solver so future IMK/BRB/wall material banks can be introduced without changing Newton/Woodbury.

## General 3D linear MPC transformation

The 3D compiler now constructs a sparse transformation

`u_full = T q`

and recursively resolves constrained DOFs into retained generalized coordinates. Cycles and conflicting constraints are rejected during compilation.

A horizontal rigid diaphragm uses

`ux_s = ux_m - dy*rz_m`

`uy_s = uy_m + dx*rz_m`

`rz_s = rz_m`.

This transformation is applied consistently to stiffness, mass, nonlinear deformation operators, responses, and base excitation.

### Generalized mass

Rigid-diaphragm condensation produces a non-diagonal generalized mass matrix

`M_r = T^T M T`.

The runtime therefore uses a general `mass_multiply()` operation rather than assuming diagonal lumped mass after condensation. Eccentric floor mass correctly creates translation-torsion coupling and torsional inertial loading.

## Same-pattern SuperLU baseline

`SuperLUSamePatternSolver` uses the SuperLU expert driver with a retained column ordering, row permutation, elimination tree, and factor storage. Subsequent tangent matrices with the same CSC pattern use `SamePattern_SameRowPerm` rather than repeating the complete symbolic/order path.

Equilibration is disabled for this retained-factor backend so direct triangular solves remain algebraically consistent with the stored factors. A regression test verifies repeated same-pattern refactorization against a fresh sparse solve.

This is a substantially fairer conventional CPU baseline than the previous fresh-factorization comparison.

## Robust Newmark driver

`run_newmark_robust()` adds a separate robustness layer while leaving the simple fast path intact.

Implemented features:

- committed/trial state rollback;
- residual-based backtracking line search;
- recursive half-step subdivision after a failed step;
- linearly interpolated ground acceleration for internal substeps;
- solver contexts cached by subdivision level/time step;
- results returned on the original record time grid;
- Woodbury direct fallback when instantaneous reduced rank exceeds a configurable fraction of global DOF or a reduced solve fails.

A regression case intentionally limits Newton to two iterations. The ordinary solver fails, while the robust driver completes the record by subdivision. This case is retained to prevent future robustness regressions.

## Constitutive-state boundary

The transient solver no longer stores `BilinearState` or asks individual springs for their initial stiffness. The model exposes:

- `initial_nonlinear_tangents()`;
- `nonlinear_state_size()`;
- material-neutral flat committed/trial state vectors;
- internal force and tangent evaluation.

Current compiled frame models still use bilinear springs internally, but the solver is now independent of that choice. The intended next step is to compile multiple material types into contiguous state banks rather than introducing virtual material dispatch inside the hot loop.

## Arbitrary generalized hinge directions

The 3D builder can define a nonlinear generalized deformation as an arbitrary linear combination of nodal DOFs. A rotational-vector helper projects relative nodal rotation onto an arbitrary global axis. This supports rotated/skew framing and provides the same mathematical operator needed for future axial BRB, panel-zone, coupling-beam, and wall-hinge components.

## Validation status

Release and ASan/UBSan test suites pass after these changes. Tests now cover:

- 3D rigid-diaphragm kinematics;
- translation-rotation generalized mass coupling;
- torsional base excitation from eccentric mass;
- arbitrary rotational-vector hinge operators;
- same-pattern direct refactorization equivalence;
- robust time-step subdivision;
- full, same-pattern, Woodbury, and robust NRHA response equivalence in convergent reference cases.


## Sparse nonlinear update basis

The compiler no longer stores generalized hinge operators in a dense `n x m` matrix. `SparseUpdateBasis` stores each column's retained DOF indices and coefficients. The runtime uses sparse column dot products for deformation extraction and sparse AXPY for nonlinear force assembly.

The lazy Woodbury solver also consumes this sparse basis directly. Only an activated influence vector `A^-1 b_j` becomes dense, which is unavoidable for a general sparse factor inverse action. This keeps candidate-hinge memory proportional to actual component connectivity and removes an `O(n*m)` constitutive/deformation-evaluation cost from each nonlinear iteration.
