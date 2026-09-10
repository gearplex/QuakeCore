# Phase 1: Compiled 2D Frame Runtime

## Purpose

Phase 0 proved the exact low-rank solve on a shear-building abstraction and a generic sparse-grid matrix. Phase 1 moves the same nonlinear solver into a structural model that more closely resembles concentrated-plasticity PBSD modeling.

## Model compiler

`Frame2DBuilder::compile()` performs work that must not occur inside the nonlinear time loop:

1. maps arbitrary node IDs into dense node slots;
2. unions equal-DOF constraints;
3. propagates fixed support conditions through constraint groups;
4. assigns reduced DOF numbers;
5. accumulates lumped mass into retained DOFs;
6. transforms and assembles elastic member stiffness;
7. builds each rotational hinge generalized deformation column `b_j`;
8. creates `K_initial = K_linear + B diag(k0) B^T`;
9. preserves explicit diagonal CSC entries for dynamic mass terms;
10. records each hinge's exact CSC value indices and `b_i b_j` coefficients.

The nonlinear tangent update therefore changes **numerical values only**. Sparse topology does not need to be searched or reconstructed.

## 2D frame element

The current element is a prismatic Euler-Bernoulli 2D frame member with six global DOFs:

`[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]`

The element computes the local elastic axial/bending stiffness and applies the standard 2D direction-cosine transformation `K_g = T^T K_l T`.

A constant compressive preload can also be supplied. The prototype subtracts the conventional beam-column geometric stiffness matrix. This allows sensitivity tests for P-Delta softening but is not yet a corotational or dynamically updated axial-force formulation.

## Concentrated plasticity

A rotational spring between two node rotations uses

`theta = b^T u`

and contributes

`f_s = b M(theta)`

`tangent = b k_t b^T`.

Because `b` has only one or two nonzero entries after constraint elimination, each hinge produces a rank-one global tangent change.

The current material is bilinear kinematic hardening. IMK deterioration is the next constitutive milestone.

## Lazy Woodbury cache

The eager Phase 0 solver computed `A^-1 B` for every possible hinge before analysis. This scales poorly when a large building contains many candidate hinges that remain elastic.

`LazyLowRankWoodburySolver` instead:

1. factorizes the baseline `A` once;
2. identifies the currently active `delta_k` columns;
3. when a column first activates, solves `A w_j = b_j` and caches `w_j`;
4. expands the cached Gram matrix only for ever-active columns;
5. solves the exact reduced Woodbury system for the currently active subset.

The approach is exact for the represented tangent decomposition; laziness changes only when influence columns are computed, not the mathematical solution.

## Prepared record suites

`PreparedWoodburyNewmark` binds a model and integration time step to a reusable low-rank solver. Independent records start with clean structural state but reuse the model/dt-dependent baseline factorization and accumulated influence cache.

This is the CPU analogue of the intended GPU model where topology and compiled matrix information remain resident while record state changes independently.

## Current robustness boundary

A deliberately weak-hinge six/ten-story frame subjected to a stronger synthetic excitation can enter a nonconvergent regime with the present plain Newton algorithm. The weaker validation motion converges and gives identical full-Newton/Woodbury histories.

The stronger case is retained conceptually as a future regression target for:

- line search;
- adaptive tangent refresh/fallback;
- time-step subdivision;
- near-singular detection;
- collapse-state classification.

It should not be "fixed" by simply loosening equilibrium tolerances.
