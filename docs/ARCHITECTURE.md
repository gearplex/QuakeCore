# Quake Engine Architecture Notes

## Product thesis

The target is a building-specific nonlinear response-history engine, not a general-purpose finite-element system. OpenSees/xara should be treated as a modeling/validation ecosystem; the runtime should be a compiled, data-oriented engine optimized around PBSD building abstractions.

## Execution layers

1. **Python/OpenSees-compatible front end**
2. **Structural intermediate representation (IR)**
3. **Model compiler**
   - DOF numbering and constraint elimination
   - fixed sparsity and scatter maps
   - element/material type grouping
   - sparse nonlinear-deformation basis `B` (column connectivity only; no dense `n x m` storage)
   - output/reduction plans
4. **CPU reference runtime**
5. **GPU runtime**
   - structure-of-arrays / tiled AoSoA state
   - batched ground motions
   - GPU material/element kernels
   - GPU residual and reductions
6. **Solver hierarchy**
   - tangent reuse
   - exact low-rank/Woodbury updates
   - nonlinear Schur complement
   - Krylov + building-specific preconditioner
   - GPU sparse direct fallback (cuDSS)
   - CPU robust direct fallback

## Critical performance KPI

Do not optimize only FLOP/s. Track **global numerical factorizations per record**. A solver that cuts factorizations by 10x can outperform a much faster factorizer while preserving an exact Newton tangent solve for localized changes.

## Phase 0 result criterion

A successful Phase 0 must show:

- Woodbury and full-factorization incremental solutions agree to near machine precision.
- NRHA histories agree within numerical tolerance for the same constitutive algorithm.
- Woodbury uses one global baseline factorization versus repeated full tangent factorizations.
- Sparse-grid microbenchmark demonstrates increasing advantage as global sparse factorization cost grows while active update rank remains small.

## GPU boundary

CUDA work should begin only after the CPU solver decomposition is verified. The intended GPU mapping is:

- nonlinear update basis: sparse column operators with compiled scatter
- material state: SoA/AoSoA
- one thread/lane per hinge/fiber integration point where practical
- model topology/scatter maps compiled once
- states remain device-resident through the record
- batches share topology/CSR indices and differ only in values/state/excitation
- asynchronous active-record compaction avoids forcing records with different Newton counts into lockstep

## Validation ladder

- SDOF elastic/bilinear analytical checks
- cyclic material loops
- small frame reference solutions
- OpenSees/xara benchmark models
- DRAIN/Perform-style concentrated-plasticity comparisons
- full-building response quantities and collapse sensitivity with time-step refinement
