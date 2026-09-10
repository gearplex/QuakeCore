# Phase 2: Compiled 3D Space-Frame Runtime

## Implemented 3D kinematics

The 3D runtime uses six DOF per node:

`[UX, UY, UZ, RX, RY, RZ]`.

The current prismatic space-frame element has 12 element DOFs and includes:

- axial stiffness `EA/L`;
- Saint-Venant torsion `GJ/L`;
- Euler-Bernoulli bending about both local principal axes;
- local-to-global coordinate transformation;
- user reference vector to define local +y orientation;
- optional constant compressive-preload geometric stiffness in both bending planes.

The local triad is formed from the member axis `e_x`, the projection of the user reference vector normal to `e_x`, and `e_z=e_x×e_y`. A near-parallel reference vector is rejected.

## Compiled constraints and nonlinear basis

The builder currently supports:

- fixed individual 3D DOFs;
- same-DOF master/slave equality constraints;
- reduced DOF numbering before analysis;
- accumulated lumped mass on retained DOFs;
- generalized bilinear spring deformation between any two retained/constrained nodal DOFs;
- automatic `B`-column generation;
- permanent CSC tangent scatter locations.

Phase 3 subsequently replaced the equality-only constraint compiler with a general sparse transformation and added true rigid-diaphragm translation/rotation coupling. See `PHASE3_MPC_ROBUSTNESS.md`.

## Verification

The test suite includes a vertical 3D cantilever with distinct `Iy` and `Iz`. Unit loads in the two global horizontal axes recover the corresponding analytical Euler-Bernoulli tip deflections.

A small 3D space-frame NRHA benchmark includes:

- vertical columns;
- X-direction beams with beam-end `RY` hinges;
- Y-direction beams with beam-end `RX` hinges;
- biaxial floor mass;
- horizontal X base excitation;
- nonlinear beam-end yielding in a genuine three-dimensional stiffness matrix.

Full-Newton and exact lazy-Woodbury histories are checked for equivalence.

## Representative nonlinear 3D benchmark

One representative run in the current CPU environment:

```text
stories=10
bays_x=2
bays_y=2
reduced DOF=780
elastic elements=210
candidate hinges=240
steps=300
synthetic excitation amplitude=2.5

full runtime≈1.4 s
full global factorizations=334
full Newton iterations=634

lazy Woodbury runtime≈0.16 s
baseline global factorizations=1
Woodbury Newton iterations=634

speedup≈8.6x
max active rank=54
cached influence columns=54
```

This result is much more relevant than the original shear-building benchmark because the matrix is a genuine 3D space-frame topology and the nonlinear basis contains hundreds of candidate hinges.

Phase 3 subsequently added an optimized same-pattern SuperLU comparison. The low-rank method remains faster in yielding reference cases, but the margin is smaller than against fresh factorization. See `BENCHMARKS.md`.
