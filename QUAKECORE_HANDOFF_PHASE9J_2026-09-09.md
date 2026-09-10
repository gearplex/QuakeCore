# QuakeCore Phase 9J handoff — September 9, 2026

Phase 9J adds in-plane MVLEM and SFI-MVLEM wall elements to the Phase 9I source. Read `docs/PHASE9J_WALLS_REVIEW.md`, `docs/WALL_ELEMENTS.md`, and `docs/JOB_FORMAT.md` first.

## Implemented

- `Wall2D` six-external-DOF MVLEM and SFI-MVLEM element.
- Concrete01-compatible, Steel01-compatible bilinear, and elastic uniaxial wall laws.
- Elastic, arbitrary layered, and fixed-angle RC plane-stress panels.
- Local `sigma_x=0` Newton solve and exact panel tangent condensation.
- Compiled full 6x6 wall CSC scatter, wall self-mass, frame/MPC integration, Rayleigh damping, modal, NRHA, and IDA paths.
- JSON wall input and full/summary wall output.
- State-aware tangent stability. Coupled walls route to exact direct sparse refactorization; incompatible scalar-bank reduction is rejected.
- Mechanical, contract, cyclic static, OpenSees dynamic, rollback, unit, Release, and sanitizer tests.

## Evidence summary

Release CTest: 6/6 pass. ASan/UBSan CTest: see evidence package; leak detection disabled because the execution environment cannot support LSAN process inspection.

OpenSeesPy 3.8.0, 740-step static paths:

- MVLEM elastic: force difference 1.82e-12 kN.
- MVLEM nonlinear: force difference 9.38e-13 kN.
- SFI elastic: force difference 9.86e-9 kN.
- SFI nonlinear two-layer: force difference 9.92e-9 kN.

Three-story, 1,200-step dynamic history maximum errors:

- MVLEM: displacement 2.41e-14 m; acceleration 4.26e-12 m/s²; wall force 1.98e-9 kN.
- SFI elastic: 1.33e-14 m; 1.76e-11 m/s²; 4.07e-9 kN.
- SFI nonlinear two-layer: 7.42e-14 m; 2.81e-11 m/s²; 4.45e-9 kN.

The fixed-angle RC panel replay agrees with OpenSees Concrete01/Steel01 stresses at QuakeCore strains, but it is not FSAM and is not experimentally validated.

## Important limits

- `fixed_angle_rc` is an explicit simplified panel, not OpenSees FSAM.
- No gravity equilibrium, evolving wall P-Delta, wall failure/collapse law, 3D wall, coupling beam, experimental specimen calibration, or ASCE/ACI acceptance provider.
- IDA examples use a synthetic waveform and its polarity reversal with a configured drift threshold.
- Wall material changes currently invoke exact same-pattern direct refactorization. A block-local wall update solver is a later performance task.

## Next bounded task

Implement and verify the complete RC panel behavior before adding more building-wall failure logic:

1. Port or independently implement a documented FSAM-compatible law, preserving source/license provenance.
2. Add material-level regression histories for ConcreteCM/SteelMPF, crack formation, reversal, aggregate interlock/friction, dowel action, and biaxial softening.
3. Reproduce PEER 2015/12 specimen RW-A15-P10-S78 or another specimen with published machine-readable loading and measured response.
4. Quantify peak strength, initial/secant stiffness, energy, residual drift, shear/flexural deformation split, and local strain histories.
5. Then add gravity-state transfer, wall geometric stiffness, axial/shear failure, and acceptance-limit tracking.
