# QuakeCore Phase 9K handoff — September 9, 2026

Phase 9K adds the first coherent steel-building component set to the Phase 9J wall source. Read `docs/PHASE9K_STEEL_REVIEW.md`, `docs/STEEL_COMPONENTS.md`, and `docs/JOB_FORMAT.md` first.

## Implemented

- Six-DOF 2D steel beam/column with elastic span and two statically condensed rotational hinges.
- Rollback-safe local Newton/backtracking and consistent condensed tangent.
- Bilinear, parameterized ASCE 41-style, and IMK-family end-hinge material selection.
- Dedicated rotational panel-zone surrogate and axial bilinear BRB component.
- Axial power-law viscous damper with explicit sublinear regularization and velocity-consistent effective tangent.
- Immutable CSC topology, global assembly, modal/NRHA/IDA integration, response histories, peak envelopes, and vertical-column story cuts.
- Mechanical, finite-difference, rollback, contract, Release, sanitizer, OpenSees component-history, and integrated-system tests.

## Evidence summary

Release CTest: 8/8 pass. Sanitizer CTest is included in the evidence package; leak detection is disabled because LSAN process inspection is unsupported in this environment.

OpenSeesPy 3.8.0 comparison maxima:

- steel-member end moment: 1.94e-10;
- steel-member hinge rotation: 1.96e-15 rad;
- panel-zone force: 6.39e-14;
- BRB force: 3.38e-14;
- nonlinear damper force: 7.11e-15;
- linear-damper SDOF displacement: 1.94e-15;
- linear-damper SDOF velocity: 2.85e-14.

The integrated strong synthetic fixture activates its column hinge, panel zone, and BRB; full and requested-Woodbury histories are identical. Halving the time step changes peak drift by 1.19%. A six-run IDA fixture completes without numerical gaps or failures and brackets only a configured drift threshold.

## Important limits

- Small-displacement member with elastic axial response; no axial yield/P-M interaction, evolving gravity/P-Delta, distributed plasticity, local/lateral-torsional buckling, fracture, or composite action.
- Panel zone is a user-parameterized relative-rotation spring, not a finite joint region or automatic AISC property model.
- BRB omits fatigue/fracture, cumulative strain, asymmetry, casing and connection response.
- Damper is memoryless; no Maxwell series spring, relief valve, gap, lock-up, stroke/force limits, or manufacturer qualification.
- Coupled member/device tangents currently use exact same-pattern direct refactorization. Woodbury is accepted as a request but correctly falls back; no Phase 9K performance claim is made.
- No physical specimen or code-acceptance validation has been completed.

## Next bounded phase: soil springs and piles

1. Add a general 2D zero-length translational/rotational foundation component with independent local axes, nonlinear spring materials, dashpots, gap/uplift/compression-only branches, rollback, energy output, and provenance.
2. Add named p-y, t-z, and q-z adapters that convert distributed resistance to nodal springs using explicit tributary length. Keep geotechnical parameter derivation outside the solver unless it is transparent, cited, and licensed for distribution.
3. Add a pile-line compiler that discretizes an existing beam-column line, attaches depth-indexed local soil springs, handles pile-head constraints, and records pile force/deformation and soil reaction by depth.
4. Verify material loops and tangents, elastic beam-on-Winkler closed forms, mesh/tributary-length convergence, local-axis rotation, rollback, and direct/fallback equivalence.
5. Compare one static and one cyclic single-pile fixture with OpenSees. Defer pile groups, liquefaction, pore-pressure generation, free-field columns, and kinematic interaction until the single-pile mechanism passes.
