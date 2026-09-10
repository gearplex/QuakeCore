# QuakeCore Phase 9L handoff, September 9, 2026

Phase 9L adds the first bounded soil-spring and pile-foundation layer to the Phase 9K source. Read `docs/PHASE9L_SOIL_PILES_REVIEW.md`, `docs/SOIL_FOUNDATIONS.md`, and `docs/JOB_FORMAT.md` first.

## Implemented

- Directional translational and relative-rotational zero-length soil springs between coincident nodes.
- Symmetric bilinear and asymmetric elastic-perfectly-plastic material options, including compression-only/uplift opening and committed permanent-gap behavior.
- Directional translational or rotational linear/power-law dashpots with a consistent velocity tangent.
- Named p-y, t-z, and q-z bilinear adapters with explicit tributary length, pile surface area, or toe area conversion.
- Mandatory parameter provenance and rejection of nonzero adapter hardening that would violate the declared ultimate force cap.
- Vertical uniform-mesh pile-line compiler using existing elastic beam-column segments, explicit node/member/spring IDs, and fixed coincident soil nodes.
- Depth-wise pile displacement, soil force/deformation, member shear/moment, and toe response histories and peaks.
- NRHA and IDA integration through the existing compiled sparse, rollback, subdivision, and output paths.

## Verification evidence

- Release CTest: 10/10 pass.
- Debug ASan/UBSan CTest: 10/10 pass with leak detection disabled because LSAN thread inspection is unsupported in this environment.
- Semi-infinite beam-on-Winkler head displacement error: 0.0833% with 200 elements.
- Semi-infinite beam-on-Winkler head rotation error: 0.1249% with 200 elements.
- Independent OpenSeesPy 3.8.0 full-history pile-head comparison: maximum absolute difference 3.75e-16, or 5.75e-15 of peak.
- Integrated pile fixture: 15 active DOFs, 11 soil springs, one soil dashpot, 40 full depth-profile records, 84 Newton iterations, no subdivision.

Run the external comparison with an OpenSeesPy environment:

```bash
PYTHONPATH=/path/to/openseespy python validation/soil/compare_opensees_pile.py ./build/quake_run phase9l_opensees.json
```

## Important limits

- The p-y/t-z/q-z adapters are documented bilinear approximations, not implementations of OpenSees `PySimple1`, `TzSimple1`, or `QzSimple1`.
- The pile is vertical, uniformly spaced, small displacement, and elastic. No nonlinear concrete pile section, inclined pile, nonuniform mesh, cap, or group is compiled.
- No p-multipliers, group shadowing, liquefaction, pore-pressure generation, cyclic degradation/pinching, free-field column, kinematic interaction, or SSI-compatible input motion.
- No site-specific parameter generation, physical-test calibration, or ASCE 41-23 acceptance evaluation.

## Recommended next bounded phase

Add nonuniform and inclined local-axis pile meshing, nonlinear pile-section response, and rigid pile-cap constraints. Introduce independently specified cyclic soil backbones and close a published physical single-pile validation before attempting three-dimensional group or kinematic-interaction models.
