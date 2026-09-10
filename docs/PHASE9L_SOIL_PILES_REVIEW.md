# QuakeCore Phase 9L soil-spring and pile implementation review

**Research release: September 9, 2026**

## Decision

Phase 9L establishes a usable first foundation mechanics layer: explicit directional zero-length soil springs and dashpots, provenance-preserving p-y/t-z/q-z bilinear adapters, and a vertical pile-line compiler with depth-wise response recorders. The implemented equations pass constitutive, closed-form, JSON-contract, sanitizer, and independent OpenSees checks.

This is software verification, not site-specific geotechnical validation. The adapters deliberately do not claim OpenSees Simple1 equivalence. They perform only explicit tributary force conversion and a documented bilinear approximation.

## What changed

- Added a rollback-safe asymmetric elastic-perfectly-plastic material with independent positive and negative capacities. A zero side capacity provides a no-resistance opening branch and committed permanent gap.
- Added arbitrary global-direction translational and relative-rotation springs between coincident nodes, compiled into the existing sparse low-rank basis.
- Added explicit-direction translational and rotational power-law dashpots between coincident nodes with a consistent velocity tangent.
- Added p-y, t-z, and q-z adapters with explicit length, surface-area, or toe-area conversion and mandatory provenance.
- Added a vertical uniform-mesh pile compiler using existing elastic beam-column elements and fixed coincident soil nodes.
- Added depth-wise pile displacement, soil reaction, member shear/moment, and toe response histories and peak profiles.
- Added malformed-input checks for IDs, profile lengths, directions, coincidence, capacity, suction, compression sign, and provenance.

## Mechanical verification

| Check | Result |
|---|---:|
| p-y nodal capacity conversion | Exact hand check |
| t-z nodal capacity conversion | Exact hand check |
| q-z nodal capacity conversion | Exact hand check |
| Nonuniform nodal tributary lengths | Exact hand check |
| Compression-only cap, unloading, and opening gap | Exact one-state return-map check |
| Directional spring deformation/force | Exact generalized-coordinate check |
| Directional linear dashpot rate/force | Exact check |
| Trial-state rollback | Passed |
| Semi-infinite Winkler pile head displacement | 0.0833% relative error with 200 elements |
| Semi-infinite Winkler pile head rotation | 0.1249% relative error with 200 elements |

The Winkler check uses `beta = (k/(4EI))^(1/4)` and the free-head, lateral-load solution `y(0) = P/(2 EI beta^3)`, `theta(0) = -beta y(0)`. The finite model extends ten characteristic lengths and uses consistent nodal tributary stiffness.

## Independent OpenSees comparison

The comparison script independently creates a four-element vertical `elasticBeamColumn` pile, five `Steel01` p-y springs with the same explicit tributary conversion, one linear `Viscous` dashpot, and Newmark average acceleration. It does not import QuakeCore forces or matrices.

| Metric | Result |
|---|---:|
| Time steps | 40 |
| QuakeCore peak absolute pile-head displacement | 0.0652099371652620 |
| OpenSees peak absolute pile-head displacement | 0.0652099371652618 |
| Maximum absolute full-history difference | 3.75e-16 |
| Difference / OpenSees peak | 5.75e-15 |
| QuakeCore Newton iterations | 84 |
| Subdivided steps | 0 |

OpenSees `Path` required one repeated trailing value to bracket its right endpoint and preserve QuakeCore's documented step-end sample convention. Without that harness correction, only the final time step differed because OpenSees returned zero beyond the unbracketed path endpoint.

## Regression and sanitizer status

The Release build passes all ten CTest targets, including the original numerical, reliability, job, wall, and steel suites plus the new soil component and JSON-contract suites. The Debug ASan/UBSan build also passes all ten targets with leak detection disabled because LeakSanitizer cannot inspect process threads in this environment. Memory-leak freedom is therefore not established by this run.

The pile NRHA example compiles 15 active DOFs, 11 soil springs, and one soil dashpot. It completes the nonlinear pulse with 84 Newton iterations and no time-step subdivision. Full history contains 40 depth-profile records.

## Important model boundaries

| Capability | Phase 9L status |
|---|---|
| Explicit bilinear p-y/t-z/q-z tributary adapters | Implemented and verified |
| OpenSees PySimple1/TzSimple1/QzSimple1 internal laws | Not implemented |
| Vertical elastic pile line | Implemented and verified for defined fixture |
| Nonlinear concrete pile section | Not implemented |
| Inclined piles, pile caps, groups, and group effects | Not implemented |
| Cyclic degradation, pinching, pore-pressure effects | Not implemented |
| Free-field/kinematic interaction and SSI input motion | Not implemented |
| Site-specific soil parameter generation | Not implemented |
| ASCE 41-23 acceptance evaluation | Not assessed |

## Recommended next foundation milestone

Before advancing to liquefaction or pile groups, add nonuniform and inclined local-axis pile meshing, nonlinear pile-section response, rigid pile-cap constraints, and independently specified cyclic soil backbones. Then validate against a published single-pile physical test with sufficient input and response data. Three-dimensional group and kinematic-interaction work should follow only after that single-pile validation closes.
