# Draft Tier 1-to-model crosswalk for the C2 pilot

This is a functional research crosswalk, not the ASCE 41 checklist. It intentionally paraphrases screening concerns and does not reproduce controlled standard text. Final identifiers, applicability conditions, and exact evaluation language require the licensed edition used for the study.

## Priority pilot variables

| Research ID | Tier 1 screening concern (paraphrased) | Model change | Expected mechanism | Pilot case |
|---|---|---|---|---|
| `D-VERT-01` | Abrupt vertical discontinuity or materially weaker/softer first story | Reduce first-story lateral strength/stiffness and redistribute upper-story strength to preserve the primary global strength anchor | Concentrated first-story drift and reduced system redundancy | `C2-SF`, `C2-SB` |
| `D-BRIT-01` | Existing wall/frame components lack deformation capacity and can lose shear or axial resistance | Activate calibrated post-peak shear deterioration and axial-support loss in designated components | Early strength degradation, loss of gravity support, or transition from distributed flexure to brittle response | `C2-RB`, `C2-SB` |

The first experiment is a two-by-two factorial. It estimates the individual effects of `D-VERT-01` and `D-BRIT-01` and whether their combination is worse than the product/additive expectation on the selected response scale.

## Later C2 screening families

| Research family | Screening concern (paraphrased) | Candidate model variable | Recommended phase |
|---|---|---|---|
| `D-TOR-01` | Plan eccentricity or torsion-sensitive lateral system | 3D mass/stiffness eccentricity and diaphragm rotation; report torsional ratio and `A_T` | Phase 2 after 2D pilot |
| `D-WALL-01` | Inadequate wall shear strength or web reinforcement | Wall shear backbone, pinching, post-peak degradation, and shear-flexure interaction | Phase 1 sensitivity or Phase 2 |
| `D-WALL-02` | Inadequate boundary confinement or compression strain capacity | Boundary concrete confinement, bar buckling, fracture, and axial-force interaction | Phase 2 |
| `D-CONN-01` | Weak transfer, collector, chord, or wall-to-diaphragm connection | Explicit force-limited connection and load-path loss | Phase 2 |
| `D-DIA-01` | Flexible, weak, or discontinuous diaphragm | Semi-rigid diaphragm and collector/chord components | Phase 2/3 |
| `D-OPEN-01` | Large or poorly detailed wall openings/coupling elements | Explicit pier/spandrel geometry and coupling-beam deterioration | Separate archetype family |
| `D-RED-01` | Limited redundancy or single critical wall line | Vary wall-line count and alternate load paths at fixed global strength | Separate redundancy study |
| `D-ADJ-01` | Pounding or incompatible adjacent-building response | Contact/gap model with independently varying neighbor periods | Building-context study |

## Mapping rules

1. A checklist finding becomes a model variable only when there is a defensible physical mechanism and a traceable parameter change.
2. A single finding may map to more than one uncertain parameter; preserve that uncertainty rather than selecting a favorable deterministic value.
3. The reference case represents the absence of the studied deficiencies, not a fully code-conforming modern building unless independently established.
4. Hold non-target variables fixed for the primary comparison and list every deliberate exception.
5. Do not label numerical noncompletion as brittle failure. A constitutive event, response threshold, or instability criterion must support the physical classification.
6. Preserve the Tier 1 observation, engineering interpretation, model parameter, analysis response, and fragility consequence as separate fields in the eventual dataset.

## Data row proposed for each evaluated building or archetype

| Field group | Minimum fields |
|---|---|
| Screening | edition, checklist identifier, applicability, result, observation, confidence |
| Interpretation | expected mechanism, affected components/stories, alternative explanations |
| Model | variable ID, parameter distribution, provenance, correlation, case IDs |
| Analysis | record, scale, solver, status, limit state, peak/residual EDPs |
| Fragility | censoring interval, fitted `IM50`, `beta_RTR`, reserve ratio, `beta_eff`, confidence interval |

