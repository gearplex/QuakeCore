# QuakeCore Phase 8.1 RC1 Release Notes

**Release:** `v0.8.1-rc1`  
**Date:** 2026-09-05  
**Validation status:** `PRE-VALIDATION_PHASE8_1_RC1_PROXY_INPUT`

## Purpose

Phase 8.1 is a mechanism-validation correction to Phase 8 RC1. It resolves the apparent ~41 kip nonlinear base-shear discrepancy by separating inertial reaction from structural restoring shear and by decoupling physical post-yield hardening from the numerical zero-length-hinge penalty stiffness.

The UC Berkeley three-story RC frame remains the primary external benchmark. The recorded DT1 table acceleration has not yet been recovered, so response-history results remain proxy-input pre-validation.

## Major corrections

### 1. Base-shear semantics

Previous validation dashboards labeled

`sum(m_i * a_abs_i)`

as "base shear." That quantity is an inertial reaction and can contain the dynamic force needed to balance viscous damping as well as structural restoring forces. It is not the same quantity as the first-story restoring shear used in a plastic-mechanism calculation.

Phase 8.1 records four quantities separately:

- `peak_inertial_base_reaction_kip`
- `peak_reduced_generalized_restoring_force_kip` (diagnostic generalized-coordinate force)
- `peak_first_story_column_restoring_shear_kip` (end-moment-only diagnostic, `sum[(Mbot+Mtop)/39]`)
- `peak_first_story_component_restoring_shear_kip` (A1/B1 explicit series-shear spring forces plus C1/D1 moment-derived shear)

The benchmark-facing dynamic strength comparison uses the fourth quantity. The static pushover load factor remains the cleanest global mechanism-capacity ordinate.

### 2. Physical hinge-hardening semantics

Concentrated plastic hinges use a 100x elastic penalty stiffness to keep their pre-yield deformation small. Previously the ASCE 41 B-C hardening ratio was multiplied by this artificial stiffness, unintentionally magnifying physical post-yield hardening.

`ASCE41HingeParams` now includes optional `hardening_stiffness` in physical force/deformation units. When this value is positive, the B-C tangent is independent of the penalty stiffness. When omitted, legacy `hardening_ratio * Ke` behavior remains unchanged for backward compatibility.

For the Berkeley model, beam and ductile-column post-yield slopes are referenced to the physical member rotational stiffness rather than the 100x hinge penalty stiffness.

### 3. Static displacement-controlled pushover

A validation-oriented displacement-control solver was added to the Berkeley validation executable. The augmented Newton system is solved directly so a singular structural tangent does not automatically terminate the step through a `K^-1` Schur complement.

The pushover is controlled by Story-1 displacement using a first-mode inertial load pattern. It tracks first-story column yielding and classifies the terminal limit point separately from numerical collapse.

## Phase 8.1 validation results

| Quantity | Phase 8.1 RC1 | Reference |
|---|---:|---:|
| First-mode period | 0.480 s | Perform3D 0.480 s |
| Static pushover peak | **23.555 kip** | NIST/FEMA P-2018 `Vy = 23.3 kip` |
| Pushover error vs NIST `Vy` | **+1.10%** | |
| Story-1 drift at pushover peak | 2.125% | |
| All first-story columns yielded by | 1.45% drift | |
| Dynamic peak first-story **component** restoring shear | **28.798 kip** | Perform3D 27.9 kip; measured 29.8 kip |
| Dynamic component-shear error vs Perform3D | **+3.22%** | |
| Dynamic component-shear error vs measured | **-3.36%** | |
| End-moment-only first-story shear diagnostic | **25.737 kip** | secondary diagnostic |
| Peak inertial reaction | **38.526 kip** | intentionally not used as mechanism capacity |
| Peak S1 / S2 / S3 drift | 3.278% / 4.248% / 2.176% | Perform3D 6.07% / 3.92% / 1.91% |
| Story-1 residual drift | 2.078% | Perform3D 3.00%; measured 0.30% |
| Roof residual drift | 1.351% | Perform3D 1.90%; measured 0.24% |
| Proxy 5% `Sa(0.48 s)` | **1.900 g** | NIST/FEMA P-2018 `Sa(Te=0.50s) = 1.80 g` |
| Proxy spectrum difference | **+5.56%** | |

The static mechanism result is the main Phase 8.1 strength-validation gate. The 23.555 kip peak is within approximately 1.1% of the independent NIST/FEMA P-2018 23.3 kip yield-strength estimate, and all four first-story columns have yielded before the terminal limit point.


## Mechanism hierarchy cross-check

The Phase 8.1 model properties give a simple first-story double-curvature column mechanism of **25.641 kip** from `2*sum(Mn)/39`. The individual gravity-load nominal moments are approximately A1/B1/C1/D1 = **149.52 / 158.57 / 103.68 / 88.23 kip-in**. The beam nominal moment is about **203.50 kip-in**, while the panel-zone generalized strengths are approximately **258.12 kip-in exterior** and **344.16 kip-in interior**. The initial hierarchy therefore favors a first-story column mechanism rather than a beam or joint mechanism.

The displacement-controlled QuakeCore pushover peaks at **23.555 kip**, only **+1.10%** above the independent NIST/FEMA P-2018 `Vy = 23.3 kip`. The difference between the simple 25.641-kip hand mechanism and the 23.555-kip nonlinear pushover is consistent in scale with P-Delta softening and redistribution as the mechanism develops.

## Dynamic force interpretation

The legacy inertial reaction is **38.526 kip**, while the benchmark-facing component-recovered first-story restoring shear is **28.798 kip**. The latter uses the explicit A1/B1 series shear-spring forces and C1/D1 moment equilibrium, which is more appropriate for the Phase-8 series topology than the end-moment-only 25.737-kip diagnostic. This resolves the apparent ~40-kip structural mechanism without redefining or fitting component strengths.

The dynamic component-recovered restoring shear is within about **+3.22%** of Perform3D (27.9 kip) and **-3.36%** of the measured test (29.8 kip). Because the actual recorded DT1 table motion is still unavailable, this agreement is encouraging but is not claimed as strict response-history validation.

## FSC behavior after hardening correction

The physical hardening correction changes the accepted dynamic state path. In the current proxy run:

- B1 FSC initiation: 9.44 s, 2.804% Story-1 drift, 10.53 kip shear, 21.38 kip compression
- A1 FSC initiation: 22.32 s, 2.875% Story-1 drift
- B1 initiates before A1 in every post-failure sensitivity case
- no F / gravity-capacity redistribution is claimed

The previous Phase 8 RC1 B1 initiation state should therefore not be treated as the current release result.

## Sensitivity and solver gates

Five post-shear-failure parameter combinations were rerun. All completed, and B1 initiated before A1 in every run. The peak component-recovered first-story restoring shear remained **28.798 kip** across the bracket (the end-moment-only diagnostic remained 25.737 kip), confirming that peak mechanism strength is governed primarily by first-story flexural yielding rather than the selected post-failure FSC slope/residual ratio.

FullFactorization and SamePattern direct refactorization remain exactly invariant for the roof, inertial-reaction, end-moment shear, and component-recovered first-story shear histories in the release run.

## Remaining blockers

Phase 8.1 is not strict external DT1 validation. Remaining blockers include:

1. actual recorded Dynamic Test 1 shake-table acceleration history;
2. production moving-surface P-M-M return mapping for axial-flexural coupling;
3. static gravity-load redistribution after axial-capacity F loss;
4. exact Perform3D component definitions and joint/hinge parameters.

## Preferred dashboard interpretation

For benchmark dynamic strength comparison use:

`peak_first_story_component_restoring_shear_kip`

Retain `peak_first_story_column_restoring_shear_kip` as the end-moment-only diagnostic.

Do not substitute:

`peak_inertial_base_reaction_kip`

The dashboard shows both quantities explicitly and includes the Phase 8.1 pushover curve with the NIST/FEMA P-2018 23.3-kip reference line.
