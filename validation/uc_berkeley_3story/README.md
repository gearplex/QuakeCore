# UC Berkeley 3-story RC frame validation

## Purpose

This is QuakeCore's first formal external nonlinear response-history validation target.
It uses the three-story, three-bay, one-third-scale planar RC shake-table specimen in
NIST GCR 22-917-50, Chapter 3.

The validation intentionally has two references:

1. **Numerical parity target:** the NIST Perform3D nonlinear dynamic model.
2. **Physical target:** the measured UC Berkeley shake-table response.

Matching Perform3D tests implementation/model parity. Comparison with the experiment tests
model-form accuracy and should not be used to tune QuakeCore to the measured EDPs.

## Published global targets

| Metric | Measured | Perform3D |
|---|---:|---:|
| First-mode period | 0.34 s | 0.48 s |
| Peak base shear | 29.8 kip | 27.9 kip |
| Peak base shear / W | 0.51 | 0.48 |
| Story 1 peak drift | 5.18% | 6.07% |
| Story 2 peak drift | 4.70% | 3.92% |
| Story 3 peak drift | 2.61% | 1.91% |
| Story 1 residual drift | 0.30% | 3.0% |
| Roof residual drift | 0.24% | 1.9% |

NIST reports that the Perform3D model predicted the correct first-story-dominated mechanism
and similar peak base shear, but overpredicted permanent drift. This makes residual response
an especially useful non-calibration validation metric.

## Model definition known from the reference

- 2D, fixed-base frame.
- Three stories, three bays.
- Bay spacing = 5 ft 10 in center-to-center.
- Floor height = 4 ft.
- Columns = 6 x 6 in.
- Beams = 6 x 9 in.
- Floor weights = 19.6, 19.6, 19.3 kip from first floor through roof.
- Measured concrete compressive strength = 3.57 ksi.
- Perform3D columns use P-M-M hinges at both ends.
- Perform3D beams use lumped flexural hinges at both ends.
- Beam-column joints use nonlinear panel-zone components.
- P-delta is included.
- Total damping = 3%, implemented in the NIST model as 2.5% modal plus 0.5% Rayleigh.

## Input motion provenance rule

The strict parity run must use the **acceleration actually recorded on the shake table during
Dynamic Test 1**. The nominal command was Llolleo Component 100 scaled by 4.06, but the
recorded table motion is not identical to the original record multiplied by 4.06. Therefore,
using the original COSMOS record is permitted only for a development/pre-validation run and
must not be labeled Perform3D parity.

## Validation gates

### Gate A: model/modal parity

Before NRHA:

- QuakeCore `T1` within 1% of the Perform3D 0.48 s value.
- Confirm total weight and mass placement.
- Confirm P-delta and damping configuration.
- Compare static lateral stiffness/base-shear response if reconstructable from the reference.

### Gate B: nonlinear numerical parity

Using the recorded DT1 table motion and materially equivalent component models:

- Peak base shear target: 27.9 kip.
- Peak story drift targets: 6.07%, 3.92%, 1.91%.
- Residual drift targets: Story 1 = 3.0%; roof = 1.9%.
- Compare time histories, not only scalar peaks, when digitized/reference histories become available.

Target tolerances for an actually equivalent model are initially:

- period: <= 1%
- peak base shear: <= 2%
- peak story drift: <= 3%
- residual drift: <= 5% or an absolute 0.1 percentage-point floor, whichever is larger

These are engineering validation tolerances, not calibration objectives. If the inputs/model are
not identical, differences must be diagnosed rather than tuned away.

### Gate C: solver invariance

Run the identical QuakeCore model using fresh full factorization, same-pattern direct, and
adaptive/Woodbury. Demand parameters and response histories must agree to numerical tolerance.
Only after this check may runtime ratios be reported for the benchmark.

### Gate D: physical comparison

Compare QuakeCore against measured response separately. QuakeCore is not expected to match
experiment more closely than Perform3D without changing the model assumptions. Model-to-test
error is a model-form question, not a linear-solver verification question.

## Current blockers to strict apples-to-apples parity

1. The recorded Dynamic Test 1 shake-table acceleration file has not yet been acquired locally.
2. The current validation model uses a two-pass axial-demand-dependent P-M strength surrogate rather than the exact Perform3D coupled P-M-M hinge surface and hysteretic implementation.
3. Finite RC panel-zone compliance is now modeled, but the exact Perform3D joint property definition is not available in machine-readable form.
4. The exact ASCE 41 component parameter choices and axial-force iteration rules used by the NIST team must still be reconstructed independently.
5. The reported 2.5% modal + 0.5% Rayleigh damping is approximated with equivalent mass/stiffness-proportional terms at mode 1.

The published NIST Perform3D DT1 benchmark intentionally starts from a **virgin numerical state**; prior experimental motions were not replayed in that parity model. A separate pre-DT1 history study may be useful for specimen-to-test research, but it is not a requirement for reproducing the published Perform3D benchmark.

Until the blockers above are closed, the dynamic runs remain **pre-validation models**, not final parity benchmarks.

## Story-drift definition

The Ghannoum/Moehle experimental reporting used the **39 in clear column height** for the interstory drift ratios used in this validation comparison. Phase 8 and later therefore use 39 in for both published-drift comparison and the flexible first-story column length, with 4.5 in rigid offsets to the 48 in joint-center elevations. Older README text that described 48 in as the benchmark drift denominator was superseded after the source definition was rechecked.

## Phase 8.1 RC1 mechanism-validation correction

Phase 8.1 resolves two independent issues that made the nonlinear frame appear much stronger than its actual plastic mechanism:

1. **Base-shear semantics.** The former dashboard quantity `sum(m*a_abs)` is an inertial reaction and is not the same quantity as the first-story structural restoring shear used in a plastic-mechanism calculation. Phase 8.1 records these separately. The benchmark-facing dynamic quantity is the component-recovered first-story shear: A1/B1 use their explicit series shear-spring forces and C1/D1 use end-moment equilibrium. The end-moment-only sum is retained as a secondary diagnostic.
2. **Concentrated-hinge hardening semantics.** Zero-length hinges use a 100x elastic penalty stiffness to keep their pre-yield deformation small. The previous B-C hardening ratio multiplied this artificial penalty stiffness, unintentionally magnifying physical post-yield hardening. ASCE41 hinges can now specify an explicit physical B-C tangent independent of the penalty stiffness; legacy behavior is preserved when the explicit tangent is omitted.

A displacement-controlled first-mode pushover was also added. With the Phase 8.1 model it reaches a first-story mechanism/limit point after all four first-story columns yield.

| Quantity | Phase 8.1 RC1 | Published reference |
|---|---:|---:|
| First-mode period | 0.480 s | Perform3D 0.480 s |
| Static pushover peak | **23.56 kip** | NIST/FEMA P-2018 `Vy = 23.3 kip` |
| Story-1 drift at pushover peak | 2.125% | mechanism comparison, not a published EDP target |
| Dynamic peak first-story **component** restoring shear | **28.80 kip** | Perform3D peak base shear 27.9 kip; measured 29.8 kip |
| End-moment-only first-story shear diagnostic | **25.74 kip** | secondary diagnostic |
| Dynamic peak inertial reaction `sum(m*a_abs)` | **38.53 kip** | not used as the mechanism-capacity comparison |
| Story 1 / 2 / 3 peak drift | 3.28% / 4.25% / 2.18% | Perform3D 6.07% / 3.92% / 1.91% |
| Story 1 residual drift | 2.08% | Perform3D 3.00%; measured 0.30% |
| Roof residual drift | 1.35% | Perform3D 1.90%; measured 0.24% |
| Proxy 5% Sa at 0.48 s | 1.90 g | NIST/FEMA P-2018 reports 1.80 g at Te = 0.50 s |

The static strength result is particularly important: 23.56 kip is within about 1.1% of the independent NIST/FEMA P-2018 23.3-kip yield-strength estimate. The five-point FSC post-failure sensitivity bracket does not change the 28.80-kip component-recovered peak in the proxy run, so the global mechanism strength is controlled primarily by first-story flexural yielding rather than the selected post-shear-failure slope/residual strength.

The Phase 8.1 dynamic run still uses the deterministic proxy motion, so it remains **pre-validation**, not strict DT1 parity.

## Baseline proxy-input QuakeCore run

`run_validation.cpp` is the original transparent pre-validation representation. It uses the published geometry, floor weights, measured material strengths, constant gravity P-delta preload, scalar ASCE41-style member-end hinges, and rigid joint cores. The only fitted response quantity is the initial Perform3D period of 0.48 s. The deterministic 70 s input has PGA = 1.52 g but is **not** the measured DT1 table record.

Corrected baseline results:

| Metric | QuakeCore baseline | Perform3D | Measured |
|---|---:|---:|---:|
| First-mode period | 0.480 s | 0.480 s | 0.340 s |
| Legacy inertial reaction (previously labeled peak base shear) | 47.09 kip | 27.9 kip | 29.8 kip |
| Story 1 peak drift | 2.65% | 6.07% | 5.18% |
| Story 2 peak drift | 3.47% | 3.92% | 4.70% |
| Story 3 peak drift | 1.89% | 1.91% | 2.61% |
| Story 1 residual drift | 0.097% | 3.0% | 0.30% |
| Roof residual drift | 0.139% | 1.9% | 0.24% |

Fresh full factorization and same-pattern direct factorization still produce identical response histories. Active Woodbury nonlinear branch-path invariance remains open for this small, highly nonlinear nonsmooth model even though isolated Woodbury linear solves agree with direct solution to roundoff.

## Refined proxy-input validation run

`run_validation_refined.cpp` adds three independently defined modeling refinements without fitting nonlinear EDPs:

1. updated state-dependent P-delta;
2. finite RC panel-zone rotational compliance at elevated joints using a mechanics-derived research proxy; and
3. a two-pass axial-demand-dependent column P-M strength update, where pass-1 maximum compression is used to rebuild pass-2 column flexural strengths.

The current refined result (`quakecore_refined_result.json`) is:

| Metric | QuakeCore refined | Perform3D | Measured |
|---|---:|---:|---:|
| First-mode period | 0.480 s | 0.480 s | 0.340 s |
| Legacy inertial reaction (previously labeled peak base shear) | 40.98 kip | 27.9 kip | 29.8 kip |
| Legacy inertial reaction / W | 0.700 | 0.48 | 0.51 |
| Story 1 peak drift | 2.59% | 6.07% | 5.18% |
| Story 2 peak drift | 4.45% | 3.92% | 4.70% |
| Story 3 peak drift | 1.90% | 1.91% | 2.61% |
| Story 1 residual drift | 0.631% | 3.0% | 0.30% |
| Roof residual drift | 0.764% | 1.9% | 0.24% |

The progression under the same proxy input is:

| Model state | Base shear | S1 drift | S2 drift | S3 drift | Roof residual |
|---|---:|---:|---:|---:|---:|
| Updated P-delta, rigid joints | 44.19 kip | 2.69% | 3.80% | 1.87% | 0.179% |
| + finite panel zones | 41.25 kip | 2.61% | 4.26% | 1.96% | 0.807% |
| + second-pass P-M strengths | 40.98 kip | 2.59% | 4.45% | 1.90% | 0.764% |
| Published Perform3D | 27.9 kip | 6.07% | 3.92% | 1.91% | 1.90% |

Historical pre-Phase-8.1 diagnostic: Story 3 was essentially coincident with the Perform3D drift target and Story 2 was of comparable magnitude, while Story 1 appeared too stiff/strong when the inertial reaction was still labeled as base shear. Phase 8.1 supersedes that strength interpretation. That focuses the next modeling investigation on the first-story nonductile columns, particularly flexure-shear/gravity-capacity degradation and the exact column hinge formulation, rather than further global stiffness tuning.

The current panel-zone surrogate is intentionally not considered validated. Its component rotations are recorded in the JSON so joint response can be checked against published qualitative behavior rather than accepting a global EDP match obtained through unrealistically soft joints.

Full and same-pattern solutions are exactly invariant for the final refined run. The final refined model required 24 adaptive subdivisions with the current proxy motion.

## Interactive dashboard

- `dashboard_results.html` contains the corrected baseline proxy run.
- `dashboard_refined_results.html` contains the refined run and a model-progression table comparing each modeling step against Perform3D.

Both legacy dashboards are retained for provenance. `quakecore_phase8_1_validation_dashboard.html` is the preferred current visualization and distinguishes inertial reaction, physical first-story restoring shear, and the static pushover mechanism.
