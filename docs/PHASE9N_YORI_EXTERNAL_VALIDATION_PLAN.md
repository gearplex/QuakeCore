# Phase 9N — YORi multi-record exact-law external validation plan

**Date:** September 17, 2026  
**Predecessor:** Phase 9M / Gate 4 engineering parity  
**Status:** Initiated  
**Purpose:** Broaden independent OpenSees/QuakeCore evidence beyond the single archetype-10001 record without converting software verification into a FEMA P-695 or code-qualification claim.

## 1. Frozen starting point

The Phase 9N validation program starts from the stabilized Gate 4 production-source baseline:

- source evidence baseline: `76dd01e9d671d82d2e2163ed5219bc1f355fdb67`;
- reconstructed-surrogate SHA-256: `22203092a1d125d59ced2a3777a9120ba9209d190ea32f6c71107f9c97403cf3`;
- lost received-source SHA-256: `02eca2fbdde5d2cc1400be354802a2b26d6d79ce15ac94dd981cd26dc6b23104`;
- independent reference engine: OpenSeesPy 3.8.0;
- QuakeCore material stack: promoted ConcreteCM + Pinching4 + MinMax + Parallel implementation;
- Gate 4 acceptance implementation: `validation/yori_wall_10001/assess_engineering_parity.py`.

Documentation-only commits after the source evidence baseline do not redefine the validated mechanics.

The existing Gate 4 thresholds are frozen and must not be loosened after reviewing Phase 9N results.

## 2. Claim boundary

Phase 9N is an **external software/model validation program** for the documented reconstructed surrogate.

It does not establish:

- identity to the lost received model input;
- physical validation of the YORi wall system;
- FEMA P-695 collapse capacity, CMR, SSF, ACMR, or acceptance;
- a collapse intensity measure;
- code approval or product prequalification;
- an approved or recommended YORi response-modification coefficient.

Any later FEMA P-695 program must separately establish its record set, normalization/scaling procedure, collapse classification, censoring treatment, uncertainty model, and acceptance criteria.

## 3. Record manifest requirement

Do not add records to the suite without a machine-readable manifest entry containing at least:

- record identifier;
- source/provenance;
- component identifier and direction;
- original and analysis time step;
- units;
- number of points and duration;
- processing applied before analysis;
- source-file SHA-256;
- any resampling/filtering/baseline correction;
- scale multiplier actually applied in each run;
- the semantics of that multiplier;
- whether the record is part of a formal FEMA P-695 set or only a software-validation record.

The repository currently embeds the single processed comparison history used in Phase 9M. Additional record files/provenance must be supplied or otherwise legally/operationally available before the multi-record suite can be treated as reproducible.

## 4. Pilot suite

The first Phase 9N batch should use multiple distinct records while leaving the archetype/model definition unchanged.

For each available record:

1. run the same exact-law reconstructed surrogate in QuakeCore and OpenSees;
2. begin with a low/moderate scale that both solvers can independently complete where practical;
3. include the nominal scale used for the software-validation study;
4. add higher diagnostic multipliers only when useful for exercising nonlinear path dependence;
5. record every actual scale multiplier and solver termination status.

Acceleration multipliers used in this phase are software-validation scales only unless an explicitly documented intensity-measure calculation says otherwise.

## 5. Required response quantities

### Existing frozen Gate 4 metrics

Continue to gate on:

- first modal periods;
- peak story drift ratios;
- story-drift history RMSE normalized by peak;
- MVLEM peak shear;
- MVLEM peak axial force;
- MVLEM peak bottom moment.

Existing Gate 4 thresholds remain:

| Metric | Threshold |
|---|---:|
| Modal-period relative error | <= 1e-6 |
| Peak story-drift relative error | <= 0.1% |
| Drift-history RMSE / peak | <= 0.1% |
| Wall peak shear/axial/bottom-moment EDP relative error | <= 1% |

Comparison is over independently converged common histories unless both solvers complete the full record.

### Additional Phase 9N diagnostic metrics

Add and report, but do not assign acceptance thresholds post hoc:

- residual story drift;
- floor absolute acceleration;
- base shear;
- story restoring shear;
- peak wall deformation/curvature measures available from the current model;
- hysteretic/dissipated energy by relevant material/component and total system, where consistently defined;
- solver iterations, accepted subdivisions, rejected attempts, and termination time.

Thresholds for any new required acceptance metric must be defined and frozen **before** using it as a batch pass/fail criterion.

## 6. Run-status taxonomy

Every solver/record/scale combination must end in exactly one primary status:

- `completed` — target record duration completed;
- `response_limit_crossed` — an explicitly defined response threshold terminated the analysis;
- `numerical_noncompletion` — solver could not advance under the declared algorithm/recovery policy;
- `classified_collapse` — reserved for a future phase with an explicit collapse definition and independent qualification.

Do not infer collapse merely from numerical noncompletion.

Also record:

- last committed time;
- last nominal step;
- whether substepping occurred;
- maximum subdivision depth;
- final convergence norm/iteration count where available;
- applied scale multiplier.

## 7. Common-window comparison policy

When one solver stops before the other:

- compare EDPs only through the last independently committed common time;
- report both solvers' completion status separately;
- retain the longer QuakeCore history for diagnostic purposes, but do not use post-reference response to claim two-solver parity;
- keep full-record comparison available whenever both solvers complete.

A QuakeCore-only continuation is not evidence of OpenSees parity beyond the reference solver's last committed state.

## 8. Runtime benchmarking protocol

Use the Phase 9M paired benchmark design as the minimum standard:

- same record prefix and number of steps;
- equivalent output mode;
- one warm-up per solver;
- seven paired/interleaved timed repetitions unless a later benchmark protocol explicitly changes this;
- alternate solver execution order;
- report individual paired times, medians, paired speedup factors, and range;
- record runner/CPU/software environment sufficiently to interpret variability.

Do not combine the earlier simple-law timing result with exact-law timing as though they were one population.

No universal speedup factor should be claimed until the suite spans multiple records and more than one model size/topology.

## 9. Phase 9N deliverables

1. machine-readable record manifest;
2. reproducible QuakeCore job generation for each record/scale;
3. independent OpenSees comparison runner;
4. per-run status ledger;
5. per-record/per-scale EDP comparison JSON;
6. suite-level summary with frozen Gate 4 pass/fail metrics;
7. expanded diagnostic metrics;
8. paired runtime report;
9. documented unexplained divergences, if any;
10. explicit recommendation on whether the evidence is sufficient to proceed to a FEMA P-695-specific phase.

## 10. Phase 9N exit criteria

Phase 9N is ready to close when:

- the selected multi-record software-validation set is reproducible from its manifest;
- the current Gate 4 gated metrics remain within their frozen thresholds over independently converged windows, or any exceptions are explicitly explained and dispositioned;
- no unresolved QuakeCore-only constitutive/state-management defect is indicated by the comparisons;
- completion and numerical-noncompletion boundaries are explicitly reported rather than hidden;
- runtime results are reported using controlled paired measurements;
- new diagnostic metrics are summarized without retroactively chosen pass thresholds; and
- the team has a documented decision on whether to begin the separate FEMA P-695 collapse/statistical workflow.

## 11. Record-source status

The FEMA P-695 / ATC-63 far-field source set has now been acquired from the public SP3 / Haselton Baker Risk Group ground-motion mirror and archived in the project's private Google Drive data store.

Permanent public-repository provenance is recorded in:

`validation/yori_wall_10001/phase9n_record_source_manifest.json`

The archived package contains:

- the ATC-63/FEMA P-695 far-field summary workbook;
- the unscaled original PEER-NGA archive;
- the sorted 44-horizontal-component text archive, including `SortedEQFile_(120111).txt`, which is the Phase 9M comparison record;
- SHA-256 checksums; and
- an archive inventory.

Raw ground-motion bytes are intentionally not retained in the public QuakeCore repository. The manifest/parser and batch runner have been implemented and the initial 44-component suite has run.

## 12. September 19 reporting correction and current disposition

The archived run at `c7b0d6a23bbdb80e877f23a5ccd42b950bb1b3e8` has been
reaggregated without rerunning mechanics or modifying any response metric.
The completion taxonomy is exhaustive: 28 both completed, 5 both stopped at
the same nominal step, 2 both stopped at different steps, 9 OpenSees-reference
limited, and 0 QuakeCore-only limited. The two asymmetric dual stops are
`121011` (2374 / 1981 steps) and `121211` (4567 / 7283 steps), QuakeCore /
OpenSees respectively. Solver-specific step counts remain in each result.

Strict acceptance remains **43/44**, with `120521` failing the unchanged
minimum 1000-common-step requirement (798 available). The aggregator continues
to return a failing exit status. All original strict-window EDPs are unchanged.

The separate adaptive OpenSees diagnostic completes `120521`, but its archived
peak drifts imply relative differences of **0.312055%, 0.301707%, 0.268293%**,
using the same symmetric peak denominator as the strict comparator. Thus the
full-record diagnostic exceeds BOTH the frozen 0.1% peak-drift criterion and
the 0.1% history-RMSE criterion (maximum 0.367166%). Its maximum wall peak
force-EDP difference is 0.360549%, within 1%. The diagnostic now explicitly
reports the previously omitted peak-drift comparison.

Disposition: **strict reference-window failure with adaptive numerical-path
sensitivity; unresolved for full-record frozen parity**. This is not a pass,
not physical collapse, and not proof of a production mechanics defect. Do not
use the favorable force comparison to imply all engineering EDP gates pass.
The original 120521 diagnostic and strict result remain separate evidence.

Before using this mechanics path for collapse IDA, investigate the first
failed nominal interval (7.98–7.99 s), committed-state/branch events, and a
predeclared time-step refinement study. Keep interpolation, input phase,
constitutive laws, tolerances and thresholds fixed across comparison runs.
Report results even if refinement worsens agreement. No source-mechanics
change or relaxation of frozen criteria is justified by this audit alone.

The four original full-output artifacts have been hash-verified and preserved
in the existing private project archive. Shards 2 and 3 are stored as ordered
byte parts due to the transfer size limit; concatenation restores the original
artifact ZIP and SHA-256. The private archive ledger records file identities,
part hashes, byte counts, and reconstruction instructions. Future full-output
artifact retention is increased from 1 to 14 days.
