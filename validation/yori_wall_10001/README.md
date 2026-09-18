# YORi wall archetype 10001 validation

This directory records QuakeCore/OpenSees verification for the three-story YORi special nonbearing reinforced-concrete wall archetype 10001.

## Current status — September 17, 2026

Gate 4 is accepted for **engineering parity of the documented reconstructed three-story surrogate**. The accepted evidence covers:

- independently implemented ConcreteCM and Pinching4 laws;
- MinMax and Parallel composition with committed/trial state persistence;
- native `WallUniaxial` and `Wall2D`/MVLEM integration;
- modal properties;
- story drift histories and peak drifts;
- MVLEM shear, axial force, and bottom-moment engineering demand parameters;
- material replay against OpenSeesPy 3.8.0 over the independently converged comparison window;
- multi-intensity diagnostic comparisons; and
- paired exact-law runtime benchmarking.

The strict full-trace comparator remains in the workflow as a **non-blocking diagnostic**. It is not the engineering acceptance gate because the independent OpenSees nominal 1.0x system analysis stops at approximately 11.46 s under the current Newton settings, while QuakeCore completes the full 2,999-step record.

These results establish software/model verification for the reconstructed surrogate. They do **not** establish physical wall-system validation, FEMA P-695 acceptance, collapse qualification, code approval, or a justified response-modification coefficient.

## Provenance and exact-system boundary

The original `AllStoriesOpenseesModelINPUT.txt` used to generate the received model is not available in the repository. Its frozen expected SHA-256 is:

`02eca2fbdde5d2cc1400be354802a2b26d6d79ce15ac94dd981cd26dc6b23104`

The documented three-story reconstruction used for current system-level engineering parity has SHA-256:

`22203092a1d125d59ced2a3777a9120ba9209d190ea32f6c71107f9c97403cf3`

The reconstruction is a **surrogate** and must not be described as the recovered, exact, or received source.

The exact received MinMax limits remain unavailable. The reconstruction uses the documented protocol-scoped assumption of +/-0.075. Story 2 uses the same active row as Story 1 based on the surviving inputs; Story 3 uses the documented reconstructed roof properties. These assumptions are recorded in `reconstructed_input_provenance.json`.

If the original source is recovered, `materialize_received_job.py` remains the SHA-gated bridge to an exact received-law comparison.

## Component constitutive evidence

For the frozen OpenSees 3.8.0 reference protocols:

- unconfined and confined ConcreteCM complete the 241-point frozen cyclic protocol;
- reinforcing-steel Pinching4 completes the full 281-point protocol through ten reversals;
- the sampled MinMax + Parallel reinforcing-steel stack reproduces the frozen lifecycle through permanent MinMax failure;
- shear Pinching4 completes the frozen 221-point shear protocol;
- additional ConcreteCM reversal paths required by the reconstructed three-story NRHA were independently admitted against OpenSees behavior;
- production ConcreteCM and Pinching4 committed-deformation replay is idempotent and leaves the committed state unchanged.

The validated ConcreteCM/state-persistence implementation is now in normal production source. The reconstructed validation workflow applies **no runtime source patches**.

## Engineering-parity gate

`assess_engineering_parity.py` requires:

| Criterion | Required threshold |
|---|---:|
| QuakeCore nominal record completion | complete |
| Independently converged common OpenSees history | >= 1,000 nominal steps |
| Modal-period relative error | <= 1e-6 |
| Peak story-drift relative error | <= 0.1% |
| Drift-history RMSE / peak | <= 0.1% |
| Wall peak shear/axial/bottom-moment EDP relative error | <= 1% |

The stabilized reconstructed workflow passes every required criterion.

### Nominal 1.0x result

- QuakeCore: 2,999 nominal steps completed.
- OpenSees: 1,145 nominal steps completed; Newton noncompletion at approximately 11.46 s.
- Maximum relative difference in first three periods: `1.52235e-12`.
- Maximum relative peak drift difference: `3.25121e-13`.
- Worst drift-history RMSE: `0.0252463%` of peak.
- Maximum relative peak wall EDP difference in the dedicated five-MVLEM force comparison: `1.47183e-12`.
- Maximum absolute six-component wall-force history difference: `1.40012e-08`.

The first ConcreteCM material-replay discrepancy greater than `1e-6` occurs at 11.46 s, coincident with the independent OpenSees system noncompletion.

## Multi-intensity diagnostic

Scale factors below are dimensionless acceleration/PGA multipliers. They are **not** Sa(T1), collapse capacities, FEMA P-695 intensity measures, or code qualification results.

| Multiplier | Solver completion | Peak-drift relative difference | Worst drift RMSE / peak | Max wall peak EDP relative difference |
|---|---|---:|---:|---:|
| 0.5x | both complete 2,999 steps | 0.00768% | 0.00409% | 0.0790% |
| 1.0x | QC 2,999; OS 1,145 | ~3.54e-13 | 0.02525% | 0.117% |
| 2.0x | QC 2,999; OS 1,018 | 0.00166% | 0.000446% | 0.0207% |
| 4.0x | both reach 553 nominal steps, then numerical noncompletion | 0.00383% | 0.00131% | 0.0767% |

The 0.5x case is the cleanest full-record independent two-solver exact-law comparison currently available.

A separate OpenSees recursive-subdivision diagnostic reduced the nominal step to 1/64 near the 1.0x stopping point and still did not advance materially beyond approximately 11.453 s. This supports treating the nominal OpenSees stop as an algorithm/path limitation for this reconstructed model rather than evidence of a QuakeCore-only constitutive defect.

## Runtime evidence

Keep the two benchmark generations separate:

- Earlier matched Concrete01/Steel01 benchmark: QuakeCore `0.1897 s`, OpenSees `0.1823 s` median; QuakeCore approximately 4% slower.
- Current paired/interleaved 1,000-step exact-law benchmark: QuakeCore `0.212844 s`, OpenSees `0.284012 s`; median-ratio speedup `1.334x`, median paired speedup `1.3319x`, paired range `1.302x–1.343x`.

The current result supports a benchmark-specific roughly 33% QuakeCore speed advantage on this application path. It is not a general performance factor.

## Important files

- `reconstructed_input_provenance.json`: source reconstruction basis and explicit assumptions.
- `AllStoriesOpenseesModelINPUT.reconstructed_10001.txt`: documented surrogate input.
- `materialize_received_job.py`: SHA-gated exact-source materializer if the lost received input is recovered.
- `material_oracle_story1.json`: frozen OpenSees 3.8.0 Story-1 material oracle.
- `diagnose_material_history_parity.py`: material-history replay diagnostic.
- `compare_wall_force_edps.py`: common-window MVLEM force-demand comparison.
- `assess_engineering_parity.py`: required engineering acceptance gate.
- `compare_reconstructed_10001.py`: strict full-trace diagnostic.
- `benchmark_reconstructed_prefix.py`: equal-workload exact-law runtime benchmark.
- `compare_multi_intensity.py`: 0.5x/1x/2x/4x engineering-parity diagnostic.

## Final CI baseline

The stabilized source-only Gate 4 baseline passed:

- numerical verification with sanitizers OFF and ON;
- full CTest and smoke benchmarks;
- full QuakeCore nominal reconstructed record;
- material replay;
- common-window wall-force EDP comparison;
- required engineering-parity assessment;
- strict comparator diagnostic execution;
- robust OpenSees diagnostic;
- multi-intensity diagnostic; and
- paired exact-law runtime benchmark.

Final stabilized branch head before documentation-only commits: `76dd01e9d671d82d2e2163ed5219bc1f355fdb67`.

Final push evidence artifact: ID `10473991066`, SHA-256 `9630de38af30a44ebf3fab9c28d1a05f5fe83bfa09da80151e6469c7910f7611`.

## Next phase

Freeze this Gate 4 baseline and broaden external validation across multiple records and model demand levels. The next suite should preserve explicit record/scaling provenance and compare drift, residual drift, acceleration, base/story shear, wall force/moment/axial demand, and dissipated energy over independently converged windows. Numerical noncompletion, response-limit crossing, successful completion, and collapse classification must remain separate statuses.

FEMA P-695 collapse statistics should only follow after the multi-record external-validation evidence and intensity-measure/censoring workflow are stable.
