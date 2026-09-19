# QuakeCore Phase 9M / Gate 4 final handoff — September 17, 2026

## Identity and final status

- Base: verified Phase 9L source imported to `gearplex/QuakeCore` at commit `9ac3f8342535e27df3f7a2c6e0a876ecce479b70`.
- Current candidate branch: `phase9m-yori-wall-10001`.
- Final stabilized Gate 4 head before this documentation update: `76dd01e9d671d82d2e2163ed5219bc1f355fdb67`.
- Candidate scope: three-story YORi special nonbearing reinforced-concrete wall archetype 10001.
- Gates 1–3: passed for the corrected/matched benchmark fixture.
- Gate 4: **passed for engineering parity of the documented reconstructed three-story surrogate** using the promoted ConcreteCM + Pinching4 + MinMax + Parallel implementation.
- Strict full-trace parity remains a diagnostic rather than the engineering acceptance criterion because the independent OpenSees system analysis does not complete the nominal 1.0x record under the current Newton settings.
- FEMA P-695 CMR/SSF/ACMR, collapse qualification, physical validation, code approval, and an R-factor recommendation are **not** established.

## Source and provenance boundary

The original received `AllStoriesOpenseesModelINPUT.txt` remains unavailable. Its expected SHA-256 is:

`02eca2fbdde5d2cc1400be354802a2b26d6d79ce15ac94dd981cd26dc6b23104`

The documented reconstruction used for current three-story engineering-parity work has SHA-256:

`22203092a1d125d59ced2a3777a9120ba9209d190ea32f6c71107f9c97403cf3`

The reconstruction is a **surrogate**, not the recovered/received source. In particular, the exact received MinMax limits remain unavailable; the reconstruction uses the documented protocol-scoped assumption of +/-0.075. Do not describe current three-story exact-law results as identity to the lost source.

## Production implementation now established

Gate 4 is no longer dependent on runtime source patches.

- The validated ConcreteCM/state-persistence implementation was promoted into production source in commit `008b6dc17853d547f9d8bf7c1f525fac304a57f8`.
- Production ConcreteCM and Pinching4 committed-deformation replay is idempotent: exact re-evaluation of a committed deformation returns the committed stress/tangent without advancing path-dependent state.
- ConcreteCM persisted state is now 37 slots; wrapper/MVLEM state offsets and regression tests use the promoted layout.
- Obsolete runtime patch scripts, the rule-77 prototype source, the dead `trial_legacy` shim, and the temporary test-migration workflow were removed.
- Final normal-source CI passes with sanitizers ON and OFF.

## Engineering-parity acceptance criteria

The required engineering gate in `validation/yori_wall_10001/assess_engineering_parity.py` requires:

- QuakeCore completion of the nominal reconstructed record;
- at least 1,000 independently converged common OpenSees steps;
- modal-period relative error <= `1e-6`;
- peak story-drift relative error <= `0.1%`;
- story-drift history RMSE <= `0.1%` of peak;
- wall peak shear/axial/bottom-moment EDP relative error <= `1%`.

The stabilized reconstructed workflow passes all required checks.

## Nominal 1.0x evidence

For the nominal reconstructed acceleration history:

- QuakeCore completes all 2,999 nominal steps.
- The independent OpenSees Newton analysis completes 1,145 nominal steps and stops at approximately 11.46 s after 100 iterations.
- Maximum relative difference in the first three modal periods: `1.52235e-12`.
- Maximum relative peak story-drift difference in the common response: `3.25121e-13`.
- Worst story-drift history RMSE: `0.0252463%` of peak.
- Five-MVLEM common-window force comparison covers wall shear, axial force, and bottom moment.
- Maximum relative peak wall-force EDP difference in that comparison: `1.47183e-12`.
- Maximum absolute six-component wall-force history difference: `1.40012e-08`.
- The first ConcreteCM material-replay difference greater than `1e-6` occurs at 11.46 s, coincident with the independent OpenSees system noncompletion.

The former 9.57 s discrepancy was traced to response extraction re-trialing an already committed path-dependent material and is resolved in production source.

## Multi-intensity diagnostic evidence

The same reconstructed model and record were exercised at dimensionless acceleration multipliers of 0.5x, 1.0x, 2.0x, and 4.0x. These are acceleration/PGA multipliers only; they are **not** Sa(T1), collapse capacities, FEMA P-695 intensity measures, or code qualification results.

| Multiplier | Completion | Peak-drift relative difference | Worst drift RMSE / peak | Max wall peak EDP relative difference |
|---|---|---:|---:|---:|
| 0.5x | QuakeCore and OpenSees both complete 2,999 steps | 0.00768% | 0.00409% | 0.0790% |
| 1.0x | QuakeCore 2,999; OpenSees 1,145 | ~3.54e-13 | 0.02525% | 0.117% |
| 2.0x | QuakeCore 2,999; OpenSees 1,018 | 0.00166% | 0.000446% | 0.0207% |
| 4.0x | Both reach the same 553 nominal steps before numerical noncompletion | 0.00383% | 0.00131% | 0.0767% |

The 0.5x case is the cleanest full-record independent two-solver exact-law comparison in the reconstructed suite. The 4.0x result is numerical noncompletion under the current settings and is not classified as collapse.

A separate OpenSees recursive-subdivision diagnostic reduced the nominal 0.01 s step to 1/64 at the 1.0x convergence boundary and still did not advance materially beyond approximately 11.453 s. That supports treating the nominal OpenSees stop as an algorithm/path limitation for this reconstructed model rather than evidence of a QuakeCore-only constitutive error.

## Runtime evidence

Two distinct benchmark generations should remain separated:

1. Earlier matched Concrete01/Steel01 benchmark: QuakeCore median `0.1897 s`, OpenSees median `0.1823 s`; QuakeCore was about 4% slower on that small/simple fixture.
2. Current exact-law 1,000-step paired/interleaved benchmark: QuakeCore median `0.212844 s`, OpenSees median `0.284012 s`; median-ratio speedup `1.334x`, median paired speedup `1.3319x`, paired range `1.302x–1.343x`.

The exact-law result supports a benchmark-specific roughly 33% advantage on this three-story application path. It is **not** a general QuakeCore performance factor; broader model-size, record-count, and runner-controlled testing remains required.

## Final CI evidence

At stabilized head `76dd01e9d671d82d2e2163ed5219bc1f355fdb67`:

- push numerical-verification run `35162881764`: sanitizer ON/OFF success; full CTest and smoke benchmarks pass;
- push reconstructed-yori-parity run `35162881785`: all engineering validation stages pass;
- PR numerical-verification run `35162884468`: success;
- PR reconstructed-yori-parity run `35162884498`: success;
- final push evidence artifact ID `10473991066`, SHA-256 `9630de38af30a44ebf3fab9c28d1a05f5fe83bfa09da80151e6469c7910f7611`.

## Gate 4 acceptance statement

Gate 4 is accepted for **software/model engineering parity of the documented reconstructed surrogate**. The required evidence includes modal response, global drift response, MVLEM shear/axial/moment demands, constitutive replay through the independently converged OpenSees window, normal-source regression coverage, and clean CI.

This does **not** establish:

- identity to the lost received input;
- physical wall/system validation;
- FEMA P-695 collapse statistics or acceptance;
- a qualified collapse intensity;
- code approval or prequalification;
- an approved or recommended YORi R factor.

## Next phase

The next development phase should broaden external validation rather than continue tuning this single nominal history:

1. freeze the current Gate 4 reconstructed-surrogate baseline and do not change its acceptance thresholds post hoc;
2. run a multi-record exact-law OpenSees/QuakeCore comparison suite with explicit record provenance and scaling;
3. compare modal response, drift, residual drift, acceleration, base/story shear, wall force/moment/axial demand, and dissipated-energy metrics over independently converged windows;
4. preserve distinct statuses for completed analysis, response-limit crossing, numerical noncompletion, and collapse classification;
5. benchmark runtime across records/model sizes using paired/interleaved runs and controlled output settings;
6. only after the record-suite evidence is stable, begin the FEMA P-695-specific statistical/collapse workflow with explicit censoring, nonmonotonic-response handling, and traceable intensity measures.

