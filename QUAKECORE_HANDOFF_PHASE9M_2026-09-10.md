# QuakeCore Phase 9M candidate handoff — September 10, 2026

## Identity and status

- Base: verified Phase 9L source imported to `gearplex/QuakeCore` at commit `9ac3f8342535e27df3f7a2c6e0a876ecce479b70`.
- Candidate scope: YORi special nonbearing reinforced-concrete wall archetype 10001.
- Gates passed: corrected OpenSees baseline audit; matched-law solver/history comparison; gravity-state transfer and leaning-column P-Delta comparison.
- Open gate: exact ConcreteCM + Pinching4 + MinMax + Parallel constitutive parity.
- FEMA P-695 CMR/SSF/ACMR and system acceptance are not assessed.

## Source changes

- Added nonlinear static load stepping in `static_analysis`.
- Added Newmark initial displacement/velocity/acceleration, committed-state, and persistent-load inputs.
- Added gravity loads, controls, provenance, state transfer, and gravity-state modal analysis to `quake_run`.
- Added an optional OpenSees-compatible member `PDelta` sway transformation.
- Added a LAPACK dense factorization fallback when SuperLU is unavailable. This is a small-system verification path, not the production large sparse path.
- Added gravity-equilibrium/state-transfer regression coverage and updated the JSON job documentation.

## Quantitative evidence

Matched Concrete01/Steel01/bilinear shear, record 120111 at scale 0.03:

- Gravity residual infinity norm: `7.03e-11`.
- Maximum floor-displacement history difference: `5.35e-11 in`.
- Maximum story-drift history difference: `1.92e-13`.
- Period differences: approximately `2e-13 s` or less.
- Peak story drifts agreed to approximately `1e-14`.
- Seven-run median integration time: QuakeCore `0.1897 s`; OpenSees `0.1823 s`; ratio `1.041`.
- Release tests: 11 of 11 passed in the available environment, including the frozen material-oracle integrity check.

## Important audit findings in the received OpenSees workflow

1. Negative reinforcing-steel fourth envelope stress omitted multiplication by yield stress.
2. Rayleigh damping call omitted the computed committed-stiffness coefficient; the mass coefficient omitted a factor of two.
3. The mechanically corrected damping formulation gives about 14.7% first-mode damping with the selected anchors; the intended basis must be resolved before collapse analysis.
4. IDA refinement can alter scale without appending the actual scale to its output ledger.
5. Collapse extraction reconstructs IM from folder count rather than reading the actual applied scale.
6. Numerical failure, response-limit crossing, and successful completion are not kept distinct throughout the workflow.
7. Fragility dispersion and additional uncertainty terms are overwritten with fixed values without a traceable archetype decision record.

The normalized-spectrum column initially suspected as an error was checked and cleared: it is the median of the 44 record spectra and appears intentional for collective suite scaling.

## Reproduction

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
build/quake_run validation/yori_wall_10001/quakecore_job_10001_gravity.json result.json
python validation/yori_wall_10001/compare_gravity_matched_10001.py \
  --job validation/yori_wall_10001/quakecore_job_10001_gravity.json \
  --quake result.json --out comparison
```

The comparison script requires NumPy and OpenSeesPy 3.8.0. The repository does not redistribute the YORi model workbook or FEMA record files. The checked-in job embeds the single processed comparison acceleration history and records its provenance.

## Highest-priority next milestone

Create a material-level OpenSees oracle and independent QuakeCore implementations for ConcreteCM, Pinching4, MinMax, and Parallel. Exit criteria are force, tangent, branch-state, energy, commit/revert, and failure-state agreement under monotonic and cyclic protocols before these laws are admitted into the MVLEM archetype.
