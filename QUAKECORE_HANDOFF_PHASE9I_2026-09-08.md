# QuakeCore handoff — Phase 9I, 2026-09-08

## Current state

Updated source: `quakecore-phase9i.zip`, expanding to `quakecore-phase9i/`. Evidence: `quakecore-phase9i-evidence.zip`. Start with `docs/PHASE9I_REVIEW.md`, `docs/JOB_FORMAT.md`, and `docs/VALIDATION_PROTOCOL.md`.

The original Phase 9F2/9G archive was verified against all six supplied SHA-256 values and reproduced before changes. The new research version fixes NaN convergence, local-integration rollback, EPP/Mroz axial-domain violations, Mroz length-unit dependence, mixed-unit axial return tolerances, IDA first-bracket corruption, stale adaptive factors, failed SuperLU factors, lost substep demand peaks and an extra sample in the benchmark window. It caches accepted line-search evaluations, skips elastic Mroz tangent probes, reuses IDA preparation and includes a plane-frame JSON runner. All four Release and ASan/UBSan CTest suites pass. Leak detection is blocked by sandbox process/thread inspection and was disabled in the successful sanitizer run. Hosted CI is configured but unexecuted.

## What is established

Independent OpenSees 3.8.0 comparison passes for 3-story/39-DOF and 10-story/130-DOF bilinear frame fixtures. Full displacement/drift/acceleration histories and periods agree near floating-point precision; the hinges yield. Seven-repetition median integration speedups without optional recorders are 13.46× and 15.34×. With recording included they are 9.73× and 10.67×. Model creation/eigenanalysis and final QuakeCore JSON serialization are excluded. These small planar tests establish neither an end-to-end commercial speed claim nor degrading P–M/IMK/3D parity.

The synthetic IDA example returns threshold brackets [2.0, 2.1810154653] for its waveform and polarity reversal. Its 1.5% drift threshold is only a bracket-workflow demonstration, not a collapse or ASCE acceptance definition.

## Berkeley validation: keep these distinctions

- Input is a deterministic 70-second, 1.52 g **proxy**, not recorded DT1. SHA-256: `d0a01ac7694a708f232b6ee03c0b090e2849468243cd65b1f0ee657db02e8413`.
- Earlier “15 s” code advanced 15.01 s. New code advances exactly the requested endpoint count.
- The first 15 seconds include only 34.2522% of the proxy's ∫a²dt using trapezoidal integration. Do not compare that window to full-record experimental peaks as validation.
- Archived zero initial state omits the input's small nonzero initial acceleration (-1.1994e-5 g). New metadata exposes this assumption. A general measured-motion workflow needs explicit initialization/state transfer.
- Fixed physical settings: 39-inch clear-column drift denominator (48-inch floor spacing); existing mass, section/FSC, constant-preload and damping definitions; stiffness multiplier 0.934616. No EDP refitting occurred.
- Preferred *diagnostic* Newton starting guess is `previous_displacement`. The earlier kinematic predictor converges poorly in these nonlinear cases. The choice does not constitute an accuracy guarantee.
- EPP full and same-pattern at 10 ms complete 15 s with identical outputs: drift [5.592322, 6.215401, 3.211241]%.
- Mroz/Rayleigh: 10 ms fails at 13.45 s; 5 ms completes 15 s with [5.439878, 6.146505, 3.187661]% (use exact JSON); 2.5 ms fails at 7.83 s.
- Hybrid 2.5% modal + 0.5% initial-stiffness Rayleigh: 10/5 ms complete 15 s; 2.5 ms fails at 13.0075 s. At 5 ms the peaks are [5.535698, 6.194293, 3.139610]% (exact JSON is authoritative).
- Pure 3% modal remains a numerical failure at 5.15 s for 10 and 5 ms.
- A requested 70-second hybrid run at 5 ms fails at 19.71 s, with partial peaks [7.128153, 7.027405, 3.139610]%. This is a numerical failure, not a proven physical collapse.
- Time-step convergence and experimental validation remain open. Failed cases must remain in summaries.

## Build and execute

Prerequisites: C++20, CMake >=3.20, SuperLU headers/library, BLAS/LAPACK, threads; optional Python/NumPy/OpenSeesPy for cross-checks.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
build/quake_run examples/frame3_nrha.json frame3_result.json
build/quake_run examples/frame3_ida.json frame3_ida_result.json
python tools/check_validation_input.py validation/uc_berkeley_3story/input_manifest.json
build/ucb_phase9i validation/uc_berkeley_3story/proxy_dt1_motion.csv hybrid.json hybrid 15 same_pattern 1 2 previous_displacement
python tools/compare_opensees.py examples/frame10_nrha.json ./build/quake_run parity10 --repeats 7 --no-timing-recorders
```

The Berkeley trailing arguments are strategy, exact-modal-tangent switch, time-step refinement factor, and Newton guess. Full reproducible commands are preserved in each result's timing JSON. The JSON runner is currently 2D and has no gravity equilibrium stage or compliance assessment. IDA reports scale/PGA, not automatically Sa(T1). Raw accepted material-state arrays are for forensic inspection; they are not a complete restart state.

## Next bounded engineering task

Capture a complete rollback-safe checkpoint immediately before the 5 ms hybrid 19.71 s failure and before the 2.5 ms failures: displacement, velocity, acceleration, committed component states, original/substep time and excitation, tangent, residual by DOF, branch events and constraint map. Replay one failed step under direct Newton with consistent local/global tolerances. Separate massless-DOF equilibrium, softening branch transitions and approximate Mroz tangent effects. Test a consistent component Jacobian and bounded material substepping before making physical modeling changes. Do not loosen acceptance, introduce damping or fit strengths merely to force completion.

The other prerequisite is to obtain the exact recorded DT1 table acceleration, measured channel histories and Perform-3D model/property/damping/gravity export. No public digital DT1/model package was recovered in the present search. The full NIST report download was blocked/too large, but the NIST companion paper and official PEER/CSI sources were inspected. See source URLs in the review. The NIST benchmark evaluated ASCE 41-17; current ASCE 41-23 / ACI 369.1-22 / AISC 342-22 compliance requires separately traceable provisions.

Avoid beginning with a GUI or GPU rewrite. The highest-value next work is trustworthy full-record behavior, component parity, gravity-state transfer and a larger matched 3D benchmark. Preserve the independent OpenSees fixture and numerical-failure tests as regression gates.
