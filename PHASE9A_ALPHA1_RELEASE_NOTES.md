# QuakeCore Phase 9A Alpha 1 — Excitation / Parity Infrastructure

**Status:** infrastructure complete; exact Llolleo source bytes not yet ingested in this environment.

## Scope

Phase 9 freezes the Phase 8.1 RC1 structural model and changes excitation/provenance only. No hinge, stiffness, panel-zone, P-Delta, FSC, damping, or solver parameters are calibrated to response-history EDPs in this phase.

## Verified external input identity

The DT1 source is the 3 March 1985 Valparaiso (Central Chile) earthquake, Llolleo station, S80E / Component 100. The historical NOAA record is `chi01.055` / `chi01.055.100`; a research copy is commonly named `LLOLLEO100.th`.

Independent metadata used as ingestion gates:

- source PGA: 436.9 cm/s^2 (~0.4455 g)
- source sample interval: 0.005 s
- source duration: 116.40 s
- DT1 nominal amplitude scale: 4.06
- model time scale: divide time by sqrt(3)

The nominal source-motion reconstruction therefore has expected PGA ~1.809 g, dt ~0.002887 s, and duration ~67.20 s. This is **not** the acceleration actually recorded on the shake table, whose published PGA is approximately 1.52 g.

## New Phase 9 tooling

- `prepare_llolleo100.py`: strict ingestion/provenance tool. It will not guess units and rejects mismatched source PGA, dt, or duration by default.
- `compare_spectra.py`: independent 5%-damped pseudo-Sa calculator for QuakeCore CSV motions.
- `test_phase9_tools.py`: transformation/metadata regression test using a synthetic series solely as a software test fixture.
- `proxy_spectrum_5pct.csv`: frozen Phase 8.1 proxy spectrum from 0.05–2.0 s.
- `generated/phase9_proxy_frozen_replay.json`: exact frozen-model replay baseline.

## Frozen proxy replay

The Phase 8.1 structural model was rebuilt and replayed without source changes:

- PGA = 1.52 g
- Sa(0.48 s, 5%) = 1.900001 g
- Sa(0.34 s, 5%) = 2.422924 g
- peak first-story column restoring shear = 25.7371 kip
- peak inertial base reaction = 38.5263 kip
- pushover peak = 23.5553 kip
- FullFactorization / SamePattern recorded-history difference = 0

This is the control state for the exact-source reconstruction run.

## External data blocker

The exact file is listed in NOAA Volume 3 and has been preserved in two independent 2025 rescue archives:

1. Stanford/Data Rescue Project PURL `dx674py2928`
2. PANGAEA dataset 981245, Volume 3 archive `eq_strong_motion_v3.tar.gz`

The current execution environment can enumerate those archives but cannot fetch the large binary/archive payload. Phase 9 deliberately does not digitize or synthesize a replacement waveform.

## Next acceptance gate

When `chi01.055` or a byte-equivalent `LLOLLEO100.th` is available:

1. pass the metadata gate;
2. create source x4.06 / time/sqrt(3) command reconstruction;
3. compute 5% spectra and compare with proxy;
4. run frozen Phase 8.1 model;
5. report structural restoring shear, inertial reaction, story/residual drifts, B1/A1 initiation, solver invariance;
6. do not change model parameters in response to the result.

Only a later run using the **recorded shake-table acceleration** should be labeled full DT1 input parity.
