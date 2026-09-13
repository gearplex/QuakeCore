# YORi wall archetype 10001 validation

This directory records the first QuakeCore comparison for the three-story YORi special nonbearing reinforced-concrete wall archetype 10001.

## What is established

The checked-in QuakeCore and independently assembled OpenSeesPy matched-law models use Concrete01, Steel01, bilinear shear, gravity loads, and a leaning-column P-Delta idealization. They agree to near machine precision in gravity displacement, modal periods, full displacement/drift histories, and peak story drift.

Gate 4 now also contains independently implemented, stateful ConcreteCM and Pinching4 laws plus MinMax and Parallel composition in the native wall-material abstraction. For the frozen story-1 OpenSees 3.8.0 oracle:

- unconfined and confined ConcreteCM complete the 241-point cyclic protocol;
- reinforcing-steel Pinching4 completes the full 281-point protocol through ten reversals;
- the sampled MinMax + Parallel reinforcing-steel stack reproduces the frozen lifecycle through permanent MinMax failure and the protocol endpoint;
- shear Pinching4 uses the recovered four-point envelope `(0.00010787,171.319)`, `(0.288,285.531)`, `(0.72,299.808)`, `(1.44,57.1062)` and completes the frozen 221-point shear protocol;
- ConcreteCM and Pinching4 are available through `WallUniaxial`, including explicit trial/commit state transfer;
- the story-1 received-law material stack has been admitted into native `Wall2D`/MVLEM component checks;
- the normal JSON runner can instantiate `concrete_cm`, `pinching4`, recursive `minmax`, and recursive `parallel` laws, including Pinching4 shear.

These results establish software/model implementation evidence for the stated protocols and component assembly. They do **not** establish physical wall-system validation, FEMA P-695 acceptance, collapse qualification, code approval, or a justified response-modification coefficient.

## Remaining exact-system boundary

The original `AllStoriesOpenseesModelINPUT.txt` used to generate the oracle is not checked into this repository. Its frozen SHA-256 is:

`02eca2fbdde5d2cc1400be354802a2b26d6d79ce15ac94dd981cd26dc6b23104`

That file contains the exact story-by-story ConcreteCM, Pinching4, and MinMax parameters. In particular, the frozen story-1 wrapper history proves that the positive MinMax limit lies in the sampled interval `(0.073, 0.080]`, but it does not uniquely recover the exact source value. Therefore the current sampled wrapper test uses a protocol-scoped limit only to verify the observed failure lifecycle; it is not evidence that the exact received MinMax limit has been recovered.

`materialize_received_job.py` is the deterministic bridge to the final exact three-story comparison. Given the original source file, it:

1. verifies the source SHA-256;
2. reads the three 41-value story rows;
3. preserves the checked-in archetype geometry, masses, gravity loading, leaning-column P-Delta system, record, and solver settings;
4. replaces every MVLEM concrete, reinforcing-steel, MinMax/Parallel, and shear law with the exact received row parameters; and
5. writes a normal QuakeCore job consumable by `quake_run`.

Until that SHA-matched source file is restored, do not label a three-story nonlinear run as exact received-material parity. Copying story-1 parameters into stories 2 and 3 would be an unsupported reconstruction.

## Evidence

- `matched_gravity_summary.json`: primary full-history matched-law comparison metrics.
- `runtime_repeats.json`: seven integration-only timings.
- `matched_no_gravity_summary.json`: Newton-path isolation benchmark.
- `quakecore_job_10001_gravity.json`: reproducible matched-law structural job with embedded record 120111.
- `compare_gravity_matched_10001.py`: independent matched-law OpenSeesPy comparison.
- `make_quakecore_gravity_job.py`: historical matched-law job generator; its original auxiliary package is not retained here.
- `opensees_yori_10001_audit.py`: exact received-material OpenSees audit when supplied the original input folder.
- `material_oracle.py`: Gate 4 OpenSees material-protocol generator.
- `material_oracle_story1.json`: frozen OpenSees 3.8.0 story-1 material oracle.
- `materialize_received_job.py`: SHA-gated exact received-law QuakeCore job materializer.

## Gate 4 sequence

1. Freeze YORi story-level ConcreteCM and Pinching4 parameters. **Story 1 complete; all-story source restoration still required.**
2. Generate virgin, monotonic, reversal, nested-cycle, limit-crossing, and repeated-cycle reference paths from OpenSeesPy. **Story-1 oracle frozen.**
3. Implement MinMax and Parallel composition with exact commit/revert semantics. **Implemented; sampled story-1 lifecycle verified.**
4. Implement ConcreteCM and Pinching4 independently from their published formulations and documented behavior. **Implemented for the frozen admitted protocols.**
5. Compare stress/force, consistent tangent, cumulative hysteretic energy, branch transitions, and permanent failure state. **Story-1 component protocols verified.**
6. Admit the materials into MVLEM only after component-level tolerances pass. **Story-1 component admission implemented; exact all-story system comparison awaits source restoration.**

## Current numerical tolerances

- Linear/gravity/modal/full-history matched-law comparison: near floating-point agreement is expected.
- Smooth material branches: relative force and tangent error target `1e-8`.
- Nonsmooth transition points: compare one-sided values and event location with explicit absolute tolerances.
- Hysteretic energy over a complete protocol: relative difference target `1e-6` before system integration.

Focused Gate 4 regression tests currently use tighter protocol-specific tolerances where the frozen oracle supports them. Tolerances may be tightened with better evidence, but they may not be loosened post hoc merely to pass an implementation.
