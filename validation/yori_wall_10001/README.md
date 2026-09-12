# YORi wall archetype 10001 validation

This directory records the first QuakeCore comparison for the three-story YORi special nonbearing reinforced-concrete wall archetype 10001.

## What is established

The checked-in QuakeCore and independently assembled OpenSeesPy models use matched Concrete01, Steel01, bilinear shear, gravity loads, and a leaning-column P-Delta idealization. They agree to near machine precision in gravity displacement, modal periods, full displacement/drift histories, and peak story drift.

This establishes the equation assembly and numerical implementation for the stated idealization. It does not establish exact equivalence to the received YORi collapse model, physical validation, FEMA P-695 acceptance, or a justified response-modification coefficient.

## Evidence

- `matched_gravity_summary.json`: primary full-history comparison metrics.
- `runtime_repeats.json`: seven integration-only timings.
- `matched_no_gravity_summary.json`: Newton-path isolation benchmark.
- `quakecore_job_10001_gravity.json`: reproducible matched-law job with embedded record 120111.
- `compare_gravity_matched_10001.py`: independent OpenSeesPy comparison.
- `make_quakecore_gravity_job.py`: auditable job generator.
- `opensees_yori_10001_audit.py`: exact received-material OpenSees audit and cyclic steel check.
- `material_oracle.py`: first Gate 4 material protocol generator.

## Gate 4 sequence

1. Freeze YORi story-level ConcreteCM and Pinching4 parameters.
2. Generate virgin, monotonic, reversal, nested-cycle, limit-crossing, and repeated-cycle reference paths from OpenSeesPy.
3. Implement MinMax and Parallel composition with exact commit/revert semantics.
4. Implement ConcreteCM and Pinching4 independently from their published formulations and documented behavior.
5. Compare stress/force, consistent tangent, cumulative hysteretic energy, branch transitions, and permanent failure state.
6. Admit the materials into MVLEM only after component-level tolerances pass.

## Current numerical tolerances

- Linear/gravity/modal/full-history matched-law comparison: near floating-point agreement is expected.
- Smooth material branches: relative force and tangent error target `1e-8`.
- Nonsmooth transition points: compare one-sided values and event location with explicit absolute tolerances.
- Hysteretic energy over a complete protocol: relative difference target `1e-6` before system integration.

The exact acceptance tolerances may be tightened after the OpenSees oracle is frozen, but they may not be loosened post hoc merely to pass an implementation.
