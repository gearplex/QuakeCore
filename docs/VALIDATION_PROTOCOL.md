# Validation protocol and commercial release gates

Evidence has four different purposes: (1) code verification against equations and solver identities; (2) numerical cross-checks against an independent engine; (3) experimental validation using the same excitation and specimen; and (4) an edition-specific building assessment. Passing an earlier level does not establish a later one.

## Berkeley sequence

1. Verify `validation/uc_berkeley_3story/input_manifest.json` using `tools/check_validation_input.py`. Preserve raw waveform bytes and hash, physical units, component orientation, table versus command motion, any similitude time scaling, processing, and initial state. The current manifest identifies a synthetic proxy. It cannot establish measured DT1 validation.
2. Freeze the archived model and reproduce its reported values before changing mechanics. Preserve the 39-inch clear-column versus 48-inch floor-to-floor drift convention. Compare *simultaneous* story profiles as well as independent peak envelopes.
3. Test component elasticity, monotonic backbone, cyclic unloading/reloading, force/state compatibility, axial intercepts, flow direction, unit changes and consistent tangents independently. In particular, the research Mroz implementation is not independently equivalent to Perform-3D.
4. Freeze physical parameters while testing solver strategy, initial Newton guess, local/global tolerances and time steps. Complete the entire input duration and report failed/partial histories explicitly. The 70-second proxy contains substantial excitation after 15 seconds; a 15-second window is not a complete record.
5. Compare measured and predicted displacement/drift histories, accelerations, restoring/base reactions with consistent definitions, residual drift, energy, and critical-column P–M/axial/rotation histories. Trace preloads and gravity equilibrium. Check whether axial degradation and gravity redistribution reproduce the observed mechanism.
6. Compare the exact Perform input/properties and recorded histories independently of the experimental comparison. Equivalent EDP peaks from different hysteresis, damping or excitation definitions do not demonstrate equivalence.
7. Only after these comparisons, apply independently traceable ASCE 41-23 / ACI 369.1-22 / AISC 342-22 rules and assess the intended hazard/performance objective. Store applicability decisions and resolved parameters with their edition, clause/table provenance, units and overrides.

Example reproducible proxy command:

```bash
python tools/check_validation_input.py validation/uc_berkeley_3story/input_manifest.json
build/ucb_phase9i validation/uc_berkeley_3story/proxy_dt1_motion.csv hybrid.json hybrid 15 same_pattern 1 2 previous_displacement
```

The final three arguments select exact modal damping tangent, integer time-step refinement, and Newton guess. No strength, stiffness, damping ratio or degradation parameter is fitted by this command. The result labels proxy status and the archived omission of the small nonzero t=0 acceleration. A proper measured-motion interface will require initial equilibrium and state transfer rather than inheriting that assumption.

## Commercial milestones and measurable exit criteria

| Priority | Deliverable | Exit evidence |
|---|---|---|
| P0 | Reproducible numerical release | Passing release and instrumented tests, deterministic input/output hashes, rollback and singular-factor tests, explicit failure taxonomy |
| P0 | Verified gravity-to-dynamic state | Static equilibrium solver, continuation, reactions, retained material/geometric state, tests of load redistribution and loss of support |
| P0 | Component reference library | Independent cyclic histories and energy comparisons for bilinear, degrading RC/steel, P–M and P–M–M, shear/axial loss, joints and walls; analytical or automatic consistent tangents |
| P0 | Experimental Berkeley validation | Identified recorded DT1 and measured channels, full-duration response, time-step/tolerance convergence, independently reviewed mechanism and parameter mapping |
| P1 | Independent NRHA benchmark suite | Matched model/record/precision/output requirements for low-, mid- and high-rise 3D models; separate integration, setup, output, memory and multicore measurements |
| P1 | More than 10× speed target | Median QuakeCore time <0.1× matched OpenSees across the declared representative workload, with accuracy thresholds and failed runs counted; repeat on controlled machines |
| P1 | Current-edition assessment service | Licensed/traceable parameter providers, applicability and units checks, edition migration, suite/hazard decisions, deformation/force-controlled assessment, reviewable calculation reports |
| P1 | Robust IDA and fragility | Appropriate intensity measure (e.g. Sa), documented collapse definition, record provenance, interval/right censoring, unresolved numerical cases, statistically correct censored fitting |
| P2 | User-facing model editor | 3D model/constraint validation, import diagnostics, units, provenance, reviewable changes, result visualization and independent verification views |
| P2 | Distribution and support | Stable versioned API/state/results, reproducible installers, dependency inventory, regression corpus, profiling/telemetry, job cancellation/checkpointing and support process |

Before GPU work, profile constitutive integration and sparse factor reuse. The fixed-modal tangent currently constructs a structural factor and low-rank workspace for each Newton correction. Reuse structural symbolic work, solve multiple modal right-hand sides together, and share immutable topology across independently owned record states. Replace finite-difference Mroz/corotational derivatives only after force/tangent and energy parity is established. Sparse modal extraction is needed for large models because the current LAPACK eigenanalysis densifies the system. Reduced-order or GPU implementations must retain the declared accuracy and failure behavior.

This sequence is an engineering development proposal, not a claim that the listed code provisions have already been implemented or verified.
