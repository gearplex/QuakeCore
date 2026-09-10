# QuakeCore Phase 8 RC1 Release Notes

**Release:** `v0.8.0-rc1`  
**Benchmark:** UC Berkeley 3-story nonductile RC shake-table frame, Dynamic Test 1  
**Validation status:** **PRE-VALIDATION / PROXY INPUT**

Phase 8 completes the first mechanics-based flexure-shear-critical (FSC) column release path for QuakeCore. The release is numerically regression-tested and mechanism-tested against published component states, but it is intentionally **not** labeled final external DT1 parity because the recorded shake-table acceleration history has not yet been recovered.

## Accepted Phase 8 configuration

- 39 in clear flexible column length for the Berkeley specimen, with 4.5 in rigid joint offsets at each end of the 48 in joint-center story.
- Constant initial-stress P-Delta from the actual gravity preloads for the benchmark component-demand path.
- Dynamic `P(u)` geometric updating disabled in RC1 because the inherited research formulation altered the vertical equilibrium path and made the geometric-reference axial force unsuitable as a component axial resultant.
- Separate FSC zero-length shear springs in series with the flexural columns.
- Shear-degradation initiation based on local column end rotation with current axial compression and shear demand entering the limit.
- Cyclic shear-strength deterioration and explicit residual strength.
- E (lateral-resistance loss) tracked separately from F (effective/gravity-capacity loss).
- No automatic gravity redistribution after F; full static gravity re-equilibration remains future work.
- Experimental moving-`P` ASCE 41 flexural-capacity wrapper disabled in the release configuration pending a true moving-surface/return-mapping formulation.

## Phase 8 RC1 proxy-motion result

| Quantity | Phase 8 RC1 | Perform3D published | Measured DT1 |
|---|---:|---:|---:|
| First-mode period, s | 0.480 | 0.480 | 0.340 |
| Peak base shear, kip | 41.00 | 27.9 | 29.8 |
| V/W | 0.701 | 0.480 | 0.510 |
| Story 1 peak drift, % | 3.271 | 6.07 | 5.18 |
| Story 2 peak drift, % | 4.588 | 3.92 | 4.70 |
| Story 3 peak drift, % | 2.365 | 1.91 | 2.61 |
| Story 1 residual drift, % | 0.859 | 3.00 | 0.30 |
| Roof residual drift, % | 0.684 | 1.90 | 0.24 |

The Phase 8 result therefore remains too strong globally and underpredicts first-story peak drift relative to both external references. This is retained as a negative validation result rather than tuned away.

## B1 local failure-state validation

The strongest Phase 8 result is at the component level. Published experimental Test-1 B1 shear degradation began at approximately 3.15% first-story drift, 9.89 kip shear, and 24.7 kip axial compression. Phase 8 RC1 predicts:

| B1 initiation quantity | Phase 8 RC1 | Published Test 1 |
|---|---:|---:|
| First-story drift, % | 3.175 | 3.15 |
| Shear, kip | 9.792 | 9.89 |
| Axial compression, kip | 26.15 | 24.7 |
| Time in current proxy, s | 5.16 | not comparable without recorded table motion |

A1 uses the same nonductile-column properties but does not initiate until 19.43 s in the current proxy. B1-before-A1 is preserved across the Phase 8 post-failure sensitivity bracket. This is important because the experiment demonstrated different behavior in otherwise identically detailed A1 and B1 columns due to their loading and boundary conditions.

## Proxy response spectrum

The deterministic development motion has PGA = 1.52 g. Using a 5%-damped linear SDOF Newmark average-acceleration calculation, its pseudo-spectral accelerations are:

- `Sa(T = 0.48 s) = 1.9000 g`
- `Sa(T = 0.34 s) = 2.4229 g`

These values belong only to the current proxy. They must not be substituted for the spectrum of the recorded DT1 shake-table motion.

For context, the Phase 8 nonlinear peak `V/W = 0.701` corresponds to about 37% of the proxy elastic spectral coefficient at 0.48 s. If the same proxy demand were used illustratively with the published Perform3D `V/W = 0.48`, the implied force-reduction ratio would be about 3.96 versus about 2.71 for QuakeCore. This is not a formal R-factor calculation; it simply shows that the global base-shear discrepancy is not explained by elastic period alone and that QuakeCore remains less force-limited under this proxy.

## Numerical validation

- Existing QuakeCore test suite: PASS.
- New FSC shear-spring test covers pre-failure elasticity, E initiation, post-failure tangent finite differences, reversal degradation, and rollback/no-double-count behavior.
- SamePattern direct refactorization vs FullFactorization:
  - maximum roof-history difference: 0.0 in
  - relative roof-history difference: 0.0
  - maximum sampled base-shear-history difference: 0.0 kip
  - status: PASS
- Story-1 axial-force equilibrium in the accepted benchmark formulation:
  - minimum sum = 58.5 kip
  - maximum sum = 58.5 kip
  - reported maximum error = 0.0 kip

## Post-failure sensitivity bracket

Five parameter combinations were evaluated: the RC1 center point plus four corners spanning:

- post-failure stiffness / elastic shear stiffness: `-0.0025` to `-0.010`
- residual shear strength ratio: `0.10` to `0.30`

All analyses completed and B1 initiated before A1 in every case. Across the bracket:

- peak base shear = 40.48 to 41.08 kip
- Story 1 peak drift = 3.231 to 3.271%
- Story 2 peak drift = 4.588% in all cases to shown precision
- Story 3 peak drift = 2.357 to 2.381%
- Story 1 residual drift = 0.803 to 0.859%

The key global mismatch is therefore not sensitive to reasonable changes in the assumed post-failure shear slope or residual strength within this bracket.

## Important code correction found during Phase 8

The coupled axial-flexural tangent scratch storage in `CompiledFrame3D` was sized to the base ASCE 41 state count even though the experimental coupled wrapper stores one extra committed capacity scalar. Phase 8 replaces that fixed scratch array with correctly sized storage, removing a one-double out-of-bounds write risk. The experimental moving-surface wrapper nevertheless remains disabled because its constitutive formulation still requires a proper return-mapping treatment.

## Remaining blockers for final external validation

1. Recover the **recorded DT1 shake-table acceleration history** and replace the deterministic proxy.
2. Re-run the exact Phase 8 RC1 model without retuning global nonlinear EDPs.
3. Resolve remaining Perform3D component-definition differences where the NIST documentation is not machine-complete.
4. Implement production moving-surface P-M/P-M-M plasticity with consistent hysteretic return mapping.
5. Implement static gravity re-equilibration/redistribution after vertical-capacity F events before claiming progressive gravity-collapse simulation.

## Release artifacts

- `validation/uc_berkeley_3story/quakecore_phase8_release_result.json`
- `validation/uc_berkeley_3story/quakecore_phase8_validation_dashboard.html`
- `validation/uc_berkeley_3story/run_validation_phase8.cpp`
- `include/quake/fsc_shear_spring.hpp`
- `include/quake/fsc_shear_damage.hpp`
