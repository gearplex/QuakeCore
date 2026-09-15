# QuakeCore–OpenSees Wall Archetype 10001 Validation Report

**Date:** September 10, 2026  
**Scope:** Three-story special nonbearing reinforced-concrete wall archetype 10001, treated as the first building-frame-system comparison  
**Status:** Validation Gates 1–3 passed; Gate 4 is open; FEMA P-695 collapse/R-factor conclusions are not yet authorized

## Executive conclusion

The work establishes that QuakeCore can reproduce the three-story archetype's assembled structural equations, gravity equilibrium, gravity-state transfer, modal properties, leaning-column P-Delta response, and nonlinear time-history solution to essentially machine precision when both programs use matched constitutive laws.

For the matched gravity/P-Delta benchmark, the maximum difference between the full QuakeCore and OpenSees histories was:

- **Floor displacement:** 5.35×10^-11 in
- **Story drift ratio:** 1.92×10^-13
- **Peak story drift:** agreement to approximately 1×10^-14
- **Periods:** agreement to approximately 2×10^-13 s

Across seven timed repetitions of that benchmark, median integration time was **0.1897 s in QuakeCore** and **0.1823 s in OpenSees**, making QuakeCore **1.04× the OpenSees time** on this small model. This is functional speed parity, not yet evidence of superior large-suite throughput. The available environment used QuakeCore's dense LAPACK verification fallback because SuperLU was unavailable; the production sparse path should be benchmarked before making scaling claims.

The exact YORi collapse model is **not yet reproduced in QuakeCore**. It depends on ConcreteCM, Pinching4, MinMax, and Parallel material behavior. ConcreteCM and Pinching4 are path-dependent cyclic models whose unloading, reloading, pinching, and degradation rules can directly control collapse. Substituting Concrete01, Steel01, and bilinear shear is suitable for solver verification, but not for an R-factor determination. Accordingly, no CMR, SSF, ACMR, or FEMA P-695 pass/fail result is reported here.

## Archetype frozen for comparison

| Parameter | Archetype 10001 |
|---|---:|
| Stories | 3 |
| Story height | 144 in each |
| Wall length | 72 in |
| Wall thickness | 12 in |
| MVLEM fibers by story | 2, 2, 1 |
| Boundary length by story | 12.24, 12.24, 10.80 in |
| Boundary reinforcement ratio | 0.033347, 0.033347, 0.018868 |
| Web reinforcement ratio | 0.002564 |
| Test record | FEMA P-695 record 120111 |
| Test scale factor | 0.03 for matched-law comparisons |

The received model, the mechanically corrected OpenSees model, and the QuakeCore matched-law model are kept separate. The received source was not overwritten.

## Gate 1 — OpenSees audit and corrected baseline

### Received model response

| Result | Received OpenSees model |
|---|---:|
| T1 | 1.218047 s |
| T2 | 0.182694 s |
| Effective modal damping, mode 1 | 7.143% |
| Effective modal damping, mode 2 | 1.071% |
| Wall vertical base reaction after gravity | 555.975 kip |
| Leaning-column vertical base reaction after gravity | 1,650.713 kip |
| Pushover maximum base shear | 116.103 kip at 5.8 in roof displacement |
| First descending 80%-of-peak point | approximately 15.9 in |
| Record 120111 peak story drifts at scale 0.0785533 | 0.001367, 0.002851, 0.003509 |
| NRHA integration time | 1.069 s |

### Damping implementation correction

The received Tcl calculated a committed-stiffness coefficient but called:

`rayleigh $alphaM 0 $betaKinit 0`

This omitted the calculated committed-stiffness term. The mass coefficient also omitted the factor of 2 from the standard two-frequency Rayleigh expression. The mechanically corrected call is:

`rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm`

with the factor of 2 restored in `alphaM`.

This correction materially changed the one-record response: peak drifts became **0.001044, 0.002126, and 0.002539**, reductions of roughly 24%–28% from the received model.

However, this does not close the damping review. With the model's selected anchors (`TcoeffI = 0.2` applied to mode 1 and `TcoeffJ = 1.0` applied to mode 2), the corrected coefficients produce about **14.7% damping in mode 1** and **5.0% in mode 2**. That is mathematically consistent with those anchors, but the 14.7% first-mode damping is high enough that the intended damping basis must be confirmed before collapse analyses. The corrected file therefore implements the apparent formula; it is not an endorsement of the anchor selection.

### Reinforcing-steel cyclic branch correction

The generated negative steel envelope omitted multiplication by yield stress in the fourth negative stress point. Under a symmetric strain protocol, the received material produced approximately:

- +0.055 strain: +67.17 ksi
- -0.055 strain: **-0.67 ksi**

The corrected branch produces -67.17 ksi. The maximum stress difference between received and corrected negative branches was approximately **73.4 ksi** near -0.0529 strain. The associated fourth negative deformation point was also made symmetric with the regularized positive point. This is a critical model-generation error because it creates a strongly artificial tension/compression asymmetry in nominally symmetric reinforcing steel.

### Other corrected OpenSees workflow defects

- Removed an import that executed all initial pushovers as a side effect.
- Corrected the single-archetype generator to use the special-nonbearing design workbook.
- Added an environment-controlled archetype ID rather than an embedded ID.
- Corrected story-array indexing from `[counter-1, nST-1]` to `[counter, nST]`.
- Made output-folder creation explicit and copied the selected master model deterministically.
- Wrote only the number of periods actually calculated; the received code requested two eigenvalues but attempted to write six periods.
- Corrected pushover post-peak extraction to start at the peak and interpolate the first descending crossing of 0.8Vmax.
- Included the peak point in trapezoidal work integration.
- Replaced the brittle fixed-row gravity-weight read with the final base vertical reactions.
- Changed the `C0` selection logic to mutually exclusive branches.
- Removed a hard-coded six-story total height and interpolated the ultimate displacement.

## Gate 2 — Solver-history and speed reconciliation

The earlier full-history difference was traced to the nonlinear search path rather than ground-motion timing, response extraction, or an equilibrium error. QuakeCore's line search selected a slightly different iterative path than OpenSees full Newton. When QuakeCore line search was disabled to match the OpenSees algorithm:

| Matched no-gravity metric | Result |
|---|---:|
| Maximum full-history drift difference | 1.69×10^-5 |
| Story drift RMSE / peak | 0.107%, 0.066%, 0.067% |
| Correlation by story | greater than 0.9999966 |
| Peak drift difference | approximately 1×10^-11 |

The remaining tiny trajectory difference is ordinary nonlinear path/numerical sensitivity. It is not an engineering-significant mismatch for this fixture.

One QuakeCore statistic was also misleading: `direct_fallbacks` increased whenever state-dependent wall tangents required an ordinary direct refactorization. It did not mean thousands of solver failures. The label should be revised before using telemetry for performance or reliability conclusions.

## Gate 3 — Gravity-state transfer and P-Delta implementation

QuakeCore was extended with:

- Incremental nonlinear static gravity analysis.
- Transfer of converged displacement and committed material state into NRHA.
- Persistence of the gravity load in the dynamic residual.
- Modal analysis from the gravity-state tangent.
- An OpenSees-compatible plane-frame PDelta sway transformation option for elastic members.
- JSON inputs for nodal gravity loads, load steps, tolerance, maximum iterations, provenance, and member PDelta flags.
- A gravity-equilibrium regression test.
- A dense LAPACK factorization fallback for verification environments without SuperLU, including stale-factor invalidation testing.

The initial full beam-column geometric-stiffness trial produced a false instability in the intentionally near-zero-flexural-stiffness leaning column. It was replaced with the sway term that matches OpenSees' `PDelta` transformation for this use. The generic initial SPD screening was disabled for this particular fixture because its near-zero EI, massless rotational DOFs fall below the generic pivot threshold. The model nonetheless passed gravity equilibrium, modal, and full dynamic comparisons. This is a modeling/diagnostic caveat that should remain explicit.

### Matched-law gravity/P-Delta comparison

| Quantity | QuakeCore | OpenSees |
|---|---:|---:|
| Gravity residual infinity norm | 7.03×10^-11 | converged |
| T1 | 1.255015243746 s | 1.255015243746 s |
| T2 | 0.184586584969 s | 0.184586584969 s |
| T3 | 0.105104604056 s | 0.105104604056 s |
| Peak drift, story 1 | 0.0009702078290 | 0.0009702078290 |
| Peak drift, story 2 | 0.0020461091341 | 0.0020461091341 |
| Peak drift, story 3 | 0.0024632822729 | 0.0024632822728 |

The gravity displacements also matched directly, including wall-floor vertical displacements of -0.0161493, -0.0267132, and -0.0315603 in and leaning-column displacements of -2.3770×10^-5, -3.9522×10^-5, and -4.7257×10^-5 in.

All ten QuakeCore test targets passed after the changes.

## Analysis speed

Seven integration-only repetitions of the final matched-law benchmark produced:

| Engine | Median | Range |
|---|---:|---:|
| QuakeCore, dense LAPACK verification path | 0.1897 s | 0.1814–0.2077 s |
| OpenSeesPy 3.8.0 | 0.1823 s | 0.1749–0.1852 s |

QuakeCore was 4.1% slower at the median on this three-story fixture. That is effectively parity for the immediate validation question, but it is too small a model to assess sparse-solver scaling or parallel IDA throughput. The earlier packaged executable was about 11.6× slower than OpenSees in a different environment; the current result shows that the large deficit is not inherent to the governing equations, but a production SuperLU benchmark and multi-record wall suite are still needed.

## Confirmed and unresolved errors in the received evaluation

| Priority | Finding | Consequence / required action |
|---|---|---|
| Critical | Negative reinforcing-steel envelope omitted the yield-stress multiplier and used an inconsistent fourth deformation point | Artificially removes negative-branch strength; corrected in the supplied generator and material file |
| Critical | Rayleigh command did not apply the computed committed-stiffness coefficient; mass coefficient omitted a factor of 2 | Changes modal damping and one-record drifts by roughly 24%–28%; implementation corrected, but damping anchors still require an engineering decision |
| Critical | IDA “refinement” increases intensity after drift exceeds twice the limit and does not append the altered intensity to `IM_values` | Response and recorded IM can become inconsistent; replace with a true lower/upper collapse bracket and persist every actual scale |
| Critical | Residual/collapse extraction reconstructs incipient-collapse IM as `IMincr*(len(IMList)-2)` rather than reading the actual stored scale | Invalid when velocity normalization or any refined scale is used; extract the recorded IM directly |
| High | Missing drift files are used as the principal success test | A numerically failed or prematurely terminated analysis can be treated as successful if files exist; require explicit completion and convergence status |
| High | Drift-limit termination and numerical convergence are conflated in the dynamic Tcl workflow | Preserve distinct statuses for completed, threshold crossing, dynamic instability, and unresolved numerical failure |
| High | Recovery procedures change algorithms, tests, tolerances, and time steps without a fully auditable accepted-step ledger | Collapse capacity may become solver-policy dependent; log every recovery branch and restore the original analysis controls exactly |
| High | Fragility code fits dispersion and then overrides collapse dispersion with 0.4; additional uncertainty components are hard-coded as 0.2, 0.2, 0.2 | FEMA P-695 acceptance can change materially; document and derive each uncertainty term rather than silently overriding the fit |
| High | The fragility script sets `R = 5.0` while the selected design workbook is labeled `R=6` | Latent configuration defect even where current algebra happens to cancel; use one traceable archetype parameter source |
| Medium | Received pushover extraction could select ascending points, omit the peak from work, read weight from a fixed row, overwrite `C0`, and use a six-story height | Corrected in supplied pushover scripts; rerun all archetype summaries |
| Medium | Story drift tracking assumes equal story heights | Acceptable for archetype 10001, unsafe for generalization; calculate each drift from node elevations |
| Low | Several validation executables compile with “control reaches end of non-void function” warnings | Existing research-code hygiene issue; correct before treating the full repository as release-quality |

### Scaling issue investigated and cleared

`P695Far FieldNormalizedResponseSpect.txt` has 46 columns: period, the median normalized spectrum, and 44 individual record spectra. At T = 0.41 s, column 2 is 0.82851156 and equals the median of the 44 individual values to the displayed precision. Therefore, use of `NormalizedResposeSpect[idx,1]` as a common suite denominator appears intentional for collective record-set scaling and is **not listed as an error**. This should nevertheless be documented in code because the variable name `GM_SaT1` suggests a record-specific value.

## Gate 4 — Exact constitutive parity remains open

The received wall fibers use ConcreteCM in Parallel with reinforcing steel that includes Pinching4 and MinMax behavior. A collapse comparison cannot be considered analogous until the following are verified against OpenSees under monotonic and reversing protocols:

1. ConcreteCM compression and tension envelopes.
2. Concrete unloading/reloading, crack closure, and cyclic degradation.
3. Pinching4 positive and negative envelopes, pinched reload paths, and all selected degradation modes.
4. MinMax failure-state behavior and state persistence.
5. Parallel force summation and tangent/state commit/revert behavior.
6. MVLEM fiber strain mapping and shear-flexure interaction using those exact laws.

The current Concrete01/Steel01/bilinear-shear model is intentionally labeled **matched-law verification**, not the final YORi archetype. A clean-room implementation with point-by-point OpenSees material tests is the preferred path. A calibrated surrogate could be used only if it matches cyclic energy, residual deformation, strength loss, unloading stiffness, and collapse-driving response over representative wall protocols and the substitution is accepted as part of the FEMA P-695 methodology.

## Gate status and next sequence

| Gate | Status |
|---|---|
| 1. Freeze/correct OpenSees 10001 | Passed, with damping-anchor decision still flagged |
| 2. Resolve solver/history delta and speed | Passed for matched laws |
| 3. Gravity-state transfer and P-Delta | Passed for the benchmark fixture |
| 4. Exact cyclic material parity | Open; current technical stop |
| 5. Full exact-law QuakeCore/OpenSees comparison | Not run pending Gate 4 |
| 6. FEMA P-695 suite, CMR, SSF, ACMR, acceptance | Not run pending Gates 4 and 5 |

The next defensible work item is a material-level validation harness, followed by one exact-law record at multiple intensities, then the 44-record suite. FEMA P-695 statistics should only be generated after collapse status, record scaling, censoring, nonmonotonic response, and numerical failures are stored explicitly for every run.

## Deliverable contents

- `OpenSees_corrected_10001/`: corrected three-story model and the single comparison record.
- `QuakeCore_10001_dev_source.zip`: QuakeCore source with gravity-state transfer, PDelta member option, dense verification fallback, documentation, and tests; build products excluded.
- `benchmarks/gate1/`: exact received-versus-corrected OpenSees audit.
- `benchmarks/gate2/`: matched no-gravity QuakeCore/OpenSees histories and job.
- `benchmarks/gate3/`: final gravity/P-Delta histories, summary, and seven-run timing data.
- `scripts/`: OpenSees translation, QuakeCore job generator, and comparison harness.

## Source references for constitutive implementation review

- OpenSees `ConcreteCM` source: https://github.com/OpenSees/OpenSees/blob/master/SRC/material/uniaxial/ConcreteCM.cpp
- OpenSees `Pinching4` documentation: https://opensees.berkeley.edu/wiki/index.php/Pinching4_Material
- OpenSees Rayleigh command: https://opensees.berkeley.edu/wiki/index.php/Rayleigh_Damping_Command

