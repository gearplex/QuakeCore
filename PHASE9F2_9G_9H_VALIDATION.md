# QuakeCore Phase 9F.2 / 9G / 9H — P-M Hardening, Damping, and Remaining Story-2 Audit

**Status:** mechanics/damping validation increment; deterministic proxy input; not final external validation  
**Date:** 2026-09-06

## 1. Purpose

This increment follows Phase 9F, which showed that a true associated P-M hinge with a smooth Perform-type concrete interaction surface is a first-order correction for the UC Berkeley shake-table benchmark. The objectives here were:

1. add the documented Perform-style two-surface Mroz hardening option and determine whether it materially changes the Phase-9F result;
2. replace the Berkeley benchmark's prior Rayleigh approximation with the NIST-reported 2.5% modal + 0.5% Rayleigh damping formulation and isolate its effect;
3. investigate the remaining excess Story-2 response after the P-M correction;
4. preserve no-EDP-tuning and proxy-input caveats.

## 2. Phase 9F.2 — Mroz two-surface hardening

The P-M hinge now supports a translating inner Y surface and fixed homothetic U surface. The Y surface translates toward the conjugate point on U and transitions to EPP response once tangent to U. This is kept separate from the P-M surface definition and from ASCE/ACI parameter generation.

Over the same 15-s validation window:

| Case | S1 (%) | S2 (%) | S3 (%) | B1 max plastic rotation (rad) | B1 max |plastic axial| (in) |
|---|---:|---:|---:|---:|---:|
| P-M EPP | 5.1038 | 5.6998 | 2.9919 | 0.01828 | 0.16963 |
| P-M + Mroz | 5.0894 | 5.7159 | 2.9899 | 0.01842 | 0.17159 |

**Finding:** Mroz is a second-order correction for this benchmark over the validated window. The first-order mechanism change remains the associated P-M interaction itself.

## 3. Phase 9G — NIST damping parity

The NIST benchmark description reports **2.5% modal damping + 0.5% Rayleigh damping** for the Berkeley Perform3D model. The previous QuakeCore benchmark path represented the 2.5% part as mass-proportional Rayleigh damping. Phase 9G adds fixed-elastic-mode modal damping on the resisting-force side while retaining the small stiffness-proportional Rayleigh term.

| Case | S1 (%) | S2 (%) | S3 (%) | B1 max plastic rotation (rad) |
|---|---:|---:|---:|---:|
| Current Rayleigh approximation | 5.0894 | 5.7159 | 2.9899 | 0.01842 |
| NIST-style hybrid | 5.1152 | 5.7284 | 2.9223 | 0.02078 |

Global drift changes are small, but B1 plastic rotation increases by **12.8%** under the NIST-style hybrid. Modal damping is therefore a meaningful component-demand modeling choice even though it is not the first-order cause of the vertical drift discrepancy.

### 3.1 Pure 3% modal diagnostic

A 3% modal-only control remains numerically nonviable for this massless-DOF/concentrated-plasticity model, terminating near 2.45 s even after increasing the Newton/subdivision budget. QuakeCore also implemented an **optional exact low-rank modal-damping Newton Jacobian**:

\[
K_{d,modal} = a_1\sum_i 2\zeta_i\omega_i (M\phi_i)(M\phi_i)^T
\]

This term has rank equal to the number of damped modes and can be applied by Woodbury without forming a dense damping matrix. The modal-only control still terminates at the same state with this exact Jacobian, demonstrating that its failure is not caused simply by omitting modal damping from the Newton tangent. The exact low-rank tangent remains opt-in; the default modal implementation stays force-side-only to preserve Perform/OpenSees parity.

The practical benchmark formulation remains the NIST hybrid, which completes.

## 4. Phase 9H — remaining Story-2 state-path audit

A detailed committed-state audit was run with the Phase-9F EPP P-M mechanics and the prior Rayleigh damping formulation because those settings are computationally much cheaper. Phase 9F.2 and 9G show only small changes in the relevant global drift envelope, so the audit is used as a **mechanics surrogate**, not claimed as the exact Mroz+hybrid state snapshot.

At the Story-2 peak (about 12.32 s), the simultaneous drift profile is:

\[
[4.884,\ 5.688,\ 2.913]\%.
\]

The condensed committed unloading/reloading story tangent diagonal is approximately:

\[
[K_1,K_2,K_3]=[52.64,\ 39.17,\ 43.80]\ 	ext{kip/in}.
\]

Story 2 is therefore the softest of the three stories in this state.

The first positive tangent mode has roof-normalized floor shape approximately:

\[
\phi_1=[0.315,\ 0.768,\ 1.000]
\]

and interstory modal components:

\[
\Delta\phi_1=[0.315,\ 0.453,\ 0.232].
\]

Thus the modal Story-2/Story-1 interstory ratio is about **1.44**, compared with roughly 1.14 for the virgin Phase-9C first mode. The P-M state path has therefore evolved in a direction that **strengthens the Story-2 modal predisposition**, even while Story-1 component damage has increased substantially.

### 4.1 Associated-flow axial deformation by story

At the same Story-2 peak, the largest absolute hinge plastic axial deformation by story is:

| Story | Max |plastic rotation| (rad) | Max |plastic axial| (in) |
|---|---:|---:|
| 1 | 0.0348 | 0.2249 |
| 2 | 0.0283 | 0.1812 |
| 3 | 0.0194 | 0.1690 |

In particular, upper-story A/B hinges accumulate axial plastic deformation comparable to first-story hinges. Examples at the audited peak include approximately **0.18 in at B2 top** and **0.17 in at B3 top**, versus approximately **0.15 in at B1 bottom**.

This is consistent with Perform's own caution that associated P-M/P-M-M plasticity can generate unrealistic cyclic axial growth in reinforced-concrete columns. It is now a leading explanation for why the P-M correction simultaneously improves Story-1 response yet leaves too much Story-2 deformation.

## 5. Current interpretation

The evidence now ranks the mechanisms as follows:

1. **True associated P-M interaction:** first-order and essential. It moved Story-1 drift from about 3.3% to about 5.1%.
2. **Mroz hardening:** second-order for this benchmark window.
3. **Modal vs Rayleigh damping:** second-order globally, but important for local hinge demand and therefore worth preserving in production QuakeCore.
4. **Remaining Story-2 excess:** most strongly associated with the evolved vertical tangent distribution and substantial upper-story axial plastic flow generated by the P-M surface normal/cyclic path.

The next fidelity gate should therefore **not** be response tuning. It should be a surface-normal/cyclic-flow parity audit against the exact Perform/ASCE column P-M-M inputs where recoverable. Only after benchmark parity should QuakeCore consider a separate RC-enhanced, potentially non-associated or axial-flow-limited interaction law.

## 6. Verification

- `quake_tests`: PASS
- CTest: 100% PASS
- `ucb_validation_phase8`: builds successfully
- Phase 9F.2 Mroz executable: builds
- Phase 9G hybrid/modal executables: build
- Phase 9H state-path audit executable: builds

## 7. Validation limitations

- Ground motion remains the deterministic 1.52-g proxy, not the recovered recorded DT1 table acceleration.
- The current 15-s mechanics window contains the Phase-9F/9G peak cluster but is not a full-record validation.
- The detailed 9H component-state audit uses the computationally cheaper EPP/Rayleigh surrogate and is interpreted only for state-path diagnosis.
- Exact ASCE 41-17 / Perform3D proprietary component parameter sets are not embedded or inferred beyond publicly documented mechanics.
