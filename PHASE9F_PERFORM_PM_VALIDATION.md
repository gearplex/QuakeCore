# QuakeCore Phase 9F — Perform-Type P–M Mechanics Validation

**Status:** Mechanics validation increment; not a final external-validation release  
**Benchmark:** UC Berkeley 3-story RC shake-table frame  
**Input:** deterministic DT1 proxy, not recovered recorded shake-table motion  
**Validated nonlinear window:** 0–15 s

## 1. Purpose

Phase 9F tests whether the missing Story-1 drift localization in QuakeCore is materially affected by the column axial-force/flexure interaction mechanics used by Perform3D. The experiment changes the column hinge from a scalar moment hinge with external P-dependent strength to a true two-generalized-force zero-length P–M hinge.

The mechanics are intentionally isolated from code-table calibration and EDP tuning.

## 2. Perform-Compatible Mechanics Implemented

The hinge has generalized deformation and force vectors

\[
\mathbf d = [\delta,\theta]^T,\qquad \mathbf q=[N,M]^T,
\]

with compression-positive physical column force \(P=P_0-N\).

Plastic flow is associated / normal to the interaction surface, so flexural yielding can generate both plastic rotation and plastic axial deformation.

For the final Phase-9F experiment, the P–M yield surface follows Perform's documented concrete-type form in a P–M plane:

\[
\left|\frac{P-P_B}{P_{Y0}-P_B}\right|^{\alpha}
+
\left|\frac{M}{M_{YB}}\right|^{\beta}=1.
\]

Separate axial intercepts and \(\alpha\) values are used on the tension and compression branches. Perform's suggested \(\beta=1.1\) is used.

The surface parameters are fitted only to the transparent Whitney-block + elastic-perfectly-plastic-steel section mechanics already used in the Berkeley validation source. **No response-history EDP enters the fit.**

### A/B nonductile columns

- \(P_{YT}=-27.44\) kip
- \(P_B=45.0472\) kip
- \(M_{YB}=125.599\) kip-in
- \(P_{YC}=136.682\) kip
- \(\alpha_T=1.504\)
- \(\alpha_C=1.711\)
- \(\beta=1.1\)
- section-fit RMSE ≈ 1.78 kip-in tension branch / 2.69 kip-in compression branch

### C/D ductile columns

- \(P_{YT}=-56.32\) kip
- \(P_B=47.2249\) kip
- \(M_{YB}=177.486\) kip-in
- \(P_{YC}=165.562\) kip
- \(\alpha_T=1.524\)
- \(\alpha_C=1.386\)
- \(\beta=1.1\)
- section-fit RMSE ≈ 3.77 kip-in tension branch / 2.98 kip-in compression branch

## 3. Important Numerical Correction Found During Phase 9F

The PM hinge uses a high axial penalty stiffness to behave as a rigid-plastic zero-length hinge in series with the physical elastic column. Initially this numerical axial penalty was inadvertently present in the stiffness-proportional Rayleigh damping matrix. This produced nonphysical vertical damping forces and violated axial force parity between the PM hinge and elastic column.

Phase 9F therefore separates:

- structural / Newton initial stiffness, which includes the PM penalty, and
- physical Rayleigh damping-reference stiffness, which excludes that artificial axial penalty.

This is a modeling correction, not an EDP calibration.

## 4. Controlled Validation Result

The same geometry, mass, elastic calibration, P-Delta treatment, beams, panel zones, FSC springs, proxy excitation, solver controls, and damping definition were retained.

| Model | S1 peak drift | S2 peak drift | S3 peak drift | S1/S2 |
|---|---:|---:|---:|---:|
| QuakeCore research | 3.278% | 4.248% | 2.176% | 0.772 |
| Phase 9F Perform-type P–M, 0–15 s | **5.104%** | **5.700%** | **2.992%** | **0.895** |
| Published Perform3D | 6.07% | 3.92% | 1.91% | 1.548 |
| Measured experiment | 5.18% | 4.70% | 2.61% | 1.102 |

The Phase-9F peaks occur at approximately:

- Story 1: 12.365 s
- Story 2: 12.315 s
- Story 3: 12.305 s

At the Story-2 peak, the simultaneous drift shape is approximately:

\[
[4.853,\ 5.700,\ 2.956]\%.
\]

The Phase-9F Story-1 peak is only about **1.5% below the measured 5.18% drift**, while Story 2 is about **21% above** the measured value and Story 3 about **15% above**.

The result therefore strongly supports P–M state-path mechanics as a major missing ingredient, but does not by itself reproduce the Perform3D vertical distribution.

## 5. B1 Response

For the Phase-9F 15-second window:

- maximum B1 bottom total hinge rotation ≈ 0.01838 rad
- maximum B1 bottom plastic rotation ≈ 0.01828 rad
- maximum accumulated B1 plastic axial deformation ≈ **0.1696 in**
- at maximum plastic rotation: \(P\approx22.82\) kip compression and \(M\approx106.15\) kip-in

The magnitude and continuing growth of the plastic axial deformation are important. Perform's own P–M plasticity documentation warns that associated plasticity can overpredict cyclic axial growth/shortening in reinforced concrete columns. Phase 9F reproduces that qualitative concern directly.

## 6. Tangent-State Diagnostic

At the Phase-9F Story-2 peak, the condensed interstory tangent diagonal is approximately:

\[
[52.62,\ 39.16,\ 43.80]\ \text{kip/in}.
\]

Unlike the earlier straight-C→E scalar-hinge sensitivity, the P–M model does not create the same negative condensed Story-1/Story-2 tangents at this event. The response redistribution is therefore arising through coupled axial-flexural state evolution, not simply through a prematurely negative scalar flexural tangent.

## 7. Numerical Robustness Boundary

With the standard robust solver setting of six allowed subdivisions, the Perform-type P–M model completes through the validated 15-second window. Continuation toward 20 seconds becomes computationally expensive because the analysis enters repeated subdivision as axial-flexural plastic state evolves.

A diagnostic run limited to three subdivisions fails at about 3.59 s, while the six-subdivision formulation proceeds through 15 s. This demonstrates that the response is numerically demanding rather than proving a physical collapse at the early lower-subdivision boundary.

The full 70-second proxy record has **not** yet been completed with this formulation. Phase-9F results must therefore not be described as the final full-record envelope.

## 8. What Phase 9F Establishes

1. A true associated P–M hinge materially changes the Berkeley response.
2. The change is in the direction needed to recover Story-1 localization.
3. A smooth Perform-type concrete surface is much better behaved and more physically defensible than using the legacy section-capacity helper directly as a plasticity surface.
4. The response cannot be reproduced by scalar \(M_y(P_{max})\) adjustment alone.
5. Associated RC P–M flow generates substantial plastic axial deformation and introduces a significant state-path effect.
6. P–M interaction is therefore promoted from a hypothesis to a **first-order model ingredient** for the Berkeley benchmark and for production RC-column NLRH modeling.

## 9. Remaining Difference from Exact Perform3D

This is still labeled **Perform-type / Perform-compatible mechanics**, not an exact Perform3D reproduction, because:

- full Perform trilinear Mroz translating-surface hardening is not yet implemented;
- ASCE 41-17 proprietary component parameter tables have not been reproduced;
- the exact recorded DT1 shake-table motion has not been recovered;
- the damping implementation has intentionally remained frozen during this mechanics experiment.

## 10. Damping Finding and Recommended Next Experiment

The damping issue identified by the user is real and worth testing next, with one correction: modern OpenSees does support modal damping, and its documentation states that its modal damping implementation follows Perform3D. OpenSees also warns that Rayleigh damping can produce incorrect results in nonlinear concentrated-plasticity analysis.

The NIST benchmarking paper states that most models used either 3% Rayleigh damping or a combination of Rayleigh and modal damping, and specifically that the UC Berkeley model used a 3% damping ratio based on experimental measurements. The paper does not uniquely identify the exact modal/Rayleigh split for this Berkeley Perform3D model in the available public description.

Current QuakeCore Berkeley damping is Rayleigh with approximately:

- 3.00% at T1 = 0.48 s
- 2.36% at T2 ≈ 0.154 s
- 3.07% at T3 ≈ 0.093 s

and the PM numerical penalty stiffness is now excluded from the damping-reference matrix.

A dedicated next sensitivity should hold the Phase-9F P–M mechanics fixed and compare:

1. current physical-reference Rayleigh damping;
2. 3% modal damping in the first several elastic modes;
3. 3% modal damping plus a very small high-frequency Rayleigh term, consistent with the Powell/NIST modeling philosophy;
4. optionally an updated/instantaneous modal-damping research path, kept separate from the fixed-elastic-mode benchmark formulation.

This should be treated as a controlled **damping-model sensitivity**, not another calibration knob.

