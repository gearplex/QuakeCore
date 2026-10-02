# C2 deficiency-to-fragility pilot

**Status:** protocol frozen for the first pilot; executable nonlinear archetype not yet frozen.

This study tests whether ASCE 41 Tier 1 screening deficiencies can be mapped to explicit changes in structural mechanism and then to measurable changes in collapse fragility. The first pilot is a six-story C2-type reinforced-concrete frame-wall building patterned after the FEMA P-2018 Appendix D six-story, six-bay, single-wall studies.

The work has two simultaneous purposes:

1. answer the engineering question with OpenSees as the primary research engine; and
2. create a blinded, versioned benchmark that QuakeCore must reproduce before QuakeCore results are used as primary evidence.

This is a research protocol, not an ASCE 41 evaluation, FEMA P-695 qualification, or code-compliance claim.

## 1. Research question

At a common global strength-floor anchor, how much fragility reserve is consumed by:

- a concentrated first-story mechanism;
- nonductile shear/axial degradation in the wall-frame system; and
- their interaction?

The initial study deliberately excludes torsion. A later 3D phase will add torsional response and test the framework's `A_T` modifier only after the 2D mechanism effects are understood.

## 2. Four controlled cases

| Case | Vertical strength distribution | Deformation/degradation behavior | Purpose |
|---|---|---|---|
| `C2-RF` | Regular | Flexure-governed, distributed ductile mechanism | Reference |
| `C2-RB` | Regular | Shear/axial degradation enabled in designated nonductile components | Brittle-mechanism main effect |
| `C2-SF` | First-story weakness sufficient to concentrate drift | Flexure-governed | Story-mechanism main effect |
| `C2-SB` | Same first-story weakness | Same shear/axial degradation as `C2-RB` | Interaction effect |

The primary, mechanism-isolation comparison holds the global strength floor approximately constant across all four cases. The first-story variants redistribute strength rather than simply reducing every lateral-resisting component. A secondary, practice-realistic comparison may retain the unnormalized strength loss, but its result must be labeled as the combined effect of strength and fragility reserve.

## 3. Framework quantities

The deficiency-based framework uses a lognormal collapse fragility

\[
P(C\mid IM)=\Phi\left(\frac{\ln IM-\ln IM_{50}}{\beta}\right)
\]

and defines the fragility reserve ratio at probability `p` as

\[
FRR_p=\frac{IM_{50}}{IM_p}=\exp\left(|z_p|\,\beta_\mathrm{eff}\right).
\]

At `p = 0.01`, `|z_p| = 2.3263`, so the framework's reference `beta_eff = 0.60` corresponds to `FRR_1% = 4.04`.

For a deficiency demand-amplification factor `A`, the framework predicts

\[
\beta_{\mathrm{eff,new}}
=\beta_{\mathrm{eff,ref}}-\frac{\ln A}{|z_p|}.
\]

With a common low-probability strength anchor, the amplification inferred from simulation is

\[
A_\mathrm{obs}
=\frac{FRR_\mathrm{ref}}{FRR_\mathrm{case}}
\approx\frac{IM_{50,\mathrm{ref}}}{IM_{50,\mathrm{case}}}.
\]

### Do not conflate the two betas

- `beta_RTR` is the fitted record-to-record logarithmic dispersion from the IDA population.
- `beta_eff` is the framework's reserve-ratio index derived from the distance between the low-probability anchor and median collapse.

Raw `beta_RTR` must not be compared directly with the framework's `0.60`. If a total uncertainty is needed, modeling, material, construction-quality, and record-to-record uncertainties must be declared and combined under an explicit probabilistic model. Square-root-sum-of-squares combination is permitted only when independence is justified.

## 4. Archetype basis

The starting geometry is a documented reconstruction patterned after the FEMA P-2018 Appendix D archetypes:

- six stories and six bays;
- one planar reinforced-concrete wall coupled to a concrete frame/gravity system;
- reference wall length of 20 ft and thickness of 10 in;
- nonconforming wall reinforcement and component behaviors representative of older construction; and
- one reference configuration with shear-critical columns, corresponding to the mechanism illustrated in FEMA P-2018 Figure D-3.

These facts define the research family, not a claim of identity to FEMA's model. The exact member sizes, gravity loads, material properties, hinge/shear parameters, damping, and mass distribution must be entered in `archetype.json` with provenance before the OpenSees baseline is frozen.

## 5. Analysis variables

### Primary intensity measure

Use `Sa(T1_ref)` at 5% damping, where `T1_ref` is frozen from `C2-RF`. Holding the period definition fixed avoids moving the intensity-measure axis when a deficiency changes stiffness or period.

Report a spectral-shape sensitivity using an explicitly documented average spectral acceleration. Its period band and weighting must be frozen before reviewing the IDA results.

### Response quantities

Every run records at least:

- maximum and residual story drift by story;
- wall shear, moment, axial force, curvature, and shear deformation;
- column shear, axial force, chord rotation, and axial-degradation state;
- floor acceleration, base shear, and overturning moment;
- component IO/LS/CP or corresponding research limit-state crossings;
- story restoring-force degradation;
- energy terms available consistently in both solvers; and
- convergence, subdivision, and termination diagnostics.

### Limit states

| ID | Limit state | Use |
|---|---|---|
| `L0` | First significant inelasticity / strength-floor anchor | Calibration only; not collapse |
| `L1` | Declared component CP threshold exceeded | Damage/performance result; not automatically collapse |
| `L2` | Any story lateral resistance degrades to 50% of its peak | Primary FEMA P-2018-style global collapse proxy |
| `L3` | Dynamic instability or explicitly modeled gravity-load/axial-support loss | Physical-collapse classification |

Numerical noncompletion is a separate run status and is never silently converted to `L2` or `L3`.

## 6. IDA and statistical protocol

1. Use a traceable far-field record set with immutable raw-file hashes, units, components, processing, and scale factors.
2. Perform adaptive IDA with a declared initial scale grid and a logarithmic refinement tolerance.
3. Store the last noncollapse intensity and first collapse intensity for each record.
4. Treat records that do not collapse at the maximum analyzed intensity as right-censored.
5. Treat collapse below the first analyzed intensity as left-censored or rerun at lower intensity.
6. Fit the lognormal fragility with the interval/right/left-censored likelihood in `fragility.py`; do not replace brackets with arbitrary midpoint collapse intensities.
7. Bootstrap by record for confidence intervals on `IM50`, `beta_RTR`, and pairwise median-capacity ratios.
8. Report unresolved solver failures separately and perform bounded sensitivity analyses rather than assigning them a preferred physical outcome.

The four cases use the same record components and scaling rules. Pairwise record-level comparisons are retained even when the suite-level fragilities are fit independently.

## 7. OpenSees–QuakeCore decision gate

OpenSees is the primary engine for the first scientific result. QuakeCore is promoted from shadow solver only after it meets all of the following on the frozen archetype:

1. modal periods, gravity equilibrium, and pushover curves satisfy predeclared tolerances;
2. component cyclic histories and degradation events match the intended constitutive definitions;
3. at least four pilot records span elastic, nonlinear, and near-collapse response;
4. response histories and peak EDPs meet frozen cross-solver tolerances over the independently converged common window;
5. collapse brackets and classified mechanisms agree, or every disagreement is adjudicated without changing thresholds after the fact; and
6. the unresolved Phase 9N record `120521` issue is either resolved or demonstrated not to affect the C2 model path.

If QuakeCore passes, it may be used for throughput while OpenSees remains the independent reference on a stratified audit sample. If it does not pass, all scientific conclusions remain OpenSees-based and the mismatch becomes a software-validation work item.

## 8. Staged execution

### Gate A — frozen definitions

- Complete licensed ASCE 41 Tier 1 checklist mapping without reproducing controlled text.
- Freeze the archetype parameter table and provenance.
- Freeze record set, intensity measure, scaling, damping, and limit states.
- Register expected model/input/output hashes.

### Gate B — monotonic and cyclic mechanics

- Gravity and modal checks.
- Pushover and story-mechanism classification.
- Component cyclic tests for flexure, shear degradation, and axial loss.
- Demonstrate that the four cases change only their declared variables.

### Gate C — pilot NRHA/IDA

- Four records selected before viewing variant results.
- Elastic, nonlinear, and near-collapse scales.
- OpenSees first; QuakeCore comparison blind to the final fragility fit.

### Gate D — production suite

- Complete record set and censored fragility fits.
- Pairwise main-effect and interaction estimates.
- Compare observed reserve loss with the deficiency-based predictions.
- Predeclare whether the evidence warrants a later 3D torsion study.

## 9. Immediate next implementation task

Build and freeze the OpenSees `C2-RF` baseline. The minimum acceptable baseline package contains:

- a machine-readable model generator input;
- an OpenSeesPy model generator and analysis driver;
- gravity/modal/pushover/cyclic verification outputs;
- record manifest schema and one non-proprietary smoke-test motion;
- deterministic hashes for inputs and key outputs; and
- a QuakeCore translation manifest that lists supported, approximated, and unsupported mechanics.

No production IDA should begin until Gate A and Gate B are signed off in the repository.

## 10. Reproducible analysis environment

The current independent reference version is OpenSeesPy 3.8.0. Install the analysis-only dependencies separately from the QuakeCore C++ build:

```bash
python3 -m venv --system-site-packages .venv-c2
.venv-c2/bin/python -m pip install -r research/c2_deficiency_fragility/requirements-analysis.txt
```

In restricted Linux containers, the bundled MPICH initialization can require an explicit loopback hostname:

```bash
MPICH_INTERFACE_HOSTNAME=127.0.0.1 .venv-c2/bin/python path/to/driver.py
```

That variable is an execution-environment workaround, not a modeling parameter. Record it in the run environment but do not embed it in analysis results or model hashes.
