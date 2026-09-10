# Phase 9D.2 — RC section rules boundary and axial-demand iteration

## Status

Architecture increment / prototype. This is not a release claim and does not
embed ASCE 41-23 or ACI 369.1-22 table coefficients.

## Objective

Phase 9D.2 turns the Phase 9D provider seam into a usable analysis workflow:

`physical RC column -> authorized rules resolver -> auditable hinge parameters -> NLRH/NDP -> observed P/V demand -> regenerated hinge`

The nonlinear solver remains independent of the code-rule implementation.

## New physical input object

`RCColumnSectionInput` carries factual member information rather than code-table
outputs, including:

- component identifier;
- section width/depth and clear length;
- expected concrete and reinforcing-steel strengths;
- longitudinal reinforcement;
- transverse reinforcement geometry;
- gravity axial compression and initial shear demand;
- expected shear capacity;
- factual detailing/condition flags.

The effect of those inputs on modeling parameters is intentionally not coded in
the solver. That responsibility belongs to an edition-specific rules resolver.

## Authorized rules boundary

`RCColumnRulesResolver` is an abstract interface:

```text
resolve(section, demand) -> RCColumnRuleResolution
```

The resolution contains:

- `RCColumnResolvedParameters` (`My`, `Mc`, `a`, `b`, `c`, IO/LS/CP, etc.);
- the exact demand state used to resolve them;
- source/provenance text;
- per-parameter audit entries containing source reference and controlling condition.

`CallbackRCColumnRulesResolver` allows a licensed/firm rules engine to be
connected without linking the QuakeCore solver to that implementation.

## Axial-demand iteration

`RCColumnAxialIteration` implements a solver-agnostic fixed-point workflow:

1. Initialize with gravity compression and initial shear demand.
2. Resolve ASCE/ACI modeling parameters at that demand state.
3. Construct the QuakeCore hinge model.
4. Run NLRH/NDP through an application-supplied response callback.
5. Record maximum compression, tension, and absolute shear.
6. Re-resolve the hinge at the observed demand.
7. Repeat until both demand and generated parameters stabilize or the iteration cap is reached.

Convergence is checked independently on:

- relative demand-envelope change; and
- maximum relative change among generated backbone/acceptance parameters.

Optional under-relaxation is provided for oscillatory fixed-point behavior.

## Two-pass convenience workflow

`run_two_pass_asce41_23_aci369_1_22()` implements the practical two-analysis
workflow directly:

- Pass 1: generate from gravity/preanalysis demand and run response history.
- Pass 2: regenerate at the maximum compression observed in Pass 1 and rerun.

The returned object still reports whether Pass 2 is self-consistent with its own
observed demand, but the convenience function always stops after two analyses.

## Provenance / audit philosophy

A generated hinge should be reviewable as engineering work product. A model
record can now retain, for example:

```text
parameter: pos_a
value: 0.012 rad
source_reference: <authorized rule/table reference>
controlling_condition: <axial ratio / shear condition / detailing case>
```

QuakeCore does not need to know the copyrighted lookup coefficient that produced
the value in order to preserve an auditable record of the result.

## Tests added

The regression suite now checks:

- section-to-rules callback boundary;
- explicit provenance retention;
- per-parameter audit retention;
- axial-sensitive parameter regeneration;
- convergence after regeneration/reanalysis;
- two-pass convenience behavior;
- under-relaxation behavior;
- continued absence of silently embedded code coefficients.

`quake_tests` result after the change:

```text
All quake-core tests passed.
```

The `ucb_validation_phase8` executable also builds successfully against the new
API. A no-argument execution uses its internal development proxy and is not an
external-validation gate; Phase 9D.2 does not change that status.

## Important limitations

1. Actual ASCE 41-23 / ACI 369.1-22 coefficients are not included.
2. The rules resolver must be supplied by an authorized implementation.
3. The iteration manager currently operates on one component/model callback at a
   time. A whole-building coordinator that regenerates all axial-sensitive
   column hinges together is the next integration step.
4. This is still an iterative scalar-demand treatment, not a moving-surface P-M
   or P-M-M return-mapping formulation.
5. Shear/FSC coupling remains separate and should ultimately consume the same
   evolving demand envelope rather than being silently folded into flexure.

## Recommended next gate

Phase 9D.3 should add a model-level `RCColumnCodeModelManager` that:

- registers every code-generated RC column hinge;
- runs one complete building NLRH/NDP per iteration;
- extracts each column's P/V envelope from the same analysis;
- regenerates all affected hinges simultaneously;
- records per-column convergence and parameter provenance;
- stops when the building-wide parameter set stabilizes;
- preserves a frozen benchmark mode for NIST ASCE 41-17 comparisons.

After that manager exists, QuakeCore can perform the Berkeley three-way run
(legacy research / NIST benchmark / production resolved-code mode) without
manual per-column orchestration.
