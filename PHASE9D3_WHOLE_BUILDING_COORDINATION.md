# QuakeCore Phase 9D.3 — Whole-Building ASCE 41 RC Column Coordination

**Status:** architecture prototype / regression-clean increment  
**Basis:** Phase 9D.2 section→rules→hinge→analysis demand iteration  
**Purpose:** coordinate demand-sensitive RC-column code models across a complete building analysis without update-order dependence.

## 1. What Phase 9D.3 adds

Phase 9D.2 could regenerate one RC-column hinge from an observed axial/shear envelope. Phase 9D.3 adds a building-level coordinator that operates on the full RC-column parameter field at once:

`all section/detailing inputs -> all code-resolved hinges -> one global NLRH/NDP -> all P/V envelopes -> simultaneous regeneration -> repeat`

New files:

- `include/quake/rc_building_asce41.hpp`
- `src/rc_building_asce41.cpp`

The library target includes the new source automatically.

## 2. Deterministic synchronous update rule

The coordinator deliberately uses a Jacobi-style iteration. At iteration `k`:

1. Canonicalize the coordinated columns by unique `component_id`.
2. Resolve **every** column model from the demand field committed at the beginning of iteration `k`.
3. Pass the complete model vector to one building response callback.
4. Require exactly one P/V demand envelope for every coordinated column, matched by `component_id` rather than array order.
5. Resolve candidate parameters for **all** columns from the newly observed demand field.
6. Compute component and global convergence metrics.
7. Only after all candidate updates exist, commit the next demand field simultaneously (with optional under-relaxation).

No column can therefore see another column's `k+1` demand while it is still being generated for iteration `k`. This avoids component-order-dependent models.

## 3. Public API

### Column definition

`RCBuildingColumnDefinition`

- physical `RCColumnSectionInput`
- selected `ASCE41BackboneShape`

### Global analysis input

The building analysis callback receives:

`std::vector<RCBuildingColumnModel>`

Each entry carries the canonical `component_id` and the fully resolved `RCColumnModelSpec` for that iteration.

### Global analysis output

`RCBuildingAnalysisObservation`

- `analysis_succeeded`
- free-form analysis `status`
- `RCBuildingColumnDemandObservation` for every coordinated column

Each demand retains:

- maximum compression
- maximum tension
- maximum absolute shear

The coordinator rejects duplicate, missing, unknown, or negative demand records rather than silently accepting an incomplete building envelope.

## 4. Convergence and failure semantics

The convergence norm is the maximum over **all columns** of:

- relative demand-field change, and
- relative resolved-parameter change.

A building therefore cannot be marked converged because most columns are stable while one controlling column continues to change.

Termination is explicit:

- `Converged`
- `IterationLimit`
- `AnalysisFailure`

A failed/nonconverged NRHA is not mislabeled as failure of the code-parameter fixed-point iteration.

## 5. Two-pass workflow

`run_two_pass_asce41_23_aci369_1_22(...)` performs exactly two global analyses:

1. initial gravity/preanalysis demand field -> all hinges -> building analysis;
2. pass-one P/V envelopes -> simultaneous regeneration of all hinges -> second building analysis.

The general iterative workflow can continue beyond two passes until both demand and parameter fields stabilize.

## 6. Regression gates added

Tests now verify:

1. whole-building synchronous regeneration;
2. input-order invariance;
3. component-ID matching even when observations are returned in another order;
4. global convergence controlled by the worst-changing column;
5. whole-building two-pass regeneration;
6. explicit building-analysis failure classification.

Regression status:

- `quake_tests`: PASS (`All quake-core tests passed.`)
- `ctest`: PASS, 1/1
- `ucb_validation_phase8` build: PASS

## 7. Important scope boundary

Phase 9D.3 is now a real **whole-building coordination engine**, but the application analysis driver still owns the mapping from QuakeCore frame elements to column demand recorders and returns the P/V envelopes through `RCBuildingResponseRunner`.

This is intentional separation of concerns, but it is the remaining integration seam before a user can simply mark RC columns in a Frame2D/Frame3D model and ask QuakeCore to perform the entire ASCE iteration automatically.

No ASCE 41-23 / ACI 369.1-22 copyrighted table coefficients are embedded or guessed in this increment. The existing authorized-rules boundary remains intact.

## 8. Recommended Phase 9D.4

The next increment should integrate the coordinator with the structural model itself:

1. add an RC-column component registry/binding from physical column IDs to Frame2D/Frame3D elements and end hinges;
2. add native axial/shear envelope recorders that survive substepping and rollback correctly;
3. expose a building-analysis adapter that automatically harvests these envelopes after each NLRH/NDP;
4. connect regenerated `RCColumnModelSpec` objects back to the corresponding hinge materials simultaneously;
5. exercise the complete path on the Berkeley 3-story frame with a benchmark resolver first;
6. then compare legacy-research, NIST/ASCE41-17 benchmark, and authorized ASCE41-23/ACI369 production parameter sets without EDP tuning.

That will turn the Phase 9D.3 orchestration contract into a one-call building workflow.
