# QuakeCore Phase 9D.4 — Native RC Column Binding + Committed Demand Recording

**Status:** architecture prototype / regression-clean increment  
**Basis:** Phase 9D.3 whole-building synchronous ASCE 41 coordination  
**Purpose:** bind physical RC columns directly to compiled frame elements/hinges and harvest axial/shear demand envelopes from committed nonlinear analysis states without application-specific probe code.

## 1. Why this phase exists

Phase 9D.3 could regenerate all RC-column code models synchronously, but the global analysis callback still had to supply each column's `P/V` envelope. Phase 9D.4 moves that responsibility into QuakeCore.

The architecture is now:

`RC section/detailing -> rules resolver -> code hinge field -> native frame factory -> compiled element/hinge bindings -> robust NLRH -> committed P/V envelopes -> synchronous regeneration`

The application still owns building geometry/model construction because QuakeCore does not yet have a universal declarative building schema. It no longer owns response probing, rollback filtering, substep bookkeeping, or demand-envelope assembly.

## 2. Native physical-column binding

New `RCColumnElementBinding` records map a physical column `component_id` to:

- one or more elastic frame element IDs,
- optional nonlinear hinge/component IDs.

Bindings are validated before analysis:

- component IDs must be unique,
- at least one elastic element is required,
- every element ID must exist,
- an elastic element cannot belong to two physical-column bindings,
- supplied hinge IDs must exist,
- the binding set must exactly cover the RC column model field supplied by the Phase 9D.3 coordinator.

Multiple elastic elements are permitted for a segmented physical column; demand envelopes are taken over all bound segments.

## 3. Native committed member-force recovery

`CompiledFrame2D` and `CompiledFrame3D` now expose element-ID-based member response recovery.

### 2D

`ElasticFrame2DResponse` reports:

- compression-positive axial force,
- local end shears,
- local end moments.

### 3D

`ElasticFrame3DResponse` reports:

- compression-positive axial force,
- local y/z shears at both ends,
- torsion,
- local y/z end moments.

For RC demand recording, 3D shear demand is the maximum local vector shear magnitude:

`V = hypot(Vy, Vz)`.

For the small-displacement / updated-P-Delta formulation, axial compression is recovered as

`P = P_preload - EA/L * delta_axial`.

This is the same axial-force convention previously hand-coded in the UC Berkeley validation driver.

### Current limitation

Native committed member-force recovery intentionally throws for the finite-rotation corotational 3D reference path. An objective current-configuration section-force recovery should be implemented before native ASCE demand iteration is enabled for corotational models. Silent use of an initial-axis approximation was deliberately avoided.

## 4. Rollback-safe accepted-substep observer

The robust Newmark driver now provides `accepted_substep_state_observer`.

This differs materially from the previous output-step observers. During adaptive subdivision, legitimate force peaks can occur at internal integration substeps between original ground-motion samples. Those states should contribute to code demand envelopes.

However, an internal half-step can later be abandoned if the remainder of its parent output step cannot be completed. Such a state must not contaminate an envelope.

The new implementation therefore:

1. solves internal substeps normally,
2. buffers successful substep states locally,
3. exposes nothing while the parent original output step remains unresolved,
4. publishes all buffered committed substeps only after the complete original step succeeds,
5. discards the buffer if the parent step fails/rolls back.

Thus rejected Newton trials and rolled-back subdivision branches are never visible to the native demand recorder.

## 5. Native RC demand recorders

New classes:

- `RCFrame2DColumnDemandRecorder`
- `RCFrame3DColumnDemandRecorder`

For every bound physical column they retain:

- `max_compression_kip`,
- `max_tension_kip`,
- `max_abs_shear_kip`.

The zero-displacement preload state is sampled at construction so gravity compression is retained even if an analysis terminates before its first accepted dynamic step.

## 6. Native whole-building ASCE adapters

New `RCBuildingASCE41NativeFrame` entry points:

- `run_frame2d_two_pass_asce41_23_aci369_1_22`
- `run_frame2d_asce41_23_aci369_1_22`
- `run_frame3d_two_pass_asce41_23_aci369_1_22`
- `run_frame3d_asce41_23_aci369_1_22`

The model factory receives the complete synchronous `RCBuildingColumnModel` field and rebuilds the structural model with those hinges. QuakeCore then runs Robust Newmark and constructs the building demand observation natively.

No application-level response callback is required.

## 7. 2D material-bank upgrade

`CompiledFrame2D` previously hard-coded `BilinearSpring` state layout. That would have made a nominal native ASCE adapter misleading because the generated ASCE hinge could not actually be inserted into the 2D frame.

Phase 9D.4 upgrades Frame2D to the material-neutral `NonlinearMaterial` bank already used by Frame3D.

It now supports:

- legacy bilinear rotational springs,
- ASCE 41 hinges,
- IMK peak-oriented hinges,
- variable material state sizes,
- material diagnostics for IO/LS/CP, E/F, reversal and deterioration tracking.

Existing bilinear regression behavior remains intact.

## 8. Collapse versus numerical failure

The building-level iteration termination enum now distinguishes:

- `Converged`,
- `IterationLimit`,
- `AnalysisFailure`,
- `PhysicalCollapse`.

A physical-collapse NLRH is a valid structural result, but it does not supply a complete-record demand envelope for a code-parameter fixed-point iteration. Native adapters therefore stop the regeneration loop and preserve the physical-collapse classification rather than relabeling it numerical nonconvergence.

## 9. Berkeley validation integration

The Phase-8 Berkeley source's manual axial probe

`P = Pg - k_axial * (u_top - u_bottom)`

has been replaced by

`model.elastic_element_response(element_id, u).axial_compression`.

This removes one application-specific demand probe and exercises the native recovery path in the principal validation model source.

The uploaded Phase-9A-derived source snapshot has a pre-existing validation-run limitation: the Phase-8.1 release pushover terminates with SuperLU `info=83` after a 23.0929-kip peak at 1.625% drift. The untouched Phase 9D.3 package reproduces the same stop and values. Phase 9D.4 therefore preserves that source-snapshot behavior; this is not a new 9D.4 regression.

## 10. Regression gates

Passed:

- `quake_tests`: all tests pass,
- CTest: 100% pass,
- `ucb_validation_phase8`: builds successfully,
- committed-substep callback exercises actual adaptive subdivision,
- rolled-back failing output steps are not emitted,
- 2D and 3D element axial/shear recovery tests,
- native 2D ASCE material-bank two-pass test,
- native 3D two-pass demand/regeneration test,
- physical-collapse versus analysis-failure semantics.

## 11. Next recommended architecture gate

Phase 9D.5 should focus on **declarative hinge replacement / model cloning and objective corotational force recovery**.

The current native factory is a much narrower seam than the old arbitrary response callback, but the application still rebuilds geometry each ASCE iteration. A compiled-model cloning/rebinding API could preserve invariant geometry, mass, MPC topology and elastic matrices while replacing only the code-dependent material bank.

That would make the desired production call closer to:

`run_asce41_building_iteration(base_model, rc_bindings, motion, rules)`

with no custom model factory and with only the ASCE-dependent hinge field changing between passes.
