# QuakeCore Phase 9D.5 — Immutable Compiled Topology + Objective Corotational RC Demand Recovery

**Status:** architecture prototype / regression-clean increment  
**Basis:** Phase 9D.4 native RC-column binding and rollback-safe committed demand recording  
**Purpose:** stop rebuilding invariant frame topology during ASCE/ACI parameter iterations and enable native RC demand recording on the finite-rotation corotational 3D reference path.

## 1. Architecture change

Phase 9D.4 still used a model factory on every ASCE iteration. Phase 9D.5 adds a production path in which a frame is compiled once and only its nonlinear material field changes between analyses.

The iteration is now:

`compile geometry / MPC / mass / sparse topology / nonlinear basis once`

then, for iteration `k`:

`resolved RC code parameters -> replace hinge material field -> reset analysis state -> NLRH -> committed P/V envelopes -> synchronous rules regeneration`

The invariant data are not reconstructed between passes:

- reduced DOF map,
- nodes and elastic elements,
- MPC/diaphragm transformation,
- mass matrix,
- elastic stiffness topology,
- sparse effective-system pattern,
- nonlinear update-basis directions,
- element/component ID maps,
- nonlinear state offsets/layout.

## 2. State-layout-preserving material field

`CompiledFrame2D` and `CompiledFrame3D` now provide:

- `replace_nonlinear_material(...)`
- `replace_nonlinear_material_field(...)`

A replacement is accepted only if the new material has the same state size as the existing material. This is deliberate: the accepted-state vector, rollback layout, component offsets and compiled nonlinear topology remain immutable.

For an ASCE 41 parameter iteration, the topology should therefore be seeded with ASCE 41 hinge materials and subsequent iterations replace those ASCE 41 parameter sets in-place.

When `Ke` changes, QuakeCore updates the affected `b b^T` contribution to `K_initial` in-place and updates the component initial tangent. It does not rebuild the sparse matrix pattern or nonlinear basis.

The 3D implementation stores a second nonlinear scatter map indexed directly into `K_initial`, separate from the existing union-system scatter map used during state tangent assembly.

## 3. Compound-model safety boundary

Scalar material-field replacement is intentionally rejected for:

- flexure-shear-critical compound springs, and
- the experimental moving/coupled P-M hinge wrapper.

Those models have additional constitutive state or interaction data beyond a scalar hinge law. Allowing only the base material to change would leave stale coupled interaction data and would be physically misleading.

A future moving-surface P-M/P-M-M implementation should expose a dedicated interaction-field replacement API rather than bypass this guard.

## 4. Compiled whole-building ASCE adapters

New persistent-topology containers:

- `RCNativeCompiledFrame2D`
- `RCNativeCompiledFrame3D`

New entry points:

- `run_compiled_frame2d_two_pass_asce41_23_aci369_1_22`
- `run_compiled_frame2d_asce41_23_aci369_1_22`
- `run_compiled_frame3d_two_pass_asce41_23_aci369_1_22`
- `run_compiled_frame3d_asce41_23_aci369_1_22`

Before each global analysis QuakeCore maps the synchronous `RCBuildingColumnModel` field to the hinge IDs in each `RCColumnElementBinding`, replaces all hinge materials, and then runs the same compiled frame.

The supplied compiled frame is left carrying the final material field after the iteration, which makes the final model directly inspectable.

The Phase 9D.4 factory-based APIs remain available for backward compatibility and for cases where the structural topology itself genuinely changes between iterations.

## 5. Exact replacement equivalence gate

A 3D regression constructs the same frame in two ways:

1. compile with ASCE hinge field A, then replace it in-place with field B;
2. freshly compile the identical geometry directly with field B.

The two models are required to have:

- identical `K_initial` CSC pattern,
- `K_initial` values equal to numerical roundoff,
- identical nonlinear force response,
- identical nonlinear tangent response.

The test passes. This is the principal correctness gate for the immutable-topology architecture.

## 6. Objective corotational section-force recovery

Phase 9D.4 intentionally refused native force recovery when `CompiledFrame3D` used the finite-rotation corotational reference element. Phase 9D.5 removes that limitation.

A new `corotational3d_section_response(...)` recovers local physical member forces from the objective co-rotated basic deformation system:

- compression-positive axial force,
- local y/z end shear,
- torsion,
- local y/z end moments,
- current chord length.

For basic axial extension `q0`, the physical compression is

`P = P0 - EA/L0 * q0`.

Bending end moments are the analytic derivatives of the same co-rotated basic strain-energy function used by the reference element. With no distributed element load, local end shear follows from end-moment equilibrium on the **current** chord length.

This avoids projecting forces on the initial axes and preserves objectivity under large rigid-body motion.

## 7. Corotational objectivity regression

The section-force recorder is tested under:

- a large rigid-body rotation with prescribed compression preload,
- the same axial extension with and without that rigid-body rotation,
- the compiled `Frame3D` corotational element-ID response path.

The rigid rotation produces no spurious bending/shear and retains the preload. The rotated and unrotated axial-extension cases recover the same physical axial force within numerical tolerance.

## 8. Whole-building compiled-topology regression

Both 2D and 3D native two-pass tests now exercise the persistent compiled topology path.

The regression verifies that:

- the same frame instance is used for both global analyses,
- the sparse topology/nonlinear basis is unchanged,
- committed member-demand recording remains active,
- the second pass receives the demand-sensitive regenerated hinge,
- the compiled frame retains the final generated material field.

## 9. Verification gates

Passed after Phase 9D.5 changes:

- `quake_tests`: **All quake-core tests passed**,
- CTest: **100% tests passed, 0 failed**,
- `ucb_validation_phase8`: builds successfully.

The supplied Phase-9A-derived Berkeley source snapshot still has its pre-existing release-pushover limitation. Running the Phase-8 release mode terminates at:

- peak base shear: **23.0929 kip**,
- roof drift at peak: **1.625%**,
- SuperLU failure: **info=83**.

Phase 9D.5 does not change or resolve that source-snapshot issue.

## 10. What Phase 9D.5 does not claim

- It does not populate proprietary ASCE 41-23 / ACI 369.1-22 coefficient tables.
- It does not convert the experimental coupled P-M wrapper into a production moving-surface plasticity formulation.
- It does not yet reuse a persistent factorization across separate complete NLRH analyses; each analysis still constructs its solver/prepared solve state from the unchanged compiled matrices.
- It does not resolve the exact Berkeley DT1 recorded-input blocker.

## 11. Recommended next gate

The next useful step is **Phase 9D.6: prepared-analysis reuse + code-model comparison harness**.

Two related tasks should be combined:

1. cache/reuse analysis preparation that remains invariant when only the hinge material field changes, where mathematically valid;
2. add a single comparison harness that runs the Berkeley frame through:
   - legacy QuakeCore research hinge field,
   - NIST / ASCE 41-17 benchmark field,
   - production ASCE 41-23 / ACI 369-resolved field,
   with identical geometry, excitation, mass, damping definition, recorder definitions and solver controls.

This would turn the architectural work of Phases 9D.1-9D.5 into a direct component/system-level modeling comparison rather than another EDP-fitting exercise.
