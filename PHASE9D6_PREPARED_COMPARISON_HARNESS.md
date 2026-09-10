# QuakeCore Phase 9D.6 — Prepared Analysis Reuse and Unified Code-Model Comparison Harness

**Status:** architecture / diagnostic increment  
**Purpose:** compare constitutive model alternatives on one compiled structural topology without conflating model-definition changes with geometry, solver, excitation, or recorder changes.

## 1. Objective

Phase 9D.6 turns the Phase 9D.1–9D.5 ASCE 41 architecture into a controlled comparison environment. A building is compiled once and then analyzed repeatedly with named nonlinear-material-field alternatives, for example:

1. existing QuakeCore research hinge field,
2. NIST / ASCE 41-17 benchmark hinge field, and
3. ASCE 41-23 / ACI 369.1-22 production hinge field when authorized resolved parameters are supplied.

The comparison harness intentionally does **not** embed or infer proprietary ASCE 41-23 / ACI 369.1-22 table coefficients.

The intended experimental control is:

```text
same compiled topology
same mass and damping
same geometric nonlinearity
same excitation
same timestep / Newton / subdivision controls
same recorders
              |
              +-- research material field
              +-- NIST / ASCE 41-17 material field
              +-- production ASCE 41-23 / ACI 369 resolved field
```

Only the explicitly replaced nonlinear material field differs between variants.

## 2. Prepared robust analysis reuse

`PreparedRobustNewmark` was added as a reusable robust NRHA driver.

```cpp
PreparedRobustNewmark prepared(frame, dt, strategy, max_subdivisions);
auto r1 = prepared.run(record, options);
auto r2 = prepared.run(record, options);
```

The prepared driver can retain linear solver preparation across repeated analyses. Reuse is permitted only when the exact effective-initial matrix previously prepared at each visited subdivision depth remains unchanged.

For Newmark depth `d`, the safety invariant is evaluated against the actual matrix

\[
K_{\mathrm{eff},0}^{(d)}
 = K_0 + a_0^{(d)}M + a_1^{(d)}C.
\]

Therefore:

- changing strength, plastic rotations, acceptance limits, or deterioration while retaining the same initial tangent can reuse preparation;
- changing hinge `Ke`, or anything else that changes `K_initial`, invalidates and rebuilds the affected preparation;
- comparison correctness never depends on cache reuse.

Subdivision-depth preparations are **lazy**. QuakeCore does not factor every potential half-step depth pre-emptively; a depth is prepared only if the robust integrator actually visits it.

The driver exposes diagnostics:

- `preparation_count()`
- `setup_factorizations()`
- `last_run_reused_preparation()`

These support performance auditing without changing solver semantics.

## 3. Variant isolation and execution-order invariance

A comparison variant is a named material-field delta:

```cpp
struct ModelComparisonVariant {
    std::string name;
    std::string provenance;
    std::vector<std::pair<int, NonlinearMaterial>> replacements;
};
```

Before the first run, the harness captures the baseline material for every component touched by any variant. Before **every** variant:

1. all touched components are restored to the captured baseline;
2. that variant's replacements are applied;
3. the building is analyzed from a fresh dynamic/constitutive state.

This prevents hidden carry-over from the variant executed previously.

Regression testing runs the same named variants in reversed order and verifies that the named physical responses are unchanged.

## 4. RC-column model-field adapter

`make_rc_column_model_comparison_variant(...)` converts a Phase 9D whole-building `RCBuildingColumnModel` field into the hinge replacements required by the comparison harness using the native `RCColumnElementBinding` map.

Coverage is exact:

- duplicate RC model IDs are rejected;
- a binding with no corresponding model is rejected;
- an RC model with no physical binding is rejected;
- a hinge ID assigned to more than one physical column is rejected;
- a variant with no hinge replacements is rejected.

Unbound non-column components, such as beams or joints, deliberately retain their captured baseline materials unless a separate variant explicitly replaces them.

## 5. Aligned story response histories

For every output step, a variant can retain:

- floor response coordinates,
- interstory drift ratios,
- roof response,
- peak and residual story drifts,
- peak and residual roof response.

Because all variants share the same excitation and nominal output `dt`, the stored histories are directly time-aligned.

Peak envelopes can additionally observe rollback-safe accepted subdivision states. This means a legitimate peak occurring between two record samples is included, while a Newton branch that is subsequently rejected/rolled back cannot contaminate the envelope.

## 6. Component-level M–theta–tangent histories

Selected nonlinear components can be monitored by stable component ID. At every committed output state the recorder stores:

- generalized deformation,
- generalized force,
- current tangent,
- material diagnostics.

For a rotational hinge this is the native

\[
M(t)\;\text{vs.}\;\theta(t)
\]

history requested for the Berkeley/NIST comparison.

Rollback-safe accepted substeps also update exact component envelopes and first-event times:

- maximum absolute deformation,
- maximum absolute force,
- minimum tangent,
- first IO,
- first LS,
- first CP,
- first lateral-resistance-loss event,
- first effective/gravity-failure event.

This retains the QuakeCore semantic rule that CP exceedance is a component acceptance event, **not** structural collapse.

## 7. Native RC-column P/V histories

The same Phase 9D.4 `RCColumnElementBinding` used for ASCE demand regeneration can be supplied to the comparison recorder.

For each bound physical column, every committed output state records:

- maximum compression across bound member segments,
- maximum tension across bound member segments,
- maximum absolute shear.

For 3D members the shear measure is the resultant

\[
V = \sqrt{V_y^2 + V_z^2}.
\]

Accepted subdivision states participate in the demand envelope under the same rollback-safe policy.

This puts story drift, hinge response, and column axial/shear history on one common committed-state timeline.

## 8. Instantaneous tangent modes at peak response states

Phase 9D.6 adds `modal_analysis_from_stiffness(...)`, which accepts an explicitly supplied stiffness matrix instead of always using the virgin `K_initial`.

While the record is running, the comparison recorder retains the exact committed displacement and nonlinear-state vectors that establish each story's peak absolute drift. After the run it reassembles the structural tangent at those states and computes instantaneous tangent modes.

For each requested snapshot the result contains:

- time / output step,
- trigger description,
- eigenvalue and period,
- reduced DOF mode shape,
- mode shape projected to configured floor/story response coordinates,
- roof-normalized story shape when normalization is well-defined.

Multiple story peaks occurring at the same committed state are merged rather than redundantly eigensolved.

Additional output-step indices may also be requested explicitly.

## 9. Condensed story tangent stiffness

At the same committed tangent snapshots, Phase 9D.6 statically condenses the full tangent onto the configured floor-response coordinates.

If `c_i` is the generalized response vector for floor coordinate `i`, the floor flexibility is constructed from

\[
F_{ij}=c_i^T K_t^{-1}c_j.
\]

The floor tangent is

\[
K_{floor}=F^{-1}.
\]

With cumulative floor displacement related to interstory deformation by

\[
q_{floor}=T\,\delta_{story},
\]

the interstory-coordinate tangent is

\[
K_{story}=T^T K_{floor}T.
\]

The diagnostic stores:

- full row-major `floor_tangent_matrix`,
- full row-major `interstory_tangent_matrix`,
- the interstory tangent diagonal.

If the current tangent is singular or cannot be condensed near a collapse state, the snapshot is retained with `story_tangent_available=false` instead of converting a diagnostic failure into an analysis failure.

This allows a Berkeley comparison to answer both:

1. **How did the instantaneous mode shape redistribute?**
2. **Which story tangent softened/stiffened when that redistribution occurred?**

## 10. 2D/3D parity additions

For the comparison layer, `CompiledFrame2D` now exposes the same story-response helpers needed by the 3D diagnostics:

- story response values,
- story elevations,
- maximum drift ratio,
- stable nonlinear component IDs.

Both 2D and 3D compiled frames expose read-only nonlinear component snapshots from arbitrary committed state vectors.

Compound/coupled components that cannot be represented by one scalar force/deformation channel remain intentionally excluded from the generic scalar snapshot API and require dedicated recorders.

## 11. Verification gates

Phase 9D.6 regression coverage includes:

- prepared robust-run reuse for an unchanged model;
- safe reuse for strength-only material changes with unchanged initial tangent;
- mandatory preparation invalidation after `Ke` changes;
- identical-field variant parity;
- time-aligned story histories;
- scalar component history recording;
- native column P/V history recording;
- tangent modes at story peaks;
- condensed story tangent availability;
- comparison execution-order invariance;
- exact RC model/binding coverage;
- 3D comparison smoke coverage.

At packaging time:

```text
All quake-core tests passed.
100% CTest tests passed.
ucb_validation_phase8 builds successfully.
```

## 12. Important limits / provenance

### 12.1 No production code-table coefficients are inferred

The production path remains a resolver boundary. Phase 9D.6 can compare a resolved ASCE 41-23 / ACI 369.1-22 field, but it does not populate proprietary coefficients from memory or secondary sources.

### 12.2 The three-way Berkeley comparison is enabled, not yet claimed complete

The harness now provides the required common execution/diagnostic path, but this increment does **not** claim final Berkeley results for all three model definitions. A true production comparison still requires authorized resolved ASCE 41-23 / ACI 369.1-22 column parameters, and the current source lineage remains based on the earlier Phase-9A source snapshot with later Phase-9B/9C validation artifacts maintained separately.

### 12.3 Exact recorded Berkeley DT1 motion remains a separate validation requirement

The comparison harness does not change the earlier provenance rule: proxy-input comparisons are development diagnostics and must not be relabeled as final external validation.

## 13. Suggested immediate engineering use

The next Berkeley run can now be structured as:

```cpp
std::vector<ModelComparisonVariant> variants{
    research_variant,
    nist_asce41_17_variant,
    production_41_23_aci369_variant
};

ModelComparisonOptions opt;
opt.column_force_bindings = berkeley_columns;
opt.tangent_modes_at_story_peaks = true;
opt.tangent_mode_count = 3;

auto comparison = RCFrameModelComparison::run_frame2d(
    compiled_berkeley_frame,
    variants,
    monitored_hinges,
    common_analysis_settings,
    opt);
```

The comparison result is then sufficient to construct a dashboard with:

- global EDP comparison,
- story-drift history overlays,
- representative column `M-theta` overlays,
- axial/shear history overlays,
- IO/LS/CP/E/F event timing,
- instantaneous mode-shape comparison,
- story tangent-stiffness comparison.

That is the intended Phase 9D.6 handoff into the Berkeley constitutive/state-path audit.
