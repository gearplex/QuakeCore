# Phase 9D — Code-generated RC column hinge architecture

## Purpose

Phase 9D separates the nonlinear material algorithm from the engineering rules
that generate reinforced-concrete column modeling parameters.  This prevents a
research calibration from being mistaken for an ASCE 41 model and allows
benchmark and production code editions to coexist.

## Architectural split

1. `ASCE41HingeMaterial` is the hysteretic/material engine.
2. `ASCE41BackboneShape` selects envelope topology without changing the solver.
3. `RCColumnASCE41Provider` converts already-resolved engineering/code values
   into a material definition while retaining edition and provenance metadata.
4. Actual ASCE/ACI table coefficients are intentionally not embedded.  The
   production provider requires values resolved from an authorized source.

## Backward compatibility

`ASCE41BackboneShape::ResearchExtendedCDE` is the default. Existing QuakeCore
models therefore retain the prior C-D residual-E and optional E-F behavior
unless explicitly changed.

## NIST benchmark mode

`RCColumnASCE41Provider::nist_asce41_17_benchmark()` sets:

- explicit B (`My`) and C (`Mc`) strengths;
- a straight C-to-E degrading branch;
- F equal to E for the benchmark component-failure convention;
- user-supplied IO/LS/CP plastic-rotation thresholds;
- no hidden code-table coefficients.

The cyclic degradation controls remain independent. This is intentional because
backbone provenance and hysteretic path rules are different modeling decisions.

## Production ASCE 41-23 / ACI 369.1-22 mode

`RCColumnASCE41Provider::asce41_23_aci369_1_22()` requires explicit provenance
and resolved modeling/acceptance parameters. The API is ready for a future
authorized table/rules implementation without reconstructing or guessing
copyrighted coefficients.

## Berkeley fidelity correction

The working Berkeley Phase-8 validation source is corrected to use the
published longitudinal reinforcement assignment:

- A/B nonductile columns: 0.049 in² bar area, fy = 70 ksi;
- C/D ductile columns: 0.11 in² bar area, fy = 64 ksi.

This correction is treated as model fidelity, not EDP calibration.

## Next implementation gate

1. Add a section-level RC column input object and authorized rules lookup layer.
2. Add a two-pass maximum-compression parameter update provider.
3. Add full moving-surface P-M return mapping as a separate constitutive layer.
4. Run three-way Berkeley comparisons: legacy research backbone, NIST benchmark
   backbone, and production ASCE 41-23/ACI 369 resolved backbone.
5. Compare tangent modes and story stiffness evolution without retuning strength.
