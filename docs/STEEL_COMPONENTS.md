# Steel building components

QuakeCore Phase 9K adds four small-displacement, two-dimensional component families to the compiled `frame2d` model: concentrated-plasticity steel beams/columns, rotational panel-zone surrogates, axial bilinear BRBs, and axial power-law viscous dampers. They share the existing constraint, modal, Rayleigh damping, NRHA, IDA, rollback, subdivision, and response-output paths.

These are explicit analysis mechanisms, not section-design or acceptance-checking objects. QuakeCore does not select section properties, derive AISC/ASCE backbones, run a gravity analysis, or judge component acceptance. Every JSON definition therefore requires parameter provenance.

## Concentrated-plasticity steel member

`steel_members` use the ordinary external order `[UXi, UYi, RZi, UXj, UYj, RZj]`. The element contains an elastic Euler-Bernoulli span plus two nonlinear rotational end hinges. The hinge rotations are internal variables solved by a local two-equation Newton iteration and statically condensed, so users do not need coincident end nodes solely to represent member hinges.

Axial response is elastic. An optional positive `constant_compression` adds a tangent-only geometric stiffness, consistent with the existing frame-member baseline. It is not an equilibrated gravity load and does not evolve with response. The initial and trial six-by-six element blocks are assembled into an immutable global sparse pattern; nonlinear member changes currently use exact same-pattern sparse refactorization.

```json
{
  "id": 101,
  "type": "concentrated_plasticity",
  "role": "column",
  "i": 1,
  "j": 3,
  "E": 200000000.0,
  "A": 0.025,
  "I": 0.0015,
  "constant_compression": 0.0,
  "hinge_i": {"type":"bilinear", "k":10000000.0, "fy":450.0, "hardening_ratio":0.02},
  "hinge_j": {"type":"bilinear", "k":10000000.0, "fy":450.0, "hardening_ratio":0.02},
  "provenance": "engineer-defined verification parameters"
}
```

`role` is `beam` or `column`; it labels output and permits vertical columns in `story_cut_members`. It does not change mechanics. Hinge types are:

| Type | Required mechanical fields | Scope |
|---|---|---|
| `bilinear` | `k`, `fy`, `hardening_ratio` | Symmetric kinematic bilinear moment-rotation law |
| `asce41_parameterized` | `k`, `fy`, `a`, `b`; optional `f`, `c`, `io`, `ls`, `cp`, hardening fields | Existing symmetric research A-B-C-D-E-F hinge; values are user-supplied |
| `imk_peak_oriented` | `k`, `fy`, `up`, `upc`, `uu`, `fcap_fy`, `fres_fy`; optional deterioration fields | Existing independent IMK-family research material |

Member history includes nodal restoring force, axial deformation/force, chord rotation, end rotations, hinge rotations, end moments, and local iteration count. Peak output includes axial force, hinge rotation, and end moment at both ends.

## Panel zone

A `panel_zones` object is a relative-`RZ` zero-length spring between coincident nodes. Translational equal-DOF constraints must be supplied explicitly. Its `material` accepts the same three rotational material types as a steel-member hinge.

```json
{
  "id": 201,
  "i": 3,
  "j": 5,
  "material": {"type":"bilinear", "k":200000.0, "fy":350.0, "hardening_ratio":0.03},
  "provenance": "engineer-defined joint-shear backbone"
}
```

This is a transparent rotational joint-shear surrogate. It is not a finite eight-node joint region, does not infer web/flange/doubler-plate properties, and does not automatically convert panel shear distortion to joint rotation. Its recorder reports relative rotation, moment, and tangent.

## Buckling-restrained brace

A `brbs` object is an axial nonlinear element oriented along the line between its nodes. Its material is a symmetric bilinear force-deformation law.

```json
{
  "id": 301,
  "i": 1,
  "j": 4,
  "k": 150000.0,
  "fy": 500.0,
  "hardening_ratio": 0.01,
  "provenance": "engineer-defined BRB qualification backbone"
}
```

History and peak output contain axial deformation and force. The component does not model fatigue/fracture, cumulative plastic strain limits, tension/compression asymmetry, connection deformation, casing restraint, local buckling, or qualification-test acceptance.

## Viscous damper

A `viscous_dampers` object is an axial, memoryless velocity device. For `alpha = 1`,

`F = C v`.

For nonlinear exponents QuakeCore uses the smooth law

`F = C v (v^2 + v_r^2)^((alpha - 1)/2)`

and assembles its exact derivative with respect to velocity into the Newmark effective tangent. `regularization_velocity` is mandatory and positive for `alpha < 1`; it defines the small interval in which the law intentionally differs from the singular unregularized power law.

```json
{
  "id": 401,
  "i": 2,
  "j": 3,
  "coefficient": 300.0,
  "alpha": 0.5,
  "regularization_velocity": 0.0001,
  "provenance": "engineer-defined device law"
}
```

History and peak output contain axial deformation rate, force, and velocity tangent. This is not a Maxwell spring-dashpot model and has no relief valve, gap, lock-up, stroke/force limit, temperature/frequency dependence, or manufacturer qualification logic.

## Verification and present limits

The Phase 9K suite checks rigid-body behavior, tangent symmetry, finite-difference Jacobians, rollback, the rigid-hinge elastic limit, nonlinear activation, invalid input, a linear-damper/proportional-damping equivalence, and full versus requested-Woodbury response equality. Independent OpenSeesPy comparisons are reported in [PHASE9K_STEEL_REVIEW.md](PHASE9K_STEEL_REVIEW.md).

The steel member remains a small-displacement concentrated-plasticity model. Axial yielding/P-M interaction, gravity redistribution, evolving P-Delta, distributed plasticity, lateral-torsional/local buckling, fracture, connection/slab composite action, and automated ASCE 41-23/AISC 342-22 acceptance criteria are not implemented.

## Primary numerical references

- OpenSees Steel01 material: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel01.html
- OpenSees Viscous material: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Viscous.html
- OpenSees ViscousDamper/Maxwell material (not implemented here): https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/ViscousDamper.html
- OpenSees Joint2D element (finite joint-region comparison point): https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/Joint2D.html
