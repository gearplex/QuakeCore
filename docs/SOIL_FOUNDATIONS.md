# Soil springs, dashpots, and vertical pile lines

QuakeCore Phase 9L adds an explicit small-displacement foundation layer to the compiled `frame2d` model. It is intended for transparent mechanics verification and engineer-controlled building models. It is not a site-response or geotechnical parameter-generation program.

## Sign and topology contract

A soil spring or dashpot connects coincident nodes `i` and `j` and uses either a user-supplied global translation direction `[nx, ny]` or `"dof": "RZ"` for relative rotation. QuakeCore normalizes a translation vector and defines positive deformation as

`q = n dot (u_j - u_i)`.

For an `RZ` component, `q = RZ_j - RZ_i`. The usual topology is a fixed soil node at `i` and a structural or pile node at `j`. Coincidence is enforced. A soil spring is compiled as one sparse generalized-deformation column, so it retains the existing low-rank nonlinear solver path. A directional or rotational dashpot uses the same kinematics and contributes its consistent velocity derivative to the Newmark tangent. Rotational soil components accept only the direct `bilinear` or `asymmetric_elastic_perfectly_plastic` material types because p-y/t-z/q-z are translational adapters.

Every JSON soil object requires nonempty `provenance`. QuakeCore does not select p-y, t-z, or q-z parameters.

## Direct spring laws

`soil_springs` supports:

| Type | Behavior |
|---|---|
| `bilinear` | Symmetric kinematic bilinear law with `k`, `fy`, and optional `hardening_ratio` |
| `asymmetric_elastic_perfectly_plastic` | Independent `positive_capacity` and `negative_capacity`; a zero capacity creates an opening/no-resistance branch and committed permanent gap |
| `py_bilinear` | Explicit distributed lateral resistance to nodal force conversion |
| `tz_bilinear` | Explicit interface stress and pile surface tributary area to nodal force conversion |
| `qz_bilinear` | Explicit toe bearing stress and toe area to an asymmetric compression/suction spring |

Example:

```json
{
  "id": 41,
  "i": 101,
  "j": 1,
  "direction": [0.0, 1.0],
  "type": "qz_bilinear",
  "ultimate_bearing_stress": 2500.0,
  "toe_area": 0.20,
  "displacement_50": 0.010,
  "suction_ratio": 0.0,
  "compression_sign": -1,
  "provenance": "Geotechnical report section and engineer calculation reference"
}
```

For the named adapters,

`Pult = pult_per_length * tributary_length`

`Tult = tult_stress * pile_perimeter * tributary_length`

`Qult = qult_stress * toe_area`

and the transparent initial stiffness approximation is

`k0 = 0.5 * Fult / displacement_50`.

Thus the elastic branch reaches half the supplied ultimate nodal resistance at the supplied 50-percent displacement and reaches the force cap at twice that displacement. The p-y and t-z adapters require zero post-yield hardening so `Fult` remains a true cap. The q-z adapter uses independent compression and suction caps. `compression_sign = -1` is appropriate for a vertical spring whose positive global direction is upward.

These adapters are not implementations of OpenSees `PySimple1`, `TzSimple1`, or `QzSimple1`. Those materials have additional series components, gap rules, and radiation-damping behavior. Use `soil_dashpots` when a separately defined local radiation dashpot is intended.

## Directional dashpots

`soil_dashpots` uses the existing memoryless power-law device in a translation direction or relative `RZ`:

`F = C v` for `alpha = 1`,

or the documented regularized power law for `alpha < 1`. The explicit direction permits coincident nodes:

```json
{
  "id": 51,
  "i": 101,
  "j": 1,
  "direction": [1.0, 0.0],
  "coefficient": 8.0,
  "alpha": 1.0,
  "regularization_velocity": 0.0,
  "provenance": "Engineer-defined radiation damping calculation"
}
```

## Vertical pile-line compiler

`pile_lines` expands a vertical pile into ordinary elastic 2D beam-column segments, a fixed coincident soil node at each pile node, optional p-y and t-z springs at each depth, and an optional q-z toe spring. The head node must already exist. All remaining pile nodes, soil nodes, member IDs, and spring IDs are explicit, which prevents hidden identifier generation and makes output traceable.

The Phase 9L compiler uses uniform vertical spacing. Nodal tributary lengths are one-half segment at the head and toe and one full segment at interior nodes. The complete executable example is [`examples/soil/single_pile_nrha.json`](../examples/soil/single_pile_nrha.json).

Pile history output includes, by depth:

- lateral and vertical pile-node displacement;
- lateral and shaft spring deformation and force;
- member end moment and shear for every pile segment;
- toe spring deformation and force.

Summary output contains the corresponding peak absolute profiles. The ordinary `soil_springs` and `soil_dashpots` maps also report component-level peaks.

## Validation and limits

The native tests cover tributary conversions, asymmetric contact/gap state, rollback, directional kinematics, dashpot force, malformed topology, and a 200-element beam-on-Winkler closed form. An independently assembled OpenSeesPy model reproduces the complete nonlinear pile-head transient history for the synthetic fixture. Results are reported in [PHASE9L_SOIL_PILES_REVIEW.md](PHASE9L_SOIL_PILES_REVIEW.md).

Phase 9L does not include nonlinear pile-section response, nonvertical or nonuniform pile meshes, pile caps or groups, p-multipliers/group shadowing, liquefaction, pore-pressure generation, gapping rules beyond the defined one-state asymmetric law, cyclic degradation/pinching, free-field soil columns, kinematic interaction, or SSI-compatible input motion. ASCE 41-23 acceptance and site-specific geotechnical adequacy are not assessed.

## Primary comparison references

- OpenSees zeroLength element: https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/zeroLength.html
- OpenSees PySimple1: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/PySimple1.html
- OpenSees TzSimple1: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/TzSimple1.html
- OpenSees QzSimple1: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/QzSimple1.html
- OpenSees laterally loaded pile example: https://opensees.berkeley.edu/wiki/index.php/Laterally-Loaded_Pile_Foundation
