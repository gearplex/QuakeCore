# MVLEM and SFI-MVLEM wall elements

QuakeCore Phase 9J adds two-node, in-plane `mvlem` and `sfi_mvlem` wall components to the compiled `frame2d` model. Each element has the six external degrees of freedom `[UXi, UYi, RZi, UXj, UYj, RZj]`. This is the standard line-element topology used by the original OpenSees MVLEM and SFI-MVLEM formulations. It serves the same modeling role as a story-height shear-wall component with rigid floor interfaces. It is not the separate four-node, 24-DOF `MVLEM_3D` or `SFI_MVLEM_3D` formulation.

## Kinematics and forces

For height `h`, center-of-rotation ratio `c`, and macro-fiber coordinate `x_j` measured from the section centroid, the generalized strains are

```
epsilon_y,j = (-UY_i - x_j RZ_i + UY_j + x_j RZ_j) / h
gamma_xy    = (-UX_i + c h RZ_i + UX_j + (1-c) h RZ_j) / h
curvature   = (RZ_j - RZ_i) / h
```

The sign convention produces positive `gamma_xy` for positive top translation with fixed end rotations. Output `shear_deformation = h gamma_xy`.

MVLEM assigns a uniaxial concrete and steel law to each fiber. The fiber stress is `(1-rho) sigma_c + rho sigma_s`. Flexural and axial forces follow virtual work with fiber volume `h b_j t_j`. An independent uniaxial shear force-versus-shear-displacement spring acts at `c h`; its initial and trial stiffness therefore have force/length units. Flexure and shear are uncoupled except through global equilibrium.

SFI-MVLEM replaces each uniaxial fiber and the common shear spring with a plane-stress panel. The panel strain is `[epsilon_x,j, epsilon_y,j, gamma_xy]`. The internal transverse deformation is solved locally so `sigma_x,j = 0`. QuakeCore then statically condenses the full 3x3 panel tangent:

```
D_condensed = D_aa - D_ax * inv(D_xx) * D_xa
```

where the retained strain coordinates are `[epsilon_y, gamma_xy]`. This preserves all axial-shear tangent coupling supplied by the panel law without placing massless internal DOFs in the global matrix. Local Newton iterations use backtracking, a relative stress residual, and transactional trial states. A failed local solve raises a recoverable `ConstitutiveIntegrationError`, allowing the transient driver to roll back and subdivide the step.

Both elements support lumped self-mass from `density * h * sum(b_j t_j)`, divided equally between end nodes in each translational direction. Fixed-end mass contributes reactions but is eliminated from the active equations. Element gravity loads and evolving wall geometric stiffness are not implemented in this release.

## Material laws

MVLEM supports these explicit uniaxial laws:

| Type | Meaning | Required fields |
|---|---|---|
| `elastic` | Linear stress-strain | `E` |
| `steel_bilinear` | OpenSees Steel01-compatible kinematic bilinear law | `E`, `fy`, `b` |
| `concrete01` | OpenSees Concrete01-compatible compression envelope and cyclic unloading/reloading, zero tension | `fc`, `epsc`, `fcu`, `epsu`, all compression values negative except `fcu` may be zero |

SFI-MVLEM supports:

| Type | Meaning |
|---|---|
| `elastic_plane_stress` | Isotropic elastic plane stress with `E`, `nu` |
| `layered_plane_stress` | Elastic isotropic background plus any number of oriented uniaxial layers |
| `fixed_angle_rc` | Two fixed orthogonal concrete struts plus horizontal and vertical steel layers and an optional elastic background |

`fixed_angle_rc` is a transparent, history-dependent panel approximation. It does not implement the complete OpenSees FSAM model. In particular, the Phase 9J law does not reproduce FSAM crack-angle history, aggregate-interlock/friction, dowel-action, biaxial concrete softening, or the ConcreteCM/SteelMPF calibrations used in the PEER validation work. An input declaring `type: "FSAM"` is rejected so the software cannot silently substitute the simpler law.

## JSON model definition

Walls are listed under `model.walls`. Node `i` must be below node `j` at the same x coordinate. IDs share one namespace with frame members and hinges.

A compact MVLEM definition is:

```json
{
  "id": 1,
  "type": "mvlem",
  "i": 1,
  "j": 2,
  "c": 0.4,
  "density": 0.0,
  "fibers": [
    {
      "width": 0.5,
      "thickness": 0.2,
      "rho": 0.02,
      "concrete": {"type":"concrete01","fc":-30000,"epsc":-0.002,"fcu":-6000,"epsu":-0.006},
      "steel": {"type":"steel_bilinear","E":200000000,"fy":400000,"b":0.01}
    }
  ],
  "shear": {"type":"bilinear","stiffness":100000,"yield_force":200,"hardening_ratio":0.02}
}
```

At least two fibers are required. Every width and thickness is positive; `0 <= rho < 1`; `0 <= c <= 1`; and density is nonnegative.

An SFI-MVLEM panel uses the same geometry with `panels` in place of `fibers` and no separate shear material:

```json
{
  "id": 1,
  "type": "sfi_mvlem",
  "i": 1,
  "j": 2,
  "c": 0.4,
  "panels": [
    {
      "width": 0.5,
      "thickness": 0.2,
      "material": {
        "type": "fixed_angle_rc",
        "background_E": 1000000,
        "background_nu": 0.2,
        "angle_deg": 35,
        "concrete": {"type":"concrete01","fc":-30000,"epsc":-0.002,"fcu":-6000,"epsu":-0.006},
        "steel_x": {"type":"steel_bilinear","E":200000000,"fy":400000,"b":0.01},
        "steel_y": {"type":"steel_bilinear","E":200000000,"fy":400000,"b":0.01},
        "rho_x": 0.01,
        "rho_y": 0.02
      }
    }
  ]
}
```

Optional `local_max_iterations` and `local_relative_tolerance` control transverse panel equilibrium. Defaults are 60 and `1e-10`. Do not relax these controls to make an invalid panel law converge; diagnose the panel branch and tangent.

## Output and solver behavior

Full NRHA history rows include a `walls` object keyed by wall ID. Each wall reports nodal restoring force, curvature, shear deformation, fiber or panel strains and stresses, and the SFI transverse-stress residual. Summary output reports peak absolute shear, axial force, bottom moment, curvature, and shear deformation. IDA uses the same wall behavior and existing drift-threshold taxonomy.

A coupled wall tangent is not diagonal low rank in the hinge basis. Requested Woodbury analyses therefore route each evolving wall Newton iteration through exact same-pattern sparse refactorization. This preserves correctness but does not yet give the hinge-only speed advantage. Craig-Bampton's scalar constitutive-bank reduction is rejected for models containing coupled wall tangents.

## Current engineering limits

Phase 9J is suitable for numerical development and controlled comparison fixtures. It is not yet a general concrete-wall evaluation implementation because it lacks:

- full FSAM or another independently validated nonlinear RC panel law;
- gravity equilibrium and transfer of axial-load state into NRHA;
- wall P-Delta, buckling, axial failure, shear failure, bar buckling/fracture, and acceptance-limit state tracking;
- coupling beams, wall openings, flanges, nonrectangular 3D walls, and out-of-plane behavior;
- calibration and blind prediction against measured wall cyclic tests;
- code parameter providers and traceability for ACI 369.1-22 / ASCE 41-23 decisions.

## Primary references

- Kolozvari, Orakcal, and Wallace, *Shear-Flexure Interaction Modeling for Reinforced Concrete Structural Walls and Columns under Reversed Cyclic Loading*, PEER 2015/12: https://peer.berkeley.edu/sites/default/files/webpeer-2015-12-kristijan_kolozvari_kutay_orakcal_john_wallace.pdf
- OpenSees MVLEM description: https://opensees.berkeley.edu/wiki/index.php/MVLEM_-_Multiple-Vertical-Line-Element-Model_for_RC_Walls
- OpenSees SFI-MVLEM description: https://opensees.berkeley.edu/wiki/index.php/SFI_MVLEM_-_Cyclic_Shear-Flexure_Interaction_Model_for_RC_Walls
- OpenSees FSAM description: https://opensees.berkeley.edu/wiki/index.php/FSAM_-_2D_RC_Panel_Constitutive_Behavior
