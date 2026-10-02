# Archetype provenance and unresolved inputs

The pilot uses the public FEMA P-2018 report as a traceable archetype family, not as a claim that the original analytical model has been recovered. Appendix D documents enough to choose a meaningful target, but not enough to reproduce the model without additional assumptions or source files.

## Directly documented in FEMA P-2018 Appendix D

| Item | Pilot value | Source location | Use |
|---|---:|---|---|
| Stories | 6 | Table D-1 | Archetype selection |
| Walls | 1 | Table D-1 | Archetype selection |
| Wall length | 20 ft | Table D-1 | Reference wall geometry |
| Wall thickness | 10 in | Table D-1 | Reference wall geometry |
| Wall index, `WI` | 0.0004 | Table D-1 | Cross-check only |
| Effective period, `Te` | 1.20 sec | Table D-1 | Pushover-derived cross-check, not `T1_ref` |
| Reference column strength ratio | `Vp/Vn = 0.6` | Table D-1 | Flexure-governed reference family |
| Shear-critical comparison ratio | `Vp/Vn = 1.2` | Table D-1 family and footnote | Target for degrading-column variant |
| Wall detailing | Unconfined, #4 bars at 12 in on center | Section D.2.1 | Detailing family; orientation/face application still requires confirmation |
| Frame elements | Distributed plasticity beams and columns | Section D.2.2 | Modeling approach |
| Joints | Elastic springs | Section D.2.2 and Figure D-1 | Modeling approach |
| Shear-critical columns | Zero-length shear and axial springs in series, Elwood limit-state relations | Section D.2.2 | `C2-RB` and `C2-SB` mechanism |
| Walls | Distributed plasticity element with rigid frame link | Section D.2.2 and Figure D-2 | Modeling approach |
| Concrete | Modified Kent-Park relation | Section D.2.2 | Constitutive family |
| Reinforcing steel | Giuffré-Menegotto-Pinto relation with buckling/post-buckling | Section D.2.2 | Constitutive family |
| Flexural elements | One force-based fiber element with five-point regularized integration | Section D.2.2 | Discretization target |
| Ground motions | FEMA P-695 far-field set | Section D.2.3 | Production record-set candidate |
| Collapse proxy | 50% loss of peak lateral resistance in any story; drift threshold identified by pushover | Section D.2.4 | `L2` definition |

Figure D-3 illustrates a six-story, six-bay, one-20-ft-wall system with shear-critical columns. The report's Table D-1 does not list that complete combination as a uniquely recoverable source model. The pilot therefore treats the figure as a mechanism example and Table D-1 as the parameter-family basis.

## Not sufficiently specified by Appendix D

The following must be recovered from cited source work, an archived input model, another documented appendix, or a predeclared reconstruction assumption:

- bay widths and story heights for the selected six-story model;
- beam, column, and joint dimensions and reinforcement;
- floor masses, gravity loads, and tributary areas;
- concrete and steel strengths;
- wall reinforcement direction, face distribution, and boundary steel;
- exact shear and axial limit-state material parameters;
- joint spring stiffnesses;
- damping formulation and coefficients;
- diaphragm/rigid-link details and leaning-column gravity system;
- regularized plastic-hinge length implementation details;
- geometric-nonlinearity transformation;
- solver, convergence, recovery, and time-step policies; and
- original input files and their hashes.

## Reconstruction rule

Every filled parameter in `archetype.json` must have one of these provenance classes:

1. `reported` — directly stated in a cited source;
2. `derived` — calculated from reported values with the equation stored;
3. `adopted` — a predeclared engineering assumption with sensitivity bounds; or
4. `calibrated` — fit to an independent target, with the objective and data stored.

No value may be labeled `reported` merely because it is typical. If the original FEMA P-2018 model files become available, they should be archived by hash and compared against—not silently substituted for—the reconstruction.

