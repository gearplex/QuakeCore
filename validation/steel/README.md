# Phase 9K steel-component verification

Build `steel_probe` and `quake_run`, then run:

```bash
PYTHONPATH=deps/python python3 validation/steel/compare_opensees.py \
  --probe build/steel_probe --runner build/quake_run --out evidence/steel

python3 validation/steel/verify_system.py \
  --runner build/quake_run --example examples/steel/steel_frame_nrha.json \
  --out evidence/steel-system

python3 validation/steel/render_figure.py \
  --comparison evidence/steel --output docs/figures/phase9k_steel_verification.png
```

The comparison intentionally separates mechanisms:

- the condensed steel member is compared with an explicit OpenSees assembly of an `elasticBeamColumn` and two `Steel01` zero-length end hinges;
- panel-zone and BRB cyclic force-deformation histories are compared with `Steel01`;
- the nonlinear power-law dashpot force is compared with OpenSees `Viscous` away from the zero-velocity regularization interval;
- a linear-damper SDOF transient is independently integrated in both engines.

`verify_system.py` exercises an integrated steel frame using full and requested-Woodbury strategies, repeats the analysis at half the time step, confirms nonlinear component activation, and runs a synthetic IDA drift-threshold fixture. Coupled member and velocity-dependent tangent changes intentionally route the Woodbury request to exact same-pattern direct refactorization.

These are software-verification cases. They do not establish experimental validation, AISC 342 acceptance criteria, BRB fatigue/fracture behavior, panel-zone geometry rules, or manufacturer-specific damper qualification.
