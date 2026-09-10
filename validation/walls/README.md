# Wall verification harness

Build `wall_probe`, `quake_run`, and tests, then run the three Python adapters with OpenSeesPy and NumPy available:

```bash
python validation/walls/compare_opensees.py --probe build/wall_probe --runner build/quake_run --out evidence/comparison
python validation/walls/verify_rc_panel.py --probe build/wall_probe --out evidence/rc_panel
python validation/walls/verify_dynamic.py --runner build/quake_run --out evidence/dynamic --examples examples/walls
```

`compare_opensees.py` applies 740 cyclic static increments to elastic and nonlinear MVLEM/SFI fixtures. `verify_dynamic.py` compares three 1,200-step wall-frame NRHAs and runs two synthetic IDA examples. `verify_rc_panel.py` replays the fixed-angle panel layer strains through OpenSees uniaxial laws; that replay is not a global SFI or FSAM validation.
