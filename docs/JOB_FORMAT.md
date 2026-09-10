# QuakeCore Phase 9L job interface

`quake_run input.json result.json` executes a consistent-unit plane-frame NRHA or IDA job. The examples are numerical fixtures, not code-designed buildings. A successful numerical run is not an ASCE 41 assessment.

## Build and run

Prerequisites: a C++20 compiler, CMake >=3.20, SuperLU development headers/library, BLAS/LAPACK, and threads. On a Debian/Ubuntu development machine these are supplied by `build-essential cmake libsuperlu-dev liblapack-dev`. No internet access is required at runtime.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
build/quake_run examples/frame3_nrha.json frame3_result.json
build/quake_run examples/frame3_ida.json frame3_ida_result.json
```

The frame examples cover 3-story/10-story NRHA and a two-record IDA demonstration. `examples/walls` adds three-story MVLEM and SFI-MVLEM NRHA and IDA fixtures. `examples/steel/steel_frame_nrha.json` integrates steel beams, columns, a panel zone, BRB, and viscous damper. `examples/soil/single_pile_nrha.json` expands a vertical elastic pile with p-y, t-z, q-z, and dashpot components. The IDA records are synthetic, so they do not constitute an independent hazard suite.

Exit codes: **0** completed analysis or configured threshold reached; **2** invalid input/model or an uncaught setup error; **3** numerical failure or initial instability, with the result written. Always inspect every run's `termination`, not just the process status. Result files are local JSON and may contain partial histories for a failed run.

## Schema `quakecore.job.v1`

Top-level fields are `schema`, `name`, `units`, `model`, `analysis`, `records`, and optional `provenance`. Unknown fields are rejected. IDs must be 32-bit integers. Units are explicit labels; the runner does **not** convert between unit systems. Use one consistent system throughout: forces, lengths, seconds, nodal masses in force·s²/length, stiffness, moments, and ground acceleration in length/s². Thus acceleration supplied in g must first be multiplied by the acceleration of gravity in the chosen units. Rotations are radians. Force and moment convergence tolerances still use the raw assembled residual norm; dimensional residual scaling remains a development item.

`model.type` is `frame2d`. Supported arrays:

| Field | Row definition |
|---|---|
| `nodes` | `[id, x, y, mass_x, mass_y, rotational_mass]` |
| `fixities` | `[node_id, fix_ux, fix_uy, fix_rz]`, actual JSON booleans |
| `equal_dofs` | `[retained_node, constrained_node, "UX" / "UY" / "RZ"]` |
| `members` | `[id, node_i, node_j, E, A, I, constant_compressive_preload]` |
| `rayleigh` | `[alpha_mass, beta_initial_stiffness]` |
| `story_nodes` | One output node per floor, bottom to top, with increasing positive elevations |
| `story_cut_members` | Optional arrays of upward vertical member IDs whose signed end shears form each story cut |
| `walls` | Optional objects defining vertical two-node `mvlem` or `sfi_mvlem` components; see [WALL_ELEMENTS.md](WALL_ELEMENTS.md) |
| `steel_members` | Optional two-node concentrated-plasticity steel beams or columns with two internal end hinges |
| `panel_zones` | Optional zero-length rotational joint-shear springs between coincident nodes |
| `brbs` | Optional two-node axial bilinear buckling-restrained-brace surrogates |
| `viscous_dampers` | Optional two-node axial power-law velocity devices |
| `soil_springs` | Optional directional translational or relative-`RZ` zero-length springs between coincident nodes |
| `soil_dashpots` | Optional directional translational or relative-`RZ` zero-length power-law velocity devices |
| `pile_lines` | Optional vertical elastic pile meshes with explicit nodes, members, soil springs, and depth recorders |

`response_node` selects the roof/primary horizontal response. The supported drift output assumes the lowest reference elevation is zero and the base horizontal displacement is fixed at zero. Stories with moving foundations or multiple differential floor reference points require the C++ API and an explicit response reducer.

`hinges` contains coincident-node rotational springs. Their translations must be tied by the supplied constraints as appropriate; the runner does not silently add constraints. Bilinear objects require `id`, `i`, `j`, `type: "bilinear"`, `k`, `fy`, and `hardening_ratio` (0 <= ratio < 1). Symmetric `asce41_parameterized` hinges accept `k`, `fy`, `a`, `b`, optional `f`, `c`, `io`, `ls`, `cp`, `hardening_ratio`, and `hardening_stiffness`, with mandatory nonempty `provenance`. They use the existing straight C–E research backbone. Values are resolved by the engineer; the interface supplies no licensed code tables and does not establish code compliance.

Walls use the small-displacement, two-node in-plane formulation described in [WALL_ELEMENTS.md](WALL_ELEMENTS.md). Their six external DOFs join the same frame and constraint system. Full NRHA output includes wall force, strain, stress, curvature and shear-deformation histories; summary output includes wall peak envelopes. A requested Woodbury strategy falls back to exact same-pattern sparse refactorization for changing coupled wall tangents. The `fixed_angle_rc` panel is explicitly not FSAM and the runner rejects an FSAM label.

Steel component definitions and limitations are detailed in [STEEL_COMPONENTS.md](STEEL_COMPONENTS.md). All steel components require nonempty parameter provenance. `steel_members` use six external frame DOFs and statically condensed local end-hinge rotations; `role` must be `beam` or `column`. Each end hinge may use `bilinear`, `asce41_parameterized`, or `imk_peak_oriented`. Panel zones connect coincident nodes and use a nonlinear rotational `material`. BRBs use `k`, `fy`, and `hardening_ratio` in axial force-deformation space. Viscous dampers use `coefficient`, optional `alpha` (default 1), and `regularization_velocity`; a positive regularization is mandatory when `alpha < 1`.

Full NRHA history rows report steel-member nodal force, axial force/deformation, chord and hinge rotations, end moments, and local iterations; panel-zone rotation/moment; BRB axial deformation/force; and damper deformation rate/force/velocity tangent. Summary output reports component peak envelopes. These are response quantities, not automated acceptance determinations. Coupled steel-member and velocity-dependent device tangents currently route through exact same-pattern sparse refactorization even when `strategy: "woodbury"` is requested.

Soil and pile definitions are detailed in [SOIL_FOUNDATIONS.md](SOIL_FOUNDATIONS.md). Translational components use `q = n dot (u_j-u_i)`; rotational components use `"dof":"RZ"` and `q = RZ_j-RZ_i`. Named p-y/t-z/q-z adapters convert explicit engineer-supplied distributed capacities to nodal force and use `k0 = 0.5 Fult/u50`; they are not OpenSees Simple1 implementations. A pile line generates a uniform vertical beam mesh and fixed coincident soil nodes from explicit ID arrays. Full output includes depth-wise displacement, soil reaction, pile shear/moment, and toe response; summary output includes peak absolute profiles.

Members use the existing small-displacement formulation with optional constant-preload geometric stiffness. The preload influences the tangent but is **not** a gravity equilibrium analysis. Loads, geometric nonlinearity, material state and gravity redistribution are not interchangeable. The result explicitly reports `gravity_analysis_performed: false`. The 3D, P–M and corotational research APIs remain in C++; the JSON steel-member path exposes IMK-family rotational hinges but not axial-flexural interaction.

## Analysis and records

`analysis.type` is `nrha` or `ida`. `strategy` is `woodbury`, `same_pattern`, or `full`. Controls include `tolerance`, `relative_tolerance` (default true), `max_iterations`, `max_subdivisions`, `line_search`, and `initial_guess` (`kinematic`, default, or `previous_displacement`). With relative tolerance the accepted infinity norm is tolerance × max(1, norm of the applied force); otherwise the norm is absolute. Subdivision halves the time step and linearly interpolates the input; it is convergence recovery, not proof of integration accuracy.

For NRHA, `output_mode` is `full` (default), `summary` (substep envelopes without full histories), or `minimal` (core diagnostics and roof envelope only). `history_stride` is a positive integer. Full histories are recorded at original input endpoints, while reported demand envelopes include all accepted subdivisions and exclude rolled-back branches. Member story shears are elastic restoring-force cuts and do not include damping-force or inertia contributions to support reactions. Acceleration histories are absolute, obtained by adding the supplied base acceleration to relative horizontal acceleration.

Each record requires `name`, positive finite `dt`, a nonempty finite `acceleration` array, and `sample_convention: "step_end"`. Samples represent **a(dt), a(2dt), …**, with a zero displacement/velocity/acceleration initial state and zero assumed a(0). Do not include an extra t=0 sample. Nonzero initial excitation needs an explicitly equilibrated initial-state extension; simply discarding it is not a general recorded-motion importer. `provenance` should identify the waveform, units, component, scaling, processing, original time interval, source and checksum.

For IDA, add positive `scales`, nonnegative `refinements`, positive `workers`, `stop_after_first_collapse`, and a positive `drift_limit`. This runner uses that limit as a **configured response threshold**, not a verified dynamic-collapse definition or ASCE 41 acceptance limit. Document its meaning in job provenance. PGA is reported in the supplied acceleration units; the intensity is not automatically Sa(T1). The first crossing is bracketed and refined geometrically. Numerical failures remain unresolved and are identified separately. Brackets also report unresolved numerical gaps, nonmonotonic responses, and right censoring. The C++ suite's simple uncensored lognormal statistics are not a full censored fragility fit and are intentionally omitted from this interface.

## Independent numerical comparison

Install OpenSeesPy and NumPy in a separate verification environment, then run:

```bash
python tools/compare_opensees.py examples/frame3_nrha.json ./build/quake_run comparison3 --repeats 7 --no-timing-recorders
python tools/compare_opensees.py examples/frame10_nrha.json ./build/quake_run comparison10 --repeats 7 --no-timing-recorders
PYTHONPATH=/path/to/openseespy python validation/soil/compare_opensees_pile.py ./build/quake_run phase9l_opensees.json
```

Omit `--no-timing-recorders` for the recording-inclusive comparison. The adapter requires one fixed-step record, bilinear hinges, zero geometric preload, absolute tolerance and the previous-displacement Newton guess. OpenSees uses Steel01, transformation constraints, UmfPack, full Newton, Newmark β=0.25/γ=0.5, and initial-stiffness Rayleigh damping including zeroLength hinges. It advances the record with native `analyze(n, dt)`, not a Python loop. A separate recorded pass checks periods, displacement/drift/acceleration histories and actual hinge yielding before timing is interpreted. Model construction and modal analysis are excluded from both integration timers. Results include all repetitions, scope and engine version. This is not a degrading RC/P–M validation.
