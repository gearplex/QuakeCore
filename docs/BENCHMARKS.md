# Benchmarks

Measured in the current ChatGPT Linux CPU environment on 2026-09-04/05. These are prototype measurements and vary modestly run-to-run. They are **not claims about OpenSees, xara, Perform-3D, or GPU performance**.

## Benchmark fairness changes

Two accounting corrections are important:

1. Woodbury single-record timing includes baseline factorization/setup.
2. The compiled frame full-Newton path no longer rebuilds triplets/topology every iteration. Hinge tangent changes are scattered directly into precompiled CSC value locations.

These changes lowered the apparent frame speedup and make the comparison substantially more credible.

The current full-Newton SuperLU backend still performs a fresh factorization/order analysis for every tangent solve. A same-pattern numerical-refactorization backend remains required before comparing the low-rank algorithm against an optimized production CPU baseline.

## 1. Phase 0 nonlinear shear benchmark

Representative run:

```text
stories=120 nonlinear_springs=20 steps=1200
full_seconds≈0.15
full_factorizations=1274
woodbury_seconds≈0.016
woodbury_factorizations=1
speedup≈9x
```

This matrix is extremely sparse and easy; it primarily validates exact NRHA path equivalence.

## 2. Sparse-grid low-rank microbenchmark

Representative 10,000-DOF grid, 200 candidate updates, rank 5 active, 300 solves:

```text
full factorization total≈7.7 s
Woodbury setup≈0.50 s
Woodbury repeated solves≈0.19 s
Woodbury total≈0.70 s
speedup including setup≈11x
relative solution error≈roundoff
```

This remains useful for a topology with more sparse fill than the shear model.

## 3. Compiled 2D nonlinear frame benchmark

Representative 20-story, 2-bay model:

```text
reduced DOF=220
elastic elements=100
candidate rotational hinges=80
steps=600

full_seconds≈0.16-0.17
full_factorizations=786
full_newton_iters=1386

lazy_woodbury_seconds≈0.039-0.041
woodbury_factorizations=1
woodbury_newton_iters=1386

speedup≈4.2x
max_active_rank=50
mean_active_rank≈9.30
ever-cached influence columns=56 of 80
```

The roof displacement and maximum story-drift histories agree within test tolerance. Full Newton and exact Woodbury use the same Newton iteration count.

### Size trend

Representative single-record runs after compiled CSC scatter and lazy influence caching:

| Model | Reduced DOF | Hinges | Full (s) | Lazy Woodbury (s) | Speedup | Max active rank |
|---|---:|---:|---:|---:|---:|---:|
| 10-story, 1-bay | 70 | 20 | ~0.044 | ~0.0048 | ~9.3x | 4 |
| 20-story, 1-bay | 140 | 40 | ~0.075 | ~0.014 | ~5.5x | 12 |
| 20-story, 2-bay | 220 | 80 | ~0.17 | ~0.040 | ~4.2x | 50 |
| 40-story, 1-bay | 280 | 80 | ~0.12 | ~0.036 | ~3.3x | 16 |
| 40-story, 2-bay | 440 | 160 | ~0.18-0.20 | ~0.068 | ~2.7x | 10 |
| 80-story, 2-bay | 880 | 320 | ~0.44 | ~0.23 | ~2.0x | 16 |

Interpretation: at these still-small 2D matrix sizes, SuperLU can refactor the global tangent cheaply, while the low-rank method pays dense/reduced-system and influence-cache costs. This is precisely why the next meaningful performance gate is a **larger 3D building** plus an optimized same-pattern direct-solver baseline.

## 4. Prepared record-suite benchmark

Representative 20-story, 2-bay model, 12 independent 600-step records:

```text
full_suite_seconds≈1.7-1.8
full_factorizations≈7900

lazy Woodbury setup≈0.0002 s
Woodbury record runtime≈0.39-0.40 s
Woodbury total≈0.39-0.40 s

suite speedup≈4.3x
```

The prepared suite shares only model/dt-dependent solver data. Each record still starts with zero committed nonlinear state.

## What Phase 1 establishes

1. The low-rank decomposition survives a real 2D frame coordinate/constraint/hinge compilation path.
2. Full-Newton and Woodbury NRHA histories can follow the same nonlinear iteration path.
3. Constraint elimination and nonlinear tangent scatter can be compiled before time stepping.
4. Lazy influence generation avoids precomputing all candidate hinge columns.
5. Record-suite solver preparation can be amortized across many independent motions.
6. Honest baseline optimization can materially reduce apparent speedup; performance claims must therefore be made only after matching production-quality solver practices.

## What remains unproven

- speedup against OpenSees/xara/Perform-3D;
- 3D concentrated-plasticity building performance;
- optimized same-pattern sparse numerical refactorization comparison;
- robustness near severe cyclic deterioration/collapse;
- IMK deterioration;
- fiber/distributed plasticity;
- GPU/cuDSS performance.

## 5. Compiled 3D nonlinear space-frame benchmark

Representative 10-story, 2×2-bay space frame with a nonlinear synthetic record:

```text
reduced DOF=780
elastic elements=210
candidate beam-end hinges=240
steps=300
excitation amplitude=2.5

full runtime≈1.4 s
full factorizations=334
full Newton iterations=634

lazy Woodbury runtime≈0.16 s
Woodbury baseline factorizations=1
Woodbury Newton iterations=634

speedup≈8.6x
max active rank=54
cached low-rank columns=54
```

The space-frame model has columns and beams in both horizontal directions, torsion and biaxial bending in the linear skeleton, and separate beam-end hinge directions for X and Y framing. This is the first benchmark in the repository using a genuinely three-dimensional global stiffness topology.

The same caveat remains: the reference full solver performs fresh SuperLU factorization rather than a production-quality same-pattern numerical refactorization.

## 6. Phase 3 optimized same-pattern and rigid-diaphragm benchmarks

Representative current runs after implementing `SamePattern_SameRowPerm`, true 3D MPC condensation, robust transient infrastructure, and the material-neutral solver-state boundary:

### 10-story, 2x2-bay space frame, yielding motion

```text
reduced DOF=780
elastic elements=210
candidate hinges=240
steps=300
excitation amplitude=4.0

fresh full SuperLU       ~=1.63 s, 363 global factorizations
same-pattern SuperLU     ~=0.53 s, 147 numerical refactorizations + setup
lazy Woodbury            ~=0.21 s, 1 baseline global factorization

same-pattern / fresh     ~=3.05x
Woodbury / fresh         ~=7.75x
Woodbury / same-pattern  ~=2.54x
max active rank=84
peak roof response=0.113700 for all three strategies
Newton iterations=663 for all three strategies
```

The **2.5x-class advantage over the same-pattern baseline**, rather than the 7-8x advantage over fresh factorization, is the relevant algorithmic result at this stage.

### Same building with true rigid-diaphragm condensation

```text
reduced DOF=540
mass matrix nnz=580
elastic elements=210
candidate hinges=240
steps=300

fresh full SuperLU       ~=0.73 s
same-pattern SuperLU     ~=0.22 s
lazy Woodbury            ~=0.12 s

same-pattern / fresh     ~=3.27x
Woodbury / fresh         ~=5.92x
Woodbury / same-pattern  ~=1.81x
max active rank=84
peak roof response=0.113657 for all three strategies
Newton iterations=664 for all three strategies
```

The lower Woodbury advantage after diaphragm condensation is expected: reducing the global system from 780 to 540 DOF makes conventional sparse numerical refactorization cheaper while leaving the nonlinear reduced rank similar.

### Interpretation

These remain internal CPU prototype measurements, not OpenSees/xara/Perform claims. They do establish that the low-rank method retains a measurable advantage after correcting two major baseline biases: fixed sparse topology and same-pattern direct refactorization. The next credibility gate is an equivalent external OpenSees/xara model on the same hardware.


## 7. Sparse nonlinear-basis checkpoint

The Phase 3 model initially retained `B` as a dense `n x m` array. That representation was mathematically correct but computationally inappropriate for concentrated plasticity: evaluating `q_j=b_j^T u` scanned every global DOF for every candidate/active hinge.

`SparseUpdateBasis` now stores each generalized deformation operator as a sparse column and the lazy Woodbury solver computes `B^T y`/Gram terms from those nonzeros. For the 10-story 2x2-bay benchmark:

```text
n = 780
m = 240 candidate hinges
dense B entries = 187,200
sparse B nnz = 480
```

Five consecutive representative runs at excitation amplitude 4.0 produced approximately:

```text
fresh full SuperLU       1.38 - 1.46 s
same-pattern SuperLU     0.39 - 0.43 s
lazy sparse-B Woodbury   0.053 - 0.062 s

Woodbury / fresh         23 - 28x
Woodbury / same-pattern  6.9 - 7.8x
max active rank           84
Newton iterations         663 for all strategies
peak roof response        0.113700 for all strategies
```

A larger 20-story, 2x2-bay case (`n=1560`, 480 candidate hinges, 960 basis nonzeros) at amplitude 6.0 produced one representative run of:

```text
fresh full SuperLU       ~=3.09 s
same-pattern SuperLU     ~=0.90 s
lazy sparse-B Woodbury   ~=0.143 s
Woodbury / same-pattern  ~=6.3x
max active rank           102
ever-cached columns       204
Newton iterations         655 for all strategies
```

These results remain internal and hardware-specific. They do, however, identify a concrete architectural lesson: **sparsity of the nonlinear deformation operators matters as much as sparsity of the global stiffness matrix**. A dense candidate-hinge basis would have hidden much of the low-rank algorithm's advantage and would scale poorly in memory for realistic towers.

## Phase 7: corotational, geometric-update, and IDA checkpoints

### Corotational Euler stability convergence

For a cantilever with analytical Euler load `Pcr = 2313.1885`, the finite-rotation energy reference converges as the member mesh is refined:

```text
1 element:   Pcr ratio 1.21615  (+21.62%)
2 elements:  Pcr ratio 1.05264  (+5.26%)
4 elements:  Pcr ratio 1.01316  (+1.32%)
8 elements:  Pcr ratio 1.00332  (+0.33%)
16 elements: Pcr ratio 1.00081  (+0.08%)
```

A small-motion dynamic cantilever at excitation amplitude 0.05 gives a corotational/linear peak-roof ratio of approximately 0.9988. At amplitude 1000 the ratio is approximately 0.974, demonstrating finite-rotation divergence from the linear model.

### Generalized geometric update crossover

The exact `A + U C U^T` representation reproduces updated-P-Delta tangent solves to numerical tolerance. It is not automatically faster. A representative 20-story/120-DOF column chain has 138 active generalized directions and gives roughly:

```text
generalized Woodbury: ~0.00235 s
direct tangent solve: ~0.000054 s
speed ratio:          ~0.023
```

This negative result is intentional: global geometric Woodbury is rejected when the update space becomes comparable to or larger than the global system. The cost model/direct fallback remains authoritative.

### IDA suite

Representative 8-record synthetic collapse suite with adaptive scale refinement:

```text
runs:                  89
physical collapses:     8
numerical failures:     0
median collapse PGA:    7.88450
descriptive beta_ln:    0.05136
1 worker wall time:    ~0.466 s
4 worker wall time:    ~0.165 s
parallel speedup:      ~2.8x
```

Collapse brackets and cause classifications are identical across worker counts. These are internal prototype timings, not external-software comparisons.
