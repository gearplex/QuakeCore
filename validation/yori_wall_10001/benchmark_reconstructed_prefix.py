#!/usr/bin/env python3
"""Paired common-prefix runtime benchmark for reconstructed YORi 10001.

Both solvers integrate the same truncated acceleration record, safely before the
current OpenSees full-record convergence stop. Timers are the solvers' existing
integration timers; model construction, gravity setup, and modal analysis are
outside the reported transient integration time in the comparator paths.

The measured repetitions are interleaved and alternate which solver runs first
to reduce runner-frequency / temporal bias. One unreported warm-up is executed
for each solver before timing samples are retained.
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import tempfile
from pathlib import Path

from compare_reconstructed_exact_10001 import run as run_opensees


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake-exe", type=Path, required=True)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--repetitions", type=int, default=7)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    if args.steps <= 0 or args.repetitions <= 0:
        raise SystemExit("steps and repetitions must be positive")

    job = json.loads(args.job.read_text())
    values = job["records"][0]["acceleration"]
    if args.steps > len(values):
        raise SystemExit(f"requested {args.steps} steps but record has {len(values)}")
    job["records"][0]["acceleration"] = values[: args.steps]
    job["name"] = f"{job.get('name', 'reconstructed')}_runtime_prefix_{args.steps}"

    quake_times: list[float] = []
    opensees_times: list[float] = []
    quake_iterations: list[int] = []
    opensees_iterations: list[int] = []
    execution_order: list[str] = []

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        prefix_job = td / "prefix_job.json"
        prefix_job.write_text(json.dumps(job))
        quake_exe = args.quake_exe.resolve()

        counter = 0
        def one_quake(retain: bool):
            nonlocal counter
            result_path = td / f"quake_{counter}.json"
            counter += 1
            proc = subprocess.run(
                [str(quake_exe), str(prefix_job), str(result_path)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            if proc.returncode != 0:
                raise SystemExit(f"QuakeCore benchmark run failed:\n{proc.stdout}")
            q = json.loads(result_path.read_text())
            qrun = q["runs"][0]
            if len(qrun.get("history", [])) != args.steps:
                raise SystemExit(
                    f"QuakeCore completed {len(qrun.get('history', []))} "
                    f"history steps; expected {args.steps}"
                )
            if retain:
                quake_times.append(float(qrun["stats"]["elapsed_seconds"]))
                quake_iterations.append(int(qrun["stats"].get("iterations", 0)))

        def one_opensees(retain: bool):
            o = run_opensees(job)
            if not o["analysis_ok"] or o["steps_completed"] != args.steps:
                raise SystemExit(
                    f"OpenSees benchmark run incomplete: ok={o['analysis_ok']} "
                    f"steps={o['steps_completed']} expected={args.steps}"
                )
            if retain:
                opensees_times.append(float(o["elapsed_seconds"]))
                opensees_iterations.append(int(o["iterations"]))

        # Warm both code paths before collecting samples.
        one_quake(False)
        one_opensees(False)

        # Alternate run order while preserving one Q/OS pair per repetition.
        for i in range(args.repetitions):
            if i % 2 == 0:
                execution_order.append("quakecore_then_opensees")
                one_quake(True)
                one_opensees(True)
            else:
                execution_order.append("opensees_then_quakecore")
                one_opensees(True)
                one_quake(True)

    qmed = statistics.median(quake_times)
    omed = statistics.median(opensees_times)
    paired_speedups = [o / q for q, o in zip(quake_times, opensees_times)]
    report = {
        "scope": (
            "Reconstructed YORi 10001 exact-law common-prefix transient benchmark; "
            "same truncated acceleration record and step count"
        ),
        "claim_boundary": (
            "Runtime comparison only for this small 3-story surrogate on the GitHub runner; "
            "not a general performance claim"
        ),
        "steps": args.steps,
        "repetitions": args.repetitions,
        "warmup_runs_per_solver": 1,
        "execution_order": execution_order,
        "quakecore_seconds": quake_times,
        "opensees_seconds": opensees_times,
        "quakecore_median_seconds": qmed,
        "opensees_median_seconds": omed,
        "quakecore_to_opensees_median_ratio": qmed / omed,
        "quakecore_speedup_factor_vs_opensees": omed / qmed,
        "paired_speedup_factors": paired_speedups,
        "paired_speedup_median": statistics.median(paired_speedups),
        "quakecore_iterations": quake_iterations,
        "opensees_iterations": opensees_iterations,
        "note": (
            "Measured repetitions are paired/interleaved and alternate solver order after one "
            "warm-up per solver. QuakeCore's job-runner transient timer includes accepted-state "
            "observer/history bookkeeping; the OpenSees comparator records floor displacement "
            "and story drift each step. This is an application-path benchmark, not a microbenchmark "
            "of only linear solves."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
