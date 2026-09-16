#!/usr/bin/env python3
"""Full-record reconstructed YORi comparator with OpenSees time-step subdivision.

The established Gate-4 engineering acceptance uses the independently converged
common window and remains unchanged.  This diagnostic asks a stronger question:
can OpenSees 3.8.0 traverse the same full record if a failed nominal Newmark step
is recursively subdivided, as QuakeCore's robust integration path does?

This is software/model verification of the documented reconstructed surrogate,
not identity to the lost received input, physical validation, FEMA P-695
acceptance, collapse qualification, or an R-factor recommendation.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import numpy as np
import openseespy.opensees as ops

from compare_reconstructed_exact_10001 import build, gravity


def setup(job):
    build(job)
    grav = gravity(job)
    if not grav["ok"]:
        raise RuntimeError("OpenSees gravity analysis did not converge")
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    rec = job["records"][0]
    values = rec["acceleration"]
    ops.timeSeries("Path", 10, "-dt", rec["dt"], "-values", 0.0, *values, "-useLast")
    ops.pattern("UniformExcitation", 2, 1, "-accel", 10)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-8, 100, 0, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    return grav


def qwall(row, wid):
    walls = row["walls"]
    return walls[str(wid)] if str(wid) in walls else walls[wid]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--max-depth", type=int, default=6)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.max_depth < 0:
        raise SystemExit("max-depth must be nonnegative")

    job = json.loads(args.job.read_text())
    q = json.loads(args.quake.read_text())
    qrun = q["runs"][0]
    qh = qrun["history"]
    rec = job["records"][0]
    values = rec["acceleration"]
    if len(qh) != len(values):
        raise SystemExit(f"QuakeCore history is not full record: {len(qh)} vs {len(values)}")

    grav = setup(job)
    nodes = job["model"]["story_nodes"]
    heights = np.asarray([ops.nodeCoord(n, 2) for n in nodes], dtype=float)
    wall_ids = [int(w["id"]) for w in job["model"]["walls"]]

    failed_attempts = 0
    substep_commits = 0
    max_depth_used = 0

    def advance(h, depth=0):
        nonlocal failed_attempts, substep_commits, max_depth_used
        ok = ops.analyze(1, h)
        if ok == 0:
            if depth:
                substep_commits += 1
                max_depth_used = max(max_depth_used, depth)
            return True
        failed_attempts += 1
        if depth >= args.max_depth:
            return False
        max_depth_used = max(max_depth_used, depth + 1)
        half = 0.5 * h
        if not advance(half, depth + 1):
            return False
        return advance(half, depth + 1)

    history = []
    wall_force_history = {wid: [] for wid in wall_ids}
    started = time.perf_counter()
    completed = 0
    for _ in values:
        if not advance(rec["dt"], 0):
            break
        u = np.asarray([ops.nodeDisp(n, 1) for n in nodes], dtype=float)
        drift = np.diff(np.r_[0.0, u]) / np.diff(np.r_[0.0, heights])
        history.append([ops.getTime(), *u.tolist(), *drift.tolist()])
        for wid in wall_ids:
            wall_force_history[wid].append(np.asarray(ops.eleForce(wid), dtype=float))
        completed += 1
    elapsed = time.perf_counter() - started

    oh = np.asarray(history, dtype=float)
    qarr = np.asarray([[r["time_s"], *r["floor_displacement"], *r["story_drift_ratio"]] for r in qh], dtype=float)
    n = min(len(qarr), len(oh))
    if n:
        max_u = float(np.max(np.abs(qarr[:n, 1:4] - oh[:n, 1:4])))
        max_d = float(np.max(np.abs(qarr[:n, 4:7] - oh[:n, 4:7])))
        drift_rmse = [
            float(100.0 * np.sqrt(np.mean((qarr[:n, 4+i] - oh[:n, 4+i])**2)) /
                  max(np.max(np.abs(qarr[:n, 4+i])), np.max(np.abs(oh[:n, 4+i])), 1e-30))
            for i in range(3)
        ]
        q_peak_common = np.max(np.abs(qarr[:n, 4:7]), axis=0)
        o_peak_common = np.max(np.abs(oh[:n, 4:7]), axis=0)
    else:
        max_u = max_d = float("inf")
        drift_rmse = [float("inf")] * 3
        q_peak_common = o_peak_common = np.asarray([float("inf")] * 3)

    per_wall = {}
    max_wall_peak_rel = 0.0
    max_wall_history_abs = 0.0
    for wid in wall_ids:
        if not wall_force_history[wid]:
            continue
        of = np.asarray(wall_force_history[wid], dtype=float)
        qf = np.asarray([qwall(r, wid)["nodal_force"] for r in qh[:len(of)]], dtype=float)
        max_hist = float(np.max(np.abs(qf - of)))
        max_wall_history_abs = max(max_wall_history_abs, max_hist)
        edps = {}
        for name, j in {"shear": 3, "axial_force": 4, "bottom_moment": 2}.items():
            qp = float(np.max(np.abs(qf[:, j])))
            op = float(np.max(np.abs(of[:, j])))
            rel = abs(qp - op) / max(abs(qp), abs(op), 1e-30)
            max_wall_peak_rel = max(max_wall_peak_rel, rel)
            edps[name] = {"quakecore_peak_abs": qp, "opensees_peak_abs": op, "relative_difference": rel}
        per_wall[str(wid)] = {"max_abs_force_history_error": max_hist, "edps": edps}

    result = {
        "scope": "full-record robust OpenSees diagnostic for reconstructed YORi 10001 surrogate",
        "claim_boundary": "diagnostic software/model comparison only; not received-source identity, physical validation, FEMA P-695 acceptance, collapse qualification, code approval, or R-factor qualification",
        "opensees_version": ops.version(),
        "gravity_ok": grav["ok"],
        "steps_requested": len(values),
        "steps_completed": completed,
        "analysis_completed": completed == len(values),
        "end_time_s": float(ops.getTime()),
        "elapsed_seconds": elapsed,
        "failed_nominal_or_substep_attempts": failed_attempts,
        "substep_commits": substep_commits,
        "max_subdivision_depth_used": max_depth_used,
        "max_floor_displacement_difference_in": max_u,
        "max_story_drift_difference": max_d,
        "story_drift_rmse_percent_of_peak": drift_rmse,
        "quakecore_peak_story_drift_common": q_peak_common.tolist(),
        "opensees_peak_story_drift_common": o_peak_common.tolist(),
        "max_wall_peak_edp_relative_difference": max_wall_peak_rel,
        "max_abs_wall_force_history_error": max_wall_history_abs,
        "walls": per_wall,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
