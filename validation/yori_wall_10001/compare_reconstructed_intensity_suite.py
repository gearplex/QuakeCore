#!/usr/bin/env python3
"""Multi-intensity QuakeCore/OpenSees engineering comparison for YORi 10001.

The supplied acceleration record is multiplied by dimensionless scale factors.
These factors are ground-motion/PGA multipliers only; they are NOT Sa(T1),
collapse capacities, FEMA P-695 intensity measures, or code qualification.

For each scale, QuakeCore runs its documented robust transient path while the
independent OpenSeesPy 3.8.0 model runs the same nominal Newmark/Newton setup.
Engineering response is compared only through the nominal steps that both
programs complete. OpenSees nonconvergence is reported, not reclassified as
physical collapse or as a QuakeCore error.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import numpy as np
import openseespy.opensees as ops

from compare_reconstructed_exact_10001 import build, gravity


def qwall(row, wid):
    walls = row["walls"]
    return walls[str(wid)] if str(wid) in walls else walls[wid]


def setup_opensees(job):
    build(job)
    grav = gravity(job)
    if not grav["ok"]:
        raise RuntimeError("OpenSees gravity analysis did not converge")
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    rec = job["records"][0]
    ops.timeSeries("Path", 10, "-dt", rec["dt"], "-values", 0.0, *rec["acceleration"], "-useLast")
    ops.pattern("UniformExcitation", 2, 1, "-accel", 10)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-8, 100, 0, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    return grav


def max_relative(a, b):
    out = 0.0
    for x, y in zip(a, b):
        out = max(out, abs(float(x) - float(y)) / max(abs(float(x)), abs(float(y)), 1e-30))
    return out


def run_scale(base_job, quake_exe: Path, scale: float, temp_dir: Path):
    job = copy.deepcopy(base_job)
    job["name"] = f"{base_job.get('name','YORi_10001')}_scale_{scale:g}"
    rec = job["records"][0]
    rec["acceleration"] = [scale * float(v) for v in rec["acceleration"]]
    rec["name"] = f"{rec.get('name','record')}_x{scale:g}"

    job_path = temp_dir / f"job_scale_{scale:g}.json"
    q_path = temp_dir / f"quake_scale_{scale:g}.json"
    job_path.write_text(json.dumps(job))
    started = time.perf_counter()
    proc = subprocess.run([str(quake_exe), str(job_path), str(q_path)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    q_wall = time.perf_counter() - started
    if not q_path.exists():
        return {
            "scale_factor": scale,
            "quakecore_process_returncode": proc.returncode,
            "quakecore_process_wall_seconds": q_wall,
            "quakecore_output_missing": True,
            "quakecore_stdout_tail": proc.stdout[-4000:],
        }
    q = json.loads(q_path.read_text())
    qrun = q["runs"][0]
    qh = qrun.get("history", [])

    grav = setup_opensees(job)
    nodes = job["model"]["story_nodes"]
    heights = np.asarray([ops.nodeCoord(n, 2) for n in nodes], dtype=float)
    wall_ids = [int(w["id"]) for w in job["model"]["walls"]]
    oh = []
    oforce = {wid: [] for wid in wall_ids}
    os_iterations = 0
    os_started = time.perf_counter()
    os_code = 0
    for _ in rec["acceleration"]:
        os_code = ops.analyze(1, rec["dt"])
        if os_code:
            break
        u = np.asarray([ops.nodeDisp(n, 1) for n in nodes], dtype=float)
        drift = np.diff(np.r_[0.0, u]) / np.diff(np.r_[0.0, heights])
        oh.append([ops.getTime(), *u.tolist(), *drift.tolist()])
        os_iterations += ops.testIter()
        for wid in wall_ids:
            oforce[wid].append(np.asarray(ops.eleForce(wid), dtype=float))
    os_elapsed = time.perf_counter() - os_started

    qa = np.asarray([[r["time_s"], *r["floor_displacement"], *r["story_drift_ratio"]] for r in qh], dtype=float)
    oa = np.asarray(oh, dtype=float)
    n = min(len(qa), len(oa))
    if n:
        qcommon = qa[:n]
        ocommon = oa[:n]
        max_u = float(np.max(np.abs(qcommon[:, 1:4] - ocommon[:, 1:4])))
        max_d = float(np.max(np.abs(qcommon[:, 4:7] - ocommon[:, 4:7])))
        rmse = [
            float(100.0 * np.sqrt(np.mean((qcommon[:, 4+i] - ocommon[:, 4+i])**2)) /
                  max(np.max(np.abs(qcommon[:, 4+i])), np.max(np.abs(ocommon[:, 4+i])), 1e-30))
            for i in range(3)
        ]
        qp = np.max(np.abs(qcommon[:, 4:7]), axis=0)
        op = np.max(np.abs(ocommon[:, 4:7]), axis=0)
        peak_drift_rel = max_relative(qp, op)
    else:
        max_u = max_d = float("inf")
        rmse = [float("inf")] * 3
        qp = op = np.asarray([float("inf")] * 3)
        peak_drift_rel = float("inf")

    per_wall = {}
    max_wall_peak_rel = 0.0
    max_wall_hist_abs = 0.0
    for wid in wall_ids:
        m = min(n, len(oforce[wid]))
        if m == 0:
            continue
        qf = np.asarray([qwall(r, wid)["nodal_force"] for r in qh[:m]], dtype=float)
        of = np.asarray(oforce[wid][:m], dtype=float)
        max_hist = float(np.max(np.abs(qf - of)))
        max_wall_hist_abs = max(max_wall_hist_abs, max_hist)
        edps = {}
        for name, j in {"shear": 3, "axial_force": 4, "bottom_moment": 2}.items():
            qpk = float(np.max(np.abs(qf[:, j])))
            opk = float(np.max(np.abs(of[:, j])))
            rel = abs(qpk - opk) / max(abs(qpk), abs(opk), 1e-30)
            max_wall_peak_rel = max(max_wall_peak_rel, rel)
            edps[name] = {"quakecore_peak_abs": qpk, "opensees_peak_abs": opk, "relative_difference": rel}
        per_wall[str(wid)] = {"max_abs_force_history_error": max_hist, "edps": edps}

    q_stats = qrun.get("stats", {})
    return {
        "scale_factor": scale,
        "intensity_interpretation": "dimensionless acceleration/PGA multiplier; not Sa(T1)",
        "quakecore_process_returncode": proc.returncode,
        "quakecore_steps": len(qh),
        "quakecore_full_record": len(qh) == len(rec["acceleration"]),
        "quakecore_termination": qrun.get("termination"),
        "quakecore_termination_reason": qrun.get("termination_reason", qrun.get("reason")),
        "quakecore_integration_seconds": q_stats.get("elapsed_seconds"),
        "quakecore_process_wall_seconds": q_wall,
        "opensees_analysis_code": os_code,
        "opensees_steps": len(oh),
        "opensees_full_record": len(oh) == len(rec["acceleration"]),
        "opensees_end_time_s": float(ops.getTime()),
        "opensees_integration_seconds": os_elapsed,
        "opensees_iterations": os_iterations,
        "common_steps": n,
        "common_end_time_s": float(oh[n-1][0]) if n else None,
        "quakecore_peak_story_drift_common": qp.tolist(),
        "opensees_peak_story_drift_common": op.tolist(),
        "peak_story_drift_relative_difference": peak_drift_rel,
        "max_floor_displacement_difference_in": max_u,
        "max_story_drift_difference": max_d,
        "story_drift_rmse_percent_of_peak": rmse,
        "max_wall_peak_edp_relative_difference": max_wall_peak_rel,
        "max_abs_wall_force_history_error": max_wall_hist_abs,
        "walls": per_wall,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake-exe", type=Path, required=True)
    ap.add_argument("--scales", type=float, nargs="+", default=[0.5, 1.0, 2.0, 4.0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if any((not np.isfinite(s) or s <= 0.0) for s in args.scales):
        raise SystemExit("all scale factors must be finite and positive")

    base_job = json.loads(args.job.read_text())
    results = []
    with tempfile.TemporaryDirectory(prefix="yori_intensity_") as td:
        temp = Path(td)
        for scale in args.scales:
            results.append(run_scale(base_job, args.quake_exe.resolve(), float(scale), temp))

    valid = [r for r in results if not r.get("quakecore_output_missing")]
    summary = {
        "scope": "multi-intensity engineering-response comparison for the documented reconstructed YORi 10001 surrogate",
        "claim_boundary": "software/model verification only; scale factors are acceleration/PGA multipliers, not Sa(T1), collapse capacity, FEMA P-695 acceptance, code approval, or R-factor qualification",
        "opensees_version": ops.version(),
        "scales": [float(s) for s in args.scales],
        "results": results,
        "suite_summary": {
            "quakecore_full_record_scales": [r["scale_factor"] for r in valid if r["quakecore_full_record"]],
            "opensees_full_record_scales": [r["scale_factor"] for r in valid if r["opensees_full_record"]],
            "minimum_common_steps": min((r["common_steps"] for r in valid), default=0),
            "maximum_common_peak_drift_relative_difference": max((r["peak_story_drift_relative_difference"] for r in valid), default=float("inf")),
            "maximum_common_drift_rmse_percent_of_peak": max((max(r["story_drift_rmse_percent_of_peak"]) for r in valid), default=float("inf")),
            "maximum_common_wall_peak_edp_relative_difference": max((r["max_wall_peak_edp_relative_difference"] for r in valid), default=float("inf")),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
