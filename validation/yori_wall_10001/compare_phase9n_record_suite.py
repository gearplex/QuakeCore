#!/usr/bin/env python3
"""Phase 9N multi-record QuakeCore/OpenSees exact-law comparison.

This is a software/model external-validation harness for the documented YORi
10001 reconstructed surrogate.  It does not calculate FEMA P-695 collapse
statistics or establish physical/code qualification.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import time

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import numpy as np
import openseespy.opensees as ops

from compare_reconstructed_exact_10001 import build, gravity
from phase9n_materialize_record_job import load_record, materialize


def qwall(row, wid):
    walls = row["walls"]
    return walls[str(wid)] if str(wid) in walls else walls[wid]


def max_relative(a, b):
    out = 0.0
    for x, y in zip(a, b):
        out = max(out, abs(float(x) - float(y)) / max(abs(float(x)), abs(float(y)), 1e-30))
    return out


def setup_opensees(job):
    build(job)
    grav = gravity(job)
    if not grav["ok"]:
        raise RuntimeError("OpenSees gravity analysis did not converge")
    lam = np.asarray(ops.eigen(3), dtype=float)
    periods = (2.0 * np.pi / np.sqrt(lam)).tolist()
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
    return grav, periods


def classify_quake(qrun, requested):
    steps = len(qrun.get("history", []))
    if steps == requested and qrun.get("termination") not in {"numerical_failure", "initial_instability"}:
        return "completed"
    return "numerical_noncompletion"


def run_record(base_job, source_archive, manifest, record_id, scale, quake_exe, output_dir):
    meta, source_g = load_record(source_archive, manifest, record_id)
    job = materialize(base_job, meta, source_g, scale)
    record_dir = output_dir / record_id
    record_dir.mkdir(parents=True, exist_ok=True)
    job_path = record_dir / "job.json"
    q_path = record_dir / "quakecore.json"
    job_path.write_text(json.dumps(job, indent=2) + "\n")

    q_started = time.perf_counter()
    proc = subprocess.run(
        [str(quake_exe), str(job_path), str(q_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    q_process_wall = time.perf_counter() - q_started
    (record_dir / "quakecore_stdout.txt").write_text(proc.stdout)
    if not q_path.exists():
        return {
            "record_id": record_id,
            "pair_id": meta["pair_id"],
            "scale": scale,
            "quakecore_process_returncode": proc.returncode,
            "quakecore_process_wall_seconds": q_process_wall,
            "quakecore_output_missing": True,
            "passed_frozen_gate4_metrics": False,
        }

    q = json.loads(q_path.read_text())
    qrun = q["runs"][0]
    qh = qrun.get("history", [])
    rec = job["records"][0]
    requested = len(rec["acceleration"])

    grav, o_periods = setup_opensees(job)
    nodes = job["model"]["story_nodes"]
    heights = np.asarray([ops.nodeCoord(n, 2) for n in nodes], dtype=float)
    wall_ids = [int(w["id"]) for w in job["model"]["walls"]]
    oh = []
    oforce = {wid: [] for wid in wall_ids}
    os_iterations = 0
    os_started = time.perf_counter()
    os_code = 0
    for i, ag in enumerate(rec["acceleration"]):
        os_code = ops.analyze(1, rec["dt"])
        if os_code:
            break
        u = np.asarray([ops.nodeDisp(n, 1) for n in nodes], dtype=float)
        drift = np.diff(np.r_[0.0, u]) / np.diff(np.r_[0.0, heights])
        # UniformExcitation nodal acceleration is relative; add support motion to
        # form the absolute floor acceleration used by QuakeCore output.
        a_abs = np.asarray([ops.nodeAccel(n, 1) + ag for n in nodes], dtype=float)
        oh.append([ops.getTime(), *u.tolist(), *drift.tolist(), *a_abs.tolist()])
        os_iterations += ops.testIter()
        for wid in wall_ids:
            oforce[wid].append(np.asarray(ops.eleForce(wid), dtype=float))
    os_elapsed = time.perf_counter() - os_started

    qa = np.asarray([
        [r["time_s"], *r["floor_displacement"], *r["story_drift_ratio"], *r["floor_absolute_acceleration"]]
        for r in qh
    ], dtype=float)
    oa = np.asarray(oh, dtype=float)
    n = min(len(qa), len(oa))
    if n:
        qc, oc = qa[:n], oa[:n]
        max_u = float(np.max(np.abs(qc[:, 1:4] - oc[:, 1:4])))
        max_d = float(np.max(np.abs(qc[:, 4:7] - oc[:, 4:7])))
        drift_rmse = [
            float(100.0 * np.sqrt(np.mean((qc[:, 4+i] - oc[:, 4+i]) ** 2)) /
                  max(np.max(np.abs(qc[:, 4+i])), np.max(np.abs(oc[:, 4+i])), 1e-30))
            for i in range(3)
        ]
        q_peak_drift = np.max(np.abs(qc[:, 4:7]), axis=0)
        o_peak_drift = np.max(np.abs(oc[:, 4:7]), axis=0)
        peak_drift_rel = max_relative(q_peak_drift, o_peak_drift)
        q_resid = qc[-1, 4:7]
        o_resid = oc[-1, 4:7]
        residual_drift_abs_diff = np.abs(q_resid - o_resid)
        q_peak_accel = np.max(np.abs(qc[:, 7:10]), axis=0)
        o_peak_accel = np.max(np.abs(oc[:, 7:10]), axis=0)
        peak_accel_rel = max_relative(q_peak_accel, o_peak_accel)
        accel_rmse = [
            float(100.0 * np.sqrt(np.mean((qc[:, 7+i] - oc[:, 7+i]) ** 2)) /
                  max(np.max(np.abs(qc[:, 7+i])), np.max(np.abs(oc[:, 7+i])), 1e-30))
            for i in range(3)
        ]
    else:
        max_u = max_d = peak_drift_rel = peak_accel_rel = float("inf")
        drift_rmse = accel_rmse = [float("inf")] * 3
        q_peak_drift = o_peak_drift = np.asarray([float("inf")] * 3)
        q_resid = o_resid = residual_drift_abs_diff = np.asarray([float("inf")] * 3)
        q_peak_accel = o_peak_accel = np.asarray([float("inf")] * 3)

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

    q_periods = [float(x["period_s"]) for x in q["modes"][:3]]
    period_rel = max_relative(q_periods, o_periods)
    q_status = classify_quake(qrun, requested)
    o_status = "completed" if len(oh) == requested else "numerical_noncompletion"

    checks = {
        "quakecore_completed_full_record": q_status == "completed",
        "minimum_common_converged_steps": n >= 1000,
        "period_relative_difference_le_1e-6": period_rel <= 1e-6,
        "peak_story_drift_relative_difference_le_0.1_percent": peak_drift_rel <= 1e-3,
        "story_drift_rmse_le_0.1_percent_of_peak": max(drift_rmse) <= 0.1,
        "wall_peak_edp_relative_difference_le_1_percent": max_wall_peak_rel <= 0.01,
    }
    passed = all(checks.values())
    q_stats = qrun.get("stats", {})
    result = {
        "record_id": record_id,
        "pair_id": meta["pair_id"],
        "component_index": meta["component_index"],
        "peer_header": meta["peer_header"],
        "source_file": meta["original_file"],
        "source_npts": meta["npts"],
        "source_dt_s": meta["dt_s"],
        "source_pga_g": meta["pga_g"],
        "dimensionless_acceleration_scale": scale,
        "analysis_peak_input_acceleration_in_per_s2": max(abs(v) for v in rec["acceleration"]),
        "quakecore_process_returncode": proc.returncode,
        "quakecore_status": q_status,
        "quakecore_steps": len(qh),
        "quakecore_termination": qrun.get("termination"),
        "quakecore_reason": qrun.get("reason"),
        "quakecore_integration_seconds": q_stats.get("elapsed_seconds"),
        "quakecore_process_wall_seconds": q_process_wall,
        "opensees_status": o_status,
        "opensees_analysis_code": os_code,
        "opensees_steps": len(oh),
        "opensees_end_time_s": float(ops.getTime()),
        "opensees_integration_seconds": os_elapsed,
        "opensees_iterations": os_iterations,
        "common_steps": n,
        "common_end_time_s": float(oh[n-1][0]) if n else None,
        "quakecore_periods_seconds": q_periods,
        "opensees_periods_seconds": o_periods,
        "period_relative_difference": period_rel,
        "quakecore_peak_story_drift_common": q_peak_drift.tolist(),
        "opensees_peak_story_drift_common": o_peak_drift.tolist(),
        "peak_story_drift_relative_difference": peak_drift_rel,
        "max_floor_displacement_difference_in": max_u,
        "max_story_drift_difference": max_d,
        "story_drift_rmse_percent_of_peak": drift_rmse,
        "quakecore_residual_story_drift_common": q_resid.tolist(),
        "opensees_residual_story_drift_common": o_resid.tolist(),
        "residual_story_drift_absolute_difference": residual_drift_abs_diff.tolist(),
        "quakecore_peak_floor_absolute_acceleration_common": q_peak_accel.tolist(),
        "opensees_peak_floor_absolute_acceleration_common": o_peak_accel.tolist(),
        "peak_floor_absolute_acceleration_relative_difference": peak_accel_rel,
        "floor_absolute_acceleration_rmse_percent_of_peak": accel_rmse,
        "floor_acceleration_is_diagnostic_not_gate": True,
        "max_wall_peak_edp_relative_difference": max_wall_peak_rel,
        "max_abs_wall_force_history_error": max_wall_hist_abs,
        "walls": per_wall,
        "frozen_gate4_checks": checks,
        "passed_frozen_gate4_metrics": passed,
    }
    (record_dir / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-job", type=Path, required=True)
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--quake-exe", type=Path, required=True)
    ap.add_argument("--record-ids", nargs="+", required=True)
    ap.add_argument("--scale", type=float, default=0.03)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--assert-engineering-parity", action="store_true")
    args = ap.parse_args()
    if not math.isfinite(args.scale) or args.scale <= 0.0:
        raise SystemExit("scale must be finite and positive")

    base_job = json.loads(args.base_job.read_text())
    manifest = json.loads(args.manifest.read_text())
    known = {str(r["record_id"]) for r in manifest["records"]}
    unknown = [rid for rid in args.record_ids if rid not in known]
    if unknown:
        raise SystemExit(f"unknown record IDs: {unknown}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for record_id in args.record_ids:
        print(f"=== Phase 9N record {record_id} @ scale {args.scale:g} ===", flush=True)
        results.append(run_record(
            base_job,
            args.source_archive,
            manifest,
            record_id,
            args.scale,
            args.quake_exe.resolve(),
            args.out_dir,
        ))

    valid = [r for r in results if not r.get("quakecore_output_missing")]
    summary = {
        "scope": "Phase 9N multi-record exact-law external software/model validation of the documented YORi 10001 reconstructed surrogate",
        "claim_boundary": "not source-model identity, physical validation, FEMA P-695 collapse qualification/acceptance, code approval, or R-factor recommendation",
        "opensees_version": ops.version(),
        "dimensionless_acceleration_scale": args.scale,
        "record_ids": args.record_ids,
        "frozen_gate4_thresholds": {
            "minimum_common_steps": 1000,
            "maximum_period_relative_difference": 1e-6,
            "maximum_peak_story_drift_relative_difference": 1e-3,
            "maximum_story_drift_rmse_percent_of_peak": 0.1,
            "maximum_wall_peak_edp_relative_difference": 0.01,
        },
        "new_diagnostics_are_not_acceptance_gates": [
            "residual_story_drift_absolute_difference",
            "peak_floor_absolute_acceleration_relative_difference",
            "floor_absolute_acceleration_rmse_percent_of_peak",
        ],
        "results": results,
        "suite_summary": {
            "records_requested": len(args.record_ids),
            "records_with_quakecore_output": len(valid),
            "quakecore_completed_records": [r["record_id"] for r in valid if r["quakecore_status"] == "completed"],
            "opensees_completed_records": [r["record_id"] for r in valid if r["opensees_status"] == "completed"],
            "records_passing_frozen_gate4_metrics": [r["record_id"] for r in valid if r["passed_frozen_gate4_metrics"]],
            "minimum_common_steps": min((r["common_steps"] for r in valid), default=0),
            "maximum_period_relative_difference": max((r["period_relative_difference"] for r in valid), default=float("inf")),
            "maximum_peak_story_drift_relative_difference": max((r["peak_story_drift_relative_difference"] for r in valid), default=float("inf")),
            "maximum_story_drift_rmse_percent_of_peak": max((max(r["story_drift_rmse_percent_of_peak"]) for r in valid), default=float("inf")),
            "maximum_wall_peak_edp_relative_difference": max((r["max_wall_peak_edp_relative_difference"] for r in valid), default=float("inf")),
            "maximum_peak_floor_acceleration_relative_difference_diagnostic": max((r["peak_floor_absolute_acceleration_relative_difference"] for r in valid), default=float("inf")),
        },
    }
    summary["passed_frozen_gate4_metrics_for_all_records"] = (
        len(valid) == len(args.record_ids) and all(r["passed_frozen_gate4_metrics"] for r in valid)
    )
    (args.out_dir / "suite_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["suite_summary"], indent=2))

    if args.assert_engineering_parity and not summary["passed_frozen_gate4_metrics_for_all_records"]:
        raise SystemExit("Phase 9N pilot did not pass frozen Gate 4 engineering metrics for all requested records")


if __name__ == "__main__":
    main()
