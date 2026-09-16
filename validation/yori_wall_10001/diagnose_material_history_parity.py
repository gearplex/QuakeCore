#!/usr/bin/env python3
"""Replay QuakeCore wall-fiber histories through OpenSeesPy uniaxial materials."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import openseespy.opensees as ops

from compare_reconstructed_exact_10001 import Tags, material


def replay_response(spec, strains, seed_steps=100):
    ops.wipe()
    tags = Tags()
    tag = material(spec, tags)
    ops.testUniaxialMaterial(tag)
    if strains:
        first = strains[0]
        for i in range(1, seed_steps + 1):
            ops.setStrain(first * i / seed_steps)
    stress = []
    tangent = []
    for e in strains:
        ops.setStrain(e)
        stress.append(float(ops.getStress()))
        tangent.append(float(ops.getTangent()))
    return stress, tangent


def first_over(errors, threshold):
    for i, e in enumerate(errors):
        if e > threshold:
            return i
    return None


def trace_tangent_report(trace_path, spec, history):
    if not trace_path.exists():
        return {"available": False, "path": str(trace_path)}

    rows = list(csv.DictReader(trace_path.open()))
    strains = [float(r["strain"]) for r in rows]
    _, otangent = replay_response(spec, strains)
    qtangent = [float(r["tangent"]) for r in rows]
    errors = [abs(a - b) for a, b in zip(qtangent, otangent)]
    imax = max(range(len(errors)), key=errors.__getitem__)

    thresholds = [1e-8, 1e-6, 1e-4, 1e-2, 1.0]
    onset = {}
    for threshold in thresholds:
        idx = first_over(errors, threshold)
        onset[f"first_gt_{threshold:g}"] = None if idx is None else {
            "index": idx,
            "time_s": history[idx]["time_s"] if idx < len(history) else None,
            "strain": strains[idx],
            "quake_tangent": qtangent[idx],
            "opensees_tangent": otangent[idx],
            "abs_error": errors[idx],
            "quake_rule": int(rows[idx]["rule"]),
            "quake_increment": float(rows[idx]["increment"]),
            "quake_landmarks": {
                "eunn": float(rows[idx]["eunn"]),
                "funn": float(rows[idx]["funn"]),
                "espln": float(rows[idx]["espln"]),
                "Epln": float(rows[idx]["Epln"]),
                "Te0": float(rows[idx]["Te0"]),
                "Teunp": float(rows[idx]["Teunp"]),
                "Tfunp": float(rows[idx]["Tfunp"]),
                "esplp": float(rows[idx]["esplp"]),
                "Eplp": float(rows[idx]["Eplp"]),
                "fnewn": float(rows[idx]["fnewn"]),
                "Enewn": float(rows[idx]["Enewn"]),
                "esren": float(rows[idx]["esren"]),
                "origin12": float(rows[idx]["origin12"]),
                "Tea": float(rows[idx]["Tea"]),
                "Teb": float(rows[idx]["Teb"]),
            },
        }

    return {
        "available": True,
        "path": str(trace_path),
        "wall": 5,
        "fiber": 0,
        "max_abs_tangent_error": errors[imax],
        "max_index": imax,
        "max_time_s": history[imax]["time_s"] if imax < len(history) else None,
        "quake_tangent_at_max": qtangent[imax],
        "opensees_tangent_at_max": otangent[imax],
        "onset": onset,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--trace",
        type=Path,
        default=Path("reconstructed_comparison/concrete_cm_state_trace.csv"),
        help="optional QuakeCore ConcreteCM state trace for wall 5 / fiber 0 tangent replay",
    )
    args = ap.parse_args()

    job = json.loads(args.job.read_text())
    q = json.loads(args.quake.read_text())
    history = q["runs"][0]["history"]
    walls = {str(w["id"]): w for w in job["model"]["walls"]}
    findings = []

    for wall_id, spec in walls.items():
        if not history or wall_id not in history[0]["walls"]:
            continue
        nf = len(spec["fibers"])
        for fi in range(nf):
            strains = [h["walls"][wall_id]["fiber_strain"][fi] for h in history]
            qstress = [h["walls"][wall_id]["concrete_stress"][fi] for h in history]
            ostress, _ = replay_response(spec["fibers"][fi]["concrete"], strains)
            err = [abs(a-b) for a,b in zip(qstress, ostress)]
            imax = max(range(len(err)), key=err.__getitem__)
            i1 = first_over(err, 1e-6)
            findings.append({
                "kind": "concrete", "wall": int(wall_id), "fiber": fi,
                "max_abs_stress_error": err[imax], "max_index": imax,
                "max_time_s": history[imax]["time_s"],
                "first_gt_1e-6_index": i1,
                "first_gt_1e-6_time_s": None if i1 is None else history[i1]["time_s"],
                "quake_at_first": None if i1 is None else qstress[i1],
                "opensees_at_first": None if i1 is None else ostress[i1],
            })

            qsteel = [h["walls"][wall_id]["steel_stress"][fi] for h in history]
            osteel, _ = replay_response(spec["fibers"][fi]["steel"], strains)
            err = [abs(a-b) for a,b in zip(qsteel, osteel)]
            imax = max(range(len(err)), key=err.__getitem__)
            i1 = first_over(err, 1e-6)
            findings.append({
                "kind": "steel", "wall": int(wall_id), "fiber": fi,
                "max_abs_stress_error": err[imax], "max_index": imax,
                "max_time_s": history[imax]["time_s"],
                "first_gt_1e-6_index": i1,
                "first_gt_1e-6_time_s": None if i1 is None else history[i1]["time_s"],
                "quake_at_first": None if i1 is None else qsteel[i1],
                "opensees_at_first": None if i1 is None else osteel[i1],
            })

    findings.sort(key=lambda x: (float("inf") if x["first_gt_1e-6_index"] is None else x["first_gt_1e-6_index"], -x["max_abs_stress_error"]))
    trace_report = trace_tangent_report(
        args.trace,
        walls["5"]["fibers"][0]["concrete"],
        history,
    ) if "5" in walls and walls["5"].get("fibers") else {"available": False}
    report = {
        "note": "OpenSees material replay is seeded monotonically from zero to the first recorded dynamic strain; use onset timing more strongly than absolute initial offset.",
        "trace_tangent_replay": trace_report,
        "earliest_divergences": findings[:20],
        "largest_divergences": sorted(findings, key=lambda x: x["max_abs_stress_error"], reverse=True)[:20],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
