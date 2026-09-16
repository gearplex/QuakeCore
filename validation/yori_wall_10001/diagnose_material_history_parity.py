#!/usr/bin/env python3
"""Replay QuakeCore wall-fiber histories through OpenSeesPy uniaxial materials."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import openseespy.opensees as ops

from compare_reconstructed_exact_10001 import Tags, material


def replay(spec, strains, seed_steps=100):
    ops.wipe()
    tags = Tags()
    tag = material(spec, tags)
    ops.testUniaxialMaterial(tag)
    if strains:
        first = strains[0]
        for i in range(1, seed_steps + 1):
            ops.setStrain(first * i / seed_steps)
    out = []
    for e in strains:
        ops.setStrain(e)
        out.append(float(ops.getStress()))
    return out


def first_over(errors, threshold):
    for i, e in enumerate(errors):
        if e > threshold:
            return i
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
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
            ostress = replay(spec["fibers"][fi]["concrete"], strains)
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
            osteel = replay(spec["fibers"][fi]["steel"], strains)
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
    report = {
        "note": "OpenSees material replay is seeded monotonically from zero to the first recorded dynamic strain; use onset timing more strongly than absolute initial offset.",
        "earliest_divergences": findings[:20],
        "largest_divergences": sorted(findings, key=lambda x: x["max_abs_stress_error"], reverse=True)[:20],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
