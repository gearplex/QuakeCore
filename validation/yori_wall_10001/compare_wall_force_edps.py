#!/usr/bin/env python3
"""Compare QuakeCore and OpenSees MVLEM wall force EDPs on a common converged prefix.

This is an engineering-demand comparator for the reconstructed surrogate.  It
uses identical model construction, gravity state, Newmark parameters, and the
same ground-motion prefix.  Absolute peak force/moment EDPs are emphasized so
opposite element-end sign conventions cannot masquerade as an engineering
response difference; signed force-vector differences are retained as a
convention diagnostic.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import numpy as np
import openseespy.opensees as ops

from compare_reconstructed_exact_10001 import build, gravity


def setup_transient(job, values):
    build(job)
    grav = gravity(job)
    if not grav["ok"]:
        raise RuntimeError("OpenSees gravity analysis did not converge")
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    rec = job["records"][0]
    ops.timeSeries("Path", 10, "-dt", rec["dt"], "-values", 0.0, *values, "-useLast")
    ops.pattern("UniformExcitation", 2, 1, "-accel", 10)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-8, 100, 0, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def qwall(row, wall_id):
    walls = row["walls"]
    return walls[str(wall_id)] if str(wall_id) in walls else walls[wall_id]


def relerr(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    job = json.loads(args.job.read_text())
    q = json.loads(args.quake.read_text())
    qh = q["runs"][0]["history"]
    rec = job["records"][0]
    n_requested = min(args.steps, len(qh), len(rec["acceleration"]))
    values = rec["acceleration"][:n_requested]
    setup_transient(job, values)

    wall_ids = [int(w["id"]) for w in job["model"]["walls"]]
    qforces = {wid: [] for wid in wall_ids}
    oforces = {wid: [] for wid in wall_ids}
    completed = 0
    for i in range(n_requested):
        ok = ops.analyze(1, rec["dt"])
        if ok:
            break
        for wid in wall_ids:
            qf = np.asarray(qwall(qh[i], wid)["nodal_force"], dtype=float)
            of = np.asarray(ops.eleForce(wid), dtype=float)
            if qf.shape != (6,) or of.shape != (6,):
                raise RuntimeError(f"wall {wid} force vector is not length 6: Quake={qf.shape} OpenSees={of.shape}")
            qforces[wid].append(qf)
            oforces[wid].append(of)
        completed += 1

    if completed == 0:
        raise RuntimeError("OpenSees completed no comparison steps")

    per_wall = {}
    max_peak_rel = 0.0
    max_abs_magnitude_history_error = 0.0
    edp_indices = {"shear": 3, "axial_force": 4, "bottom_moment": 2}
    for wid in wall_ids:
        qf = np.asarray(qforces[wid])
        of = np.asarray(oforces[wid])
        same = np.abs(qf - of)
        flipped = np.abs(qf + of)
        mag = np.abs(np.abs(qf) - np.abs(of))
        item = {
            "signed_force_convention": "same" if float(np.max(same)) <= float(np.max(flipped)) else "opposite",
            "max_abs_signed_force_error_same_convention": float(np.max(same)),
            "max_abs_signed_force_error_opposite_convention": float(np.max(flipped)),
            "max_abs_force_magnitude_history_error": float(np.max(mag)),
            "edps": {},
        }
        max_abs_magnitude_history_error = max(max_abs_magnitude_history_error, item["max_abs_force_magnitude_history_error"])
        for name, j in edp_indices.items():
            qp = float(np.max(np.abs(qf[:, j])))
            op = float(np.max(np.abs(of[:, j])))
            e = relerr(qp, op)
            max_peak_rel = max(max_peak_rel, e)
            item["edps"][name] = {
                "quakecore_peak_abs": qp,
                "opensees_peak_abs": op,
                "relative_difference": e,
            }
        per_wall[str(wid)] = item

    summary = {
        "scope": "reconstructed 3-story surrogate, common OpenSees-converged prefix; engineering force EDP diagnostic only",
        "claim_boundary": "does not establish identity to the lost received source, physical validation, FEMA P-695 acceptance, or R-factor qualification",
        "opensees_version": ops.version(),
        "steps_requested": n_requested,
        "steps_completed": completed,
        "end_time_s": float(ops.getTime()),
        "edp_definition": {
            "shear": "absolute peak of wall nodal force component 3 (top horizontal resisting force in QuakeCore convention)",
            "axial_force": "absolute peak of wall nodal force component 4 (top vertical resisting force in QuakeCore convention)",
            "bottom_moment": "absolute peak of wall nodal force component 2 (bottom nodal moment in QuakeCore convention)",
        },
        "max_relative_peak_edp_difference": max_peak_rel,
        "max_abs_force_magnitude_history_error": max_abs_magnitude_history_error,
        "walls": per_wall,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
