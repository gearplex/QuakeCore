#!/usr/bin/env python3
"""Independent OpenSeesPy 3.8.0 comparator for the reconstructed exact-law job."""
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


class Tags:
    def __init__(self):
        self.value = 1000

    def next(self):
        self.value += 1
        return self.value


def material(spec, tags):
    kind = spec["type"]
    tag = tags.next()
    if kind == "concrete_cm":
        ops.uniaxialMaterial(
            "ConcreteCM", tag,
            spec["fc"], spec["epsc"], spec["Ec"], spec["rc"], spec["xcrn"],
            spec["ft"], spec["et"], spec["rt"], spec["xcrp"],
            "-GapClose", int(bool(spec.get("gap_close", False))),
        )
    elif kind == "elastic":
        ops.uniaxialMaterial("Elastic", tag, spec["E"])
    elif kind == "pinching4":
        args = []
        for deformation, force in spec["positive"]:
            args.extend([force, deformation])
        for deformation, force in spec["negative"]:
            args.extend([force, deformation])
        args.extend([
            spec["r_disp_positive"], spec["r_force_positive"], spec["u_force_positive"],
            spec["r_disp_negative"], spec["r_force_negative"], spec["u_force_negative"],
        ])
        args.extend([*spec["gamma_k"], spec["gamma_k_limit"]])
        args.extend([*spec["gamma_d"], spec["gamma_d_limit"]])
        args.extend([*spec["gamma_f"], spec["gamma_f_limit"]])
        args.extend([spec["gamma_e"], spec["damage_mode"]])
        ops.uniaxialMaterial("Pinching4", tag, *args)
    elif kind == "minmax":
        child = material(spec["material"], tags)
        ops.uniaxialMaterial("MinMax", tag, child, "-min", spec["min"], "-max", spec["max"])
    elif kind == "parallel":
        children = [material(child, tags) for child in spec["materials"]]
        ops.uniaxialMaterial("Parallel", tag, *children)
    else:
        raise ValueError(f"unsupported material: {kind}")
    return tag


def build(job):
    m = job["model"]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for n in m["nodes"]:
        ops.node(n[0], n[1], n[2])
        ops.mass(n[0], *n[3:])
    for f in m["fixities"]:
        ops.fix(f[0], *[int(x) for x in f[1:]])
    for c in m["equal_dofs"]:
        dof = {"UX": 1, "UY": 2, "RZ": 3}[c[2]]
        ops.equalDOF(c[0], c[1], dof)

    ops.geomTransf("Linear", 1)
    ops.geomTransf("PDelta", 2)
    for e in m["members"]:
        transf = 2 if len(e) >= 8 and bool(e[7]) else 1
        ops.element("elasticBeamColumn", e[0], e[1], e[2], e[4], e[3], e[5], transf)

    tags = Tags()
    for w in m["walls"]:
        fs = w["fibers"]
        conc = [material(x["concrete"], tags) for x in fs]
        steel = [material(x["steel"], tags) for x in fs]
        shear = material(w["shear"], tags)
        ops.element(
            "MVLEM", w["id"], w["density"], w["i"], w["j"], len(fs), w["c"],
            "-thick", *[x["thickness"] for x in fs],
            "-width", *[x["width"] for x in fs],
            "-rho", *[x["rho"] for x in fs],
            "-matConcrete", *conc,
            "-matSteel", *steel,
            "-matShear", shear,
        )


def gravity(job):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, fx, fy, mz in job["gravity"]["loads"]:
        ops.load(n, fx, fy, mz)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.01)
    ops.analysis("Static")
    started = time.perf_counter()
    ok = ops.analyze(100)
    elapsed = time.perf_counter() - started
    ops.reactions()
    return {
        "ok": ok == 0,
        "elapsed_seconds": elapsed,
        "wall_base_reaction": [ops.nodeReaction(1, i) for i in (1, 2, 3)],
        "lean_base_reaction": [ops.nodeReaction(100, i) for i in (1, 2, 3)],
        "story_wall_vertical_displacement": [ops.nodeDisp(n, 2) for n in job["model"]["story_nodes"]],
        "story_lean_vertical_displacement": [ops.nodeDisp(101 + i, 2) for i in range(3)],
    }


def run(job):
    build(job)
    grav = gravity(job)
    lam = np.asarray(ops.eigen(3), dtype=float)
    periods = (2 * np.pi / np.sqrt(lam)).tolist()

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

    nodes = job["model"]["story_nodes"]
    heights = np.asarray([ops.nodeCoord(n, 2) for n in nodes])
    history = []
    peak = np.zeros(3)
    iterations = 0
    started = time.perf_counter()
    ok = 0
    for _ in values:
        ok = ops.analyze(1, rec["dt"])
        if ok:
            break
        u = np.asarray([ops.nodeDisp(n, 1) for n in nodes])
        d = np.diff(np.r_[0.0, u]) / np.diff(np.r_[0.0, heights])
        peak = np.maximum(peak, abs(d))
        iterations += ops.testIter()
        history.append([ops.getTime(), *u.tolist(), *d.tolist()])
    elapsed = time.perf_counter() - started
    return {
        "opensees_version": ops.version(),
        "gravity": grav,
        "periods_seconds": periods,
        "analysis_ok": ok == 0,
        "steps_completed": len(history),
        "steps_requested": len(values),
        "elapsed_seconds": elapsed,
        "iterations": iterations,
        "peak_story_drift_ratio": peak.tolist(),
        "history": history,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--assert-parity", action="store_true")
    args = ap.parse_args()

    job = json.loads(args.job.read_text())
    q = json.loads(args.quake.read_text())
    o = run(job)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "opensees_reconstructed_exact.json").write_text(json.dumps(o, indent=2) + "\n")

    qrun = q["runs"][0]
    qh = np.asarray([
        [x["time_s"], *x["floor_displacement"], *x["story_drift_ratio"]]
        for x in qrun["history"]
    ], dtype=float)
    oh = np.asarray(o["history"], dtype=float)
    n = min(len(qh), len(oh))
    if n == 0:
        max_u = float("inf")
        max_d = float("inf")
        rmse = [float("inf")] * 3
    else:
        qhc, ohc = qh[:n], oh[:n]
        max_u = float(np.max(abs(qhc[:, 1:4] - ohc[:, 1:4])))
        max_d = float(np.max(abs(qhc[:, 4:7] - ohc[:, 4:7])))
        rmse = [
            float(100 * np.sqrt(np.mean((qhc[:, 4 + i] - ohc[:, 4 + i]) ** 2)) /
                  max(np.max(abs(qhc[:, 4 + i])), np.max(abs(ohc[:, 4 + i])), 1e-30))
            for i in range(3)
        ]

    qperiods = [x["period_s"] for x in q["modes"]]
    period_diff = [abs(a - b) for a, b in zip(qperiods, o["periods_seconds"])]
    peak_diff = [abs(a - b) for a, b in zip(qrun["peak_story_drift_ratio"], o["peak_story_drift_ratio"])]
    summary = {
        "scope": "reconstructed 3-story ConcreteCM + Pinching4 + MinMax + Parallel surrogate",
        "claim_boundary": "parity of reconstructed surrogate only; not identity to lost received source or FEMA P-695 acceptance",
        "opensees_version": o["opensees_version"],
        "quakecore_gravity_residual": q["gravity"]["residual_inf_norm"],
        "opensees_gravity_ok": o["gravity"]["ok"],
        "quakecore_periods_seconds": qperiods,
        "opensees_periods_seconds": o["periods_seconds"],
        "period_absolute_difference_seconds": period_diff,
        "quakecore_peak_story_drift_ratio": qrun["peak_story_drift_ratio"],
        "opensees_peak_story_drift_ratio": o["peak_story_drift_ratio"],
        "peak_story_drift_absolute_difference": peak_diff,
        "history_steps_compared": n,
        "quakecore_history_steps": len(qh),
        "opensees_history_steps": len(oh),
        "max_floor_displacement_difference_in": max_u,
        "max_story_drift_difference": max_d,
        "story_drift_rmse_percent_of_peak": rmse,
        "quakecore_integration_seconds": qrun["stats"]["elapsed_seconds"],
        "opensees_integration_seconds": o["elapsed_seconds"],
    }
    (args.out / "reconstructed_exact_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

    if args.assert_parity:
        if not str(o["opensees_version"]).startswith("3.8.0"):
            raise SystemExit("OpenSees version is not 3.8.0")
        if not o["gravity"]["ok"] or not o["analysis_ok"]:
            raise SystemExit("OpenSees reconstructed analysis did not complete")
        if len(qh) != len(oh):
            raise SystemExit(f"history length mismatch: QuakeCore={len(qh)} OpenSees={len(oh)}")
        if max(period_diff) > 1e-7:
            raise SystemExit(f"period parity failed: {period_diff}")
        if max_u > 1e-6:
            raise SystemExit(f"floor displacement parity failed: {max_u}")
        if max_d > 1e-8:
            raise SystemExit(f"story drift parity failed: {max_d}")


if __name__ == "__main__":
    main()
