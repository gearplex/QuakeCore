#!/usr/bin/env python3
"""OpenSees comparator for the matched-law QuakeCore gravity/P-Delta job."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
import numpy as np
import openseespy.opensees as ops


class Tags:
    def __init__(self): self.value = 1000
    def next(self): self.value += 1; return self.value


def material(spec, tags):
    tag = tags.next()
    if spec["type"] == "concrete01":
        ops.uniaxialMaterial("Concrete01", tag, spec["fc"], spec["epsc"], spec["fcu"], spec["epsu"])
    elif spec["type"] == "steel_bilinear":
        ops.uniaxialMaterial("Steel01", tag, spec["fy"], spec["E"], spec["b"])
    else:
        raise ValueError(spec)
    return tag


def build(job):
    m = job["model"]
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
    for n in m["nodes"]:
        ops.node(n[0], n[1], n[2]); ops.mass(n[0], *n[3:])
    for f in m["fixities"]: ops.fix(f[0], *[int(x) for x in f[1:]])
    for c in m["equal_dofs"]:
        dof = {"UX": 1, "UY": 2, "RZ": 3}[c[2]]
        ops.equalDOF(c[0], c[1], dof)
    ops.geomTransf("PDelta", 1)
    for e in m["members"]:
        ops.element("elasticBeamColumn", e[0], e[1], e[2], e[4], e[3], e[5], 1)
    tags = Tags()
    for w in m["walls"]:
        fs = w["fibers"]
        conc = [material(x["concrete"], tags) for x in fs]
        steel = [material(x["steel"], tags) for x in fs]
        sh = w["shear"]; shear = tags.next()
        ops.uniaxialMaterial("Steel01", shear, sh["yield_force"], sh["stiffness"], sh["hardening_ratio"])
        ops.element("MVLEM", w["id"], w["density"], w["i"], w["j"], len(fs), w["c"],
                    "-thick", *[x["thickness"] for x in fs],
                    "-width", *[x["width"] for x in fs],
                    "-rho", *[x["rho"] for x in fs],
                    "-matConcrete", *conc, "-matSteel", *steel, "-matShear", shear)


def gravity(job):
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n, fx, fy, mz in job["gravity"]["loads"]: ops.load(n, fx, fy, mz)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-8, 100, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.01); ops.analysis("Static")
    started = time.perf_counter(); ok = ops.analyze(100); elapsed = time.perf_counter() - started
    ops.reactions()
    return {
        "ok": ok == 0, "elapsed_seconds": elapsed,
        "wall_base_reaction": [ops.nodeReaction(1, i) for i in (1, 2, 3)],
        "lean_base_reaction": [ops.nodeReaction(100, i) for i in (1, 2, 3)],
        "story_wall_vertical_displacement": [ops.nodeDisp(n, 2) for n in job["model"]["story_nodes"]],
        "story_lean_vertical_displacement": [ops.nodeDisp(101 + i, 2) for i in range(3)],
    }


def run(job):
    build(job); grav = gravity(job)
    lam = np.asarray(ops.eigen(3), dtype=float)
    periods = (2 * np.pi / np.sqrt(lam)).tolist()
    ops.loadConst("-time", 0.0); ops.wipeAnalysis()
    rec = job["records"][0]; values = rec["acceleration"]
    ops.timeSeries("Path", 10, "-dt", rec["dt"], "-values", 0.0, *values, "-useLast")
    ops.pattern("UniformExcitation", 2, 1, "-accel", 10)
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-8, 100, 0, 0); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    nodes = job["model"]["story_nodes"]
    heights = np.asarray([ops.nodeCoord(n, 2) for n in nodes])
    history, peak, iterations = [], np.zeros(3), 0
    started = time.perf_counter(); ok = 0
    for _ in values:
        ok = ops.analyze(1, rec["dt"])
        if ok: break
        u = np.asarray([ops.nodeDisp(n, 1) for n in nodes])
        d = np.diff(np.r_[0.0, u]) / np.diff(np.r_[0.0, heights])
        peak = np.maximum(peak, abs(d)); iterations += ops.testIter()
        history.append([ops.getTime(), *u.tolist(), *d.tolist()])
    elapsed = time.perf_counter() - started
    return {"opensees_version": ops.version(), "gravity": grav, "periods_seconds": periods,
            "analysis_ok": ok == 0, "elapsed_seconds": elapsed, "iterations": iterations,
            "peak_story_drift_ratio": peak.tolist(), "history": history}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    job = json.loads(args.job.read_text()); q = json.loads(args.quake.read_text())
    o = run(job); args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "opensees_matched_gravity.json").write_text(json.dumps(o, indent=2))
    qrun = q["runs"][0]; qh = np.asarray([[x["time_s"], *x["floor_displacement"], *x["story_drift_ratio"]] for x in qrun["history"]])
    oh = np.asarray(o["history"]); n = min(len(qh), len(oh)); qh, oh = qh[:n], oh[:n]
    summary = {
        "scope": "matched Concrete01/Steel01/bilinear shear with gravity state and leaning-column P-Delta",
        "quakecore_gravity_residual": q["gravity"]["residual_inf_norm"],
        "opensees_gravity_reactions": o["gravity"],
        "quakecore_periods_seconds": [x["period_s"] for x in q["modes"]],
        "opensees_periods_seconds": o["periods_seconds"],
        "quakecore_peak_story_drift_ratio": qrun["peak_story_drift_ratio"],
        "opensees_peak_story_drift_ratio": o["peak_story_drift_ratio"],
        "max_floor_displacement_difference_in": float(np.max(abs(qh[:, 1:4] - oh[:, 1:4]))),
        "max_story_drift_difference": float(np.max(abs(qh[:, 4:7] - oh[:, 4:7]))),
        "story_drift_rmse_percent_of_peak": [float(100*np.sqrt(np.mean((qh[:,4+i]-oh[:,4+i])**2))/max(np.max(abs(qh[:,4+i])),np.max(abs(oh[:,4+i])))) for i in range(3)],
        "quakecore_integration_seconds": qrun["stats"]["elapsed_seconds"],
        "opensees_integration_seconds": o["elapsed_seconds"],
    }
    (args.out / "matched_gravity_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__": main()
