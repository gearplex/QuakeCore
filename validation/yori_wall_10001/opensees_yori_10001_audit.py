#!/usr/bin/env python3
"""Auditable OpenSeesPy baseline for YORi wall archetype 10001.

Recreates the received Tcl model (ConcreteCM + Pinching4 + MinMax + Parallel,
MVLEM wall, gravity preload, and leaning-column P-Delta) directly in
OpenSeesPy so it can run on Linux.  It compares the received and corrected
steel backbones while keeping all other model assumptions fixed.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import time

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import numpy as np
import openseespy.opensees as ops

G = 386.09


def load_inputs(folder: Path):
    raw_g = np.loadtxt(folder / "AllStoriesOpenseesGeomINPUT.txt")
    nstory = int(raw_g[0])
    return raw_g.reshape(nstory, 9), np.loadtxt(
        folder / "AllStoriesOpenseesModelINPUT.txt"
    ).reshape(nstory, 41)


def pinching4(tag: int, p, d, np_, nd):
    args = []
    for f, x in zip(p, d):
        args += [float(f), float(x)]
    for f, x in zip(np_, nd):
        args += [float(f), float(x)]
    args += [0.6, 0.99, 0.4, 0.6, 0.99, 0.4]
    args += [0.0, 0.0, 0.0, 0.0, 0.0]
    args += [0.1, 0.0, 0.0, 0.0, 2.0]
    args += [0.0, 0.0, 0.0, 0.0, 0.0, 10000.0, "energy"]
    ops.uniaxialMaterial("Pinching4", tag, *args)


def fiber_layout(g, m):
    nele, lw, tw, lb, rho_b, rho_w = int(g[2]), *map(float, g[4:9])
    nbar = int(m[40])
    if nele > 1:
        nb = int(min(max(nbar, 2), 12))
        nw = int(min(max((lw - 2.0 * lb) / 12.0, 2), 12))
    else:
        nb = int(min(max(nbar / 2.0, 2), 12))
        nw = int(min(max((lw - 2.0 * lb) / 48.0, 2), 12))
    thick, width, rho, is_boundary = [], [], [], []
    for k in range(2 * nb + nw):
        boundary = k < nb or k >= nb + nw
        thick.append(tw)
        width.append(lb / nb if boundary else (lw - 2.0 * lb) / nw)
        rho.append(rho_b if boundary else rho_w)
        is_boundary.append(boundary)
    return thick, width, rho, is_boundary


def floor_mass(story, nstory, lw, tw, story_h):
    area = 102.0 * 102.0 / 2.0
    if story == nstory:
        slab = area * (1.05 * 0.125 + 0.25 * 0.020) / G
        self_mass = 1.05 * (lw / 12.0) * (story_h / 12.0) * 0.150 / 2.0 / G
    else:
        slab = area * (1.05 * 0.125 + 0.25 * 0.040) / G
        self_mass = 1.05 * (lw / 12.0) * (story_h / 12.0) * 0.150 / G
    return slab + self_mass * tw / 12.0


def build_model(folder: Path):
    gg, mm = load_inputs(folder)
    nstory = len(gg)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    wall_nodes = [100000]
    lean_nodes = [200000]
    ops.node(100000, 0.0, 0.0)
    ops.node(200000, float(gg[0, 4]), 0.0)
    ops.fix(100000, 1, 1, 1)
    ops.fix(200000, 1, 1, 1)
    height = 0.0
    story_wall_nodes, story_lean_nodes = [], []
    for s, g in enumerate(gg, start=1):
        for j in range(1, int(g[2]) + 1):
            height += float(g[3])
            tag = 100000 + 100 * s + j
            ops.node(tag, 0.0, height)
            wall_nodes.append(tag)
        ltag = 200000 + 100 * s + int(g[2])
        ops.node(ltag, float(g[4]), height)
        lean_nodes.append(ltag)
        story_wall_nodes.append(wall_nodes[-1])
        story_lean_nodes.append(ltag)
        mass = floor_mass(s, nstory, float(g[4]), float(g[5]), float(g[1]))
        ops.mass(wall_nodes[-1], mass, mass, 0.0)

    ops.uniaxialMaterial("Elastic", 100001, 1.0e5)
    ops.uniaxialMaterial("Elastic", 100002, 1.0e-2)
    mats = []
    for s, m in enumerate(mm, start=1):
        uncon, conf = 201000 + s, 301000 + s
        for tag, off, xcrn in ((uncon, 0, 1.030), (conf, 7, 1.015)):
            fc, eps0, ft = float(m[off]), float(m[off + 1]), float(m[off + 5])
            ec = 57.0 * math.sqrt(-fc * 1000.0)
            et = 2.0 * ft / ec
            ops.uniaxialMaterial(
                "ConcreteCM", tag, fc, eps0, ec, 7.0, xcrn,
                ft, et, 1.2, 10000.0, "-GapClose", 1,
            )
        raw_steel, limited, steel = 401000 + s, 402000 + s, 400000 + s
        pinching4(raw_steel, m[14:21:2], m[15:22:2], m[22:29:2], m[23:30:2])
        ops.uniaxialMaterial(
            "MinMax", limited, raw_steel, "-min", float(m[30]), "-max", float(m[31])
        )
        ops.uniaxialMaterial("Parallel", steel, limited, 100002)
        shear = 500000 + s
        pinching4(shear, m[32:39:2], m[33:40:2], -m[32:39:2], -m[33:40:2])
        mats.append((uncon, conf, steel, shear))

    ops.geomTransf("Linear", 1)
    ops.geomTransf("PDelta", 2)
    ops.geomTransf("Corotational", 3)
    prev_lean = 200000
    for s, (wn, ln) in enumerate(zip(story_wall_nodes, story_lean_nodes), start=1):
        ops.element("elasticBeamColumn", 200000 + 100 * s, prev_lean, ln, 1e7, 1e3, 1e-5, 2)
        ops.element("corotTruss", 300000 + 100 * s, wn, ln, 1e2, 100001)
        prev_lean = ln

    prev_wall = 100000
    for s, (g, m, mat) in enumerate(zip(gg, mm, mats), start=1):
        thick, width, rho, boundary = fiber_layout(g, m)
        conc = [mat[1] if b else mat[0] for b in boundary]
        steel = [mat[2]] * len(boundary)
        for j in range(1, int(g[2]) + 1):
            jnode = 100000 + 100 * s + j
            etag = jnode
            ops.element(
                "MVLEM", etag, 0.0, prev_wall, jnode, len(width), 0.4,
                "-thick", *thick, "-width", *width, "-rho", *rho,
                "-matConcrete", *conc, "-matSteel", *steel, "-matShear", mat[3],
            )
            prev_wall = jnode
    return gg, mm, story_wall_nodes, story_lean_nodes


def apply_gravity(gg, story_wall_nodes, story_lean_nodes):
    nstory = len(gg)
    grav_a = 21.0 * 60.0
    lean_a = (102.0 * 102.0) / 2.0 - grav_a
    typ_lean = -(1.05 * 125.0 + 0.25 * 40.0) / 1000.0 * lean_a
    roof_lean = -(1.05 * 125.0 + 0.25 * 20.0) / 1000.0 * lean_a
    typ_grav = -(1.05 * 125.0 + 0.25 * 40.0) / 1000.0 * grav_a
    roof_grav = -(1.05 * 125.0 + 0.25 * 20.0) / 1000.0 * grav_a
    typ_self = -float(gg[0, 4]) / 12.0 * float(gg[0, 1]) / 12.0 * 0.150
    roof_self = typ_self / 2.0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in range(nstory):
        roof = s == nstory - 1
        wall_load = (roof_grav if roof else typ_grav) + 1.05 * (
            roof_self if roof else typ_self
        ) * float(gg[s, 5]) / 12.0
        ops.load(story_wall_nodes[s], 0.0, wall_load, 0.0)
        ops.load(story_lean_nodes[s], 0.0, roof_lean if roof else typ_lean, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-6, 1000, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 0.01)
    ops.analysis("Static")
    started = time.perf_counter()
    ok = ops.analyze(100)
    elapsed = time.perf_counter() - started
    ops.reactions()
    reactions = {
        "wall": [ops.nodeReaction(100000, i) for i in (1, 2, 3)],
        "leaning": [ops.nodeReaction(200000, i) for i in (1, 2, 3)],
    }
    applied = []
    for n in story_wall_nodes + story_lean_nodes:
        applied.append(ops.nodeUnbalance(n, 2))
    return ok, elapsed, reactions, applied


def modal_and_damping(received: bool):
    lam = np.asarray(ops.eigen(2), dtype=float)
    periods = 2.0 * np.pi / np.sqrt(lam)
    zeta = 0.05
    wi, wj = math.sqrt(lam[0]) / 0.2, math.sqrt(lam[1])
    alpha = (1.0 if received else 2.0) * zeta * wi * wj / (wi + wj)
    beta = 0.0 if received else 2.0 * zeta / (wi + wj)
    ops.rayleigh(alpha, 0.0, 0.0, beta)
    achieved = alpha / (2.0 * np.sqrt(lam)) + beta * np.sqrt(lam) / 2.0
    return periods.tolist(), {"alphaM": alpha, "betaKcomm": beta, "modal_zeta": achieved.tolist()}


def run_pushover(folder: Path, received: bool, target_drift=0.05):
    gg, _, wn, ln = build_model(folder)
    gok, gsec, reactions, _ = apply_gravity(gg, wn, ln)
    periods, damping = modal_and_damping(received)
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for i, node in enumerate(wn, start=1):
        ops.load(node, float(i + 1), 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-6, 200, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("DisplacementControl", wn[-1], 1, 0.1)
    ops.analysis("Static")
    roof, shear = [0.0], [0.0]
    limit = target_drift * float(gg[:, 1].sum())
    started = time.perf_counter()
    ok = 0
    while roof[-1] < limit:
        ok = ops.analyze(1)
        if ok != 0:
            break
        ops.reactions()
        roof.append(float(ops.nodeDisp(wn[-1], 1)))
        shear.append(float(-(ops.nodeReaction(100000, 1) + ops.nodeReaction(200000, 1))))
    elapsed = time.perf_counter() - started
    return {
        "gravity_ok": gok == 0, "gravity_seconds": gsec, "gravity_reactions": reactions,
        "periods_seconds": periods, "damping": damping, "analysis_ok": ok == 0,
        "elapsed_seconds": elapsed, "roof_displacement_in": roof,
        "base_shear_kip": shear,
    }


def run_nrha(folder: Path, received: bool, record: Path, scale: float, dt: float):
    gg, _, wn, ln = build_model(folder)
    gok, gsec, reactions, _ = apply_gravity(gg, wn, ln)
    periods, damping = modal_and_damping(received)
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    values = np.loadtxt(record).reshape(-1)
    ops.timeSeries("Path", 10, "-dt", dt, "-values", *values.tolist(), "-factor", scale * G)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 10)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-7, 2000, 0)
    ops.algorithm("KrylovNewton", "-maxDim", 20)
    ops.integrator("TRBDF2")
    ops.analysis("Transient")
    heights = [ops.nodeCoord(n, 2) for n in wn]
    peak = np.zeros(len(wn))
    history, iterations, ok = [], 0, 0
    started = time.perf_counter()
    for step in range(len(values)):
        ok = ops.analyze(1, dt)
        if ok != 0:
            break
        disp = np.asarray([ops.nodeDisp(n, 1) for n in wn])
        drift = np.diff(np.r_[0.0, disp]) / np.diff(np.r_[0.0, heights])
        peak = np.maximum(peak, np.abs(drift))
        iterations += ops.testIter()
        history.append([ops.getTime(), *disp.tolist(), *drift.tolist()])
    elapsed = time.perf_counter() - started
    return {
        "gravity_ok": gok == 0, "gravity_seconds": gsec, "gravity_reactions": reactions,
        "periods_seconds": periods, "damping": damping, "scale_factor": scale,
        "dt_seconds": dt, "steps_requested": len(values), "steps_completed": len(history),
        "analysis_ok": ok == 0, "elapsed_seconds": elapsed, "iterations": iterations,
        "peak_story_drift_ratio": peak.tolist(), "history": history,
    }


def run_steel_cycle(folder: Path):
    _, mm = load_inputs(folder)
    m = mm[0]
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    pinching4(1, m[14:21:2], m[15:22:2], m[22:29:2], m[23:30:2])
    ops.uniaxialMaterial("MinMax", 2, 1, "-min", float(m[30]), "-max", float(m[31]))
    ops.uniaxialMaterial("Elastic", 4, 1.0e-2)
    ops.uniaxialMaterial("Parallel", 3, 2, 4)
    ops.element("zeroLength", 1, 1, 2, "-mat", 3, "-dir", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    targets = [0.003, 0.0, -0.003, 0.0, 0.015, 0.0, -0.015, 0.0, 0.06, 0.0, -0.06, 0.0]
    u, f = [0.0], [0.0]
    for target in targets:
        while abs(target - u[-1]) > 5e-8:
            du = math.copysign(min(0.0001, abs(target - u[-1])), target - u[-1])
            ops.integrator("DisplacementControl", 2, 1, du)
            if ops.analyze(1) != 0:
                return {"analysis_ok": False, "strain": u, "stress_ksi": f}
            ops.reactions()
            u.append(float(ops.nodeDisp(2, 1)))
            f.append(float(-ops.nodeReaction(1, 1)))
    return {"analysis_ok": True, "strain": u, "stress_ksi": f}


def summarize_pushover(result):
    v = np.asarray(result["base_shear_kip"])
    d = np.asarray(result["roof_displacement_in"])
    imax = int(np.argmax(v))
    out = {"vmax_kip": float(v[imax]), "d_at_vmax_in": float(d[imax])}
    ix = np.flatnonzero(v[imax:] <= 0.8 * v[imax])
    out["descending_80pct_reached"] = bool(ix.size)
    if ix.size:
        out["d_at_descending_80pct_in"] = float(d[imax + ix[0]])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--received", type=Path, required=True)
    ap.add_argument("--corrected", type=Path, required=True)
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--scale", type=float, default=0.0785532564)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    all_results = {"opensees_version": ops.version(), "archetype": 10001}
    for label, folder, received in (
        ("received", args.received, True), ("corrected", args.corrected, False)
    ):
        cyclic = run_steel_cycle(folder)
        push = run_pushover(folder, received)
        nrha = run_nrha(folder, received, args.record, args.scale, args.dt)
        all_results[label] = {
            "cyclic": cyclic, "pushover": push, "pushover_summary": summarize_pushover(push),
            "nrha": nrha,
        }
    (args.out / "opensees_baseline_10001.json").write_text(json.dumps(all_results, indent=2))
    compact = {"opensees_version": all_results["opensees_version"], "archetype": 10001}
    for label in ("received", "corrected"):
        x = all_results[label]
        compact[label] = {
            "cyclic_ok": x["cyclic"]["analysis_ok"],
            "periods_seconds": x["nrha"]["periods_seconds"],
            "damping": x["nrha"]["damping"],
            "gravity_reactions": x["nrha"]["gravity_reactions"],
            "pushover": x["pushover_summary"],
            "pushover_seconds": x["pushover"]["elapsed_seconds"],
            "nrha_ok": x["nrha"]["analysis_ok"],
            "nrha_steps": x["nrha"]["steps_completed"],
            "nrha_seconds": x["nrha"]["elapsed_seconds"],
            "peak_story_drift_ratio": x["nrha"]["peak_story_drift_ratio"],
        }
    (args.out / "opensees_baseline_10001_summary.json").write_text(json.dumps(compact, indent=2))
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
