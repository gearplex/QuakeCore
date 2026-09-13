#!/usr/bin/env python3
"""Freeze OpenSees 3.8.0 ConcreteCM rule-77 reversal checkpoints.

The path follows the admitted YORi story-1 ConcreteCM protocol through step 77
(compression envelope -> first positive-going unloading, still on rule 3), then
reverses negative before reaching the zero-stress point. OpenSees therefore
enters rule 77 and subsequently rejoins the compression envelope.
"""
from __future__ import annotations

import json
import math
import os

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import openseespy.opensees as ops


def interpolate(targets, subdivisions=20):
    out = [float(targets[0])]
    for a, b in zip(targets, targets[1:]):
        for i in range(1, subdivisions + 1):
            out.append(float(a + (b - a) * i / subdivisions))
    return out


def material(tag: int, fc: float, epsc: float, ft: float, xcrn: float) -> None:
    ec = 57.0 * math.sqrt(-fc * 1000.0)
    et = 2.0 * ft / ec
    ops.uniaxialMaterial(
        "ConcreteCM", tag, fc, epsc, ec, 7.0, xcrn,
        ft, et, 1.2, 10000.0, "-GapClose", 1,
    )


def evaluate(name: str, fc: float, epsc: float, ft: float, xcrn: float) -> dict:
    ops.wipe()
    material(1, fc, epsc, ft, xcrn)
    ops.testUniaxialMaterial(1)

    base = interpolate([0.0, -0.0005, -0.0020, -0.0040, -0.0010])[:78]
    probes = [-0.0016, -0.0020, -0.0030, -0.0040, -0.0041, -0.0042, -0.0044, -0.0046]
    path = base + probes
    history = []
    for step, strain in enumerate(path):
        ops.setStrain(strain)
        if step >= len(base) - 1:
            history.append({
                "step": step,
                "strain": strain,
                "stress": float(ops.getStress()),
                "tangent": float(ops.getTangent()),
            })
    return {"name": name, "history": history}


def main() -> None:
    result = {
        "schema": "quakecore.opensees_concretecm_rule77_oracle.v1",
        "opensees_version": ops.version(),
        "base_protocol_commit_step": 77,
        "materials": [
            evaluate("unconfined", -6.5, -0.002, 0.0604669, 1.030),
            evaluate("confined", -8.01435, -0.00432976, 0.0671422, 1.015),
        ],
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
