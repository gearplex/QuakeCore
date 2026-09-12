#!/usr/bin/env python3
"""Generate frozen OpenSees material paths for YORi archetype 10001.

This is an external oracle, not a QuakeCore constitutive implementation.  The
protocols exercise envelope points, reversals, nested cycles, degradation, and
MinMax failure.  OpenSees' uniaxial tester commits the preceding trial when a
new strain is applied, matching the material's normal accepted-step lifecycle.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")

import numpy as np
import openseespy.opensees as ops


def interpolate_targets(targets: list[float], subdivisions: int) -> list[float]:
    values = [float(targets[0])]
    for start, stop in zip(targets, targets[1:]):
        values.extend(np.linspace(start, stop, subdivisions + 1)[1:].tolist())
    return values


def pinching4(tag: int, positive_force, positive_deformation,
              negative_force, negative_deformation) -> None:
    args: list[float | str] = []
    for force, deformation in zip(positive_force, positive_deformation):
        args.extend([float(force), float(deformation)])
    for force, deformation in zip(negative_force, negative_deformation):
        args.extend([float(force), float(deformation)])
    args.extend([0.6, 0.99, 0.4, 0.6, 0.99, 0.4])
    args.extend([0.0, 0.0, 0.0, 0.0, 0.0])
    args.extend([0.1, 0.0, 0.0, 0.0, 2.0])
    args.extend([0.0, 0.0, 0.0, 0.0, 0.0, 10000.0, "energy"])
    ops.uniaxialMaterial("Pinching4", tag, *args)


def define_material(name: str, m: np.ndarray) -> int:
    ops.wipe()
    if name in {"concrete_unconfined", "concrete_confined"}:
        offset = 0 if name == "concrete_unconfined" else 7
        tag = 1
        fc, eps0, ft = map(float, (m[offset], m[offset + 1], m[offset + 5]))
        ec = 57.0 * math.sqrt(-fc * 1000.0)
        et = 2.0 * ft / ec
        xcrn = 1.030 if offset == 0 else 1.015
        ops.uniaxialMaterial(
            "ConcreteCM", tag, fc, eps0, ec, 7.0, xcrn,
            ft, et, 1.2, 10000.0, "-GapClose", 1,
        )
        return tag

    if name in {"steel_raw", "steel_minmax_parallel"}:
        raw = 1
        pinching4(raw, m[14:21:2], m[15:22:2], m[22:29:2], m[23:30:2])
        if name == "steel_raw":
            return raw
        ops.uniaxialMaterial(
            "MinMax", 2, raw, "-min", float(m[30]), "-max", float(m[31])
        )
        ops.uniaxialMaterial("Elastic", 4, 1.0e-2)
        ops.uniaxialMaterial("Parallel", 3, 2, 4)
        return 3

    if name == "shear_pinching4":
        tag = 1
        pinching4(tag, m[32:39:2], m[33:40:2], -m[32:39:2], -m[33:40:2])
        return tag
    raise ValueError(f"unknown material {name}")


def protocol(name: str) -> list[float]:
    if name.startswith("concrete_"):
        targets = [0.0, -0.0005, -0.0020, -0.0040, -0.0010, 0.0002,
                   0.0010, 0.0025, 0.0, -0.0060, -0.0120, -0.0020, 0.0]
        return interpolate_targets(targets, 20)
    if name.startswith("steel_"):
        targets = [0.0, 0.0010, 0.0040, 0.0, -0.0040, 0.0,
                   0.0120, -0.0120, 0.0350, -0.0350, 0.0600, -0.0600,
                   0.0800, -0.0800, 0.0]
        return interpolate_targets(targets, 20)
    targets = [0.0, 0.0001, 0.010, 0.10, 0.30, 0.0, -0.30,
               0.72, -0.72, 1.20, -1.20, 0.0]
    return interpolate_targets(targets, 20)


def evaluate(name: str, m: np.ndarray) -> dict:
    tag = define_material(name, m)
    ops.testUniaxialMaterial(tag)
    history = []
    previous_strain = previous_stress = energy = 0.0
    for step, strain in enumerate(protocol(name)):
        ops.setStrain(float(strain))
        stress = float(ops.getStress())
        tangent = float(ops.getTangent())
        if step:
            energy += 0.5 * (previous_stress + stress) * (strain - previous_strain)
        history.append({
            "step": step,
            "strain_or_deformation": float(strain),
            "stress_or_force": stress,
            "tangent": tangent,
            "signed_work": energy,
        })
        previous_strain, previous_stress = strain, stress
    return {"name": name, "material_tag": tag, "history": history}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-input", type=Path, required=True)
    parser.add_argument("--story", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.model_input.read_bytes()
    data = np.loadtxt(args.model_input).reshape(-1, 41)
    if args.story < 1 or args.story > len(data):
        raise ValueError("story is outside the model input")
    row = data[args.story - 1]
    names = ["concrete_unconfined", "concrete_confined", "steel_raw",
             "steel_minmax_parallel", "shear_pinching4"]
    result = {
        "schema": "quakecore.opensees_material_oracle.v1",
        "opensees_version": ops.version(),
        "source": str(args.model_input),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "story": args.story,
        "protocol_semantics": "each new setStrain commits the preceding trial state",
        "materials": [evaluate(name, row) for name in names],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
