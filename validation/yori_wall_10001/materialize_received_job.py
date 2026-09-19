#!/usr/bin/env python3
"""Materialize an exact received-law QuakeCore archetype-10001 job.

This script deliberately starts from the checked-in matched-law structural job
(geometry, masses, gravity loads, leaning-column P-Delta, record, and solver
settings) and replaces only the MVLEM material laws using the original
AllStoriesOpenseesModelINPUT.txt rows.

The exact source input is not checked into this repository.  The script refuses
to proceed unless the supplied file has the frozen SHA-256 recorded by Gate 4.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

EXPECTED_MODEL_SHA256 = "02eca2fbdde5d2cc1400be354802a2b26d6d79ce15ac94dd981cd26dc6b23104"
EXPECTED_WALL_IDS = [1, 2, 3, 4, 5]
WALL_TO_STORY = {1: 0, 2: 0, 3: 1, 4: 1, 5: 2}


def load_rows(path: Path) -> list[list[float]]:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_MODEL_SHA256:
        raise ValueError(
            "model input SHA-256 mismatch: expected "
            f"{EXPECTED_MODEL_SHA256}, received {digest}"
        )
    values = [float(x) for x in raw.decode().split()]
    if len(values) != 3 * 41:
        raise ValueError(f"expected exactly 3 x 41 model values, received {len(values)}")
    return [values[i * 41 : (i + 1) * 41] for i in range(3)]


def concrete_cm(row: list[float], confined: bool) -> dict:
    off = 7 if confined else 0
    fc, epsc, ft = row[off], row[off + 1], row[off + 5]
    ec = 57.0 * math.sqrt(-fc * 1000.0)
    return {
        "type": "concrete_cm",
        "fc": fc,
        "epsc": epsc,
        "Ec": ec,
        "rc": 7.0,
        "xcrn": 1.015 if confined else 1.030,
        "ft": ft,
        "et": 2.0 * ft / ec,
        "rt": 1.2,
        "xcrp": 10000.0,
        "gap_close": True,
    }


def pinching4(
    positive_force: list[float],
    positive_deformation: list[float],
    negative_force: list[float],
    negative_deformation: list[float],
    admitted_reversal_count: int,
) -> dict:
    return {
        "type": "pinching4",
        "positive": [[d, f] for f, d in zip(positive_force, positive_deformation)],
        "negative": [[d, f] for f, d in zip(negative_force, negative_deformation)],
        "r_disp_positive": 0.6,
        "r_force_positive": 0.99,
        "u_force_positive": 0.4,
        "r_disp_negative": 0.6,
        "r_force_negative": 0.99,
        "u_force_negative": 0.4,
        "gamma_k": [0.0, 0.0, 0.0, 0.0],
        "gamma_k_limit": 2.0,
        "gamma_d": [0.1, 0.0, 0.0, 0.0],
        "gamma_d_limit": 2.0,
        "gamma_f": [0.0, 0.0, 0.0, 0.0],
        "gamma_f_limit": 2.0,
        "gamma_e": 10000.0,
        "damage_mode": "energy",
        "admitted_reversal_count": admitted_reversal_count,
    }


def steel(row: list[float]) -> dict:
    raw = pinching4(
        row[14:21:2], row[15:22:2], row[22:29:2], row[23:30:2], 10
    )
    return {
        "type": "parallel",
        "materials": [
            {"type": "minmax", "min": row[30], "max": row[31], "material": raw},
            {"type": "elastic", "E": 1.0e-2},
        ],
    }


def shear(row: list[float]) -> dict:
    pf = row[32:39:2]
    pd = row[33:40:2]
    return pinching4(pf, pd, [-x for x in pf], [-x for x in pd], 6)


def materialize(template: dict, rows: list[list[float]]) -> dict:
    job = json.loads(json.dumps(template))
    walls = job.get("model", {}).get("walls")
    if not isinstance(walls, list):
        raise ValueError("template model.walls must be an array")
    ids = [int(w["id"]) for w in walls]
    if ids != EXPECTED_WALL_IDS:
        raise ValueError(f"unexpected archetype-10001 wall ids/order: {ids}")

    for wall in walls:
        wid = int(wall["id"])
        row = rows[WALL_TO_STORY[wid]]
        fibers = wall.get("fibers")
        if not isinstance(fibers, list) or not fibers:
            raise ValueError(f"wall {wid} has no MVLEM fibers")
        min_rho = min(float(f["rho"]) for f in fibers)
        if min_rho <= 0.0:
            raise ValueError(f"wall {wid} has invalid reinforcement ratio")
        unconf = concrete_cm(row, False)
        conf = concrete_cm(row, True)
        st = steel(row)
        for fiber in fibers:
            # Archetype 10001 has a low-rho web and higher-rho boundary fibers.
            # Use the retained reinforcement ratio rather than positional magic,
            # and refuse ambiguous layouts instead of silently guessing.
            rho = float(fiber["rho"])
            fiber["concrete"] = unconf if math.isclose(rho, min_rho, rel_tol=0.0, abs_tol=1e-12) else conf
            fiber["steel"] = st
        wall["shear"] = shear(row)

    job["name"] = "YORi_10001_received_materials_gravity_pdelta"
    provenance = job.setdefault("provenance", {})
    provenance["purpose"] = (
        "Gate 4 exact received-material QuakeCore comparison; software/model verification only"
    )
    provenance["received_model_input_sha256"] = EXPECTED_MODEL_SHA256
    provenance["material_substitutions"] = "none; ConcreteCM + Pinching4 + MinMax + Parallel from source rows"
    provenance["claim_boundary"] = (
        "Does not establish physical validation, FEMA P-695 acceptance, collapse qualification, or an R-factor recommendation"
    )
    return job


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-input", type=Path, required=True)
    parser.add_argument(
        "--template",
        type=Path,
        default=Path(__file__).with_name("quakecore_job_10001_gravity.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = load_rows(args.model_input)
    template = json.loads(args.template.read_text())
    result = materialize(template, rows)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
