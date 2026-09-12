#!/usr/bin/env python3
"""Integrity and lifecycle checks for the frozen story-1 material oracle."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def main() -> None:
    source = Path(sys.argv[1])
    data = json.loads(source.read_text())
    assert data["schema"] == "quakecore.opensees_material_oracle.v1"
    assert data["opensees_version"] == "3.8.0"
    assert data["source_sha256"] == "02eca2fbdde5d2cc1400be354802a2b26d6d79ce15ac94dd981cd26dc6b23104"
    materials = {item["name"]: item["history"] for item in data["materials"]}
    assert set(materials) == {
        "concrete_unconfined", "concrete_confined", "steel_raw",
        "steel_minmax_parallel", "shear_pinching4",
    }
    for history in materials.values():
        assert history[0]["strain_or_deformation"] == 0.0
        for index, row in enumerate(history):
            assert row["step"] == index
            assert all(math.isfinite(row[key]) for key in
                       ("strain_or_deformation", "stress_or_force", "tangent", "signed_work"))

    raw = materials["steel_raw"]
    wrapped = materials["steel_minmax_parallel"]
    assert max(row["stress_or_force"] for row in raw) > 75.0
    assert min(row["stress_or_force"] for row in raw) < -75.0
    # Once +0.08 crosses the +0.0785049 MinMax bound, the Pinching4 child
    # remains failed; only the 0.01-stiffness Parallel elastic regularizer acts.
    assert abs(wrapped[240]["strain_or_deformation"] - 0.08) < 1e-14
    assert abs(wrapped[240]["stress_or_force"] - 0.0008) < 1e-12
    assert abs(wrapped[-1]["stress_or_force"]) < 1e-14
    assert min(row["stress_or_force"] for row in materials["concrete_unconfined"]) <= -6.5
    assert max(row["stress_or_force"] for row in materials["shear_pinching4"]) > 299.0
    print("YORi story-1 OpenSees material oracle integrity: PASS")


if __name__ == "__main__":
    main()
