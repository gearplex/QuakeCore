#!/usr/bin/env python3
"""Create the matched-law archetype-10001 QuakeCore job with gravity/P-Delta."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent / "yori_package" / "deliverable_10001"))
from benchmark_10001 import make_job  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", type=Path, required=True)
    ap.add_argument("--motion", type=Path, required=True)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--scale", type=float, default=0.03)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    job = make_job(args.folder, args.motion, args.dt, args.scale, "full")
    job["name"] = "YORi_10001_matched_gravity_pdelta"
    job["provenance"]["purpose"] = "gravity-state/P-Delta solver verification; not FEMA P-695 acceptance"
    job["provenance"]["substitutions"] = (
        "Concrete01, Steel01, bilinear shear; YORi gravity loads; "
        "leaning-column constant segment preloads; rigid-floor UX constraints"
    )
    model = job["model"]
    model["nodes"].append([100, 72.0, 0.0, 0.0, 0.0, 0.0])
    model["fixities"].append([100, True, True, True])
    for story, y in enumerate((144.0, 288.0, 432.0), start=1):
        model["nodes"].append([100 + story, 72.0, y, 0.0, 0.0, 0.0])
        model["equal_dofs"].append([model["story_nodes"][story - 1], 100 + story, "UX"])

    gravity_a = 21.0 * 60.0
    leaning_a = 102.0 * 102.0 / 2.0 - gravity_a
    wall_typ = -(1.05 * 125.0 + 0.25 * 40.0) / 1000.0 * gravity_a - 1.05 * 10.8
    wall_roof = -(1.05 * 125.0 + 0.25 * 20.0) / 1000.0 * gravity_a - 1.05 * 5.4
    lean_typ = -(1.05 * 125.0 + 0.25 * 40.0) / 1000.0 * leaning_a
    lean_roof = -(1.05 * 125.0 + 0.25 * 20.0) / 1000.0 * leaning_a
    lean_loads = [lean_typ, lean_typ, lean_roof]
    segment_preloads = [-(sum(lean_loads[s:])) for s in range(3)]
    previous = 100
    for s in range(3):
        node = 101 + s
        model["members"].append([201 + s, previous, node, 1.0e3, 1.0e7, 1.0e-5, segment_preloads[s], True])
        previous = node
    loads = []
    for s, wall_node in enumerate(model["story_nodes"]):
        loads.append([wall_node, 0.0, wall_roof if s == 2 else wall_typ, 0.0])
        loads.append([101 + s, 0.0, lean_loads[s], 0.0])
    job["gravity"] = {
        "loads": loads,
        "steps": 100,
        "tolerance": 1.0e-8,
        "max_iterations": 100,
        "provenance": "YORi 2000_Gravity.tcl load combinations and tributary areas",
    }
    job["analysis"]["line_search"] = False
    # The YORi leaning column intentionally uses nearly zero EI, which creates
    # massless rotational pivots below the generic SPD certification scale.
    job["analysis"]["check_initial_stability"] = False
    args.out.write_text(json.dumps(job, indent=2))


if __name__ == "__main__":
    main()
