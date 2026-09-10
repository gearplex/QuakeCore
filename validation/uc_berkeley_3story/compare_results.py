#!/usr/bin/env python3
"""Compare a QuakeCore validation result JSON against published NIST targets."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent

def pct_err(value: float, ref: float) -> float:
    return 100.0 * (value - ref) / ref if ref else math.nan

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("result", type=Path, help="QuakeCore result JSON")
    ap.add_argument("--reference", choices=("perform3d", "measured"), default="perform3d")
    args = ap.parse_args()
    targets = json.loads((HERE / "reference_targets.json").read_text())
    result = json.loads(args.result.read_text())
    ref = targets["targets"][args.reference]

    rows = []
    def scalar(name, result_key, ref_key):
        if result_key in result:
            rows.append((name, float(result[result_key]), float(ref[ref_key])))
    scalar("T1 (s)", "first_mode_period_s", "first_mode_period_s")
    scalar("Peak base shear (kip)", "peak_base_shear_kip", "peak_base_shear_kip")
    scalar("Story 1 residual drift (%)", "residual_story1_drift_percent", "residual_story1_drift_percent")
    scalar("Roof residual drift (%)", "residual_roof_drift_percent", "residual_roof_drift_percent")
    if "peak_story_drift_percent" in result:
        for i, (v, r) in enumerate(zip(result["peak_story_drift_percent"], ref["peak_story_drift_percent"]), 1):
            rows.append((f"Story {i} peak drift (%)", float(v), float(r)))

    print(f"Reference: {args.reference}")
    print(f"{'Metric':34s} {'QuakeCore':>12s} {'Reference':>12s} {'Error':>10s}")
    print("-" * 72)
    for name, val, r in rows:
        print(f"{name:34s} {val:12.6g} {r:12.6g} {pct_err(val,r):9.3f}%")

if __name__ == "__main__":
    main()
