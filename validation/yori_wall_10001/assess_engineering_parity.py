#!/usr/bin/env python3
"""Gate-4 engineering-parity acceptance for the reconstructed YORi surrogate.

This gate is intentionally distinct from strict trace parity.  It tests whether
QuakeCore reproduces the modal and engineering-demand quantities used for model
comparison over an independently converged OpenSees window.  It is software
verification of the reconstructed surrogate only; it is not physical
validation, FEMA P-695 acceptance, code qualification, or an R-factor finding.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def max_relative(a, b):
    out = 0.0
    for x, y in zip(a, b):
        out = max(out, abs(x - y) / max(abs(x), abs(y), 1e-30))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=Path, required=True)
    ap.add_argument("--quake", type=Path, required=True)
    ap.add_argument("--response-summary", type=Path, required=True)
    ap.add_argument("--wall-force-summary", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-common-steps", type=int, default=1000)
    ap.add_argument("--max-period-relative", type=float, default=1e-6)
    ap.add_argument("--max-peak-drift-relative", type=float, default=1e-3)
    ap.add_argument("--max-drift-rmse-percent-of-peak", type=float, default=0.1)
    ap.add_argument("--max-wall-peak-edp-relative", type=float, default=0.01)
    args = ap.parse_args()

    job = json.loads(args.job.read_text())
    q = json.loads(args.quake.read_text())
    response = json.loads(args.response_summary.read_text())
    wall = json.loads(args.wall_force_summary.read_text())

    qrun = q["runs"][0]
    record_steps = len(job["records"][0]["acceleration"])
    q_steps = len(qrun.get("history", []))
    q_complete = q_steps == record_steps and qrun.get("termination") not in {"numerical_failure", "initial_instability"}

    period_relative = max_relative(response["quakecore_periods_seconds"], response["opensees_periods_seconds"])
    peak_drift_relative = max_relative(response["quakecore_peak_story_drift_ratio"], response["opensees_peak_story_drift_ratio"])
    drift_rmse = max(response["story_drift_rmse_percent_of_peak"])
    common_steps = int(response["history_steps_compared"])
    wall_steps = int(wall["steps_completed"])
    wall_peak_edp_relative = float(wall["max_relative_peak_edp_difference"])

    checks = {
        "quakecore_completed_full_record": {
            "passed": q_complete,
            "value": q_steps,
            "limit_or_target": record_steps,
        },
        "minimum_common_converged_steps": {
            "passed": common_steps >= args.min_common_steps,
            "value": common_steps,
            "limit_or_target": args.min_common_steps,
        },
        "period_relative_difference": {
            "passed": period_relative <= args.max_period_relative,
            "value": period_relative,
            "limit_or_target": args.max_period_relative,
        },
        "peak_story_drift_relative_difference": {
            "passed": peak_drift_relative <= args.max_peak_drift_relative,
            "value": peak_drift_relative,
            "limit_or_target": args.max_peak_drift_relative,
        },
        "story_drift_rmse_percent_of_peak": {
            "passed": drift_rmse <= args.max_drift_rmse_percent_of_peak,
            "value": drift_rmse,
            "limit_or_target": args.max_drift_rmse_percent_of_peak,
        },
        "wall_force_common_steps": {
            "passed": wall_steps >= args.min_common_steps,
            "value": wall_steps,
            "limit_or_target": args.min_common_steps,
        },
        "wall_peak_force_edp_relative_difference": {
            "passed": wall_peak_edp_relative <= args.max_wall_peak_edp_relative,
            "value": wall_peak_edp_relative,
            "limit_or_target": args.max_wall_peak_edp_relative,
        },
    }
    passed = all(v["passed"] for v in checks.values())
    result = {
        "scope": "Gate-4 engineering parity of the reconstructed 3-story YORi surrogate",
        "claim_boundary": "software/model verification only; not identity to the lost received model input, physical validation, FEMA P-695 acceptance, collapse qualification, code approval, or an R-factor recommendation",
        "acceptance_basis": "modal periods, peak/global drift response, common-window drift-history error, and MVLEM peak wall shear/axial/moment EDPs",
        "thresholds_are_project_software_verification_criteria": True,
        "passed": passed,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not passed:
        raise SystemExit("engineering parity gate failed")


if __name__ == "__main__":
    main()
