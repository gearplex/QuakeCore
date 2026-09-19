#!/usr/bin/env python3
"""Aggregate Phase 9N sharded suite summaries without changing acceptance rules."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def vmax(results, key, default=float("-inf")):
    vals = [float(r[key]) for r in results if key in r]
    return max(vals) if vals else default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    summary_files = sorted(args.root.glob("**/suite_summary.json"))
    if not summary_files:
        raise SystemExit("no shard suite_summary.json files found")
    shards = [json.loads(p.read_text()) for p in summary_files]
    results = []
    for shard in shards:
        results.extend(shard.get("results", []))
    results.sort(key=lambda r: str(r["record_id"]))

    ids = [str(r["record_id"]) for r in results]
    if len(ids) != len(set(ids)):
        raise SystemExit("duplicate record IDs across shards")
    if len(ids) != 44:
        raise SystemExit(f"expected 44 unique records, found {len(ids)}")

    completed_both = [r["record_id"] for r in results if r["completion_status"]["both_completed"]]
    shared_stop = [r["record_id"] for r in results if r["completion_status"]["shared_numerical_noncompletion_same_step"]]
    os_limited = [r["record_id"] for r in results if r["completion_status"]["opensees_reference_limited"]]
    qc_limited = [r["record_id"] for r in results if r["completion_status"]["quakecore_limited"]]
    response_pass = [r["record_id"] for r in results if r["passed_frozen_gate4_response_metrics"]]

    out = {
        "scope": "Phase 9N full 44-component FEMA P-695/ATC-63 far-field external software/model validation at the Phase 9M 0.03 source-acceleration multiplier",
        "claim_boundary": "software/model engineering response parity of the documented reconstructed surrogate only; not source identity, physical validation, FEMA P-695 collapse qualification/acceptance, code approval, or R-factor recommendation",
        "shard_count": len(shards),
        "record_count": len(results),
        "record_ids": ids,
        "dimensionless_acceleration_scale": shards[0]["dimensionless_acceleration_scale"],
        "frozen_gate4_thresholds": shards[0]["frozen_gate4_thresholds"],
        "completion_status_is_separate_from_response_parity": True,
        "records_passing_frozen_gate4_response_metrics": response_pass,
        "passed_frozen_gate4_response_metrics_for_all_records": len(response_pass) == 44,
        "completion_status": {
            "both_completed": completed_both,
            "shared_numerical_noncompletion_same_step": shared_stop,
            "opensees_reference_limited": os_limited,
            "quakecore_limited": qc_limited,
        },
        "suite_summary": {
            "minimum_common_steps": min(int(r["common_steps"]) for r in results),
            "maximum_period_relative_difference": vmax(results, "period_relative_difference"),
            "maximum_peak_story_drift_relative_difference": vmax(results, "peak_story_drift_relative_difference"),
            "maximum_story_drift_rmse_percent_of_peak": max(max(r["story_drift_rmse_percent_of_peak"]) for r in results),
            "maximum_wall_peak_edp_relative_difference": vmax(results, "max_wall_peak_edp_relative_difference"),
            "maximum_peak_floor_acceleration_relative_difference_diagnostic": vmax(results, "peak_floor_absolute_acceleration_relative_difference"),
            "maximum_floor_acceleration_rmse_percent_of_peak_diagnostic": max(max(r["floor_absolute_acceleration_rmse_percent_of_peak"]) for r in results),
            "maximum_residual_story_drift_absolute_difference_diagnostic": max(max(abs(float(x)) for x in r["residual_story_drift_absolute_difference"]) for r in results),
        },
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({
        "record_count": out["record_count"],
        "all_response_metrics_pass": out["passed_frozen_gate4_response_metrics_for_all_records"],
        "completion_counts": {k: len(v) for k, v in out["completion_status"].items()},
        **out["suite_summary"],
    }, indent=2))
    if not out["passed_frozen_gate4_response_metrics_for_all_records"]:
        raise SystemExit("one or more records failed frozen Gate 4 engineering response metrics")


if __name__ == "__main__":
    main()
