#!/usr/bin/env python3
"""Phase 9N multi-record exact-law validation for YORi archetype 10001.

Raw ATC-63/FEMA P-695 accelerations are read from the archived sorted source
package in g, converted to in/s^2 using 386.09, and multiplied by an explicitly
requested software-validation scale. The scale is not Sa(T1), collapse
capacity, or a FEMA P-695 qualification intensity measure.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import zipfile

from compare_reconstructed_intensity_suite import run_scale

SORTED_ZIP = "2c_ATC-63_Far-Field_GroundMotionAccelTextFiles_Unscaled_Sorted.zip"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def classify(full_record: bool, termination) -> str:
    if full_record:
        return "completed"
    if termination in {"response_limit", "response_limit_crossed"}:
        return "response_limit_crossed"
    return "numerical_noncompletion"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--base-job", type=Path, required=True)
    ap.add_argument("--quake-exe", type=Path, required=True)
    ap.add_argument("--record-ids", nargs="+", required=True)
    ap.add_argument("--record-scale", type=float, default=0.03)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-common-steps", type=int, default=1000)
    ap.add_argument("--max-peak-drift-relative", type=float, default=1e-3)
    ap.add_argument("--max-drift-rmse-percent-of-peak", type=float, default=0.1)
    ap.add_argument("--max-wall-peak-edp-relative", type=float, default=0.01)
    args = ap.parse_args()

    if args.record_scale <= 0.0:
        raise SystemExit("record-scale must be positive")

    manifest = json.loads(args.manifest.read_text())
    records = {str(r["record_id"]): r for r in manifest["records"]}
    missing = [rid for rid in args.record_ids if str(rid) not in records]
    if missing:
        raise SystemExit(f"record IDs not present in manifest: {missing}")

    outer_bytes = args.source_archive.read_bytes()
    with zipfile.ZipFile(io.BytesIO(outer_bytes)) as outer:
        sorted_bytes = outer.read(SORTED_ZIP)
    expected_sorted_sha = manifest["source_archive"]["sorted_archive_sha256"]
    actual_sorted_sha = sha256_bytes(sorted_bytes)
    if actual_sorted_sha != expected_sorted_sha:
        raise SystemExit(
            f"sorted archive SHA mismatch: {actual_sorted_sha} != {expected_sorted_sha}"
        )

    base_job = json.loads(args.base_job.read_text())
    conversion = manifest["analysis_conversion"]
    gravity = float(conversion["gravity_in_per_s2"])
    phase9m_scale = float(conversion["phase9m_reproduction_scale"])
    results = []

    with zipfile.ZipFile(io.BytesIO(sorted_bytes)) as zf, tempfile.TemporaryDirectory(
        prefix="phase9n_"
    ) as td_name:
        td = Path(td_name)
        member_by_base = {Path(name).name: name for name in zf.namelist()}

        for requested_id in args.record_ids:
            rid = str(requested_id)
            meta = records[rid]
            source_name = meta["sorted_file"]
            if source_name not in member_by_base:
                raise SystemExit(f"record {rid}: {source_name} not found in sorted archive")
            member = member_by_base[source_name]
            raw_bytes = zf.read(member)
            digest = sha256_bytes(raw_bytes)
            if digest != meta["sorted_file_sha256"]:
                raise SystemExit(f"record {rid}: source member SHA mismatch")

            raw_g = [float(x) for x in raw_bytes.decode().split()]
            if len(raw_g) != int(meta["npts"]):
                raise SystemExit(
                    f"record {rid}: NPTS mismatch {len(raw_g)} != {meta['npts']}"
                )

            acceleration = [x * gravity * args.record_scale for x in raw_g]
            job = copy.deepcopy(base_job)
            job["name"] = f"YORi_10001_phase9n_record_{rid}_scale_{args.record_scale:g}"
            job["records"] = [{
                "name": f"SortedEQFile_({rid})",
                "dt": float(meta["dt_s"]),
                "sample_convention": "step_end",
                "acceleration": acceleration,
                "source_units": "g",
                "source_record_id": rid,
                "source_scale": args.record_scale,
                "gravity_in_per_s2": gravity,
            }]
            job.setdefault("provenance", {})["phase9n_record"] = {
                "record_id": rid,
                "sorted_file_sha256": digest,
                "source_scale": args.record_scale,
                "scale_interpretation": (
                    "dimensionless multiplier on original unscaled acceleration in g "
                    "before conversion to in/s^2; not Sa(T1)"
                ),
            }

            regression = None
            if rid == "120111" and abs(args.record_scale - phase9m_scale) <= 1e-15:
                base_rec = base_job["records"][0]
                same_length = len(base_rec["acceleration"]) == len(acceleration)
                same_dt = abs(float(base_rec["dt"]) - float(meta["dt_s"])) <= 1e-15
                max_abs = (
                    max(
                        abs(float(a) - float(b))
                        for a, b in zip(base_rec["acceleration"], acceleration)
                    )
                    if same_length
                    else float("inf")
                )
                regression = {
                    "same_dt": same_dt,
                    "same_length": same_length,
                    "max_abs_acceleration_difference": max_abs,
                    "passed": same_dt and same_length and max_abs <= 1e-12,
                }
                if not regression["passed"]:
                    raise SystemExit(f"Phase 9M 120111 regression failed: {regression}")

            comparison = run_scale(job, args.quake_exe.resolve(), 1.0, td)
            qstatus = classify(
                bool(comparison.get("quakecore_full_record")),
                comparison.get("quakecore_termination"),
            )
            ostatus = (
                "completed"
                if comparison.get("opensees_full_record")
                else "numerical_noncompletion"
            )

            drift_rmse = max(
                comparison.get("story_drift_rmse_percent_of_peak", [float("inf")])
            )
            checks = {
                "quakecore_completed_full_record": bool(
                    comparison.get("quakecore_full_record")
                ),
                "minimum_common_steps": int(comparison.get("common_steps", 0))
                >= args.min_common_steps,
                "peak_story_drift_relative_difference": float(
                    comparison.get("peak_story_drift_relative_difference", float("inf"))
                )
                <= args.max_peak_drift_relative,
                "story_drift_rmse_percent_of_peak": float(drift_rmse)
                <= args.max_drift_rmse_percent_of_peak,
                "wall_peak_edp_relative_difference": float(
                    comparison.get("max_wall_peak_edp_relative_difference", float("inf"))
                )
                <= args.max_wall_peak_edp_relative,
            }

            results.append({
                "record_id": rid,
                "pair_id": meta["pair_id"],
                "component_index": meta["component_index"],
                "peer_header": meta["peer_header"],
                "original_file": meta["original_file"],
                "dt_seconds": meta["dt_s"],
                "npts": meta["npts"],
                "source_pga_g": meta["pga_g"],
                "record_scale": args.record_scale,
                "input_peak_acceleration_in_per_s2": max(abs(x) for x in acceleration),
                "quakecore_status": qstatus,
                "opensees_status": ostatus,
                "phase9m_120111_regression": regression,
                "engineering_gate_passed": all(checks.values()),
                "engineering_gate_checks": checks,
                "comparison": comparison,
            })

    report = {
        "scope": (
            "Phase 9N multi-record exact-law external software/model validation pilot "
            "for the documented reconstructed YORi 10001 surrogate"
        ),
        "claim_boundary": (
            "not exact-source recovery, physical validation, FEMA P-695 collapse "
            "qualification/acceptance, code approval, or an R-factor recommendation"
        ),
        "record_scale": args.record_scale,
        "scale_interpretation": (
            "dimensionless multiplier on original unscaled acceleration; not Sa(T1)"
        ),
        "source_sorted_archive_sha256": actual_sorted_sha,
        "frozen_gate_thresholds": {
            "min_common_steps": args.min_common_steps,
            "max_peak_drift_relative": args.max_peak_drift_relative,
            "max_drift_rmse_percent_of_peak": args.max_drift_rmse_percent_of_peak,
            "max_wall_peak_edp_relative": args.max_wall_peak_edp_relative,
        },
        "records": results,
        "suite": {
            "record_count": len(results),
            "engineering_gate_pass_count": sum(
                bool(x["engineering_gate_passed"]) for x in results
            ),
            "quakecore_completed_count": sum(
                x["quakecore_status"] == "completed" for x in results
            ),
            "opensees_completed_count": sum(
                x["opensees_status"] == "completed" for x in results
            ),
            "minimum_common_steps": min(
                (x["comparison"].get("common_steps", 0) for x in results), default=0
            ),
            "maximum_peak_drift_relative_difference": max(
                (
                    x["comparison"].get(
                        "peak_story_drift_relative_difference", float("inf")
                    )
                    for x in results
                ),
                default=float("inf"),
            ),
            "maximum_drift_rmse_percent_of_peak": max(
                (
                    max(
                        x["comparison"].get(
                            "story_drift_rmse_percent_of_peak", [float("inf")]
                        )
                    )
                    for x in results
                ),
                default=float("inf"),
            ),
            "maximum_wall_peak_edp_relative_difference": max(
                (
                    x["comparison"].get(
                        "max_wall_peak_edp_relative_difference", float("inf")
                    )
                    for x in results
                ),
                default=float("inf"),
            ),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["suite"], indent=2))

    if not all(x["engineering_gate_passed"] for x in results):
        raise SystemExit("one or more Phase 9N pilot records failed the frozen gate")


if __name__ == "__main__":
    main()
