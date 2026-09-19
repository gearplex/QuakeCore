#!/usr/bin/env python3
"""Inject one archived FEMA P-695 far-field component into a Phase 9N job.

Raw source accelerations are in g.  QuakeCore jobs use inch-second units, so
analysis acceleration is source_g * 386.09 * dimensionless_scale.  A scale of
0.03 reproduces the Phase 9M 120111 benchmark history to floating-point
roundoff.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path, PurePosixPath

G_IN_PER_S2 = 386.09
SORTED_ZIP = "2c_ATC-63_Far-Field_GroundMotionAccelTextFiles_Unscaled_Sorted.zip"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_record(source_archive: Path, manifest: dict, record_id: str) -> tuple[dict, list[float]]:
    by_id = {str(r["record_id"]): r for r in manifest["records"]}
    if record_id not in by_id:
        raise ValueError(f"record {record_id!r} is not present in the Phase 9N manifest")
    meta = by_id[record_id]

    outer_bytes = source_archive.read_bytes()
    expected_outer = manifest.get("source_archive", {}).get("sha256")
    if expected_outer and sha256_bytes(outer_bytes) != expected_outer:
        raise ValueError("source archive SHA-256 does not match the Phase 9N manifest")

    with zipfile.ZipFile(io.BytesIO(outer_bytes)) as outer:
        sorted_bytes = outer.read(SORTED_ZIP)
    expected_sorted = manifest.get("source_archive", {}).get("sorted_archive_sha256")
    if expected_sorted and sha256_bytes(sorted_bytes) != expected_sorted:
        raise ValueError("sorted-record archive SHA-256 does not match the Phase 9N manifest")

    target = meta["sorted_file"]
    with zipfile.ZipFile(io.BytesIO(sorted_bytes)) as zf:
        candidates = [name for name in zf.namelist() if PurePosixPath(name).name == target]
        if len(candidates) != 1:
            raise ValueError(f"expected one {target} entry, found {len(candidates)}")
        data = zf.read(candidates[0])
    if sha256_bytes(data) != meta["sorted_file_sha256"]:
        raise ValueError(f"record {record_id} SHA-256 does not match the manifest")
    values = [float(x) for x in data.decode("ascii", errors="strict").split()]
    if len(values) != int(meta["npts"]):
        raise ValueError(f"record {record_id} has {len(values)} values, expected {meta['npts']}")
    return meta, values


def materialize(base_job: dict, meta: dict, source_g: list[float], scale: float) -> dict:
    job = json.loads(json.dumps(base_job))
    acceleration = [float(v) * G_IN_PER_S2 * float(scale) for v in source_g]
    job["name"] = f"YORi_10001_phase9n_record_{meta['record_id']}_scale_{scale:g}"
    job["records"] = [{
        "name": meta["sorted_file"].removesuffix(".txt"),
        "dt": float(meta["dt_s"]),
        "sample_convention": "step_end",
        "acceleration": acceleration,
    }]
    provenance = job.setdefault("provenance", {})
    provenance["purpose"] = (
        "Phase 9N multi-record exact-law QuakeCore/OpenSees software/model validation; "
        "not FEMA P-695 collapse qualification"
    )
    provenance["phase9n_record"] = {
        "record_id": meta["record_id"],
        "pair_id": meta["pair_id"],
        "component_index": meta["component_index"],
        "peer_header": meta["peer_header"],
        "original_file": meta["original_file"],
        "sorted_file": meta["sorted_file"],
        "sorted_file_sha256": meta["sorted_file_sha256"],
        "source_units": "g",
        "npts": meta["npts"],
        "dt_s": meta["dt_s"],
        "pga_g": meta["pga_g"],
        "dimensionless_acceleration_scale": float(scale),
        "gravity_in_per_s2": G_IN_PER_S2,
        "analysis_acceleration_units": "in/s^2",
    }
    provenance["claim_boundary"] = (
        "software/model external validation of the documented reconstructed surrogate only; "
        "not identity to the lost received model, physical validation, FEMA P-695 acceptance, "
        "collapse qualification, code approval, or an R-factor recommendation"
    )
    return job


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-job", type=Path, required=True)
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--record-id", required=True)
    ap.add_argument("--scale", type=float, default=0.03)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.scale < 0:
        raise ValueError("scale must be nonnegative")
    base_job = json.loads(args.base_job.read_text())
    manifest = json.loads(args.manifest.read_text())
    meta, source_g = load_record(args.source_archive, manifest, str(args.record_id))
    job = materialize(base_job, meta, source_g, args.scale)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(job, indent=2) + "\n")
    print(json.dumps({
        "output": str(args.out),
        "record_id": meta["record_id"],
        "npts": meta["npts"],
        "dt_s": meta["dt_s"],
        "scale": args.scale,
        "analysis_peak_acceleration_in_per_s2": max(abs(v) for v in job["records"][0]["acceleration"]),
    }, indent=2))


if __name__ == "__main__":
    main()
