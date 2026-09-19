#!/usr/bin/env python3
"""Build a metadata-only manifest for the FEMA P-695/ATC-63 far-field set.

The input is the private/source archive assembled for Phase 9N.  No acceleration
histories are written to the manifest.  Sorted components are matched back to
original PEER .AT2 files by the parsed acceleration vector, not workbook row
ordering.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import struct
import zipfile
from pathlib import Path, PurePosixPath

ORIGINAL_ZIP = "2b_ATC-63_Far-Field_GroundMotionAccelTextFiles_Unscaled_Original.zip"
SORTED_ZIP = "2c_ATC-63_Far-Field_GroundMotionAccelTextFiles_Unscaled_Sorted.zip"
SUMMARY_XLS = "1_ReadMe_ATC-63_GroundMotionSet_FarField_SummaryFile.xls"
SORTED_RE = re.compile(r"SortedEQFile_\((\d{6})\)\.txt$")
NPTS_DT_RE = re.compile(r"([0-9]+)\s+([0-9.]+)\s+NPTS\s*,\s*DT", re.I)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def vector_fingerprint(values: list[float]) -> str:
    h = hashlib.sha256()
    for value in values:
        h.update(struct.pack("!d", float(value)))
    return h.hexdigest()


def parse_sorted(data: bytes) -> list[float]:
    return [float(x) for x in data.decode("ascii", errors="strict").split()]


def parse_at2(data: bytes) -> dict:
    lines = data.decode("ascii", errors="replace").splitlines()
    if len(lines) < 5:
        raise ValueError("PEER AT2 file has fewer than five lines")
    match = NPTS_DT_RE.search(lines[3])
    if not match:
        raise ValueError(f"could not parse NPTS/DT line: {lines[3]!r}")
    npts = int(match.group(1))
    dt = float(match.group(2))
    values: list[float] = []
    for line in lines[4:]:
        for token in line.replace("D", "E").split():
            values.append(float(token))
    if len(values) < npts:
        raise ValueError(f"AT2 data contains {len(values)} values, expected {npts}")
    values = values[:npts]
    return {
        "header_database": lines[0].strip(),
        "header_record": lines[1].strip(),
        "header_units": lines[2].strip(),
        "npts": npts,
        "dt_s": dt,
        "values": values,
    }


def canonical_candidate(record_id: str, candidates: list[dict]) -> tuple[dict, list[str]]:
    if len(candidates) == 1:
        return candidates[0], []
    paths = [c["original_file"] for c in candidates]
    # The source archive contains duplicate POE files under LANDERS and SUPERST.
    # Their headers and acceleration vectors identify them as Superstition Hills.
    if record_id in {"121221", "121222"}:
        preferred = [c for c in candidates if c["original_file"].startswith("SUPERST/")]
        if len(preferred) == 1:
            return preferred[0], sorted(p for p in paths if p != preferred[0]["original_file"])
    raise ValueError(f"record {record_id} has unresolved source candidates: {paths}")


def build_manifest(source_archive: Path) -> dict:
    outer_bytes = source_archive.read_bytes()
    with zipfile.ZipFile(io.BytesIO(outer_bytes)) as outer:
        names = set(outer.namelist())
        for required in (ORIGINAL_ZIP, SORTED_ZIP, SUMMARY_XLS):
            if required not in names:
                raise ValueError(f"source archive is missing {required}")
        original_bytes = outer.read(ORIGINAL_ZIP)
        sorted_bytes = outer.read(SORTED_ZIP)
        summary_bytes = outer.read(SUMMARY_XLS)

    originals_by_fp: dict[str, list[dict]] = {}
    with zipfile.ZipFile(io.BytesIO(original_bytes)) as zf:
        for name in sorted(zf.namelist()):
            if not name.upper().endswith(".AT2"):
                continue
            data = zf.read(name)
            parsed = parse_at2(data)
            relative = str(PurePosixPath(name).relative_to(PurePosixPath(name).parts[0]))
            item = {
                "original_file": relative,
                "original_file_sha256": sha256_bytes(data),
                **{k: parsed[k] for k in ("header_database", "header_record", "header_units", "npts", "dt_s")},
            }
            originals_by_fp.setdefault(vector_fingerprint(parsed["values"]), []).append(item)

    records = []
    with zipfile.ZipFile(io.BytesIO(sorted_bytes)) as zf:
        for name in sorted(zf.namelist()):
            base = PurePosixPath(name).name
            match = SORTED_RE.fullmatch(base)
            if not match:
                continue
            record_id = match.group(1)
            data = zf.read(name)
            values = parse_sorted(data)
            fp = vector_fingerprint(values)
            candidates = originals_by_fp.get(fp, [])
            if not candidates:
                raise ValueError(f"no PEER AT2 source vector matches {base}")
            canonical, duplicates = canonical_candidate(record_id, candidates)
            if canonical["npts"] != len(values):
                raise ValueError(f"NPTS mismatch for {record_id}")
            records.append({
                "record_id": record_id,
                "pair_id": record_id[:-1],
                "component_index": int(record_id[-1]),
                "sorted_file": base,
                "sorted_file_sha256": sha256_bytes(data),
                "original_file": canonical["original_file"],
                "duplicate_original_file_candidates": duplicates,
                "original_file_sha256": canonical["original_file_sha256"],
                "peer_header": canonical["header_record"],
                "source_units": "g",
                "npts": canonical["npts"],
                "dt_s": canonical["dt_s"],
                "source_sample_span_s": (canonical["npts"] - 1) * canonical["dt_s"],
                "analysis_step_end_duration_s": canonical["npts"] * canonical["dt_s"],
                "pga_g": max(abs(v) for v in values),
            })

    if len(records) != 44:
        raise ValueError(f"expected 44 horizontal components, found {len(records)}")
    if len({r['record_id'] for r in records}) != 44:
        raise ValueError("record IDs are not unique")

    return {
        "schema_version": 1,
        "record_set": "FEMA P-695 / ATC-63 far-field set",
        "scope": "metadata-only Phase 9N record manifest; raw acceleration bytes excluded from public repository",
        "source_archive": {
            "file_name": source_archive.name,
            "sha256": sha256_bytes(outer_bytes),
            "summary_workbook_sha256": sha256_bytes(summary_bytes),
            "original_peer_archive_sha256": sha256_bytes(original_bytes),
            "sorted_archive_sha256": sha256_bytes(sorted_bytes),
        },
        "analysis_conversion": {
            "source_acceleration_units": "g",
            "quakecore_length_units": "in",
            "gravity_in_per_s2": 386.09,
            "formula": "analysis_acceleration_in_per_s2 = source_acceleration_g * 386.09 * dimensionless_scale",
            "phase9m_reproduction_scale": 0.03,
            "sample_convention": "step_end",
        },
        "component_count": len(records),
        "pair_count": len({r['pair_id'] for r in records}),
        "records": records,
        "claim_boundary": "record provenance and software-validation input metadata only; not FEMA P-695 collapse qualification or acceptance",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    manifest = build_manifest(args.source_archive)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {len(manifest['records'])} records / {manifest['pair_count']} pairs to {args.out}")


if __name__ == "__main__":
    main()
