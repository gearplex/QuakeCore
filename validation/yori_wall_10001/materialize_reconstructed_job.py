#!/usr/bin/env python3
"""Materialize the documented reconstructed archetype-10001 exact-law job.

This intentionally does NOT claim that the reconstructed 3x41 row set is the
lost received AllStoriesOpenseesModelINPUT.txt. It reuses the Gate-4 exact-law
materializer after parsing the surrogate rows, then replaces provenance so the
result cannot be mistaken for SHA-verified received input.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from materialize_received_job import materialize


def load_surrogate_rows(path: Path) -> tuple[list[list[float]], str]:
    raw = path.read_bytes()
    values = [float(x) for x in raw.decode().split()]
    if len(values) != 3 * 41:
        raise ValueError(f"expected exactly 3 x 41 model values, received {len(values)}")
    return [values[i * 41:(i + 1) * 41] for i in range(3)], hashlib.sha256(raw).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-input", type=Path, required=True)
    parser.add_argument(
        "--template", type=Path,
        default=Path(__file__).with_name("quakecore_job_10001_gravity.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows, digest = load_surrogate_rows(args.model_input)
    template = json.loads(args.template.read_text())
    result = materialize(template, rows)
    result["name"] = "YORi_10001_reconstructed_materials_gravity_pdelta"
    provenance = result.setdefault("provenance", {})
    provenance.pop("received_model_input_sha256", None)
    provenance["reconstructed_model_input_sha256"] = digest
    provenance["purpose"] = (
        "Gate 4 surrogate all-story QuakeCore/OpenSees parity test using reconstructed "
        "active constitutive rows; software/model verification only"
    )
    provenance["material_substitutions"] = (
        "ConcreteCM + Pinching4 + MinMax + Parallel reconstructed from frozen story-1 "
        "oracle, retained matched-law parameters, ASCE modeling workbook, and recovered Tcl; "
        "MinMax limits explicitly assumed +/-0.075"
    )
    provenance["claim_boundary"] = (
        "Does not establish identity to the lost received source file, physical validation, "
        "FEMA P-695 acceptance, collapse qualification, or an R-factor recommendation"
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(args.output)
    print(digest)


if __name__ == "__main__":
    main()
