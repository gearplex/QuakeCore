#!/usr/bin/env python3
"""Render the deterministic Phase 9K steel verification histories."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt


def load(path):
    with open(path, encoding="utf-8") as stream:
        return json.load(stream)


def values(rows, key, index=None):
    result = [row[key] for row in rows]
    return result if index is None else [row[index] for row in result]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    qc = load(args.comparison / "quakecore_component_probe.json")
    ref = load(args.comparison / "opensees_component_reference.json")
    qc_dyn = load(args.comparison / "damper_sdof_quakecore.json")["runs"][0]["history"]
    ref_dyn = load(args.comparison / "damper_sdof_opensees.json")["history"]

    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6), constrained_layout=True)

    pairs = [
        ("steel_member", "rotation", "end_moment", 0, "Steel member end i", "Chord rotation (rad)", "Moment"),
        ("panel_zone", "deformation", "force", None, "Panel-zone surrogate", "Relative rotation (rad)", "Moment"),
        ("brb", "deformation", "force", None, "BRB", "Axial deformation", "Axial force"),
    ]
    for ax, (name, xkey, ykey, index, title, xlabel, ylabel) in zip(axes[0], pairs):
        ax.plot(values(ref[name], xkey), values(ref[name], ykey, index), lw=2.4, color="#E69F00", label="OpenSees")
        ax.plot(values(qc[name], xkey), values(qc[name], ykey, index), lw=1.0, ls="--", color="#0072B2", label="QuakeCore")
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)

    ax = axes[1, 0]
    ax.plot(values(ref["viscous_damper"], "velocity"), values(ref["viscous_damper"], "force"), "o", ms=6, color="#E69F00", label="OpenSees")
    ax.plot(values(qc["viscous_damper"], "velocity"), values(qc["viscous_damper"], "force"), "x", ms=6, color="#0072B2", label="QuakeCore")
    ax.set(title="Power-law viscous damper", xlabel="Axial velocity", ylabel="Force")

    ax = axes[1, 1]
    ax.plot(values(ref_dyn, "time_s"), values(ref_dyn, "displacement"), lw=2.4, color="#E69F00", label="OpenSees")
    ax.plot(values(qc_dyn, "time_s"), values(qc_dyn, "floor_displacement", 0), lw=1.0, ls="--", color="#0072B2", label="QuakeCore")
    ax.set(title="Linear-damper SDOF", xlabel="Time (s)", ylabel="Relative displacement")

    ax = axes[1, 2]
    differences = [abs(a["displacement"] - b["floor_displacement"][0]) for a, b in zip(ref_dyn, qc_dyn)]
    ax.semilogy(values(ref_dyn, "time_s"), [max(v, 1e-18) for v in differences], color="#009E73")
    ax.set(title="SDOF absolute difference", xlabel="Time (s)", ylabel="|uOS - uQC|")

    axes[0, 0].legend(loc="upper left", frameon=False)
    fig.suptitle("QuakeCore Phase 9K steel-component verification", fontsize=14)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
