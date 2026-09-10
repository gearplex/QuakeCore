#!/usr/bin/env python3
"""Phase 9K independent steel-component verification against OpenSeesPy 3.8.

The script compares four mechanisms separately so a mismatch has a physical
location: condensed steel-member end hinges, panel-zone/BRB Steel01 histories,
power-law dashpot force, and a linear-damper SDOF transient response.
"""
import argparse
import json
import math
import os
import pathlib
import subprocess

os.environ.setdefault("MPIR_CVAR_CH3_INTERFACE_HOSTNAME", "127.0.0.1")
os.environ.setdefault("MPICH_INTERFACE_HOSTNAME", "127.0.0.1")
import openseespy.opensees as ops


def write(path, value):
    pathlib.Path(path).write_text(json.dumps(value, indent=2) + "\n")


def steel01_history(items, fy, stiffness, hardening):
    ops.wipe()
    ops.uniaxialMaterial("Steel01", 1, fy, stiffness, hardening)
    ops.testUniaxialMaterial(1)
    out = []
    for item in items:
        x = item["deformation"]
        ops.setStrain(x)
        out.append({"deformation": x, "force": ops.getStress(), "tangent": ops.getTangent()})
    return out


def expanded_member_history(qc):
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
    for node, x in ((1, 0.0), (2, 0.0), (3, 6.0), (4, 6.0)):
        ops.node(node, x, 0.0)
    for node in (1, 2, 3):
        ops.fix(node, 1, 1, 0)
    ops.fix(4, 1, 1, 1)
    ops.uniaxialMaterial("Steel01", 1, 300.0, 5.0e6, 0.02)
    ops.uniaxialMaterial("Steel01", 2, 300.0, 5.0e6, 0.02)
    # QuakeCore hinge deformation is outer rotation minus elastic-member end
    # rotation, hence inner-to-outer zeroLength ordering at both ends.
    ops.element("zeroLength", 101, 2, 1, "-mat", 1, "-dir", 3)
    ops.element("zeroLength", 102, 3, 4, "-mat", 2, "-dir", 3)
    ops.geomTransf("Linear", 1)
    ops.element("elasticBeamColumn", 103, 2, 3, 0.02, 2.0e8, 8.0e-4, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, 1.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0); ops.analysis("Static")
    previous = 0.0; out = []
    for index, item in enumerate(qc):
        target = item["rotation"]; increment = target - previous
        # Small internal reference increments avoid an OpenSees displacement-
        # control reversal failure without changing recorded protocol points.
        subdivisions = max(1, int(abs(increment) / 5e-5) + 1)
        for _ in range(subdivisions):
            ops.integrator("DisplacementControl", 1, 3, increment / subdivisions)
            if ops.analyze(1) != 0:
                raise RuntimeError(f"OpenSees expanded steel member failed at point {index}")
        previous = target
        local_force = ops.eleResponse(103, "localForce")
        out.append({
            "rotation": ops.nodeDisp(1, 3),
            "end_moment": [local_force[2], local_force[5]],
            "hinge_rotation": [ops.eleResponse(101, "material", 1, "strain")[0],
                               ops.eleResponse(102, "material", 1, "strain")[0]],
        })
    return out


def damper_force_history(items):
    ops.wipe(); ops.uniaxialMaterial("Viscous", 1, 125.0, 0.5); ops.testUniaxialMaterial(1)
    out = []
    for item in items:
        ops.setStrain(0.0, item["velocity"])
        out.append({"velocity": item["velocity"], "force": ops.getStress()})
    return out


def sdof_job(n=600, dt=0.01):
    acceleration = [0.8 * math.sin(0.031 * i) + 0.15 * math.sin(0.19 * i) for i in range(1, n + 1)]
    return {
        "schema": "quakecore.job.v1", "name": "linear viscous damper SDOF verification",
        "units": {"force": "kN", "length": "m", "time": "s"},
        "model": {"type": "frame2d", "nodes": [[1, 0, 1, 0, 0, 0], [2, 2, 1, 5, 0, 0]],
                  "fixities": [[1, True, True, True], [2, False, True, True]], "equal_dofs": [],
                  "members": [[1, 1, 2, 1000, 2, 1, 0]], "hinges": [], "walls": [],
                  "steel_members": [], "panel_zones": [], "brbs": [],
                  "viscous_dampers": [{"id": 2, "i": 1, "j": 2, "coefficient": 3.5,
                                        "alpha": 1.0, "regularization_velocity": 0.0,
                                        "provenance": "linear verification device"}],
                  "rayleigh": [0, 0], "response_node": 2, "story_nodes": [2], "story_cut_members": []},
        "analysis": {"type": "nrha", "strategy": "woodbury", "tolerance": 1e-10,
                     "relative_tolerance": False, "initial_guess": "kinematic", "max_iterations": 40,
                     "max_subdivisions": 4, "line_search": True, "history_stride": 1,
                     "output_mode": "full", "drift_limit": 0},
        "records": [{"name": "synthetic", "dt": dt, "sample_convention": "step_end",
                     "acceleration": acceleration, "provenance": "synthetic verification input"}],
        "provenance": {"purpose": "software verification, not physical validation"},
    }


def opensees_sdof(acceleration, dt):
    ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.node(2, 2.0); ops.fix(1, 1); ops.mass(2, 5.0)
    ops.uniaxialMaterial("Elastic", 1, 2000.0); ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.uniaxialMaterial("Viscous", 2, 3.5, 1.0); ops.element("twoNodeLink", 2, 1, 2, "-mat", 2, "-dir", 1)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *([0.0] + acceleration))
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-10, 30); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    out = []
    for i, ag in enumerate(acceleration):
        if ops.analyze(1, dt) != 0:
            raise RuntimeError(f"OpenSees SDOF failed at step {i}")
        out.append({"time_s": (i + 1) * dt, "displacement": ops.nodeDisp(2, 1),
                    "velocity": ops.nodeVel(2, 1), "relative_acceleration": ops.nodeAccel(2, 1),
                    "ground_acceleration": ag})
    return out


def max_abs_diff(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--probe", required=True); ap.add_argument("--runner", required=True); ap.add_argument("--out", required=True)
    args = ap.parse_args(); outdir = pathlib.Path(args.out); outdir.mkdir(parents=True, exist_ok=True)
    probe_path = outdir / "quakecore_component_probe.json"
    subprocess.run([args.probe, probe_path], check=True)
    qc = json.loads(probe_path.read_text())
    reference = {
        "engine": f"OpenSeesPy {ops.version()}",
        "panel_zone": steel01_history(qc["panel_zone"], 50.0, 5000.0, 0.02),
        "brb": steel01_history(qc["brb"], 20.0, 1000.0, 0.01),
        "steel_member": expanded_member_history(qc["steel_member"]),
        "viscous_damper": damper_force_history(qc["viscous_damper"]),
    }
    write(outdir / "opensees_component_reference.json", reference)

    job = sdof_job(); write(outdir / "damper_sdof_input.json", job)
    quake_result_path = outdir / "damper_sdof_quakecore.json"
    subprocess.run([args.runner, outdir / "damper_sdof_input.json", quake_result_path], check=True)
    quake_result = json.loads(quake_result_path.read_text())
    os_dynamic = opensees_sdof(job["records"][0]["acceleration"], job["records"][0]["dt"])
    write(outdir / "damper_sdof_opensees.json", {"engine": f"OpenSeesPy {ops.version()}", "history": os_dynamic})
    qdyn = quake_result["runs"][0]["history"]

    summary = {
        "schema": "quakecore.phase9k.steel-validation.v1", "opensees_version": ops.version(),
        "component_steps": len(qc["panel_zone"]), "steel_member_steps": len(qc["steel_member"]),
        "panel_zone_max_force_error": max_abs_diff([x["force"] for x in qc["panel_zone"]], [x["force"] for x in reference["panel_zone"]]),
        "panel_zone_max_tangent_error": max_abs_diff([x["tangent"] for x in qc["panel_zone"]], [x["tangent"] for x in reference["panel_zone"]]),
        "brb_max_force_error": max_abs_diff([x["force"] for x in qc["brb"]], [x["force"] for x in reference["brb"]]),
        "brb_max_tangent_error": max_abs_diff([x["tangent"] for x in qc["brb"]], [x["tangent"] for x in reference["brb"]]),
        "steel_member_max_moment_error": max(abs(qc["steel_member"][i]["end_moment"][j] - reference["steel_member"][i]["end_moment"][j]) for i in range(len(qc["steel_member"])) for j in range(2)),
        "steel_member_max_hinge_rotation_error": max(abs(qc["steel_member"][i]["hinge_rotation"][j] - reference["steel_member"][i]["hinge_rotation"][j]) for i in range(len(qc["steel_member"])) for j in range(2)),
        "power_law_damper_max_force_error": max_abs_diff([x["force"] for x in qc["viscous_damper"]], [x["force"] for x in reference["viscous_damper"]]),
        "linear_damper_sdof_max_displacement_error": max_abs_diff([x["floor_displacement"][0] for x in qdyn], [x["displacement"] for x in os_dynamic]),
        "linear_damper_sdof_max_velocity_error": max_abs_diff([x["viscous_dampers"]["2"]["deformation_rate"] for x in qdyn], [x["velocity"] for x in os_dynamic]),
        "quakecore_dynamic_termination": quake_result["runs"][0]["termination"],
        "scope": {"steel_member": "expanded elasticBeamColumn plus two Steel01 zeroLength hinges",
                  "panel_zone_and_brb": "Steel01 material-history reference; component kinematics checked in C++ tests",
                  "damper": "OpenSees Viscous force law and independent linear SDOF transient",
                  "physical_test_validation": False, "code_acceptance_validation": False},
    }
    write(outdir / "summary.json", summary); print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
