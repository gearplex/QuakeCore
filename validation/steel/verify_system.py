#!/usr/bin/env python3
"""Integrated steel-system NRHA, time-step refinement, and IDA checks."""
import argparse
import copy
import json
import pathlib
import subprocess


def write(path, obj): pathlib.Path(path).write_text(json.dumps(obj, indent=2) + "\n")


def run(runner, job, path):
    inp = path.with_suffix(".input.json"); write(inp, job)
    p = subprocess.run([runner, inp, path], capture_output=True, text=True)
    if p.returncode not in (0,): raise RuntimeError(f"runner failed: {p.stderr}")
    return json.loads(path.read_text())


def max_history_diff(a, b, field):
    return max(abs(x[field][0] - y[field][0]) for x, y in zip(a, b))


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--runner",required=True);ap.add_argument("--example",required=True);ap.add_argument("--out",required=True)
    args=ap.parse_args();out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
    base=json.loads(pathlib.Path(args.example).read_text())

    full=copy.deepcopy(base);full["name"]="Phase 9K integrated nonlinear steel system";full["records"][0]["acceleration"]=[100*x for x in full["records"][0]["acceleration"]];full["analysis"]["strategy"]="full"
    direct=run(args.runner,full,out/"integrated_full.json")
    wood=copy.deepcopy(full);wood["analysis"]["strategy"]="woodbury";woodbury=run(args.runner,wood,out/"integrated_woodbury.json")
    a=direct["runs"][0];b=woodbury["runs"][0]

    refined=copy.deepcopy(full);refined["analysis"]["strategy"]="woodbury";old=refined["records"][0]["acceleration"];new=[];previous=0.0
    for current in old:new.extend([.5*(previous+current),current]);previous=current
    refined["records"][0]["acceleration"]=new;refined["records"][0]["dt"]*=.5
    fine=run(args.runner,refined,out/"integrated_half_dt.json")["runs"][0]

    ida=copy.deepcopy(base);ida["name"]="Phase 9K steel-system synthetic IDA";ida["analysis"]={"type":"ida","strategy":"woodbury","tolerance":1e-8,"relative_tolerance":True,"initial_guess":"kinematic","max_iterations":40,"max_subdivisions":6,"line_search":True,"scales":[20,50,100,200],"refinements":2,"workers":1,"drift_limit":.02,"stop_after_first_collapse":True}
    ida_result=run(args.runner,ida,out/"integrated_ida.json")

    summary={
      "schema":"quakecore.phase9k.integrated-steel-validation.v1",
      "full_termination":a["termination"],"woodbury_termination":b["termination"],"half_dt_termination":fine["termination"],
      "full_vs_woodbury_max_displacement_error":max_history_diff(a["history"],b["history"],"floor_displacement"),
      "full_vs_woodbury_max_drift_error":max_history_diff(a["history"],b["history"],"story_drift_ratio"),
      "peak_drift_dt":a["peak_story_drift_ratio"][0],"peak_drift_half_dt":fine["peak_story_drift_ratio"][0],
      "peak_drift_time_step_relative_change":abs(fine["peak_story_drift_ratio"][0]-a["peak_story_drift_ratio"][0])/fine["peak_story_drift_ratio"][0],
      "nonlinear_activation":{
        "brb_yielded":a["brbs"]["301"]["peak_abs_axial_force"]>=500.0,
        "panel_zone_yielded":a["panel_zones"]["201"]["peak_abs_moment"]>=350.0,
        "column_hinge_yielded":a["steel_members"]["101"]["peak_abs_end_moment_1"]>=450.0,
        "damper_peak_force":a["viscous_dampers"]["401"]["peak_abs_force"]},
      "ida_record":ida_result["records"][0],"ida_run_count":len(ida_result["runs"]),
      "scope":"synthetic numerical integration and configured drift-threshold exercise; not collapse-capacity validation"
    }
    write(out/"system_summary.json",summary);print(json.dumps(summary,indent=2))


if __name__=="__main__":main()
