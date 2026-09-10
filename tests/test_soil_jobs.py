#!/usr/bin/env python3
import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path

runner=Path(sys.argv[1]);source=Path(sys.argv[2]);example=source/'examples'/'soil'/'single_pile_nrha.json'

def run(job,expect=0):
    with tempfile.TemporaryDirectory() as td:
        inp=Path(td)/'job.json';out=Path(td)/'result.json';inp.write_text(json.dumps(job))
        p=subprocess.run([str(runner),str(inp),str(out)],capture_output=True,text=True)
        if p.returncode!=expect:raise RuntimeError(f'expected rc {expect}, got {p.returncode}: {p.stderr}')
        return json.loads(out.read_text()) if out.exists() else p.stderr

job=json.loads(example.read_text());result=run(job);one=result['runs'][0]
assert result['engine']=='QuakeCore Phase 9L soil springs and piles research'
assert one['termination']=='completed' and len(one['soil_springs'])==11 and len(one['soil_dashpots'])==1
assert len(one['pile_lines']['P1']['depth'])==5 and len(one['pile_lines']['P1']['peak_abs_member_moment_i'])==4
assert len(one['history'])==len(job['records'][0]['acceleration'])
assert one['history'][0]['pile_lines']['P1']['node_ids']==[1,2,3,4,5]
assert result['soil_scope']['acceptance_criteria']=='not_assessed'

ida=copy.deepcopy(job);ida['analysis']={"type":"ida","strategy":"woodbury","tolerance":1e-9,"relative_tolerance":False,"max_iterations":40,"max_subdivisions":4,"line_search":True,"scales":[0.25,0.5],"refinements":1,"workers":1,"drift_limit":0.5,"stop_after_first_collapse":True}
ida_result=run(ida);assert len(ida_result['runs'])==2 and ida_result['records'][0]['right_censored']

bad=copy.deepcopy(job);bad['model']['pile_lines'][0]['lateral']['provenance']=''
assert 'provenance' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['pile_lines'][0]['pile_node_ids']=[1,2]
assert 'wrong length' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['pile_lines'][0]['toe']['suction_ratio']=1.1
assert 'suction_ratio' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['pile_lines'][0]['lateral']['hardening_ratio']=0.01
assert 'ultimate force cap' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['soil_springs']=[{"id":900,"i":101,"j":2,"direction":[1,0],"type":"bilinear","k":10,"fy":1,"provenance":"negative test"}]
assert 'coincident' in run(bad,2)
print('soil spring, dashpot, pile compiler, recorder, and validation contracts passed')
