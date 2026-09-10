#!/usr/bin/env python3
import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path

runner=Path(sys.argv[1]);source=Path(sys.argv[2])
example=source/'examples'/'steel'/'steel_frame_nrha.json'

def run(job,expect=0):
    with tempfile.TemporaryDirectory() as td:
        inp=Path(td)/'job.json';out=Path(td)/'result.json'
        inp.write_text(json.dumps(job))
        p=subprocess.run([str(runner),str(inp),str(out)],capture_output=True,text=True)
        if p.returncode!=expect:
            raise RuntimeError(f'expected rc {expect}, got {p.returncode}: {p.stderr}')
        return json.loads(out.read_text()) if out.exists() else p.stderr

job=json.loads(example.read_text())
result=run(job)
assert result['engine']=='QuakeCore Phase 9L soil springs and piles research'
assert result['runs'][0]['termination']=='completed'
for key in ('steel_members','panel_zones','brbs','viscous_dampers'):
    assert result['runs'][0][key],key
    assert result['runs'][0]['history'][0][key],key+' history'
assert len(result['runs'][0]['peak_story_restoring_shear'])==1
assert result['steel_scope']['acceptance_criteria']=='not_assessed'

bad=copy.deepcopy(job);bad['model']['viscous_dampers'][0]['regularization_velocity']=0.0
assert 'regularization' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['panel_zones'][0]['j']=4
assert 'coincident' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['steel_members'][0]['provenance']=''
assert 'provenance' in run(bad,2)
bad=copy.deepcopy(job);bad['model']['brbs'][0]['hardening_ratio']=1.0
assert 'bilinear' in run(bad,2).lower()
print('steel JSON integration and validation contracts passed')
