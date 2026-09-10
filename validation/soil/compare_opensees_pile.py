#!/usr/bin/env python3
"""Independent OpenSeesPy comparison for the Phase 9L vertical pile fixture."""
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import openseespy.opensees as ops

source=Path(__file__).resolve().parents[2]
runner=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else source.parent/'build'/'quake_run'
output=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else None
job_path=source/'examples'/'soil'/'single_pile_nrha.json';job=json.loads(job_path.read_text())

with tempfile.TemporaryDirectory() as td:
    qc_path=Path(td)/'quakecore.json'
    subprocess.run([str(runner),str(job_path),str(qc_path)],check=True,capture_output=True,text=True)
    qc=json.loads(qc_path.read_text())

p=job['model']['pile_lines'][0];head=p['head_node'];n=p['segments']+1;dz=p['embedded_length']/p['segments']
depth=[i*dz for i in range(n)];trib=[dz/2]+[dz]*(n-2)+[dz/2]
ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
head_data=job['model']['nodes'][0];x0,y0=head_data[1],head_data[2]
for i,node in enumerate(p['pile_node_ids']):
    ops.node(node,x0,y0-depth[i]);ops.fix(node,0,1,0)
for i,node in enumerate(p['soil_node_ids']):
    ops.node(node,x0,y0-depth[i]);ops.fix(node,1,1,1)
ops.mass(head,head_data[3],0.0,0.0)
ops.geomTransf('Linear',1)
for tag,i,j in zip(p['member_ids'],p['pile_node_ids'][:-1],p['pile_node_ids'][1:]):
    ops.element('elasticBeamColumn',tag,i,j,p['A'],p['E'],p['I'],1)
lat=p['lateral'];b=lat['hardening_ratio']
for idx,(tag,anchor,node) in enumerate(zip(lat['spring_ids'],p['soil_node_ids'],p['pile_node_ids'])):
    pu=lat['ultimate_resistance_per_length'][idx]*trib[idx]
    k=0.5*pu/lat['displacement_50'][idx]
    mat=1000+idx;ops.uniaxialMaterial('Steel01',mat,pu,k,b)
    ops.element('zeroLength',tag,anchor,node,'-mat',mat,'-dir',1)
d=job['model']['soil_dashpots'][0];ops.uniaxialMaterial('Viscous',2000,d['coefficient'],d['alpha'])
ops.element('zeroLength',d['id'],d['i'],d['j'],'-mat',2000,'-dir',1)

# OpenSees Path returns zero at the right endpoint unless one additional value
# brackets it. QuakeCore defines every supplied sample at a step end, including
# the final sample, so repeat the last value solely to preserve that convention.
record=job['records'][0];values=[0.0]+record['acceleration']+[record['acceleration'][-1]]
ops.timeSeries('Path',1,'-dt',record['dt'],'-values',*values)
ops.pattern('UniformExcitation',1,1,'-accel',1)
ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-10,40)
ops.algorithm('Newton');ops.integrator('Newmark',0.5,0.25);ops.analysis('Transient')
opensees=[]
for _ in record['acceleration']:
    if ops.analyze(1,record['dt'])!=0:raise RuntimeError('OpenSees pile analysis failed')
    opensees.append(ops.nodeDisp(head,1))
quake=[step['pile_lines'][p['id']]['lateral_displacement'][0] for step in qc['runs'][0]['history']]
if len(quake)!=len(opensees):raise RuntimeError('history length mismatch')
diff=[abs(a-b) for a,b in zip(quake,opensees)];scale=max(max(abs(x) for x in opensees),1e-30)
imax=max(range(len(diff)),key=diff.__getitem__)
first=next((i for i,x in enumerate(diff) if x>1e-12),None)
summary={
    'schema':'quakecore.validation.phase9l.v1',
    'opensees_version':ops.version(),
    'quakecore_engine':qc['engine'],
    'fixture':'single vertical elastic pile; nodal bilinear p-y springs; linear local dashpot',
    'quakecore_termination':qc['runs'][0]['termination'],
    'steps':len(quake),
    'quakecore_peak_abs_head_displacement':max(abs(x) for x in quake),
    'opensees_peak_abs_head_displacement':max(abs(x) for x in opensees),
    'max_abs_history_difference':max(diff),
    'max_difference_step':imax+1,
    'first_difference_over_1e-12_step':None if first is None else first+1,
    'quakecore_at_max_difference':quake[imax],
    'opensees_at_max_difference':opensees[imax],
    'max_difference_over_opensees_peak':max(diff)/scale,
    'quakecore_newton_iterations':qc['runs'][0]['stats']['newton_iterations'],
    'quakecore_subdivided_steps':qc['runs'][0]['stats']['subdivided_steps'],
    'scope':'Defined bilinear tributary springs only; not PySimple1/TzSimple1/QzSimple1 and not a site-specific SSI validation'
}
text=json.dumps(summary,indent=2)+'\n'
if output:output.write_text(text)
print(text,end='')
if summary['max_difference_over_opensees_peak']>2e-8:raise SystemExit(1)
