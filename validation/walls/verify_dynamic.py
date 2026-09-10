#!/usr/bin/env python3
"""Three-story wall/frame NRHA and IDA integration checks, synthetic motions."""
import argparse,json,pathlib,subprocess,copy,sys
import numpy as np
from compare_opensees import fixture,ops,add_wall,Tags,STEEL

def job(kind,nonlinear,layered=False):
    n=3;dt=.005;t=np.arange(1,1201)*dt
    ag=(1.5*np.sin(np.pi*t/6)**2*(np.sin(2*np.pi*1.6*t)+.35*np.sin(2*np.pi*4.7*t))).tolist()
    nodes=[[1,0,0,0,0,0],[101,6,0,0,0,0]];walls=[];members=[];ties=[]
    for i in range(1,n+1):
        nodes.extend([[i+1,0,3*i,50,20,0],[101+i,6,3*i,10,5,0]])
        w=fixture(kind,nonlinear);w.update(id=i,i=i,j=i+1)
        if layered:
            for k,f in enumerate(w['panels']):
                angle=30 if k%2==0 else -30
                f['material']=dict(type='layered_plane_stress',background_E=0.,background_nu=0.,layers=[dict(angle_deg=angle,weight=.1,material=STEEL),dict(angle_deg=angle+90,weight=.15,material=STEEL)])
        walls.append(w);members.append([100+i,100+i,101+i,2.5e7,.25,.01,0]);ties.append([i+1,101+i,'UX'])
    return dict(schema='quakecore.job.v1',name=f'{kind}_three_story_wall_frame',units=dict(force='kN',length='m',time='s'),provenance=dict(purpose='synthetic implementation verification, not a building assessment',gravity='zero; no gravity equilibration in this NRHA example'),model=dict(type='frame2d',nodes=nodes,fixities=[[1,True,True,True],[101,True,True,True]],equal_dofs=ties,members=members,hinges=[],walls=walls,rayleigh=[.05,.001],response_node=4,story_nodes=[2,3,4]),analysis=dict(type='nrha',strategy='same_pattern',tolerance=1e-8,relative_tolerance=False,initial_guess='previous_displacement',max_iterations=80,max_subdivisions=4,line_search=True,output_mode='full'),records=[dict(name='synthetic',dt=dt,sample_convention='step_end',acceleration=ag)])

def reference(j):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3);m=j['model'];tags=Tags()
    for n in m['nodes']:ops.node(n[0],n[1],n[2]);ops.mass(n[0],*n[3:])
    for f in m['fixities']:ops.fix(f[0],*[int(x) for x in f[1:]])
    for a,b,d in m['equal_dofs']:ops.equalDOF(a,b,1)
    ops.geomTransf('Linear',1)
    for tag,i,k,E,A,I,P in m['members']:ops.element('elasticBeamColumn',tag,i,k,A,E,I,1)
    for w in m['walls']:add_wall(w,tags)
    ops.rayleigh(m['rayleigh'][0],0.,m['rayleigh'][1],0.)
    r=j['records'][0];ops.timeSeries('Path',1,'-dt',r['dt'],'-values',0.,*r['acceleration'],'-useLast');ops.pattern('UniformExcitation',1,1,'-accel',1)
    ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack');ops.test('NormUnbalance',1e-8,150,0,0)
    if m['walls'][0]['type']=='sfi_mvlem':ops.algorithm('KrylovNewton','-maxDim',12)
    else:ops.algorithm('Newton')
    ops.integrator('Newmark',.5,.25);ops.analysis('Transient');hist=[]
    for i in range(len(r['acceleration'])):
        if ops.analyze(1,r['dt'])!=0:raise RuntimeError(f'OpenSees transient failed at step {i}')
        hist.append(dict(time_s=ops.getTime(),floor_displacement=[ops.nodeDisp(n,1) for n in m['story_nodes']],floor_absolute_acceleration=[ops.nodeAccel(n,1)+r['acceleration'][i] for n in m['story_nodes']],wall_force=[ops.eleForce(w['id'])[:6] for w in m['walls']]))
    return hist

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--runner',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);ap.add_argument('--examples',type=pathlib.Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);summary={}
    cases=[('mvlem_rc',job('mvlem',True),True),('sfi_elastic',job('sfi_mvlem',False),True),('sfi_two_layer',job('sfi_mvlem',True,True),True),('sfi_fixed_angle_rc',job('sfi_mvlem',True),False)]
    for name,j,compare in cases:
        ip=a.out/(name+'_job.json');op=a.out/(name+'_quake.json');ip.write_text(json.dumps(j,indent=2));subprocess.run([str(a.runner),str(ip),str(op)],check=True)
        q=json.loads(op.read_text());run=q['runs'][0];metrics=dict(termination=run['termination'],max_drift=max(run['peak_story_drift_ratio']),elapsed_seconds=run['stats']['elapsed_seconds'] if 'stats' in run else None,passed=run['termination']=='completed',independent_comparison=compare)
        if compare:
            r=reference(j);(a.out/(name+'_opensees.json')).write_text(json.dumps(r));u=np.array([v['floor_displacement'] for v in run['history']]);ur=np.array([v['floor_displacement'] for v in r]);acc=np.array([v['floor_absolute_acceleration'] for v in run['history']]);ar=np.array([v['floor_absolute_acceleration'] for v in r]);fq=np.array([[v['walls'][str(w)]['nodal_force'] for w in range(1,4)] for v in run['history']]);fr=np.array([v['wall_force'] for v in r]);metrics.update(max_displacement_error=float(abs(u-ur).max()),max_acceleration_error=float(abs(acc-ar).max()),max_wall_force_error=float(abs(fq-fr).max()));metrics['passed']=metrics['passed'] and metrics['max_displacement_error']<2e-6 and metrics['max_acceleration_error']<2e-4 and metrics['max_wall_force_error']<.02
        summary[name]=metrics;print(name,metrics,flush=True);(a.out/'dynamic_summary.json').write_text(json.dumps(summary,indent=2))
        if a.examples and name in ['mvlem_rc','sfi_fixed_angle_rc']:
            a.examples.mkdir(parents=True,exist_ok=True);(a.examples/(name+'_nrha.json')).write_text(json.dumps(j,indent=2))
    for name in ['mvlem_rc','sfi_fixed_angle_rc']:
        j=copy.deepcopy(next(j for n,j,_ in cases if n==name));j['analysis'].update(type='ida',scales=[.5,1,2],refinements=1,workers=2,drift_limit=.004);j['records'].append(copy.deepcopy(j['records'][0]));j['records'][1]['name']='synthetic_reverse';j['records'][1]['acceleration']=[-x for x in j['records'][1]['acceleration']]
        ip=a.out/(name+'_ida_job.json');op=a.out/(name+'_ida_quake.json');ip.write_text(json.dumps(j,indent=2));subprocess.run([str(a.runner),str(ip),str(op)],check=True);res=json.loads(op.read_text());summary[name+'_ida']=dict(runs=len(res['runs']),passed=all(x['termination'] in ['completed','configured_collapse_criterion'] for x in res['runs']),records=res['records'])
        if a.examples:(a.examples/(name+'_ida.json')).write_text(json.dumps(j,indent=2))
    (a.out/'dynamic_summary.json').write_text(json.dumps(summary,indent=2));sys.exit(0 if all(v['passed'] for v in summary.values()) else 1)
if __name__=='__main__':main()
