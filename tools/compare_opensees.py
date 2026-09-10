"""Independent OpenSees 3.x comparison for v1 bilinear plane-frame jobs.

Both engines use identical nodes/members/masses, Steel01-compatible kinematic
hardening, initial-stiffness Rayleigh damping including zeroLength hinges,
Newmark (gamma=.5,beta=.25), and absolute infinity-norm residual convergence.
OpenSees advances the whole record in native code. JSON serialization and
history comparison are outside its analysis timer.
"""
import argparse, json, math, os, subprocess, time
from pathlib import Path
import numpy as np

def compare(job_path, quake_exe, outdir, repeats=5, timing_recorders=True):
    if repeats<1:raise ValueError("repeats must be positive")
    import openseespy.opensees as ops
    job=json.loads(Path(job_path).read_text());m=job['model'];a=job['analysis'];record=job['records'][0]
    if len(job['records'])!=1 or a['type']!='nrha' or a.get('relative_tolerance',True):
        raise ValueError('parity requires one NRHA record and an absolute tolerance')
    if a.get('max_subdivisions',0)!=0 or a.get('line_search',True) or a.get('initial_guess')!='previous_displacement':
        raise ValueError('parity requires fixed dt and full Newton without line search')
    if any(e[6]!=0 for e in m['members']) or any(h['type']!='bilinear' for h in m['hinges']):
        raise ValueError('this adapter supports bilinear hinges and zero geometric preload only')
    p=Path(outdir);p.mkdir(parents=True,exist_ok=True)
    # Preserve values at the time-step endpoints; no resampling or scaling.
    ag=record['acceleration'];dt=record['dt'];floors=m['story_nodes'];coords={n[0]:n[1:3] for n in m['nodes']}
    eq={}
    for master,slave,dof in m.get('equal_dofs',[]):eq.setdefault((master,slave),[]).append({'UX':1,'UY':2,'RZ':3}[dof])
    if any(slave==master2 for master,slave in eq for master2,slave2 in eq):
        raise ValueError('chained MPCs require a flattened reference adapter')
    times=[];periods=[]
    for repeat in range(repeats+1):
        ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
        for tag,x,y,mx,my,mr in m['nodes']:ops.node(tag,x,y);ops.mass(tag,mx,my,mr)
        for tag,x,y,r in m['fixities']:ops.fix(tag,int(x),int(y),int(r))
        for (master,slave),dofs in eq.items():ops.equalDOF(master,slave,*dofs)
        ops.geomTransf('Linear',1)
        for tag,i,j,E,A,I,preload in m['members']:ops.element('elasticBeamColumn',tag,i,j,A,E,I,1)
        for h in m['hinges']:
            ops.uniaxialMaterial('Steel01',h['id'],h['fy'],h['k'],h['hardening_ratio'])
            ops.element('zeroLength',h['id'],h['i'],h['j'],'-mat',h['id'],'-dir',3,'-doRayleigh',1)
        ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack')
        eig=ops.eigen(min(3,len(floors)));periods=[2*math.pi/math.sqrt(v) for v in eig]
        alpha,beta=m['rayleigh'];ops.rayleigh(alpha,0,beta,0)
        ops.timeSeries('Path',1,'-dt',dt,'-values',0.0,*ag);ops.pattern('UniformExcitation',1,1,'-accel',1)
        disp=p/'opensees_displacement.txt';acc=p/'opensees_acceleration.txt'
        hfile=p/'opensees_hinges.txt'
        if repeat==0 or timing_recorders:
            ops.recorder('Node','-file',str(disp),'-precision',16,'-time','-node',*floors,'-dof',1,'disp')
            ops.recorder('Node','-file',str(acc),'-precision',16,'-time','-node',*floors,'-dof',1,'accel')
            ops.recorder('Element','-file',str(hfile),'-precision',16,'-time','-ele',*[h['id'] for h in m['hinges']],'deformation')
        ops.test('NormUnbalance',a['tolerance'],a['max_iterations'],0,0)
        ops.algorithm('Newton');ops.integrator('Newmark',.5,.25);ops.analysis('Transient')
        start=time.perf_counter();rc=ops.analyze(len(ag),dt);elapsed=time.perf_counter()-start
        if repeat>0:times.append(elapsed)
        ops.remove('recorders')
        if rc!=0:raise RuntimeError(f'OpenSees failed rc={rc} at {ops.getTime()}')
    qtimes=[];r=None
    quiet_job=p/'timing_job.json';timing_job=json.loads(json.dumps(job));timing_job['analysis']['output_mode']='full' if timing_recorders else 'minimal';quiet_job.write_text(json.dumps(timing_job))
    for repeat in range(repeats+1):
        qpath=p/('quakecore.json' if repeat==0 else 'quakecore_timing.json');run=subprocess.run([str(quake_exe),str(job_path if repeat==0 else quiet_job),str(qpath)],capture_output=True,text=True)
        if run.returncode:raise RuntimeError(f'QuakeCore return code {run.returncode}; see {qpath}; {run.stderr}')
        trial=json.loads(qpath.read_text())
        if repeat==0:r=trial
        else:qtimes.append(trial['runs'][0]['stats']['elapsed_seconds'])
    u=np.loadtxt(disp);accel=np.loadtxt(acc);q=r['runs'][0];qu=np.array([h['floor_displacement'] for h in q['history']]);qt=np.array([h['time_s'] for h in q['history']]);qa=np.array([h['floor_absolute_acceleration'] for h in q['history']]);
    if len(u)!=len(qu) or not np.allclose(u[:,0],qt,rtol=0,atol=1e-9):raise RuntimeError('history grids differ')
    heights=np.diff([0]+[coords[n][1] for n in floors]);od=np.diff(np.column_stack([np.zeros(len(u)),u[:,1:]]),axis=1)/heights
    qd=np.array([h['story_drift_ratio'] for h in q['history']]);oa=accel[:,1:]+np.array(ag)[:,None]
    qp=np.array([x['period_s'] for x in r['modes']]);op=np.array(periods)
    rotations=np.loadtxt(hfile)[:,1:]
    ductility=max(float(np.max(np.abs(rotations[:,i])))/(h['fy']/h['k']) for i,h in enumerate(m['hinges']))
    result={'schema':'quakecore.external_parity.v1','job':str(job_path),'opensees_version':ops.version(),'records':1,'steps':len(ag),'dof':r['dof'],
        'max_displacement_error':float(np.max(np.abs(qu-u[:,1:]))),'relative_displacement_error':float(np.max(np.abs(qu-u[:,1:]))/max(1e-15,np.max(np.abs(u[:,1:])))),
        'max_drift_ratio_error':float(np.max(np.abs(qd-od))),'max_absolute_acceleration_error':float(np.max(np.abs(qa-oa))),
        'periods_quakecore_s':qp.tolist(),'periods_opensees_s':op.tolist(),'max_period_relative_error':float(np.max(np.abs(qp-op)/op)),
        'peak_drifts_quakecore':np.max(np.abs(qd),axis=0).tolist(),'peak_drifts_opensees':np.max(np.abs(od),axis=0).tolist(),
        'max_hinge_rotation_ductility':ductility,'quakecore_seconds':qtimes,'opensees_seconds':times,'speedup_opensees_over_quakecore':float(np.median(times)/np.median(qtimes)),
        'timing_recorders_enabled':timing_recorders,'timing_scope':('integration with recording: OpenSees native text recorders; QuakeCore EDPs and JSON-history construction, final JSON serialization excluded' if timing_recorders else 'integration without optional recorders in both engines; QuakeCore retains built-in roof history, diagnostics and energy ledger; model creation/eigenanalysis excluded'),
        'scope':'synthetic plane frame, bilinear beam hinges, fixed dt, no P-Delta, no degrading RC/P-M/IMK parity'}
    result['passed']=result['relative_displacement_error']<1e-5 and result['max_drift_ratio_error']<1e-7 and result['max_period_relative_error']<1e-5 and result['max_absolute_acceleration_error']<1e-6 and ductility>1
    (p/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('job');p.add_argument('executable');p.add_argument('output');p.add_argument('--repeats',type=int,default=5);p.add_argument('--no-timing-recorders',action='store_true');a=p.parse_args()
    result=compare(a.job,a.executable,a.output,a.repeats,not a.no_timing_recorders)
    raise SystemExit(0 if result['passed'] else 1)
