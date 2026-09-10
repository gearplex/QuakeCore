#!/usr/bin/env python3
"""Independent wall verification. Requires OpenSeesPy; never a solver dependency.
python compare_opensees.py --probe build/wall_probe --runner build/quake_run --out evidence
"""
import argparse, json, pathlib, subprocess, os, copy, math
os.environ.setdefault('MPIR_CVAR_CH3_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('MPICH_INTERFACE_HOSTNAME','127.0.0.1')
import numpy as np
import openseespy.opensees as ops

CONC=dict(type='concrete01',fc=-30000.,epsc=-.002,fcu=-6000.,epsu=-.006)
STEEL=dict(type='steel_bilinear',E=2e8,fy=4e5,b=.01)
WIDTH=[.3,.7,.7,.3]

def fixture(kind,nonlinear):
    w=dict(id=1,type=kind,i=1,j=2,c=.4,density=0.)
    if kind=='mvlem':
        w['fibers']=[dict(width=b,thickness=.2,rho=(.025 if i in (0,3) else .01),concrete=CONC if nonlinear else dict(type='elastic',E=3e7),steel=STEEL if nonlinear else dict(type='elastic',E=2e8)) for i,b in enumerate(WIDTH)]
        w['shear']=dict(type='bilinear',stiffness=1e5,yield_force=200.,hardening_ratio=.02) if nonlinear else dict(type='elastic',stiffness=1e5)
    else:
        w['panels']=[dict(width=b,thickness=.2,material=dict(type='fixed_angle_rc',background_E=1e6,background_nu=.2,angle_deg=35.,concrete=CONC,steel_x=STEEL,steel_y=STEEL,rho_x=.01,rho_y=(.025 if i in (0,3) else .01)) if nonlinear else dict(type='elastic_plane_stress',E=3e7,nu=.2)) for i,b in enumerate(WIDTH)]
    return copy.deepcopy(w)

class Tags:
    def __init__(self):self.n=100
    def next(self):self.n+=1;return self.n

def uniaxial(p,tags):
    t=tags.next()
    if p['type']=='elastic':ops.uniaxialMaterial('Elastic',t,p['E'])
    elif p['type']=='steel_bilinear':ops.uniaxialMaterial('Steel01',t,p['fy'],p['E'],p['b'])
    elif p['type']=='concrete01':ops.uniaxialMaterial('Concrete01',t,p['fc'],p['epsc'],p['fcu'],p['epsu'])
    else:raise ValueError(p)
    return t

def panel(p,tags):
    if p['type']=='elastic_plane_stress':
        t=tags.next();ops.nDMaterial('ElasticIsotropic',t,p['E'],p['nu']);return t
    layers=[]
    if p['type']=='fixed_angle_rc':
        E,nu=p['background_E'],p['background_nu']
        layers=[dict(angle_deg=p['angle_deg'],weight=1,material=p['concrete']),dict(angle_deg=p['angle_deg']+90,weight=1,material=p['concrete']),dict(angle_deg=0,weight=p['rho_x'],material=p['steel_x']),dict(angle_deg=90,weight=p['rho_y'],material=p['steel_y'])]
    else:E,nu,layers=p['background_E'],p['background_nu'],p['layers']
    if E==0 and len(layers)==2:
        l1,l2=layers;u1=uniaxial(l1['material'],tags);u2=uniaxial(l2['material'],tags);t=tags.next()
        ops.nDMaterial('SmearedSteelDoubleLayer',t,u1,u2,l1['weight'],l2['weight'],math.radians(l1['angle_deg']));return t
    parts=[]
    if E>0:
        bg=tags.next();ops.nDMaterial('ElasticIsotropic',bg,E,nu);parts.extend([bg,1.])
    for l in layers:
        u=uniaxial(l['material'],tags);n=tags.next();ops.nDMaterial('PlaneStressRebarMaterial',n,u,l['angle_deg']);parts.extend([n,l['weight']])
    t=tags.next();ops.nDMaterial('PlaneStressLayeredMaterial',t,len(parts)//2,*parts);return t

def add_wall(w,tags):
    if w['type']=='mvlem':
        fs=w['fibers'];cs=[uniaxial(f['concrete'],tags) for f in fs];ss=[uniaxial(f['steel'],tags) for f in fs]
        sh=w['shear'];s=tags.next()
        if sh['type']=='elastic':ops.uniaxialMaterial('Elastic',s,sh['stiffness'])
        else:ops.uniaxialMaterial('Steel01',s,sh['yield_force'],sh['stiffness'],sh['hardening_ratio'])
        ops.element('MVLEM',w['id'],w.get('density',0.),w['i'],w['j'],len(fs),w['c'],'-thick',*[f['thickness'] for f in fs],'-width',*[f['width'] for f in fs],'-rho',*[f['rho'] for f in fs],'-matConcrete',*cs,'-matSteel',*ss,'-matShear',s)
    else:
        fs=w['panels'];ms=[panel(f['material'],tags) for f in fs]
        ops.element('SFI_MVLEM',w['id'],w['i'],w['j'],len(fs),w['c'],'-thick',*[f['thickness'] for f in fs],'-width',*[f['width'] for f in fs],'-mat',*ms,'-Coupling',2)

def build(w):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3);ops.node(1,0.,0.);ops.node(2,0.,3.);ops.fix(1,1,1,1);add_wall(w,Tags())

def protocol():
    p=[[0.,-300.*i/20,0.] for i in range(1,21)]
    for peak in [.003,-.003,.006,-.006,.012,-.012,.024,-.024,0]:
        start=p[-1][0];p.extend([[float(u),-300.,0.] for u in np.linspace(start,peak,81)[1:]])
    return p

def static_ref(w,p):
    build(w)
    ops.timeSeries('Path',1,'-dt',1.,'-values',0.,*[v[0] for v in p],'-useLast');ops.pattern('Plain',1,1);ops.sp(2,1,1.)
    ops.timeSeries('Path',2,'-dt',1.,'-values',0.,*[v[1] for v in p],'-useLast');ops.pattern('Plain',2,2);ops.load(2,0.,1.,0.)
    ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack');ops.test('NormUnbalance',1e-8,150,0,0);ops.algorithm('KrylovNewton','-maxDim',12) if w['type']=='sfi_mvlem' else ops.algorithm('Newton');ops.integrator('LoadControl',1.);ops.analysis('Static')
    rows=[]
    for step in range(len(p)):
        if ops.analyze(1)!=0:raise RuntimeError(f'OpenSees static failed at {step} for {w["type"]}')
        rows.append(dict(u=[*ops.nodeDisp(1),*ops.nodeDisp(2)],nodal_force=ops.eleResponse(w['id'],'globalForce')[:6],shear_deformation=(ops.eleResponse(w['id'],'ShearDef')[0] if w['type']=='sfi_mvlem' else -ops.eleResponse(w['id'],'Shear_Force_Deformation')[0])))
    return rows

def run_probe(exe,out,name,inp):
    ip=out/(name+'_input.json');op=out/(name+'_quake.json');ip.write_text(json.dumps(inp,indent=2));subprocess.run([str(exe),str(ip),str(op)],check=True);return json.loads(op.read_text())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--probe',type=pathlib.Path,required=True);ap.add_argument('--runner',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    summary=dict(opensees_version=ops.version(),scope='Numerical verification of named input laws, not FSAM or experimental validation',cases={});p=protocol()
    for kind in ['mvlem','sfi_mvlem']:
      for nl in [False,True]:
        name=kind+('_nonlinear' if nl else '_elastic');w=fixture(kind,nl)
        if kind=='sfi_mvlem' and nl:
            name='sfi_mvlem_two_layer_nonlinear'
            for i,f in enumerate(w['panels']):
                angle=30 if i%2==0 else -30
                f['material']=dict(type='layered_plane_stress',background_E=0.,background_nu=0.,layers=[dict(angle_deg=angle,weight=.1,material=STEEL),dict(angle_deg=angle+90,weight=.15,material=STEEL)])
        q=run_probe(a.probe,a.out,name,dict(height=3.,wall=w,protocol=p))['history'];r=static_ref(w,p);(a.out/(name+'_opensees.json')).write_text(json.dumps(r))
        uq=np.array([v['u'] for v in q]);ur=np.array([v['u'] for v in r]);fq=np.array([v['nodal_force'] for v in q]);fr=np.array([v['nodal_force'] for v in r]);sq=np.array([v['shear_deformation'] for v in q]);sr=np.array([v['shear_deformation'] for v in r])
        metrics=dict(steps=len(p),max_displacement_error=float(abs(uq-ur).max()),max_force_error=float(abs(fq-fr).max()),force_relative_error=float(abs(fq-fr).max()/max(1,abs(fr).max())),max_shear_deformation_error=float(abs(sq-sr).max()),peak_shear=float(abs(fq[:,3]).max()),peak_drift=float(abs(uq[:,3]).max()/3))
        metrics['passed']=metrics['max_displacement_error']<2e-7 and metrics['force_relative_error']<2e-6 and metrics['max_shear_deformation_error']<2e-7
        summary['cases'][name]=metrics;print(name,metrics,flush=True)
        (a.out/'summary.json').write_text(json.dumps(summary,indent=2))
    if not all(v['passed'] for v in summary['cases'].values()):raise SystemExit(1)
if __name__=='__main__':main()
