#!/usr/bin/env python3
import pathlib,sys,json,tempfile,subprocess,copy
exe=pathlib.Path(sys.argv[1]);source=pathlib.Path(sys.argv[2])
base=json.loads((source/'examples/walls/mvlem_rc_nrha.json').read_text())
base['records'][0]['acceleration']=base['records'][0]['acceleration'][:40]
with tempfile.TemporaryDirectory() as tmp:
    root=pathlib.Path(tmp)
    def run(j,code):
        (root/'j.json').write_text(json.dumps(j));r=subprocess.run([str(exe),str(root/'j.json'),str(root/'r.json')],capture_output=True,text=True)
        assert r.returncode==code,(r.returncode,r.stderr)
        return json.loads((root/'r.json').read_text()) if code==0 else None
    r=run(base,0);assert len(r['runs'][0]['history'][0]['walls'])==3
    assert r['wall_scope']['wall_geometric_stiffness'] is False
    j=copy.deepcopy(base);j['analysis']['output_mode']='summary';r=run(j,0);assert not r['runs'][0]['history'];assert 'peak_abs_axial_force' in r['runs'][0]['walls']['1']
    for mutate in [lambda w:w.update(type='SFI_MVLEM_3D'),lambda w:w.update(c=2),lambda w:w.update(i=2,j=1),lambda w:w.update(id=101),lambda w:w['fibers'][0].update(width=0),lambda w:w['fibers'][0].update(rho=1),lambda w:w['shear'].update(unrecognized=2)]:
        j=copy.deepcopy(base);mutate(j['model']['walls'][0]);run(j,2)
    j=copy.deepcopy(base);j['model']['walls']={};run(j,2)
    j=json.loads((source/'examples/walls/sfi_fixed_angle_rc_nrha.json').read_text());j['records'][0]['acceleration']=j['records'][0]['acceleration'][:40];r=run(j,0);assert len(r['runs'][0]['history'][0]['walls']['1']['panel_stress'])==4
    for value in ['FSAM','unsupported']:
        k=copy.deepcopy(j);k['model']['walls'][0]['panels'][0]['material']['type']=value;run(k,2)
    k=copy.deepcopy(j);k['model']['walls'][0]['local_max_iterations']=2.5;run(k,2)

    # Gate 4 runner contract: instantiate ConcreteCM + Pinching4 wrapped in
    # MinMax/Parallel, with Pinching4 also used as MVLEM shear. Zero excitation
    # keeps this test focused on parsing, construction, state allocation, and
    # the normal quake_run execution path rather than constitutive calibration.
    k=copy.deepcopy(base);k['records'][0]['acceleration']=[0.0]*5
    cm={'type':'concrete_cm','fc':-6.5,'epsc':-.002,'Ec':4595.486916530173,'rc':7.0,'xcrn':1.03,'ft':.0604669,'et':2*.0604669/4595.486916530173,'rt':1.2,'xcrp':10000.0,'gap_close':True}
    pinch={'type':'pinching4','positive':[[1.0,10.0],[2.0,20.0],[3.0,30.0],[4.0,40.0]],'negative':[[-1.0,-10.0],[-2.0,-20.0],[-3.0,-30.0],[-4.0,-40.0]],'r_disp_positive':.6,'r_force_positive':.99,'u_force_positive':.4,'r_disp_negative':.6,'r_force_negative':.99,'u_force_negative':.4,'gamma_k':[0,0,0,0],'gamma_k_limit':2.0,'gamma_d':[.1,0,0,0],'gamma_d_limit':2.0,'gamma_f':[0,0,0,0],'gamma_f_limit':2.0,'gamma_e':10000.0,'damage_mode':'energy','admitted_reversal_count':2}
    wrapped={'type':'parallel','materials':[{'type':'minmax','min':-10.0,'max':10.0,'material':pinch},{'type':'elastic','E':.01}]}
    k['model']['walls'][0]['fibers'][0]['concrete']=cm;k['model']['walls'][0]['fibers'][0]['steel']=wrapped;k['model']['walls'][0]['shear']=pinch
    r=run(k,0);assert r['runs'][0]['history'][0]['walls']['1']['fiber_strain']
print('wall JSON contracts passed')
