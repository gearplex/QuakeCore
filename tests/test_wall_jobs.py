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
print('wall JSON contracts passed')
