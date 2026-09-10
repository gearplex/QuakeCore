"""End-to-end contract checks for engineering job/result boundaries."""
import copy, json, subprocess, sys, tempfile
from pathlib import Path
exe=Path(sys.argv[1]).resolve();source=Path(sys.argv[2]).resolve()
job=json.loads((source/'examples/frame3_nrha.json').read_text())
with tempfile.TemporaryDirectory() as d:
    root=Path(d)
    def run(j, name):
        p=root/(name+'.json');p.write_text(json.dumps(j))
        out=root/(name+'_result.json');r=subprocess.run([str(exe),str(p),str(out)],capture_output=True,text=True)
        return r, json.loads(out.read_text()) if out.exists() else None
    r,out=run(job,'nrha');assert r.returncode==0,r.stderr
    a=out['runs'][0];assert a['termination']=='completed' and len(a['history'])==2400
    assert abs(a['history'][-1]['time_s']-12)<1e-12
    assert max(a['peak_story_drift_ratio'])>.004
    assert out['code_compliance']=='not_assessed' and out['gravity_analysis_performed'] is False
    summary=copy.deepcopy(job);summary['analysis']['output_mode']='minimal';r,small=run(summary,'minimal')
    assert r.returncode==0 and 'history' not in small['runs'][0]
    assert small['runs'][0]['final_displacement']==a['final_displacement']
    for field,value in [('unexpected_field',1),('records',[])]:
        bad=copy.deepcopy(job);bad[field]=value;r,out=run(bad,'invalid_'+field);assert r.returncode==2
    bad=copy.deepcopy(job);bad['model']['nodes'][0][0]=1.5;r,out=run(bad,'fractional_id');assert r.returncode==2
    bad=copy.deepcopy(job);bad['records'][0]['acceleration'][3]=float('nan');r,out=run(bad,'nan');assert r.returncode==2
    ida=json.loads((source/'examples/frame3_ida.json').read_text());r,out=run(ida,'ida');assert r.returncode==0,r.stderr
    for record in out['records']:
        assert 0<record['lower_scale']<record['upper_scale']
        assert not record['numerical_failure'] and not record['bracket_has_numerical_gap']
    assert any(a['termination']=='configured_collapse_criterion' for a in out['runs'])
print('Job contracts passed: NRHA, minimal-output equivalence, malformed input, and parallel IDA.')
