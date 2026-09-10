"""Verify the waveform identity and expose validation readiness without a fitted EDP gate."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def check(path):
    path=Path(path);m=json.loads(path.read_text())
    if m['schema']!='quakecore.validation_input.v1':raise ValueError('unsupported manifest schema')
    p=path.parent/m['motion_file'];actual=hashlib.sha256(p.read_bytes()).hexdigest()
    if actual!=m['sha256']:raise ValueError('waveform checksum differs from manifest')
    if m['sample_units']!='g' or m['applied_scale']!=1:raise ValueError('normalize units/scaling explicitly before using this archived Berkeley adapter')
    if not m.get('source') or not m.get('initial_state_convention'):raise ValueError('missing waveform provenance/initial-state convention')
    with p.open() as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=['time_s','accel_g']:raise ValueError('invalid CSV header')
        rows=[(float(r['time_s']),float(r['accel_g'])) for r in reader]
    if len(rows)<3 or any(not math.isfinite(v) for row in rows for v in row):raise ValueError('invalid waveform samples')
    dt=rows[1][0]-rows[0][0]
    if dt<=0 or rows[0][0]!=0 or any(abs(t-i*dt)>1e-7*dt for i,(t,g) in enumerate(rows)):raise ValueError('nonuniform or shifted waveform time grid')
    if m['motion_kind'] not in ['synthetic_proxy','recorded_table']:raise ValueError('unsupported motion_kind')
    result={'sha256':actual,'samples_including_t0':len(rows),'driver_steps':len(rows)-1,'dt_s':dt,'duration_s':rows[-1][0],
        'pga_g':max(abs(g) for t,g in rows),'initial_acceleration_g':rows[0][1],
        'input_kind':m['motion_kind'],'status':'proxy_only' if m['motion_kind']=='synthetic_proxy' else 'recorded_input_ready_for_reference_comparison',
        'experimental_validation':'not_established','missing_reference_artifacts':m.get('missing_reference_artifacts',[])}
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('manifest');p.add_argument('--output');a=p.parse_args()
    result=check(a.manifest);text=json.dumps(result,indent=2)+'\n'
    if a.output:Path(a.output).write_text(text)
    print(text,end='')
