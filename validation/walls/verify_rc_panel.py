#!/usr/bin/env python3
"""RC panel force replay against OpenSees uniaxial laws at accepted strains.
This is a constitutive/assembly comparison, NOT an independent global RC wall solve.
"""
import argparse, pathlib,json,sys,math
import numpy as np
from compare_opensees import ops,Tags,uniaxial,fixture,protocol,run_probe

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--probe',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 w=fixture('sfi_mvlem',True);p=protocol();out=run_probe(a.probe,a.out,'sfi_fixed_angle_rc',dict(height=3.,wall=w,protocol=p));rows=out['history']
 ops.wipe();tags=Tags();bank=[];width=sum(f['width'] for f in w['panels']);left=-width/2
 for f in w['panels']:
  m=f['material'];angle=math.radians(m['angle_deg']);layers=[]
  for aa,weight,mat in [(angle,1,m['concrete']),(angle+math.pi/2,1,m['concrete']),(0,m['rho_x'],m['steel_x']),(math.pi/2,m['rho_y'],m['steel_y'])]:
   c,s=math.cos(aa),math.sin(aa);layers.append((np.array([c*c,s*s,c*s]),weight,uniaxial(mat,tags)))
  E,nu=m['background_E'],m['background_nu'];D=E/(1-nu*nu)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
  bank.append((left+f['width']/2,f['width']*f['thickness'],D,layers));left+=f['width']
 force_err=stress_err=transverse=0.;strain_max=0.;ref=[]
 for row in rows:
  u=np.array(row['u']);F=np.zeros(6);ps=[]
  gamma=(u[3]-u[0])/3+.4*u[2]+.6*u[5]
  for i,(x,A,D,layers) in enumerate(bank):
   e=np.array([row['panel_strain'][i][0],(u[4]-u[1]+x*(u[5]-u[2]))/3,gamma]);stress=D@e
   for n,weight,tag in layers:
    strain=float(n@e);strain_max=max(strain_max,abs(strain));ops.testUniaxialMaterial(tag);ops.setStrain(strain);stress+=weight*ops.getStress()*n
   by=np.array([0,-1/3,-x/3,0,1/3,x/3]);bg=np.array([-1/3,0,.4,1/3,0,.6]);F+=3*A*(stress[1]*by+stress[2]*bg);ps.append(stress.tolist());transverse=max(transverse,abs(stress[0]));stress_err=max(stress_err,float(abs(stress-np.array(row['panel_stress'][i])).max()))
  force_err=max(force_err,float(abs(F-np.array(row['nodal_force'])).max()));ref.append(dict(nodal_force=F.tolist(),panel_stress=ps))
 metrics=dict(steps=len(rows),comparison='OpenSees uniaxial stress replay at QuakeCore accepted panel transverse strains; no independent global RC solve',max_force_error=force_err,max_panel_stress_error=stress_err,max_transverse_stress=transverse,max_layer_abs_strain=strain_max,passed=bool(force_err<1e-6 and stress_err<1e-5 and transverse<1e-4))
 (a.out/'sfi_fixed_angle_rc_replay.json').write_text(json.dumps(ref));(a.out/'rc_replay_summary.json').write_text(json.dumps(metrics,indent=2));print(json.dumps(metrics,indent=2));sys.exit(0 if metrics['passed'] else 1)
if __name__=='__main__':main()
