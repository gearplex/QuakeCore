#!/usr/bin/env python3
import argparse,json,pathlib
import numpy as np
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser();ap.add_argument('evidence',type=pathlib.Path);ap.add_argument('output',type=pathlib.Path);a=ap.parse_args()
plt.style.use('seaborn-v0_8-whitegrid');fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
colors={'q':'#087e8b','o':'#ff5a5f'}
for panel,name,title in [(ax[0,0],'mvlem_nonlinear','MVLEM: nonlinear fibers and shear'),(ax[0,1],'sfi_mvlem_two_layer_nonlinear','SFI-MVLEM: nonlinear layered panel')]:
 q=json.load(open(a.evidence/'comparison'/f'{name}_quake.json'))['history'];o=json.load(open(a.evidence/'comparison'/f'{name}_opensees.json'))
 panel.plot([100*x['u'][3]/3 for x in o],[x['nodal_force'][3] for x in o],lw=2.4,color=colors['o'],label='OpenSees 3.8.0')
 panel.plot([100*x['u'][3]/3 for x in q],[x['nodal_force'][3] for x in q],lw=1.1,ls='--',color=colors['q'],label='QuakeCore 9J')
 panel.set(title=title,xlabel='Top drift (%)',ylabel='Top shear (kN)');panel.legend(frameon=True)
for panel,name,title in [(ax[1,0],'mvlem_rc','3-story MVLEM wall-frame NRHA'),(ax[1,1],'sfi_two_layer','3-story SFI-MVLEM wall-frame NRHA')]:
 q=json.load(open(a.evidence/'dynamic'/f'{name}_quake.json'))['runs'][0]['history'];o=json.load(open(a.evidence/'dynamic'/f'{name}_opensees.json'))
 panel.plot([x['time_s'] for x in o],[100*x['floor_displacement'][-1]/9 for x in o],lw=2.4,color=colors['o'],label='OpenSees 3.8.0')
 panel.plot([x['time_s'] for x in q],[100*x['floor_displacement'][-1]/9 for x in q],lw=1.1,ls='--',color=colors['q'],label='QuakeCore 9J')
 panel.set(title=title,xlabel='Time (s)',ylabel='Roof drift ratio (%)');panel.legend(frameon=True)
fig.suptitle('QuakeCore Phase 9J wall numerical verification',fontsize=16,fontweight='bold')
fig.savefig(a.output,dpi=220);print(a.output)
