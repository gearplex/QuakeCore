"""Render scientific review figures from the preserved Phase 9I evidence bundle.
Usage: python tools/render_review_figures.py /path/to/quakecore-phase9i-evidence output_dir
Requires NumPy and Matplotlib. Does not modify the underlying results.
"""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def render(evidence,output):
    evidence=Path(evidence);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','grid.alpha':.2})
    colors=['#067a9e','#d78130','#747896']
    def load(p):return json.loads(p.read_text())
    fig,axs=plt.subplots(2,2,figsize=(12.8,8.8),constrained_layout=True)
    summaries=[]
    for ax,n in zip(axs[0],(3,10)):
        p=evidence/f'final_opensees/frame{n}_integration';r=load(p/'quakecore.json');h=r['runs'][0]['history'];o=np.loadtxt(p/'opensees_displacement.txt')
        t=np.array([x['time_s'] for x in h]);u=np.array([x['floor_displacement'] for x in h]);ax.plot(t,1000*u[:,-1],color=colors[0],lw=2,label='QuakeCore');ax.plot(o[:,0],o[:,-1]*1000,'--',color=colors[1],label='OpenSees 3.8.0',lw=1)
        ax.set(title=f'{n}-story frame: roof displacement',xlabel='Time (s)',ylabel='Relative displacement (mm)');ax.legend(frameon=False);ax.grid()
        axs[1,1].plot(t,np.max(np.abs(u-o[:,1:]),axis=1)*1e12,label=f'{n}-story');summaries.append(load(p/'comparison.json'))
    ax=axs[1,0];x=np.arange(2)
    for j,(key,label) in enumerate([('quakecore_seconds','QuakeCore'),('opensees_seconds','OpenSees')]):
        med=np.array([np.median(s[key])*1000 for s in summaries]);lo=np.array([min(s[key])*1000 for s in summaries]);hi=np.array([max(s[key])*1000 for s in summaries]);ax.bar(x+(j-.5)*.3,med,.3,color=colors[j],label=label,yerr=np.vstack([med-lo,hi-med]),capsize=4)
    for i,s in enumerate(summaries):ax.text(i,max(s['opensees_seconds'])*1000+25,f"{s['speedup_opensees_over_quakecore']:.1f}×",ha='center',fontweight='bold')
    ax.set_xticks(x,['3-story / 39 DOF','10-story / 130 DOF']);ax.set(title='Timing without optional recorders',ylabel='Integration time (ms)',ylim=(0,1000));ax.legend(frameon=False);ax.grid(axis='y')
    axs[1,1].set(title='Maximum floor displacement difference',xlabel='Time (s)',ylabel='Absolute difference (10⁻¹² m)');axs[1,1].legend(frameon=False);axs[1,1].grid()
    fig.suptitle('QuakeCore Phase 9I • Independent numerical verification',fontsize=18,fontweight='bold');fig.supxlabel('Synthetic bilinear frames • 2,400 steps at Δt = 0.005 s • no P–Delta or degrading P–M law\nSeven timing repetitions; medians and observed ranges. Setup/eigenanalysis excluded.',fontsize=10)
    fig.savefig(output/'QUAKECORE_PHASE9I_OPENSEES.png',dpi=180);plt.close(fig)
    p=evidence/'phase9i_release';motion=evidence/'inputs/proxy_dt1_motion.csv';a=np.loadtxt(motion,delimiter=',',skiprows=1)
    frac=np.trapezoid(a[:1501,1]**2,a[:1501,0])/np.trapezoid(a[:,1]**2,a[:,0])
    fig,axs=plt.subplots(2,2,figsize=(12.8,9.1),constrained_layout=True);ax=axs[0,0]
    ax.plot(a[:,0],a[:,1],lw=.7,color='#64778d');ax.axvspan(0,15,color=colors[0],alpha=.15);ax.axvline(15,color=colors[0],ls='--');ax.set(title='The earlier window covers only part of the proxy',xlabel='Time (s)',ylabel='Synthetic acceleration (g)');ax.text(.97,.96,f'0–15 s: {frac:.2%} of ∫a² dt\nFull input: 70 s',ha='right',va='top',transform=ax.transAxes,bbox={'facecolor':'white','edgecolor':'none','alpha':.85});ax.grid()
    ax=axs[0,1];r=load(p/'hybrid_70s.json');h=r['history'];t=[x['time_s'] for x in h];d=np.array([x['story_drift_percent'] for x in h])
    for i in range(3):ax.plot(t,d[:,i],color=colors[i],lw=1,label=f'Story {i+1}')
    ax.axvline(r['termination_time_s'],color='#b63340',ls='--');ax.set(title='Extended hybrid run: numerical failure at 19.71 s',xlabel='Time (s)',ylabel='Drift over 39-in clear height (%)',xlim=(0,22));ax.grid();ax.legend(frameon=False,ncol=3,fontsize=9)
    ax=axs[1,0]
    for i,(name,label) in enumerate([('hybrid','Δt = 0.010 s'),('hybrid_dt2','Δt = 0.005 s')]):ax.plot(load(p/(name+'.json'))['peak_story_drift_percent'],[1,2,3],'-o' if i==0 else '--o',color=colors[i],label=label)
    ax.set(title='Completed 15-s hybrid proxy runs',xlabel='Peak story drift (%)',ylabel='Story',yticks=[1,2,3]);ax.grid();ax.legend(frameon=False);ax.text(.02,.03,'Δt = 0.0025 s fails at 13.0075 s.\nTime-step convergence is not established.',transform=ax.transAxes,fontsize=9,bbox={'facecolor':'white','edgecolor':'#ddc8c8'})
    ax=axs[1,1];names=['mroz_rayleigh','mroz_dt2','mroz_dt4','hybrid','hybrid_dt2','hybrid_dt4','hybrid_70s'];labels=['Mroz • 10 ms','Mroz • 5 ms','Mroz • 2.5 ms','Hybrid • 10 ms','Hybrid • 5 ms','Hybrid • 2.5 ms','Hybrid • 5 ms / 70 s']
    for i,name in enumerate(names):
        r=load(p/(name+'.json'));done=r['termination']=='completed';ax.barh(i,r['termination_time_s'],height=.55,color='#138578' if done else '#b63340');ax.text(r['termination_time_s']+.3,i,'complete' if done else 'numerical failure',va='center',fontsize=8)
    ax.set_yticks(range(len(names)),labels);ax.invert_yaxis();ax.set(xlim=(0,30),xlabel='Completion or failed-step time (s)',title='Failure status stays visible');ax.grid(axis='x')
    fig.suptitle('Berkeley / NIST continuation • Proxy validation remains open',fontsize=18,fontweight='bold');fig.supxlabel('Fixed physical parameters; recorded DT1 and measured/reference histories are missing.\nPartial-run peaks cover accepted history only; they are not validated building demands.',fontsize=10)
    fig.savefig(output/'QUAKECORE_PHASE9I_BERKELEY.png',dpi=180);plt.close(fig)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence');p.add_argument('output');a=p.parse_args();render(a.evidence,a.output)
