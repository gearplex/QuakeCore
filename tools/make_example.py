"""Generate a transparent plane-frame job; numerical verification, not a design archetype."""
import argparse, json, math
from pathlib import Path

def frame_job(stories=3, bays=2, steps=2400, dt=.005):
    cols=bays+1; nodes=[]; fix=[]; equal=[]; members=[]; hinges=[]; cuts=[]
    for s in range(stories+1):
        for c in range(cols):
            tag=1000*s+c+1
            nodes.append([tag,6*c,3.5*s,100/cols if s else 0,0,0])
            if not s: fix.append([tag,True,True,True])
    for s in range(1,stories+1):
        cut=[]
        for c in range(cols):
            tag=100000+100*s+c
            members.append([tag,1000*(s-1)+c+1,1000*s+c+1,25e6,.25,.005208333333333333,0])
            cut.append(tag)
        cuts.append(cut)
        for bay in range(bays):
            end_i=200000+1000*s+2*bay; end_j=end_i+1
            joint_i=1000*s+bay+1; joint_j=joint_i+1
            for end,joint,x in [(end_i,joint_i,6*bay),(end_j,joint_j,6*(bay+1))]:
                nodes.append([end,x,3.5*s,0,0,0])
                equal.extend([[joint,end,'UX'],[joint,end,'UY']])
                hinges.append({'id':end+200000,'i':joint,'j':end,'type':'bilinear','k':12e6,'fy':150,'hardening_ratio':.002})
            members.append([300000+1000*s+bay,end_i,end_j,25e6,.18,.0054,0])
    ag=[]
    for i in range(1,steps+1):
        t=i*dt; env=math.sin(math.pi*t/(steps*dt))**2
        ag.append(2*env*(.62*math.sin(2*math.pi*.8*t)+.28*math.sin(2*math.pi*2.4*t+.7)+.10*math.sin(2*math.pi*4.7*t+1.2)))
    return {'schema':'quakecore.job.v1','name':f'{stories}-story {bays}-bay beam-hinge verification frame',
        'units':{'force':'kN','length':'m','time':'s'},
        'provenance':{'purpose':'synthetic numerical verification fixture, not a code-designed building','gravity':'no gravity analysis or P-Delta in this cross-solver fixture'},
        'model':{'type':'frame2d','nodes':nodes,'fixities':fix,'equal_dofs':equal,'members':members,'hinges':hinges,'rayleigh':[.05,.0005],
        'response_node':1000*stories+1,'story_nodes':[1000*s+1 for s in range(1,stories+1)],'story_cut_members':cuts},
        'analysis':{'type':'nrha','strategy':'woodbury','tolerance':1e-8,'relative_tolerance':False,'initial_guess':'previous_displacement','max_iterations':35,'max_subdivisions':0,'line_search':False,'history_stride':1},
        'records':[{'name':'deterministic_multifrequency','dt':dt,'acceleration':ag,'sample_convention':'step_end','provenance':{'kind':'synthetic','acceleration_units':'m/s2'}}]}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--stories',type=int,default=3);p.add_argument('--steps',type=int,default=2400);a=p.parse_args()
    a.output.write_text(json.dumps(frame_job(a.stories,steps=a.steps),indent=2)+'\n')
