from pathlib import Path

path = Path("apps/quake_run.cpp")
text = path.read_text()
old = '''   o.accepted_substep_state_observer=[&](std::size_t,std::size_t,int,double time,double ag,const auto& u,const auto& vel,const auto& a,const auto& s){auto d=drift(u),acc=model.story_response_values(a),shear=shears(u,s);'''
new = '''   o.accepted_substep_state_observer=[&](std::size_t parent_step,std::size_t substep_index,int depth,double time,double ag,const auto& u,const auto& vel,const auto& a,const auto& s){if(depth>0){auto wr=model.wall_response(3,u,s);std::cerr<<std::setprecision(17)<<"Gate4 accepted-substep parent_step="<<parent_step<<" substep_index="<<substep_index<<" depth="<<depth<<" time="<<time<<" wall3_f7_strain="<<wr.fiber_strain.at(7)<<" wall3_f7_concrete_stress="<<wr.concrete_stress.at(7)<<"\\n";}auto d=drift(u),acc=model.story_response_values(a),shear=shears(u,s);'''
if text.count(old) != 1:
    raise SystemExit("expected accepted_substep_state_observer anchor not found exactly once")
path.write_text(text.replace(old, new, 1))
