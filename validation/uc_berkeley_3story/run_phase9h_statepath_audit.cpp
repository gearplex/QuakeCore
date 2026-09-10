#define main quakecore_phase8_embedded_main
#include "run_validation_phase8.cpp"
#undef main

#include "quake/modal_damping.hpp"
#include "quake/superlu_solver.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
struct SavedState { bool valid{}; double time{}; std::array<double,3> drift{}; std::vector<double> u,state; };

std::array<double,3> drift_pct(const CompiledFrame3D& m,const std::vector<double>& u){
    auto f=m.story_response_values(u);std::array<double,3>d{};double lo=0;
    for(int i=0;i<3;++i){d[i]=(f[i]-lo)/drift_story_h*100.0;lo=f[i];}return d;
}

struct StoryTangent { bool ok{}; std::array<double,9> k{}; std::array<double,3> diag{}; };
StoryTangent story_tangent(const CompiledFrame3D& frame,const SparseMatrixCSC& kt){
    StoryTangent out;constexpr int ns=3;const int n=frame.dof();
    std::array<std::vector<double>,ns> c;for(auto& v:c)v.assign(static_cast<std::size_t>(n),0.0);
    std::vector<double> e(static_cast<std::size_t>(n),0.0);
    for(int k=0;k<n;++k){std::fill(e.begin(),e.end(),0.0);e[static_cast<std::size_t>(k)]=1.0;auto y=frame.story_response_values(e);for(int i=0;i<ns;++i)c[i][static_cast<std::size_t>(k)]=y[i];}
    try{
        std::array<double,9> F{};
        for(int j=0;j<ns;++j){auto x=superlu_solve_once(kt,c[j]);for(int i=0;i<ns;++i){double v=0;for(int k=0;k<n;++k)v+=c[i][static_cast<std::size_t>(k)]*x[static_cast<std::size_t>(k)];F[static_cast<std::size_t>(i*ns+j)]=v;}}
        std::vector<Triplet> ft;for(int i=0;i<ns;++i)for(int j=0;j<ns;++j)ft.push_back({i,j,F[static_cast<std::size_t>(i*ns+j)]});
        auto Fs=SparseMatrixCSC::from_triplets(ns,ns,ft,0.0);std::array<double,9> kfloor{};
        for(int j=0;j<ns;++j){std::vector<double> rhs(ns,0.0);rhs[j]=1.0;auto x=superlu_solve_once(Fs,rhs);for(int i=0;i<ns;++i)kfloor[static_cast<std::size_t>(i*ns+j)]=x[i];}
        for(int a=0;a<ns;++a)for(int b=0;b<ns;++b){double v=0;for(int i=a;i<ns;++i)for(int j=b;j<ns;++j)v+=kfloor[static_cast<std::size_t>(i*ns+j)];out.k[static_cast<std::size_t>(a*ns+b)]=v;}
        for(int i=0;i<ns;++i)out.diag[i]=out.k[static_cast<std::size_t>(i*ns+i)];out.ok=true;
    }catch(...){out.ok=false;}
    return out;
}

void write_state(std::ostream& o,const char* name,const SavedState& ss,const CompiledFrame3D& model){
    o<<"    \""<<name<<"\": {\n";
    o<<"      \"time_s\":"<<ss.time<<",\n      \"story_drift_percent\":["<<ss.drift[0]<<','<<ss.drift[1]<<','<<ss.drift[2]<<"],\n";
    if(!ss.valid){o<<"      \"available\":false\n    }";return;}
    std::vector<double> f,t,tr;model.internal_force_and_tangent(ss.u,ss.state,f,t,tr);auto K=model.effective_state_tangent_matrix_with_state(ss.u,t,ss.state,0.0,0.0);
    auto st=story_tangent(model,K);o<<"      \"available\":true,\n      \"interstory_tangent_matrix_kip_per_in\":[";
    for(int i=0;i<9;++i){if(i)o<<',';o<<st.k[static_cast<std::size_t>(i)];}o<<"],\n      \"interstory_tangent_diagonal_kip_per_in\":["<<st.diag[0]<<','<<st.diag[1]<<','<<st.diag[2]<<"],\n";
    auto modes=modal_analysis_from_stiffness(model,K,3);o<<"      \"positive_tangent_modes\":[";
    for(std::size_t i=0;i<modes.size();++i){if(i)o<<',';auto sh=model.story_response_values(modes[i].shape);double den=(sh.size()>=3&&std::abs(sh[2])>1e-14)?sh[2]:1.0;o<<"{\"period_s\":"<<modes[i].period<<",\"roof_normalized_floor_shape\":["<<(sh.size()>0?sh[0]/den:0)<<','<<(sh.size()>1?sh[1]/den:0)<<','<<(sh.size()>2?sh[2]/den:0)<<"]}";}o<<"],\n";
    o<<"      \"column_hinges\":[\n";
    for(int id=1;id<=24;++id){auto h=model.pm_interaction_snapshot(id,ss.u,ss.state);int idx=(id-1)/2,story=idx/4+1,grid=idx%4;const char gridc=static_cast<char>('A'+grid);const char* end=(id%2)?"bottom":"top";
        o<<"        {\"spring_id\":"<<id<<",\"column\":\""<<gridc<<story<<"\",\"end\":\""<<end<<"\",\"P_kip\":"<<h.compression<<",\"M_kip_in\":"<<h.moment<<",\"plastic_rotation_rad\":"<<h.plastic_rotation<<",\"plastic_axial_in\":"<<h.plastic_axial_deformation<<",\"kappa_pos_rad\":"<<h.positive_plastic_rotation<<",\"kappa_neg_rad\":"<<h.negative_plastic_rotation<<",\"active_capacity_kip_in\":"<<h.active_capacity<<",\"events\":"<<h.events<<",\"tangent_dN_dd\":"<<h.tangent[0]<<",\"tangent_dN_dtheta\":"<<h.tangent[1]<<",\"tangent_dM_dd\":"<<h.tangent[2]<<",\"tangent_dM_dtheta\":"<<h.tangent[3]<<"}"<<(id<24?",":"")<<"\n";
    }
    o<<"      ]\n    }";
}
}

int main(int argc,char**argv){
    try{
        const std::string motion_path=argc>1?argv[1]:"validation/uc_berkeley_3story/proxy_dt1_motion.csv";
        const std::string out_json=argc>2?argv[2]:"validation/uc_berkeley_3story/phase9/phase9h_statepath_audit.json";
        auto motion=read_motion(motion_path);const std::size_t n=std::min(motion.ag.size(),static_cast<std::size_t>(std::ceil(12.5/motion.dt))+1);motion.ag.resize(n);
        const auto gravity=gravity_column_preloads();constexpr double flex_scale=.934616,target_T=.48,slope=-.005,residual=.20,acc=.55;const double w=2*M_PI/target_T,alpha=2*.025*w,beta=2*.005/w;
        auto raw=build_frame_phase8(flex_scale,gravity,true,true,false,true,false,3,alpha,beta,true,true,slope,residual,acc,-1.0,true,false);
        auto modes=modal_analysis(raw.model,3);ModelBuild mb{std::move(raw.model),flex_scale,modes.empty()?0:modes[0].period,alpha,beta,raw.probes,true};
        SavedState s518,s1peak,s2peak,s3peak;double p1=0,p2=0,p3=0;
        RobustNewmarkOptions ro;ro.tolerance=1e-7;ro.max_iterations=50;ro.line_search=true;ro.max_backtracks=10;ro.nonmonotone_line_search_factor=1.05;ro.max_subdivisions=6;ro.return_numerical_failure=true;
        auto capture=[&](SavedState& s,double time,const std::vector<double>&u,const std::vector<double>&state,const std::array<double,3>&d){s.valid=true;s.time=time;s.drift=d;s.u=u;s.state=state;};
        auto observe=[&](double time,const std::vector<double>&u,const std::vector<double>&state,bool output){auto d=drift_pct(mb.model,u);if(std::abs(d[0])>p1){p1=std::abs(d[0]);capture(s1peak,time,u,state,d);}if(std::abs(d[1])>p2){p2=std::abs(d[1]);capture(s2peak,time,u,state,d);}if(std::abs(d[2])>p3){p3=std::abs(d[2]);capture(s3peak,time,u,state,d);}if(output && std::abs(time-5.18)<=0.51*motion.dt)capture(s518,time,u,state,d);};
        // This audit needs representative committed peak states, not substep
        // envelope precision.  Record only accepted output states so the
        // diagnostic itself does not dominate the expensive PM replay.
        ro.accepted_state_observer=[&](std::size_t,double time,double,const std::vector<double>&u,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&s){observe(time,u,s,true);};
        auto ar=run_newmark_robust(mb.model,motion.ag,motion.dt,LinearStrategy::FullFactorization,ro);
        std::ofstream o(out_json);o<<std::setprecision(12);o<<"{\n  \"phase\":\"9H_REMAINING_STORY2_STATEPATH_AUDIT\",\n  \"input_status\":\"PRE-VALIDATION / deterministic proxy\",\n  \"mechanics\":\"Phase 9F.1 P-M associated flow EPP; used for detailed state-path audit because Phase 9F.2 Mroz shifted peak drifts by <0.3%\",\n  \"damping\":\"current QuakeCore 2.5% mass-proportional + 0.5% initial-stiffness Rayleigh; used for detailed mechanics audit because Phase 9G showed near-identical global drift to NIST hybrid\",\n  \"termination\":"<<(int)ar.termination<<",\n  \"termination_time_s\":"<<ar.termination_time<<",\n  \"states\":{\n";
        write_state(o,"legacy_reference_time_5_18s",s518,mb.model);o<<",\n";write_state(o,"story1_peak",s1peak,mb.model);o<<",\n";write_state(o,"story2_peak",s2peak,mb.model);o<<",\n";write_state(o,"story3_peak",s3peak,mb.model);o<<"\n  }\n}\n";
        std::cout<<out_json<<"\n";return 0;
    }catch(const std::exception&e){std::cerr<<"9H audit error: "<<e.what()<<"\n";return 2;}
}
