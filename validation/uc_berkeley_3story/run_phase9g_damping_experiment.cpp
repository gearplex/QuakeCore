#define main quakecore_phase8_embedded_main
#include "run_validation_phase8.cpp"
#undef main

#include "quake/modal_damping.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {
struct RunOut {
    AnalysisResult analysis;
    std::array<double,3> peak{};
    std::array<double,3> peak_time{};
    std::array<double,3> sim_at_s2{};
    double b1_pr{},b1_pa{},b1_P{},b1_M{};
    std::vector<double> time;
    std::vector<std::array<double,3>> drift;
    std::vector<double> pr,pa,P,M;
};
std::array<double,3> drift_pct9g(const CompiledFrame3D& m,const std::vector<double>& u){
    auto f=m.story_response_values(u);std::array<double,3>d{};double lo=0;
    for(int i=0;i<3;++i){d[i]=(f[i]-lo)/drift_story_h*100.0;lo=f[i];}return d;
}
RunOut run15(const NonlinearDynamicModel& analysis_model,const ModelBuild& physical,const Motion& motion){
    RunOut r;RobustNewmarkOptions ro;ro.tolerance=1e-7;ro.max_iterations=50;ro.line_search=true;ro.max_backtracks=10;ro.nonmonotone_line_search_factor=1.05;ro.max_subdivisions=6;ro.return_numerical_failure=true;
    auto obs=[&](double t,const std::vector<double>&u,const std::vector<double>&s,bool output){
        auto d=drift_pct9g(physical.model,u);for(int i=0;i<3;++i)if(std::abs(d[i])>r.peak[i]){r.peak[i]=std::abs(d[i]);r.peak_time[i]=t;if(i==1)r.sim_at_s2=d;}
        auto h=physical.model.pm_interaction_snapshot(3,u,s);if(std::abs(h.plastic_rotation)>r.b1_pr){r.b1_pr=std::abs(h.plastic_rotation);r.b1_P=h.compression;r.b1_M=h.moment;}r.b1_pa=std::max(r.b1_pa,std::abs(h.plastic_axial_deformation));
        if(output){r.time.push_back(t);r.drift.push_back(d);r.pr.push_back(h.plastic_rotation);r.pa.push_back(h.plastic_axial_deformation);r.P.push_back(h.compression);r.M.push_back(h.moment);}
    };
    ro.accepted_substep_state_observer=[&](std::size_t,std::size_t,int,double t,double,const std::vector<double>&u,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&s){obs(t,u,s,false);};
    ro.accepted_state_observer=[&](std::size_t,double t,double,const std::vector<double>&u,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&s){obs(t,u,s,true);};
    r.analysis=run_newmark_robust(analysis_model,motion.ag,motion.dt,LinearStrategy::FullFactorization,ro);return r;
}
void write_case(std::ostream& o,const char* name,const RunOut&r,bool comma){
    o<<"    \""<<name<<"\": {\n"
     <<"      \"termination\":"<<(int)r.analysis.termination<<",\n"
     <<"      \"termination_time_s\":"<<r.analysis.termination_time<<",\n"
     <<"      \"termination_reason\":\""<<r.analysis.termination_reason<<"\",\n"
     <<"      \"peak_story_drift_percent\":["<<r.peak[0]<<","<<r.peak[1]<<","<<r.peak[2]<<"],\n"
     <<"      \"time_at_peak_story_drift_s\":["<<r.peak_time[0]<<","<<r.peak_time[1]<<","<<r.peak_time[2]<<"],\n"
     <<"      \"simultaneous_at_s2_peak_percent\":["<<r.sim_at_s2[0]<<","<<r.sim_at_s2[1]<<","<<r.sim_at_s2[2]<<"],\n"
     <<"      \"b1_bottom_max_abs_plastic_rotation_rad\":"<<r.b1_pr<<",\n"
     <<"      \"b1_bottom_max_abs_plastic_axial_in\":"<<r.b1_pa<<",\n"
     <<"      \"b1_P_at_max_plastic_rotation_kip\":"<<r.b1_P<<",\n"
     <<"      \"b1_M_at_max_plastic_rotation_kip_in\":"<<r.b1_M<<"\n"
     <<"    }"<<(comma?",":"")<<"\n";
}
void write_hist(const std::string& path,const RunOut&r){
    std::ofstream c(path);c<<"time_s,s1_pct,s2_pct,s3_pct,b1_plastic_rotation_rad,b1_plastic_axial_in,b1_P_kip,b1_M_kip_in\n";
    for(std::size_t i=0;i<r.time.size();++i)c<<r.time[i]<<','<<r.drift[i][0]<<','<<r.drift[i][1]<<','<<r.drift[i][2]<<','<<r.pr[i]<<','<<r.pa[i]<<','<<r.P[i]<<','<<r.M[i]<<'\n';
}
}
int main(int argc,char**argv){
    try{
        const std::string motion_path=argc>1?argv[1]:"validation/uc_berkeley_3story/proxy_dt1_motion.csv";
        const std::string out_json=argc>2?argv[2]:"validation/uc_berkeley_3story/phase9/phase9g_damping_15s.json";
        const std::string out_prefix=argc>3?argv[3]:"validation/uc_berkeley_3story/phase9/phase9g";
        auto motion=read_motion(motion_path);const std::size_t n=std::min(motion.ag.size(),static_cast<std::size_t>(std::ceil(15.0/motion.dt))+1);motion.ag.resize(n);
        const auto gravity=gravity_column_preloads();constexpr double flex_scale=0.934616,target_T=.48,slope=-.005,residual=.20,acc=.55;const double w=2*M_PI/target_T;
        const double alpha_current=2*.025*w,beta_half=2*.005/w;
        auto make=[&](double alpha,double beta){
            auto raw=build_frame_phase8(flex_scale,gravity,true,true,false,true,false,3,alpha,beta,true,true,slope,residual,acc,-1.0,true,true);
            auto modes=modal_analysis(raw.model,3);return ModelBuild{std::move(raw.model),flex_scale,modes.empty()?0:modes[0].period,alpha,beta,raw.probes,true};
        };

        auto current=make(alpha_current,beta_half);
        auto r_current=run15(current.model,current,motion);
        std::cerr<<"9G current Rayleigh d="<<r_current.peak[0]<<","<<r_current.peak[1]<<","<<r_current.peak[2]<<"\n";

        auto hybrid=make(0.0,beta_half);
        FixedModalDampingModel hybrid_modal(hybrid.model,0.025,3);
        auto r_hybrid=run15(hybrid_modal,hybrid,motion);
        std::cerr<<"9G NIST hybrid d="<<r_hybrid.peak[0]<<","<<r_hybrid.peak[1]<<","<<r_hybrid.peak[2]<<"\n";

        auto modal3=make(0.0,0.0);
        FixedModalDampingModel modal3_model(modal3.model,0.030,3);
        auto r_modal3=run15(modal3_model,modal3,motion);
        std::cerr<<"9G modal3 d="<<r_modal3.peak[0]<<","<<r_modal3.peak[1]<<","<<r_modal3.peak[2]<<"\n";

        std::ofstream o(out_json);o<<std::setprecision(12);
        o<<"{\n  \"phase\":\"9G_DAMPING_SENSITIVITY\",\n  \"input_status\":\"PRE-VALIDATION / deterministic proxy\",\n  \"window_s\":15.0,\n"
         <<"  \"mechanics\":\"Phase 9F.2 Perform-type P-M associated flow + Mroz two-surface; frozen across damping cases\",\n"
         <<"  \"damping_cases\": {\n"
         <<"    \"current_rayleigh\":{\"modal_zeta\":0.0,\"alpha_m\":"<<alpha_current<<",\"beta_k\":"<<beta_half<<",\"description\":\"2.5% mass-proportional + 0.5% initial-stiffness-proportional at T1\"},\n"
         <<"    \"nist_hybrid\":{\"modal_zeta\":0.025,\"alpha_m\":0.0,\"beta_k\":"<<beta_half<<",\"description\":\"2.5% fixed-elastic modal + 0.5% initial-stiffness-proportional Rayleigh anchored at T1\"},\n"
         <<"    \"modal_3pct\":{\"modal_zeta\":0.03,\"alpha_m\":0.0,\"beta_k\":0.0,\"description\":\"3.0% fixed-elastic modal damping only\"}\n  },\n"
         <<"  \"results\": {\n";write_case(o,"current_rayleigh",r_current,true);write_case(o,"nist_hybrid",r_hybrid,true);write_case(o,"modal_3pct",r_modal3,false);o<<"  }\n}\n";
        write_hist(out_prefix+"_current_rayleigh_history.csv",r_current);
        write_hist(out_prefix+"_nist_hybrid_history.csv",r_hybrid);
        write_hist(out_prefix+"_modal3_history.csv",r_modal3);
        std::cout<<out_json<<"\n";return 0;
    }catch(const std::exception&e){std::cerr<<"9G error: "<<e.what()<<"\n";return 2;}
}
