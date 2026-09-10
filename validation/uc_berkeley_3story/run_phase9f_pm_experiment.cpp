#define main quakecore_phase8_embedded_main
#include "run_validation_phase8.cpp"
#undef main

#include "quake/superlu_solver.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {

struct Phase9FRun {
    std::string name;
    AnalysisResult analysis;
    std::array<double,3> peak_drift_pct{};
    std::array<double,3> peak_time_s{};
    std::array<std::array<double,3>,3> simultaneous_pct{};
    std::array<double,12> max_compression_kip{};
    double peak_component_shear_kip{};
    double b1_bottom_max_abs_rotation{};
    double b1_bottom_max_abs_plastic_rotation{};
    double b1_bottom_max_abs_plastic_axial{};
    double b1_bottom_time_at_max_plastic_rotation{};
    double b1_bottom_compression_at_max_plastic_rotation{};
    double b1_bottom_moment_at_max_plastic_rotation{};
    std::vector<double> s2_peak_u,s2_peak_state;
    std::array<double,3> s2_story_tangent_diag{{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()}};
    std::vector<double> s2_mode_periods;
    std::vector<std::array<double,3>> s2_mode_story_shapes;
    std::vector<double> time;
    std::vector<std::array<double,3>> drift_history_pct;
    std::vector<double> b1_rotation_history,b1_plastic_rotation_history,b1_plastic_axial_history,b1_p_history,b1_m_history;
};

std::array<double,3> drift_pct(const CompiledFrame3D& m,const std::vector<double>& u){
    const auto f=m.story_response_values(u);std::array<double,3> d{};double lo=0.0;
    for(int i=0;i<3;++i){d[static_cast<std::size_t>(i)]=(f[static_cast<std::size_t>(i)]-lo)/drift_story_h*100.0;lo=f[static_cast<std::size_t>(i)];}
    return d;
}

std::array<double,3> story_tangent_diag(const CompiledFrame3D& frame,const SparseMatrixCSC& kt){
    constexpr int ns=3;const int n=frame.dof();
    std::array<std::vector<double>,ns> c;for(auto& v:c)v.assign(static_cast<std::size_t>(n),0.0);
    std::vector<double> e(static_cast<std::size_t>(n),0.0);
    for(int k=0;k<n;++k){std::fill(e.begin(),e.end(),0.0);e[static_cast<std::size_t>(k)]=1.0;const auto y=frame.story_response_values(e);for(int i=0;i<ns;++i)c[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]=y[static_cast<std::size_t>(i)];}
    std::vector<double> F(static_cast<std::size_t>(ns*ns),0.0);
    for(int j=0;j<ns;++j){const auto x=superlu_solve_once(kt,c[static_cast<std::size_t>(j)]);for(int i=0;i<ns;++i){double v=0;for(int k=0;k<n;++k)v+=c[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*x[static_cast<std::size_t>(k)];F[static_cast<std::size_t>(i*ns+j)]=v;}}
    std::vector<Triplet> ft;for(int i=0;i<ns;++i)for(int j=0;j<ns;++j)ft.push_back({i,j,F[static_cast<std::size_t>(i*ns+j)]});
    auto Fs=SparseMatrixCSC::from_triplets(ns,ns,ft,0.0);std::vector<double> kfloor(static_cast<std::size_t>(ns*ns),0.0);
    for(int j=0;j<ns;++j){std::vector<double> rhs(ns,0.0);rhs[static_cast<std::size_t>(j)]=1.0;auto x=superlu_solve_once(Fs,rhs);for(int i=0;i<ns;++i)kfloor[static_cast<std::size_t>(i*ns+j)]=x[static_cast<std::size_t>(i)];}
    std::array<double,3> out{};for(int a=0;a<ns;++a){double v=0;for(int i=a;i<ns;++i)for(int j=a;j<ns;++j)v+=kfloor[static_cast<std::size_t>(i*ns+j)];out[static_cast<std::size_t>(a)]=v;}return out;
}

Phase9FRun run_case(const std::string& name,const ModelBuild& mb,const Motion& motion,bool pm){
    Phase9FRun r;r.name=name;const auto& model=mb.model;
    RobustNewmarkOptions ro;ro.tolerance=1e-7;ro.max_iterations=50;ro.line_search=true;ro.max_backtracks=10;ro.nonmonotone_line_search_factor=1.05;ro.max_subdivisions=6;ro.return_numerical_failure=true;ro.collapse.check_initial_stability=true;
    auto observe=[&](double time,const std::vector<double>& u,const std::vector<double>& state,bool output){
        const auto d=drift_pct(model,u);
        for(int i=0;i<3;++i)if(std::abs(d[static_cast<std::size_t>(i)])>r.peak_drift_pct[static_cast<std::size_t>(i)]){
            r.peak_drift_pct[static_cast<std::size_t>(i)]=std::abs(d[static_cast<std::size_t>(i)]);r.peak_time_s[static_cast<std::size_t>(i)]=time;r.simultaneous_pct[static_cast<std::size_t>(i)]=d;
            if(i==1){r.s2_peak_u=u;r.s2_peak_state=state;}
        }
        for(std::size_t k=0;k<mb.column_probes.size();++k){const auto e=model.elastic_element_response(mb.column_probes[k].elastic_element_id,u);r.max_compression_kip[k]=std::max(r.max_compression_kip[k],e.axial_compression);}
        const int b1=model.nonlinear_component_index(3);if(b1<0)throw std::runtime_error("B1 bottom hinge not found");
        const double q=model.nonlinear_basis().column_dot(b1,u);
        double pr=0,pa=0,P=0,M=0;
        if(pm){const auto h=model.pm_interaction_snapshot(3,u,state);pr=h.plastic_rotation;pa=h.plastic_axial_deformation;P=h.compression;M=h.moment;}
        else {const auto h=model.nonlinear_component_snapshot(3,u,state);M=h.force;P=model.elastic_element_response(mb.column_probes[1].elastic_element_id,u).axial_compression;}
        r.b1_bottom_max_abs_rotation=std::max(r.b1_bottom_max_abs_rotation,std::abs(q));
        r.b1_bottom_max_abs_plastic_axial=std::max(r.b1_bottom_max_abs_plastic_axial,std::abs(pa));
        if(std::abs(pr)>r.b1_bottom_max_abs_plastic_rotation){r.b1_bottom_max_abs_plastic_rotation=std::abs(pr);r.b1_bottom_time_at_max_plastic_rotation=time;r.b1_bottom_compression_at_max_plastic_rotation=P;r.b1_bottom_moment_at_max_plastic_rotation=M;}
        if(output){r.time.push_back(time);r.drift_history_pct.push_back(d);r.b1_rotation_history.push_back(q);r.b1_plastic_rotation_history.push_back(pr);r.b1_plastic_axial_history.push_back(pa);r.b1_p_history.push_back(P);r.b1_m_history.push_back(M);}
    };
    ro.accepted_substep_state_observer=[&](std::size_t,std::size_t,int,double time,double,const std::vector<double>& u,const std::vector<double>&,const std::vector<double>&,const std::vector<double>& state){observe(time,u,state,false);};
    ro.accepted_state_observer=[&](std::size_t,double time,double,const std::vector<double>& u,const std::vector<double>&,const std::vector<double>&,const std::vector<double>& state){observe(time,u,state,true);};
    r.analysis=run_newmark_robust(model,motion.ag,motion.dt,LinearStrategy::FullFactorization,ro);
    if(pm && r.analysis.termination!=AnalysisTermination::Completed && !r.analysis.final_displacement.empty()){
        std::cerr<<"PHASE9F PM failure t="<<r.analysis.termination_time<<" step="<<r.analysis.termination_step<<" hinge states:\n";
        for(int id=1;id<=24;++id){
            try{const auto h=model.pm_interaction_snapshot(id,r.analysis.final_displacement,r.analysis.final_state);
                std::cerr<<" id="<<id<<" P="<<h.compression<<" M="<<h.moment<<" pr="<<h.plastic_rotation<<" pa="<<h.plastic_axial_deformation<<" kp="<<h.positive_plastic_rotation<<" kn="<<h.negative_plastic_rotation<<" K="<<h.tangent[3]<<"\n";
            }catch(...){}
        }
    }
    if(!r.s2_peak_u.empty()){
        std::vector<double> f,t,tr;model.internal_force_and_tangent(r.s2_peak_u,r.s2_peak_state,f,t,tr);auto K=model.effective_state_tangent_matrix_with_state(r.s2_peak_u,t,r.s2_peak_state,0.0,0.0);
        try{r.s2_story_tangent_diag=story_tangent_diag(model,K);}catch(...){ }
        auto modes=modal_analysis_from_stiffness(model,K,3);for(auto& md:modes){r.s2_mode_periods.push_back(md.period);auto sh=model.story_response_values(md.shape);std::array<double,3> q{{0,0,0}};if(sh.size()>=3){double den=std::abs(sh[2])>1e-14?sh[2]:1.0;for(int i=0;i<3;++i)q[static_cast<std::size_t>(i)]=sh[static_cast<std::size_t>(i)]/den;}r.s2_mode_story_shapes.push_back(q);}
    }
    return r;
}

void write_history(const Phase9FRun& a,const Phase9FRun& b,const std::string& path){
    std::ofstream o(path);if(!o)throw std::runtime_error("cannot write history");o<<"case,time_s,s1_pct,s2_pct,s3_pct,b1_total_rotation_rad,b1_plastic_rotation_rad,b1_plastic_axial_in,b1_P_kip,b1_M_kip_in\n";
    auto emit=[&](const Phase9FRun& r){for(std::size_t i=0;i<r.time.size();++i)o<<r.name<<','<<r.time[i]<<','<<r.drift_history_pct[i][0]<<','<<r.drift_history_pct[i][1]<<','<<r.drift_history_pct[i][2]<<','<<r.b1_rotation_history[i]<<','<<r.b1_plastic_rotation_history[i]<<','<<r.b1_plastic_axial_history[i]<<','<<r.b1_p_history[i]<<','<<r.b1_m_history[i]<<'\n';};emit(a);emit(b);
}

void write_json(const Phase9FRun& base,const Phase9FRun& pm,const ModelBuild& mb,const std::string& path){
    std::ofstream o(path);if(!o)throw std::runtime_error("cannot write json");o<<std::setprecision(12);
    auto num=[&](double x){if(std::isfinite(x))o<<x;else o<<"null";};
    auto arr3=[&](const std::array<double,3>& a){o<<'[';num(a[0]);o<<',';num(a[1]);o<<',';num(a[2]);o<<']';};
    auto run=[&](const Phase9FRun& r){o<<"{\n      \"termination\":"<<static_cast<int>(r.analysis.termination)<<",\n      \"termination_reason\":\""<<r.analysis.termination_reason<<"\",\n      \"peak_story_drift_percent\":";arr3(r.peak_drift_pct);o<<",\n      \"time_at_peak_story_drift_s\":";arr3(r.peak_time_s);o<<",\n      \"simultaneous_at_s2_peak_percent\":";arr3(r.simultaneous_pct[1]);o<<",\n      \"b1_bottom_max_abs_total_rotation_rad\":"<<r.b1_bottom_max_abs_rotation<<",\n      \"b1_bottom_max_abs_plastic_rotation_rad\":"<<r.b1_bottom_max_abs_plastic_rotation<<",\n      \"b1_bottom_max_abs_plastic_axial_in\":"<<r.b1_bottom_max_abs_plastic_axial<<",\n      \"b1_bottom_time_at_max_plastic_rotation_s\":"<<r.b1_bottom_time_at_max_plastic_rotation<<",\n      \"b1_bottom_P_at_max_plastic_rotation_kip\":"<<r.b1_bottom_compression_at_max_plastic_rotation<<",\n      \"b1_bottom_M_at_max_plastic_rotation_kip_in\":"<<r.b1_bottom_moment_at_max_plastic_rotation<<",\n      \"s2_peak_interstory_tangent_diagonal_kip_per_in\":";arr3(r.s2_story_tangent_diag);o<<",\n      \"max_column_compression_kip\":[";for(std::size_t i=0;i<r.max_compression_kip.size();++i){if(i)o<<',';o<<r.max_compression_kip[i];}o<<"]\n    }";};
    o<<"{\n  \"phase\":\"9F_PERFORM_COMPATIBLE_PM_EXPERIMENT\",\n  \"input_status\":\"PRE-VALIDATION / deterministic proxy, not recorded DT1 table motion\",\n  \"mechanics_note\":\"True zero-length two-force P-M hinge with associative normal flow and return mapping on the documented PERFORM concrete-type smooth P-M surface. Surface parameters are fit only to transparent section mechanics; no dynamic EDP is used. Full Perform Mroz translating-surface hardening is not yet implemented, so this remains a P-M mechanics isolation rather than exact Perform3D reproduction.\",\n  \"common_flexural_scale\":"<<mb.scale<<",\n  \"common_rayleigh_alpha\":"<<mb.alpha_m<<",\n  \"common_rayleigh_beta\":"<<mb.beta_k<<",\n  \"research_baseline\":";run(base);o<<",\n  \"perform_compatible_pm\":";run(pm);o<<"\n}\n";
}

}

int main(int argc,char** argv){
    try{
        const std::string motion_path=argc>1?argv[1]:"validation/uc_berkeley_3story/proxy_dt1_motion.csv";
        const std::string out_json=argc>2?argv[2]:"validation/uc_berkeley_3story/phase9/phase9f_pm_experiment.json";
        const std::string out_csv=argc>3?argv[3]:"validation/uc_berkeley_3story/phase9/phase9f_pm_history.csv";
        const auto motion=read_motion(motion_path);const auto gravity=gravity_column_preloads();
        constexpr double slope=-0.005,residual=0.20,acc=0.55;
        // Frozen calibration recovered from this source line's Phase-8 calibration.
        // Use identical elastic calibration for the research and PM mechanics cases.
        constexpr double flex_scale=0.934616;
        constexpr double target_T=0.48;
        const double w=2.0*3.14159265358979323846/target_T;
        const double alpha=2.0*0.025*w,beta=2.0*0.005/w;
        auto rawb=build_frame_phase8(flex_scale,gravity,true,true,false,true,false,3,alpha,beta,true,true,slope,residual,acc,-1.0,false);
        const auto bm=modal_analysis(rawb.model,1);ModelBuild base{std::move(rawb.model),flex_scale,bm.empty()?0.0:bm[0].period,alpha,beta,rawb.probes,true};
        auto rawp=build_frame_phase8(flex_scale,gravity,true,true,false,true,false,3,alpha,beta,true,true,slope,residual,acc,-1.0,true);
        const auto pmodes=modal_analysis(rawp.model,1);ModelBuild pm{std::move(rawp.model),flex_scale,pmodes.empty()?0.0:pmodes[0].period,alpha,beta,rawp.probes,true};
        std::cerr<<"PHASE9F frozen_scale="<<flex_scale<<" baseT="<<base.period<<" pmT="<<pm.period<<" baseDOF="<<base.model.dof()<<" pmDOF="<<pm.model.dof()<<"\n";
        const bool pm_only=argc>4 && std::string(argv[4])=="pm-only";
        Phase9FRun rb;
        if(pm_only){
            rb.name="research_baseline_reference";rb.peak_drift_pct={3.27759621,4.2482242732,2.176073734};rb.peak_time_s={30.37,5.18,12.31};rb.simultaneous_pct[1]={2.96455287118,4.24822427323,2.11445484447};
        }else{
            rb=run_case("research_baseline",base,motion,false);std::cerr<<"PHASE9F baseline term="<<static_cast<int>(rb.analysis.termination)<<" d="<<rb.peak_drift_pct[0]<<","<<rb.peak_drift_pct[1]<<","<<rb.peak_drift_pct[2]<<"\n";
        }
        auto rp=run_case("perform_compatible_pm",pm,motion,true);std::cerr<<"PHASE9F PM term="<<static_cast<int>(rp.analysis.termination)<<" reason="<<rp.analysis.termination_reason<<" d="<<rp.peak_drift_pct[0]<<","<<rp.peak_drift_pct[1]<<","<<rp.peak_drift_pct[2]<<" B1pr="<<rp.b1_bottom_max_abs_plastic_rotation<<" B1pa="<<rp.b1_bottom_max_abs_plastic_axial<<"\n";
        const bool summary_only=argc>5 && std::string(argv[5])=="summary-only";
        if(summary_only){return 0;}
        write_json(rb,rp,base,out_json);if(!pm_only)write_history(rb,rp,out_csv);else {std::ofstream o(out_csv);o<<"case,time_s,s1_pct,s2_pct,s3_pct,b1_total_rotation_rad,b1_plastic_rotation_rad,b1_plastic_axial_in,b1_P_kip,b1_M_kip_in\n";for(std::size_t i=0;i<rp.time.size();++i)o<<rp.name<<','<<rp.time[i]<<','<<rp.drift_history_pct[i][0]<<','<<rp.drift_history_pct[i][1]<<','<<rp.drift_history_pct[i][2]<<','<<rp.b1_rotation_history[i]<<','<<rp.b1_plastic_rotation_history[i]<<','<<rp.b1_plastic_axial_history[i]<<','<<rp.b1_p_history[i]<<','<<rp.b1_m_history[i]<<'\n';}
        std::cout<<out_json<<"\n"<<out_csv<<"\n";return 0;
    }catch(const std::exception& e){std::cerr<<"Phase9F error: "<<e.what()<<"\n";return 2;}
}
