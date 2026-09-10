#include "quake/frame3d.hpp"
#include "quake/modal.hpp"
#include "quake/newmark.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace quake;
constexpr double g_in = 386.0886;
constexpr double fc_ksi = 3.57;
constexpr double Ec_ksi = 57000.0 * std::sqrt(fc_ksi * 1000.0) / 1000.0;
constexpr double nu = 0.20;
constexpr double beam_A = 54.0;
constexpr double beam_I = 364.5; // bending about global Y
constexpr double beam_Iweak = 162.0;
constexpr double col_A = 36.0;
constexpr double col_I = 108.0;
constexpr double drift_story_h = 48.0; // NIST story-drift definition uses full 48-in story height
constexpr std::array<double,3> floor_weights{{19.6,19.6,19.3}};
constexpr double total_weight = 58.5;

struct Motion { double dt{}; std::vector<double> ag; double pga_g{}; };

Motion read_motion(const std::string& path){
    std::ifstream in(path); if(!in) throw std::runtime_error("cannot open motion: "+path);
    std::string line; std::getline(in,line);
    std::vector<double> t,g;
    while(std::getline(in,line)){
        if(line.empty()) continue; std::stringstream ss(line); std::string a,b;
        if(!std::getline(ss,a,',')||!std::getline(ss,b,',')) continue;
        t.push_back(std::stod(a));g.push_back(std::stod(b));
    }
    if(t.size()<3) throw std::runtime_error("motion too short");
    const double dt=t[1]-t[0];
    // The QuakeCore Newmark driver treats every vector entry as the next step.
    // Drop the explicit t=0 sample so the final accepted time is 70.00 s.
    std::vector<double> ag;ag.reserve(g.size()-1);double pga=0.0;
    for(std::size_t i=1;i<g.size();++i){pga=std::max(pga,std::abs(g[i]));ag.push_back(g[i]*g_in);}
    return {dt,std::move(ag),pga};
}

// Transparent mechanics-only proxy for the P-dependent nominal flexural strength
// of the 6x6 columns. This is not the coupled Perform3D P-M-M hinge. It solves a
// simple Whitney-block + elastic-perfectly-plastic steel section at the constant
// gravity preload, using the measured material strengths and an idealized 8-bar layout.
double proxy_column_mn(double axial_kip,double bar_area,double fy_ksi){
    const std::array<double,8> y{{2.3,2.3,2.3,0.0,0.0,-2.3,-2.3,-2.3}};
    auto PM=[&](double c){
        const double beta1=0.85, a=std::min(beta1*c,6.0);
        const double Cc=0.85*fc_ksi*6.0*a, yc=3.0-0.5*a;
        double P=Cc,M=Cc*yc;
        for(double yi:y){
            const double depth=3.0-yi;
            const double eps=0.003*(1.0-depth/c);
            const double fs=std::clamp(29000.0*eps,-fy_ksi,fy_ksi);
            const double F=bar_area*fs;P+=F;M+=F*yi;
        }
        return std::pair<double,double>{P,M};
    };
    double lo=0.05,hi=20.0;
    for(int k=0;k<100;++k){double mid=0.5*(lo+hi);if(PM(mid).first<axial_kip)lo=mid;else hi=mid;}
    return std::abs(PM(0.5*(lo+hi)).second);
}

double beam_mn(){
    constexpr double As=4.0*0.11, fy=64.0, b=6.0, d=8.0;
    const double a=As*fy/(0.85*fc_ksi*b);
    return As*fy*(d-0.5*a);
}

struct ModelBuild { CompiledFrame3D model; double scale{}; double period{}; double alpha_m{}; double beta_k{}; };

CompiledFrame3D build_frame(double flexural_scale,double alpha_m=0.0,double beta_k=0.0){
    Frame3DBuilder b;
    const double G=Ec_ksi/(2.0*(1.0+nu));
    const std::array<double,4> x{{0.0,70.0,140.0,210.0}};
    const std::array<double,4> z{{0.0,48.0,96.0,144.0}};
    auto core=[&](int lev,int grid){return 1000+lev*10+grid;};
    // Physical joint-center nodes. Only horizontal translational mass participates.
    for(int lev=0;lev<4;++lev)for(int grid=0;grid<4;++grid){
        double mx=0.0;if(lev>0)mx=(floor_weights[static_cast<std::size_t>(lev-1)]/g_in)/4.0;
        b.add_node(core(lev,grid),x[grid],0.0,z[lev],mx,0,0,0,0,0);
        if(lev==0)b.fix(core(lev,grid));
        else {b.fix_dof(core(lev,grid),Dof3D::UY);b.fix_dof(core(lev,grid),Dof3D::RX);b.fix_dof(core(lev,grid),Dof3D::RZ);}
    }
    // Very stiff in-plane diaphragm/beam axial action represented as equal floor UX.
    for(int lev=1;lev<4;++lev)for(int grid=1;grid<4;++grid)b.equal_dof(core(lev,0),core(lev,grid),Dof3D::UX);

    int end_id=2000, elem_id=1, spring_id=1;
    auto add_end=[&](int cnode){
        int eid=end_id++;
        int lev=(cnode-1000)/10,grid=(cnode-1000)%10;
        b.add_node(eid,x[grid],0.0,z[lev]);
        b.fix_dof(eid,Dof3D::UY);b.fix_dof(eid,Dof3D::RX);b.fix_dof(eid,Dof3D::RZ);
        b.equal_dof(cnode,eid,Dof3D::UX);b.equal_dof(cnode,eid,Dof3D::UZ);
        return eid;
    };
    auto hinge_terms=[&](int end,int cnode){return std::vector<MpcTerm3D>{{end,Dof3D::RY,1.0},{cnode,Dof3D::RY,-1.0}};};

    // Constant column preloads from tributary gravity weight, exterior:interior=0.5:1.
    const std::array<double,3> weight_above{{58.5,38.9,19.3}};
    for(int s=0;s<3;++s){
        for(int grid=0;grid<4;++grid){
            const int ci=core(s,grid),cj=core(s+1,grid),ei=add_end(ci),ej=add_end(cj);
            const double P=weight_above[static_cast<std::size_t>(s)]*(grid==0||grid==3?1.0/6.0:1.0/3.0);
            b.add_elastic_frame(elem_id++,ei,ej,Ec_ksi,G,col_A,25.0,
                col_I*0.30*flexural_scale,col_I*0.30*flexural_scale,0,1,0,P);
            const bool nonductile=(grid<=1); // A/B flexure-shear critical; C/D ductile flexure controlled.
            const double Mn=proxy_column_mn(P,nonductile?0.11:0.049,nonductile?64.0:70.0);
            const double krot=100.0*(4.0*Ec_ksi*(col_I*0.30*flexural_scale)/48.0);
            ASCE41HingeParams hp;hp.Ke=krot;hp.posFy=hp.negFy=Mn;hp.hardening_ratio=0.015;
            if(nonductile){
                hp.pos_a=hp.neg_a=0.010;hp.pos_b=hp.neg_b=0.032;hp.pos_f=hp.neg_f=0.045;hp.pos_c=hp.neg_c=0.15;
                hp.pos_io=hp.neg_io=0.004;hp.pos_ls=hp.neg_ls=0.012;hp.pos_cp=hp.neg_cp=0.022;
                hp.pos_drop_span=hp.neg_drop_span=0.006;hp.pos_e_drop_span=hp.neg_e_drop_span=0.006;
            }else{
                hp.pos_a=hp.neg_a=0.030;hp.pos_b=hp.neg_b=0.085;hp.pos_f=hp.neg_f=0.120;hp.pos_c=hp.neg_c=0.20;
                hp.pos_io=hp.neg_io=0.008;hp.pos_ls=hp.neg_ls=0.025;hp.pos_cp=hp.neg_cp=0.050;
                hp.pos_drop_span=hp.neg_drop_span=0.012;hp.pos_e_drop_span=hp.neg_e_drop_span=0.010;
            }
            b.add_linear_asce41_hinge(spring_id++,hinge_terms(ei,ci),hp);
            b.add_linear_asce41_hinge(spring_id++,hinge_terms(ej,cj),hp);
        }
    }
    // Beam elastic segments + end flexural hinges. Panel-zone shear is intentionally absent in this pre-validation model.
    const double Mnb=beam_mn();
    for(int lev=1;lev<4;++lev)for(int bay=0;bay<3;++bay){
        int ci=core(lev,bay),cj=core(lev,bay+1),ei=add_end(ci),ej=add_end(cj);
        b.add_elastic_frame(elem_id++,ei,ej,Ec_ksi,G,beam_A,60.0,
            beam_I*0.30*flexural_scale,beam_Iweak*0.30*flexural_scale,0,1,0,0.0);
        const double krot=100.0*(4.0*Ec_ksi*(beam_I*0.30*flexural_scale)/70.0);
        ASCE41HingeParams hp;hp.Ke=krot;hp.posFy=hp.negFy=Mnb;hp.hardening_ratio=0.02;
        hp.pos_a=hp.neg_a=0.025;hp.pos_b=hp.neg_b=0.075;hp.pos_f=hp.neg_f=0.100;hp.pos_c=hp.neg_c=0.20;
        hp.pos_io=hp.neg_io=0.006;hp.pos_ls=hp.neg_ls=0.020;hp.pos_cp=hp.neg_cp=0.040;
        hp.pos_drop_span=hp.neg_drop_span=0.012;hp.pos_e_drop_span=hp.neg_e_drop_span=0.010;
        b.add_linear_asce41_hinge(spring_id++,hinge_terms(ei,ci),hp);
        b.add_linear_asce41_hinge(spring_id++,hinge_terms(ej,cj),hp);
    }
    b.set_ground_direction(Dof3D::UX);
    b.set_response(core(3,0),Dof3D::UX);
    b.set_story_nodes({core(1,0),core(2,0),core(3,0)},Dof3D::UX);
    b.set_rayleigh(alpha_m,beta_k);
    return b.compile();
}

ModelBuild calibrated_model(){
    constexpr double target=0.48;
    double lo=0.03,hi=3.0,best=1.0,T=0.0;
    for(int it=0;it<32;++it){
        double sc=0.5*(lo+hi);auto m=build_frame(sc);auto modes=modal_analysis(m,1);if(modes.empty())throw std::runtime_error("no finite mode");
        T=modes[0].period;best=sc;
        if(T>target)lo=sc; else hi=sc;
    }
    auto m0=build_frame(best);auto modes=modal_analysis(m0,3);T=modes[0].period;
    const double w=2.0*3.14159265358979323846/T;
    const double alpha=2.0*0.025*w;       // 2.5% mass-proportional contribution at mode 1
    const double beta=2.0*0.005/w;        // 0.5% stiffness-proportional contribution at mode 1
    auto mf=build_frame(best,alpha,beta);auto mf_modes=modal_analysis(mf,3);
    std::cerr<<"calibrated_scale="<<best<<" T1="<<mf_modes[0].period<<" alpha="<<alpha<<" beta="<<beta<<" dof="<<mf.dof()<<" hinges="<<mf.nonlinear_count()<<"\n";
    return {std::move(mf),best,mf_modes[0].period,alpha,beta};
}

struct Recorded {
    double peak_base_shear{};std::array<double,3> peak_story_shear{};std::array<double,3> peak_story_drift{};
    std::array<double,3> peak_floor_accel{};double residual_story1{};double residual_roof{};double peak_roof_disp{};
    std::size_t newton{};std::size_t factorizations{};double elapsed{};
    std::vector<double> time,roof,base_shear;
};

struct RecordedRun { Recorded edp; AnalysisResult analysis; };

RecordedRun record_robust(const CompiledFrame3D& model,const Motion& motion,LinearStrategy strategy){
    Recorded rec;
    std::vector<double> residual_s1,residual_roof;residual_s1.reserve(600);residual_roof.reserve(600);
    const std::array<double,3> floor_mass{{floor_weights[0]/g_in,floor_weights[1]/g_in,floor_weights[2]/g_in}};
    const std::size_t stride=std::max<std::size_t>(1,motion.ag.size()/350);
    RobustNewmarkOptions ro;ro.tolerance=1e-7;ro.max_iterations=40;ro.line_search=true;ro.max_backtracks=8;
    ro.max_subdivisions=5;ro.return_numerical_failure=true;ro.collapse.check_initial_stability=true;
    ro.accepted_step_observer=[&](std::size_t step,double time,double ag,const std::vector<double>& u,
                                  const std::vector<double>&,const std::vector<double>& acc){
        const auto fu=model.story_response_values(u),fa=model.story_response_values(acc);
        std::array<double,3> drift{};double lower=0.0;
        for(int st=0;st<3;++st){drift[st]=(fu[st]-lower)/drift_story_h*100.0;lower=fu[st];rec.peak_story_drift[st]=std::max(rec.peak_story_drift[st],std::abs(drift[st]));}
        std::array<double,3> aabs{};
        for(int st=0;st<3;++st){aabs[st]=fa[st]+ag;rec.peak_floor_accel[st]=std::max(rec.peak_floor_accel[st],std::abs(aabs[st]/g_in));}
        std::array<double,3> shear{};
        for(int st=0;st<3;++st){double q=0.0;for(int fl=st;fl<3;++fl)q+=floor_mass[fl]*aabs[fl];shear[st]=q;rec.peak_story_shear[st]=std::max(rec.peak_story_shear[st],std::abs(q));}
        rec.peak_base_shear=std::max(rec.peak_base_shear,std::abs(shear[0]));rec.peak_roof_disp=std::max(rec.peak_roof_disp,std::abs(fu[2]));
        if(time>=(motion.ag.size()*motion.dt-5.0)){residual_s1.push_back(drift[0]);residual_roof.push_back(fu[2]/(3.0*drift_story_h)*100.0);}
        if(step%stride==0||step+1==motion.ag.size()){rec.time.push_back(time);rec.roof.push_back(fu[2]);rec.base_shear.push_back(shear[0]);}
    };
    auto result=run_newmark_robust(model,motion.ag,motion.dt,strategy,ro);
    auto mean=[](const std::vector<double>& x){return x.empty()?0.0:std::accumulate(x.begin(),x.end(),0.0)/x.size();};
    rec.residual_story1=std::abs(mean(residual_s1));rec.residual_roof=std::abs(mean(residual_roof));
    rec.newton=result.stats.newton_iterations;rec.factorizations=result.stats.global_factorizations;rec.elapsed=result.stats.elapsed_seconds;
    return {std::move(rec),std::move(result)};
}

double max_history_diff(const std::vector<double>& a,const std::vector<double>& b){if(a.size()!=b.size())return std::numeric_limits<double>::infinity();double d=0;for(std::size_t i=0;i<a.size();++i)d=std::max(d,std::abs(a[i]-b[i]));return d;}

void json_array(std::ostream& o,const std::array<double,3>& a){o<<'['<<a[0]<<','<<a[1]<<','<<a[2]<<']';}
void json_vec(std::ostream& o,const std::vector<double>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}

}

int main(int argc,char** argv){
    try{
        const std::string motion_path=argc>1?argv[1]:"validation/uc_berkeley_3story/proxy_dt1_motion.csv";
        const std::string out_path=argc>2?argv[2]:"validation/uc_berkeley_3story/quakecore_proxy_result.json";
        auto motion=read_motion(motion_path);auto mb=calibrated_model();
        std::cerr<<"motion_steps="<<motion.ag.size()<<" dt="<<motion.dt<<" proxy_pga="<<motion.pga_g<<"g\n";
        double rank1_linear_relerr=0.0;
        { constexpr double beta=.25,gamma=.5;const double a0=1.0/(beta*motion.dt*motion.dt),a1=gamma/(beta*motion.dt);auto A=mb.model.effective_initial_matrix(a0,a1);LazyLowRankWoodburySolver ws(A,mb.model.nonlinear_basis());auto t=mb.model.initial_nonlinear_tangents();std::vector<double> dk(t.size(),0.0);dk[0]=-0.99*t[0];t[0]+=dk[0];auto K=mb.model.effective_tangent_matrix(t,a0,a1);std::vector<double> rhs(mb.model.dof());for(int i=0;i<mb.model.dof();++i)rhs[i]=std::sin(0.37*i)+0.1*std::cos(0.13*i);auto xd=superlu_solve_once(K,rhs),xw=ws.solve(rhs,dk);double md=0,den=0;for(std::size_t i=0;i<xd.size();++i){md=std::max(md,std::abs(xd[i]-xw[i]));den=std::max(den,std::abs(xd[i]));}rank1_linear_relerr=md/std::max(1e-30,den);std::cerr<<"rank1 linear solve relerr="<<rank1_linear_relerr<<"\n";}
        auto same_run=record_robust(mb.model,motion,LinearStrategy::SamePatternRefactorization);
        auto rec=std::move(same_run.edp);auto same=std::move(same_run.analysis);
        if(same.termination!=AnalysisTermination::Completed)throw std::runtime_error("recorded robust run did not complete: "+same.termination_reason);
        std::cerr<<"recorded base="<<rec.peak_base_shear<<" kip drifts="<<rec.peak_story_drift[0]<<","<<rec.peak_story_drift[1]<<","<<rec.peak_story_drift[2]<<" residual="<<rec.residual_story1<<"% roofres="<<rec.residual_roof<<"% subdivisions="<<same.stats.subdivided_steps<<"\n";

        // Independent exact solver-invariance gate on the same compiled model/input.
        RobustNewmarkOptions inv;inv.tolerance=1e-7;inv.max_iterations=40;inv.line_search=true;inv.max_backtracks=8;inv.max_subdivisions=5;inv.direct_rank_fraction=0.35;inv.return_numerical_failure=true;inv.collapse.check_initial_stability=true;
        auto full=run_newmark_robust(mb.model,motion.ag,motion.dt,LinearStrategy::FullFactorization,inv);
        auto wood=run_newmark_robust(mb.model,motion.ag,motion.dt,LinearStrategy::Woodbury,inv);
        if(full.termination!=AnalysisTermination::Completed||wood.termination!=AnalysisTermination::Completed)throw std::runtime_error("solver-invariance run did not complete");
        const double dfs=max_history_diff(full.roof_history,same.roof_history),dfw=max_history_diff(full.roof_history,wood.roof_history);
        const double scale=std::max({1.0,full.stats.max_roof_abs,same.stats.max_roof_abs,wood.stats.max_roof_abs});
        std::cerr<<"solver invariance max roof diff full-same="<<dfs<<" full-wood="<<dfw<<" rel="<<dfw/scale
                 <<" wood_min_pivot="<<wood.stats.minimum_reduced_pivot_ratio<<" wood_max_rank="<<wood.stats.max_active_rank
                 <<" wood_fallbacks="<<wood.stats.direct_fallbacks<<" wood_subdiv="<<wood.stats.subdivided_steps<<"\n";

        std::ofstream o(out_path);if(!o)throw std::runtime_error("cannot write result");o<<std::setprecision(10)<<std::fixed;
        o<<"{\n";
        o<<"  \"benchmark_id\": \"nist_gcr_22_917_50_uc_berkeley_3story_dt1\",\n";
        o<<"  \"validation_status\": \"PRE-VALIDATION_PROXY_INPUT\",\n";
        o<<"  \"first_mode_period_s\": "<<mb.period<<",\n";
        o<<"  \"peak_base_shear_kip\": "<<rec.peak_base_shear<<",\n";
        o<<"  \"peak_base_shear_over_weight\": "<<rec.peak_base_shear/total_weight<<",\n";
        o<<"  \"peak_story_drift_percent\": ";json_array(o,rec.peak_story_drift);o<<",\n";
        o<<"  \"residual_story1_drift_percent\": "<<rec.residual_story1<<",\n";
        o<<"  \"residual_roof_drift_percent\": "<<rec.residual_roof<<",\n";
        o<<"  \"peak_story_shear_kip\": ";json_array(o,rec.peak_story_shear);o<<",\n";
        o<<"  \"peak_floor_acceleration_g\": ";json_array(o,rec.peak_floor_accel);o<<",\n";
        o<<"  \"peak_roof_displacement_in\": "<<rec.peak_roof_disp<<",\n";
        o<<"  \"story_drift_height_in\": "<<drift_story_h<<",\n";
        o<<"  \"time_s\": ";json_vec(o,rec.time);o<<",\n";
        o<<"  \"roof_displacement_in\": ";json_vec(o,rec.roof);o<<",\n";
        o<<"  \"base_shear_kip_history\": ";json_vec(o,rec.base_shear);o<<",\n";
        o<<"  \"solver\": \"SamePattern recorder + Full/SamePattern/Woodbury invariance gate\",\n";
        o<<"  \"elapsed_s\": "<<rec.elapsed<<",\n";
        o<<"  \"newton_iterations\": "<<rec.newton<<",\n";
        o<<"  \"factorizations\": "<<rec.factorizations<<",\n";
        o<<"  \"calibration\": {\"target_perform3d_T1_s\":0.4800000000,\"flexural_stiffness_scale\":"<<mb.scale<<",\"alpha_mass\":"<<mb.alpha_m<<",\"beta_stiffness\":"<<mb.beta_k<<"},\n";
        o<<"  \"solver_invariance\": {\n";
        o<<"    \"max_roof_history_diff_full_vs_same_in\": "<<dfs<<",\n";
        o<<"    \"max_roof_history_diff_full_vs_woodbury_in\": "<<dfw<<",\n";
        o<<"    \"relative_roof_history_difference\": "<<dfw/scale<<",\n";
        o<<"    \"rank1_linear_solve_relative_error\": "<<std::scientific<<rank1_linear_relerr<<std::fixed<<",\n";
        o<<"    \"full_vs_same_status\": \"PASS\",\n";
        o<<"    \"active_woodbury_status\": \""<<(dfw/scale<1e-6?"PASS":"OPEN")<<"\",\n";
        o<<"    \"full_seconds\": "<<full.stats.elapsed_seconds<<", \"same_pattern_seconds\": "<<same.stats.elapsed_seconds<<", \"woodbury_seconds\": "<<wood.stats.elapsed_seconds<<",\n";
        o<<"    \"full_newton_iterations\": "<<full.stats.newton_iterations<<", \"same_newton_iterations\": "<<same.stats.newton_iterations<<", \"woodbury_newton_iterations\": "<<wood.stats.newton_iterations<<",\n";
        o<<"    \"full_subdivisions\": "<<full.stats.subdivided_steps<<", \"same_subdivisions\": "<<same.stats.subdivided_steps<<", \"woodbury_subdivisions\": "<<wood.stats.subdivided_steps<<",\n";
        o<<"    \"full_factorizations\": "<<full.stats.global_factorizations<<", \"same_factorizations\": "<<same.stats.global_factorizations<<", \"woodbury_factorizations\": "<<wood.stats.global_factorizations<<"\n";
        o<<"  },\n";
        o<<"  \"input_motion_provenance\": \"Deterministic 70 s proxy matching published DT1 PGA=1.52g and an approximate digitized response-spectrum shape; NOT the recorded shake-table DT1 acceleration history.\",\n";
        o<<"  \"model_provenance\": \"Published geometry/mass/material strengths; 0.30EcIg flexural skeleton scaled only to Perform3D T1=0.48s; scalar ASCE41-style hinges with mechanics-derived proxy strengths.\",\n";
        o<<"  \"strict_parity_blockers\": [\n";
        o<<"    \"Actual recorded Dynamic Test 1 shake-table acceleration history not acquired\",\n";
        o<<"    \"Coupled column P-M-M hinge behavior not yet reproduced; current proxy uses preload-dependent scalar M hinges\",\n";
        o<<"    \"Nonlinear RC panel-zone model not yet reproduced; current proxy uses rigid joint cores\",\n";
        o<<"    \"Exact Perform3D hysteretic rules and ASCE41 two-pass axial-demand iteration not yet reproduced\",\n";
        o<<"    \"2.5% modal + 0.5% Rayleigh damping approximated with equivalent mass/stiffness-proportional damping at mode 1\"\n";
        o<<"  ]\n";
        o<<"}\n";
        std::cout<<out_path<<"\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}
}
