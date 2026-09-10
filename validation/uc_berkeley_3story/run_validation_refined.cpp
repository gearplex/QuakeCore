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

struct ColumnProbe {
    int bottom_node{}, top_node{};
    double gravity_preload{};
    double axial_k{};
    int story{}, grid{};
};

using ColumnP = std::array<double,12>;

ColumnP gravity_column_preloads(){
    ColumnP p{};
    const std::array<double,3> weight_above{{58.5,38.9,19.3}};
    for(int st=0;st<3;++st)for(int grid=0;grid<4;++grid)
        p[static_cast<std::size_t>(st*4+grid)] = weight_above[static_cast<std::size_t>(st)]*(grid==0||grid==3?1.0/6.0:1.0/3.0);
    return p;
}

struct ModelBuild {
    CompiledFrame3D model;
    double scale{}; double period{}; double alpha_m{}; double beta_k{};
    std::array<ColumnProbe,12> column_probes{};
    bool panel_zones{};
};

ASCE41HingeParams panel_zone_params(bool interior){
    // Mechanics-based joint shear proxy: gamma is the beam/column relative
    // rotation, Vj = G*A_j*gamma, and the generalized joint moment is Vj*d_eff.
    // Strength uses a conventional sqrt(f'c) joint-shear form with a larger
    // coefficient for interior joints. These are transparent research values,
    // not hidden calibration to the published building EDPs.
    const double G=Ec_ksi/(2.0*(1.0+nu));
    const double joint_area=6.0*6.0;
    const double lever=8.0; // approximate beam effective depth
    const double ktheta=G*joint_area*lever;
    const double coeff=interior?20.0:15.0;
    const double vj_y=coeff*std::sqrt(fc_ksi*1000.0)*joint_area/1000.0; // psi -> kip
    const double my=vj_y*lever;
    ASCE41HingeParams hp; hp.Ke=ktheta; hp.posFy=hp.negFy=my; hp.hardening_ratio=0.02;
    hp.pos_a=hp.neg_a=0.018; hp.pos_b=hp.neg_b=0.060; hp.pos_f=hp.neg_f=0.090; hp.pos_c=hp.neg_c=0.70;
    hp.pos_io=hp.neg_io=0.004; hp.pos_ls=hp.neg_ls=0.015; hp.pos_cp=hp.neg_cp=0.035;
    hp.pos_drop_span=hp.neg_drop_span=0.010; hp.pos_e_drop_span=hp.neg_e_drop_span=0.010;
    return hp;
}

struct BuiltRaw { CompiledFrame3D model; std::array<ColumnProbe,12> probes{}; };

BuiltRaw build_frame_refined(double flexural_scale,const ColumnP& design_p,bool panel_zones,
                             double alpha_m=0.0,double beta_k=0.0){
    Frame3DBuilder b;
    const double G=Ec_ksi/(2.0*(1.0+nu));
    const std::array<double,4> x{{0.0,70.0,140.0,210.0}};
    const std::array<double,4> z{{0.0,48.0,96.0,144.0}};
    auto ccore=[&](int lev,int grid){return 1000+lev*10+grid;};
    auto bcore=[&](int lev,int grid){return 1500+lev*10+grid;};

    for(int lev=0;lev<4;++lev)for(int grid=0;grid<4;++grid){
        double mx=0.0;if(lev>0)mx=(floor_weights[static_cast<std::size_t>(lev-1)]/g_in)/4.0;
        b.add_node(ccore(lev,grid),x[grid],0.0,z[lev],mx,0,0,0,0,0);
        if(lev==0)b.fix(ccore(lev,grid));
        else {b.fix_dof(ccore(lev,grid),Dof3D::UY);b.fix_dof(ccore(lev,grid),Dof3D::RX);b.fix_dof(ccore(lev,grid),Dof3D::RZ);}
        if(lev>0 && panel_zones){
            b.add_node(bcore(lev,grid),x[grid],0.0,z[lev]);
            b.fix_dof(bcore(lev,grid),Dof3D::UY);b.fix_dof(bcore(lev,grid),Dof3D::RX);b.fix_dof(bcore(lev,grid),Dof3D::RZ);
            b.equal_dof(ccore(lev,grid),bcore(lev,grid),Dof3D::UX);
            b.equal_dof(ccore(lev,grid),bcore(lev,grid),Dof3D::UZ);
        }
    }
    for(int lev=1;lev<4;++lev)for(int grid=1;grid<4;++grid)b.equal_dof(ccore(lev,0),ccore(lev,grid),Dof3D::UX);

    int end_id=2000, elem_id=1, spring_id=1;
    auto add_end=[&](int parent,int lev,int grid){
        int eid=end_id++;b.add_node(eid,x[grid],0.0,z[lev]);
        b.fix_dof(eid,Dof3D::UY);b.fix_dof(eid,Dof3D::RX);b.fix_dof(eid,Dof3D::RZ);
        b.equal_dof(parent,eid,Dof3D::UX);b.equal_dof(parent,eid,Dof3D::UZ);return eid;
    };
    auto hinge_terms=[&](int end,int core){return std::vector<MpcTerm3D>{{end,Dof3D::RY,1.0},{core,Dof3D::RY,-1.0}};};
    const auto gravity=gravity_column_preloads();
    std::array<ColumnProbe,12> probes{};

    // Columns. The scalar hinge capacity is evaluated on a P-M interaction
    // section at the supplied design axial force. A second pass updates those
    // design axial forces from the first nonlinear run, mirroring the NIST
    // axial-demand iteration without fitting any response EDP.
    for(int st=0;st<3;++st){
        for(int grid=0;grid<4;++grid){
            const int ci=ccore(st,grid),cj=ccore(st+1,grid),ei=add_end(ci,st,grid),ej=add_end(cj,st+1,grid);
            const double Pg=gravity[static_cast<std::size_t>(st*4+grid)];
            b.add_elastic_frame(elem_id++,ei,ej,Ec_ksi,G,col_A,25.0,
                col_I*0.30*flexural_scale,col_I*0.30*flexural_scale,0,1,0,Pg);
            const bool nonductile=(grid<=1);
            const double Pd=std::max(0.0,design_p[static_cast<std::size_t>(st*4+grid)]);
            const double Mn=proxy_column_mn(Pd,nonductile?0.11:0.049,nonductile?64.0:70.0);
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
            probes[static_cast<std::size_t>(st*4+grid)]={ci,cj,Pg,Ec_ksi*col_A/48.0,st,grid};
        }
    }

    const double Mnb=beam_mn();
    for(int lev=1;lev<4;++lev)for(int bay=0;bay<3;++bay){
        const int li=panel_zones?bcore(lev,bay):ccore(lev,bay);
        const int lj=panel_zones?bcore(lev,bay+1):ccore(lev,bay+1);
        int ei=add_end(li,lev,bay),ej=add_end(lj,lev,bay+1);
        b.add_elastic_frame(elem_id++,ei,ej,Ec_ksi,G,beam_A,60.0,
            beam_I*0.30*flexural_scale,beam_Iweak*0.30*flexural_scale,0,1,0,0.0);
        const double krot=100.0*(4.0*Ec_ksi*(beam_I*0.30*flexural_scale)/70.0);
        ASCE41HingeParams hp;hp.Ke=krot;hp.posFy=hp.negFy=Mnb;hp.hardening_ratio=0.02;
        hp.pos_a=hp.neg_a=0.025;hp.pos_b=hp.neg_b=0.075;hp.pos_f=hp.neg_f=0.100;hp.pos_c=hp.neg_c=0.20;
        hp.pos_io=hp.neg_io=0.006;hp.pos_ls=hp.neg_ls=0.020;hp.pos_cp=hp.neg_cp=0.040;
        hp.pos_drop_span=hp.neg_drop_span=0.012;hp.pos_e_drop_span=hp.neg_e_drop_span=0.010;
        b.add_linear_asce41_hinge(spring_id++,hinge_terms(ei,li),hp);
        b.add_linear_asce41_hinge(spring_id++,hinge_terms(ej,lj),hp);
    }
    if(panel_zones){
        for(int lev=1;lev<4;++lev)for(int grid=0;grid<4;++grid){
            const auto pz=panel_zone_params(grid==1||grid==2);
            b.add_linear_asce41_hinge(spring_id++,{{ccore(lev,grid),Dof3D::RY,-1.0},{bcore(lev,grid),Dof3D::RY,1.0}},pz);
        }
    }
    b.set_updated_pdelta(true);
    b.set_ground_direction(Dof3D::UX);
    b.set_response(ccore(3,0),Dof3D::UX);
    b.set_story_nodes({ccore(1,0),ccore(2,0),ccore(3,0)},Dof3D::UX);
    b.set_rayleigh(alpha_m,beta_k);
    return {b.compile(),probes};
}

ModelBuild calibrated_model_refined(const ColumnP& design_p,bool panel_zones){
    constexpr double target=0.48;
    double lo=0.03,hi=3.0,best=1.0,T=0.0;
    for(int it=0;it<34;++it){
        const double sc=0.5*(lo+hi);auto raw=build_frame_refined(sc,design_p,panel_zones);auto modes=modal_analysis(raw.model,1);if(modes.empty())throw std::runtime_error("no finite mode");
        T=modes[0].period;best=sc;if(T>target)lo=sc;else hi=sc;
    }
    auto raw0=build_frame_refined(best,design_p,panel_zones);auto modes=modal_analysis(raw0.model,3);T=modes[0].period;
    const double w=2.0*3.14159265358979323846/T;
    const double alpha=2.0*0.025*w;const double beta=2.0*0.005/w;
    auto raw=build_frame_refined(best,design_p,panel_zones,alpha,beta);auto final_modes=modal_analysis(raw.model,3);
    std::cerr<<"refined calibrated panel="<<panel_zones<<" scale="<<best<<" T1="<<final_modes[0].period<<" dof="<<raw.model.dof()<<" hinges="<<raw.model.nonlinear_count()<<"\n";
    return {std::move(raw.model),best,final_modes[0].period,alpha,beta,raw.probes,panel_zones};
}

struct Recorded {
    double peak_base_shear{};std::array<double,3> peak_story_shear{};std::array<double,3> peak_story_drift{};
    std::array<double,3> peak_floor_accel{};double residual_story1{};double residual_roof{};double peak_roof_disp{};
    std::size_t newton{};std::size_t factorizations{};double elapsed{};
    std::array<double,12> max_column_compression{};
    std::array<double,12> peak_panel_zone_rotation{};
    std::vector<double> time,roof,base_shear;
};

struct RecordedRun { Recorded edp; AnalysisResult analysis; };

RecordedRun record_robust(const ModelBuild& mb,const Motion& motion,LinearStrategy strategy){
    const auto& model=mb.model;
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
        for(std::size_t k=0;k<mb.column_probes.size();++k){
            const auto& pr=mb.column_probes[k];
            const int ib=model.reduced_dof(pr.bottom_node,Dof3D::UZ),it=model.reduced_dof(pr.top_node,Dof3D::UZ);
            const double ub=ib>=0?u[static_cast<std::size_t>(ib)]:0.0,ut=it>=0?u[static_cast<std::size_t>(it)]:0.0;
            const double P=pr.gravity_preload-pr.axial_k*(ut-ub);
            rec.max_column_compression[k]=std::max(rec.max_column_compression[k],P);
        }
        if(mb.panel_zones && model.nonlinear_count()>=54){
            const auto& B=model.nonlinear_basis();
            for(int j=0;j<12;++j){
                const double q=B.column_dot(42+j,u);
                rec.peak_panel_zone_rotation[static_cast<std::size_t>(j)]=std::max(rec.peak_panel_zone_rotation[static_cast<std::size_t>(j)],std::abs(q));
            }
        }
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
void json_array12(std::ostream& o,const std::array<double,12>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}
void json_vec(std::ostream& o,const std::vector<double>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}

}

int main(int argc,char** argv){
    try{
        const std::string motion_path=argc>1?argv[1]:"validation/uc_berkeley_3story/proxy_dt1_motion.csv";
        const std::string out_path=argc>2?argv[2]:"validation/uc_berkeley_3story/quakecore_refined_result.json";
        const auto motion=read_motion(motion_path);
        const auto gravity=gravity_column_preloads();
        std::cerr<<"motion_steps="<<motion.ag.size()<<" dt="<<motion.dt<<" proxy_pga="<<motion.pga_g<<"g\n";

        // Scenario 1: current member model without panel-zone compliance, but with
        // updated P-delta and the same P-M section-capacity construction.
        auto m0=calibrated_model_refined(gravity,false);
        auto r0=record_robust(m0,motion,LinearStrategy::SamePatternRefactorization);
        if(r0.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("no-panel pass failed: "+r0.analysis.termination_reason);

        // Scenario 2: add mechanics-derived finite panel-zone compliance.
        auto m1=calibrated_model_refined(gravity,true);
        auto r1=record_robust(m1,motion,LinearStrategy::SamePatternRefactorization);
        if(r1.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("panel pass 1 failed: "+r1.analysis.termination_reason);

        // Scenario 3: NIST-style second property pass. Use the maximum compression
        // observed in pass 1 to re-evaluate each column's P-M section strength.
        ColumnP p2=gravity;
        for(std::size_t i=0;i<p2.size();++i){
            const double upper=0.80*col_A*fc_ksi;
            p2[i]=std::clamp(std::max(gravity[i],r1.edp.max_column_compression[i]),0.0,upper);
        }
        auto m2=calibrated_model_refined(p2,true);
        auto final_run=record_robust(m2,motion,LinearStrategy::SamePatternRefactorization);
        if(final_run.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("panel P-M pass 2 failed: "+final_run.analysis.termination_reason);
        auto& rec=final_run.edp; auto& same=final_run.analysis;

        // Fresh-factorization invariance check on the final refined model. Active
        // Woodbury is intentionally not a required gate for this tiny model because
        // updated geometric blocks can exceed a profitable reduced dimension.
        RobustNewmarkOptions inv;inv.tolerance=1e-7;inv.max_iterations=40;inv.line_search=true;inv.max_backtracks=8;inv.max_subdivisions=5;
        inv.direct_rank_fraction=0.35;inv.return_numerical_failure=true;inv.collapse.check_initial_stability=true;
        auto full=run_newmark_robust(m2.model,motion.ag,motion.dt,LinearStrategy::FullFactorization,inv);
        if(full.termination!=AnalysisTermination::Completed)throw std::runtime_error("full invariance run failed");
        const double dfs=max_history_diff(full.roof_history,same.roof_history);
        const double hscale=std::max({1.0,full.stats.max_roof_abs,same.stats.max_roof_abs});

        auto print=[&](const char* name,const Recorded& r){
            std::cerr<<name<<" base="<<r.peak_base_shear<<" drifts="<<r.peak_story_drift[0]<<","<<r.peak_story_drift[1]<<","<<r.peak_story_drift[2]
                     <<" residual="<<r.residual_story1<<" roofres="<<r.residual_roof<<"\n";
        };
        print("no_panel",r0.edp);print("panel_pass1",r1.edp);print("panel_pm_pass2",rec);
        std::cerr<<"full_same_rel="<<dfs/hscale<<" final_subdiv="<<same.stats.subdivided_steps<<"\n";

        std::ofstream o(out_path);if(!o)throw std::runtime_error("cannot write refined result");o<<std::setprecision(10)<<std::fixed;
        auto write_edp=[&](const Recorded& r){
            o<<"{\"peak_base_shear_kip\":"<<r.peak_base_shear<<",\"peak_base_shear_over_weight\":"<<r.peak_base_shear/total_weight<<",\"peak_story_drift_percent\":";json_array(o,r.peak_story_drift);
            o<<",\"residual_story1_drift_percent\":"<<r.residual_story1<<",\"residual_roof_drift_percent\":"<<r.residual_roof<<",\"peak_story_shear_kip\":";json_array(o,r.peak_story_shear);
            o<<",\"peak_floor_acceleration_g\":";json_array(o,r.peak_floor_accel);o<<",\"peak_roof_displacement_in\":"<<r.peak_roof_disp<<"}";
        };
        o<<"{\n";
        o<<"  \"benchmark_id\": \"nist_gcr_22_917_50_uc_berkeley_3story_dt1\",\n";
        o<<"  \"validation_status\": \"PRE-VALIDATION_REFINED_MODEL_PROXY_INPUT\",\n";
        o<<"  \"first_mode_period_s\": "<<m2.period<<",\n";
        o<<"  \"peak_base_shear_kip\": "<<rec.peak_base_shear<<",\n";
        o<<"  \"peak_base_shear_over_weight\": "<<rec.peak_base_shear/total_weight<<",\n";
        o<<"  \"peak_story_drift_percent\": ";json_array(o,rec.peak_story_drift);o<<",\n";
        o<<"  \"residual_story1_drift_percent\": "<<rec.residual_story1<<",\n";
        o<<"  \"residual_roof_drift_percent\": "<<rec.residual_roof<<",\n";
        o<<"  \"peak_story_shear_kip\": ";json_array(o,rec.peak_story_shear);o<<",\n";
        o<<"  \"peak_floor_acceleration_g\": ";json_array(o,rec.peak_floor_accel);o<<",\n";
        o<<"  \"peak_roof_displacement_in\": "<<rec.peak_roof_disp<<",\n";
        o<<"  \"story_drift_height_in\": "<<drift_story_h<<",\n";
        o<<"  \"elapsed_s\": "<<rec.elapsed<<", \"newton_iterations\": "<<rec.newton<<", \"factorizations\": "<<rec.factorizations<<",\n";
        o<<"  \"time_s\": ";json_vec(o,rec.time);o<<",\n";
        o<<"  \"roof_displacement_in\": ";json_vec(o,rec.roof);o<<",\n";
        o<<"  \"base_shear_kip_history\": ";json_vec(o,rec.base_shear);o<<",\n";
        o<<"  \"model_progression\": {\n";
        o<<"    \"updated_pdelta_no_panel_gravity_P\": ";write_edp(r0.edp);o<<",\n";
        o<<"    \"panel_zone_gravity_P\": ";write_edp(r1.edp);o<<",\n";
        o<<"    \"panel_zone_PM_second_pass\": ";write_edp(rec);o<<"\n  },\n";
        o<<"  \"column_gravity_preload_kip\": ";json_array12(o,gravity);o<<",\n";
        o<<"  \"column_pass1_max_compression_kip\": ";json_array12(o,r1.edp.max_column_compression);o<<",\n";
        o<<"  \"column_second_pass_design_compression_kip\": ";json_array12(o,p2);o<<",\n";
        o<<"  \"peak_panel_zone_rotation_rad\": ";json_array12(o,rec.peak_panel_zone_rotation);o<<",\n";
        o<<"  \"calibration\": {\"target_perform3d_T1_s\":0.4800000000,\"flexural_stiffness_scale\":"<<m2.scale<<",\"alpha_mass\":"<<m2.alpha_m<<",\"beta_stiffness\":"<<m2.beta_k<<"},\n";
        o<<"  \"solver_invariance\": {\"full_vs_same_status\":\""<<(dfs/hscale<1e-8?"PASS":"OPEN")<<"\",\"max_roof_history_diff_full_vs_same_in\":"<<dfs<<",\"relative_roof_history_difference\":"<<dfs/hscale<<",\"full_seconds\":"<<full.stats.elapsed_seconds<<",\"same_pattern_seconds\":"<<same.stats.elapsed_seconds<<"},\n";
        o<<"  \"model_refinements\": [\n";
        o<<"    \"finite RC panel-zone rotational compliance at all elevated joints using a mechanics-derived research proxy\",\n";
        o<<"    \"updated state-dependent P-delta geometric tangent\",\n";
        o<<"    \"two-pass axial-demand-dependent column P-M section strength iteration using pass-1 maximum compression\"\n";
        o<<"  ],\n";
        o<<"  \"input_motion_provenance\": \"Same deterministic 70 s, 1.52g development proxy used in Run 1. Its duration is approximately consistent with one-third-scale time compression, but it is not the recorded DT1 shake-table waveform and therefore remains unsuitable for strict parity.\",\n";
        o<<"  \"model_provenance\": \"Published geometry/mass/material strengths; stiffness calibrated only to Perform3D T1=0.48s; nonlinear EDPs not fitted. Panel-zone and P-M refinements are independently mechanics-derived research proxies.\",\n";
        o<<"  \"strict_parity_blockers\": [\n";
        o<<"    \"Actual recorded Dynamic Test 1 shake-table acceleration history not acquired\",\n";
        o<<"    \"Exact Perform3D coupled P-M-M hinge surface and hysteretic implementation not available; current two-pass P-M strength iteration is a closer but not identical surrogate\",\n";
        o<<"    \"Exact Perform3D RC panel-zone property definition not available; current joint model is a finite mechanics-derived surrogate\",\n";
        o<<"    \"Exact ASCE 41 component parameter choices and axial-force iteration rules used by the NIST analysis team are not fully disclosed in machine-readable form\",\n";
        o<<"    \"2.5% modal + 0.5% Rayleigh damping remains approximated by equivalent mass/stiffness-proportional terms at mode 1\"\n";
        o<<"  ],\n";
        o<<"  \"perform3d_initial_condition_note\": \"NIST explicitly ran the published Perform3D DT1 benchmark from a virgin numerical state; prior experimental motions were not replayed in that numerical parity model.\"\n";
        o<<"}\n";
        std::cout<<out_path<<"\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}
}
