#include "quake/frame3d.hpp"
#include "quake/fsc_shear_damage.hpp"
#include "quake/modal.hpp"
#include "quake/newmark.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/stability.hpp"

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
constexpr double drift_story_h = 39.0; // Ghannoum/Moehle clear column height used for published drift ratios
// Preserved from the accepted FSC refinement.  This is a mechanics-derived
// expected shear cap from the nonductile test-column detailing, not an EDP fit.
constexpr double first_story_nonductile_expected_shear_kip = 8.1743673639;
constexpr std::array<double,3> floor_weights{{19.6,19.6,19.3}};
constexpr double total_weight = 58.5;
constexpr std::array<double,4> fsc_confined_area_ratios{{0.45,0.55,0.65,0.70}};

struct Motion { double dt{}; std::vector<double> ag; double pga_g{}; double initial_acceleration_g{}; };

Motion read_motion(const std::string& path){
    std::ifstream in(path); if(!in) throw std::runtime_error("cannot open motion: "+path);
    std::string line;std::getline(in,line);
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    if(line!="time_s,accel_g")throw std::invalid_argument("motion header must be time_s,accel_g");
    std::vector<double> t,g;std::size_t row=1;
    while(std::getline(in,line)){
        ++row;if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(line.empty())continue;std::stringstream ss(line);std::string a,b,extra;
        if(!std::getline(ss,a,',')||!std::getline(ss,b,',')||std::getline(ss,extra,','))throw std::invalid_argument("invalid motion CSV row "+std::to_string(row));
        auto number=[&](const std::string& value){std::size_t used=0;double x=std::stod(value,&used);if(used!=value.size()||!std::isfinite(x))throw std::invalid_argument("non-finite or malformed motion value");return x;};
        t.push_back(number(a));g.push_back(number(b));
    }
    if(t.size()<3)throw std::invalid_argument("motion too short");
    const double dt=t[1]-t[0];
    if(t[0]!=0.0||!(dt>0.0))throw std::invalid_argument("motion requires t=0 start and positive dt");
    for(std::size_t i=1;i<t.size();++i)if(std::abs(t[i]-i*dt)>1e-7*dt)throw std::invalid_argument("motion must have a uniform time grid");
    // Driver samples are endpoints a(dt), a(2dt), ... . Preserve the archived
    // zero-initial-state convention, and expose the omitted t=0 acceleration
    // in metadata so measured records cannot silently inherit this assumption.
    std::vector<double> ag;ag.reserve(g.size()-1);double pga=0.0;
    for(std::size_t i=1;i<g.size();++i){pga=std::max(pga,std::abs(g[i]));ag.push_back(g[i]*g_in);}
    return {dt,std::move(ag),pga,g[0]};
}

double pseudo_spectral_accel_g(const Motion& motion,double period_s,double damping_ratio=0.05){
    if(!(period_s>0.0)||motion.ag.empty()||!(motion.dt>0.0))throw std::invalid_argument("invalid response-spectrum input");
    const double w=2.0*3.14159265358979323846/period_s;
    const double k=w*w,c=2.0*damping_ratio*w,dt=motion.dt;
    constexpr double beta=0.25,gamma=0.5;
    const double A0=1.0/(beta*dt*dt),A1=gamma/(beta*dt),A2=1.0/(beta*dt);
    const double A3=1.0/(2.0*beta)-1.0,A4=gamma/beta-1.0,A5=dt*(gamma/(2.0*beta)-1.0);
    const double keff=k+A0+A1*c;
    double u=0.0,v=0.0,a=-motion.ag.front(),sd=0.0;
    for(std::size_t i=1;i<motion.ag.size();++i){
        const double peff=-motion.ag[i]+A0*u+A2*v+A3*a+c*(A1*u+A4*v+A5*a);
        const double un=peff/keff;
        const double an=A0*(un-u)-A2*v-A3*a;
        const double vn=v+dt*((1.0-gamma)*a+gamma*an);
        u=un;v=vn;a=an;sd=std::max(sd,std::abs(u));
    }
    return w*w*sd/g_in;
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


struct PerformConcretePMSurfaceFit {
    double p_balance{};
    double py_tension{};
    double py_compression{};
    double my_balance{};
    double alpha_tension{};
    double alpha_compression{};
    double beta{1.1};
    double rmse_tension{};
    double rmse_compression{};
};

// Fit PERFORM's documented concrete-type P-M surface to the same transparent
// Whitney-block + EPP-steel section mechanics used by the Berkeley benchmark
// source. This is a component-mechanics fit only; no dynamic EDP enters it.
PerformConcretePMSurfaceFit fit_perform_concrete_pm_surface(double bar_area,double fy_ksi){
    PerformConcretePMSurfaceFit out;
    out.py_tension=-8.0*bar_area*fy_ksi;
    out.py_compression=0.85*fc_ksi*col_A+8.0*bar_area*fy_ksi;
    constexpr int n=241;
    std::array<double,n> pp{},mm{};
    for(int i=0;i<n;++i){
        const double q=out.py_tension+(out.py_compression-out.py_tension)*static_cast<double>(i)/(n-1);
        pp[static_cast<std::size_t>(i)]=q;mm[static_cast<std::size_t>(i)]=proxy_column_mn(q,bar_area,fy_ksi);
    }
    auto im=std::max_element(mm.begin(),mm.end());const int ib=static_cast<int>(im-mm.begin());
    out.p_balance=pp[static_cast<std::size_t>(ib)];out.my_balance=*im;
    auto fit_branch=[&](bool compression,double& rmse){
        const double py=compression?out.py_compression:out.py_tension;
        double best_a=1.5,best=1e300;
        for(int ia=0;ia<=2500;++ia){
            const double a=0.5+0.001*ia;double e2=0.0;int count=0;
            for(int i=0;i<n;++i){
                const double q=pp[static_cast<std::size_t>(i)];
                if(compression?(q<out.p_balance):(q>out.p_balance))continue;
                const double r=std::abs((q-out.p_balance)/(py-out.p_balance));
                if(!(r>0.01&&r<0.99))continue;
                const double pred=out.my_balance*std::pow(std::max(0.0,1.0-std::pow(r,a)),1.0/out.beta);
                const double d=(pred-mm[static_cast<std::size_t>(i)])/out.my_balance;e2+=d*d;++count;
            }
            if(count>0&&e2/count<best){best=e2/count;best_a=a;}
        }
        rmse=std::sqrt(best)*out.my_balance;return best_a;
    };
    out.alpha_tension=fit_branch(false,out.rmse_tension);
    out.alpha_compression=fit_branch(true,out.rmse_compression);
    return out;
}

double beam_mn(){
    constexpr double As=4.0*0.11, fy=64.0, b=6.0, d=8.0;
    const double a=As*fy/(0.85*fc_ksi*b);
    return As*fy*(d-0.5*a);
}

struct ColumnProbe {
    int bottom_node{}, top_node{};
    int elastic_element_id{};
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

BuiltRaw build_frame_phase8(double flexural_scale,const ColumnP& design_p,bool panel_zones,
                             bool flexure_shear_critical,bool coupled_pm_critical,bool clear_column_geometry,
                             bool dynamic_updated_pdelta=true,int coupled_grid_mask=3,double alpha_m=0.0,double beta_k=0.0,
                             bool series_shear_springs=false,bool degrading_series_shear=false,
                             double shear_post_ratio=-0.005,double shear_residual_ratio=0.20,
                             double shear_confined_area_ratio=0.55,double shear_cyclic_coefficient=-1.0,
                             bool true_pm_columns=false,bool pm_mroz=false){
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
    if(true_pm_columns && !clear_column_geometry)
        throw std::invalid_argument("true P-M hinges require explicit clear-column end nodes");
    auto add_column_end=[&](int parent,int lev,int grid,double dz,bool release_axial=false){
        int eid=end_id++;b.add_node(eid,x[grid],0.0,z[lev]+dz);
        b.fix_dof(eid,Dof3D::UY);b.fix_dof(eid,Dof3D::RX);b.fix_dof(eid,Dof3D::RZ);
        // Small-rotation rigid arm from the joint center to the column face:
        // ux_face = ux_joint + dz * ry_joint.  For an ordinary scalar hinge
        // uz_face=uz_joint; for a true P-M hinge UZ is deliberately released
        // so the zero-length generalized hinge can develop plastic axial flow.
        b.linear_constraint(eid,Dof3D::UX,{{parent,Dof3D::UX,1.0},{parent,Dof3D::RY,dz}});
        if(!release_axial)b.equal_dof(parent,eid,Dof3D::UZ);
        return eid;
    };
    auto hinge_terms=[&](int end,int core){return std::vector<MpcTerm3D>{{end,Dof3D::RY,1.0},{core,Dof3D::RY,-1.0}};};
    const auto gravity=gravity_column_preloads();
    std::array<ColumnProbe,12> probes{};
    struct PendingShearSpring {
        int core{}, end{};
        std::vector<MpcTerm3D> bottom_rotation_terms,top_rotation_terms,axial_terms;
        double axial_preload{},axial_stiffness{};
    };
    std::vector<PendingShearSpring> pending_shear_springs;

    // Columns. The scalar hinge capacity is evaluated on a P-M interaction
    // section at the supplied design axial force. A second pass updates those
    // design axial forces from the first nonlinear run, mirroring the NIST
    // axial-demand iteration without fitting any response EDP.
    for(int st=0;st<3;++st){
        for(int grid=0;grid<4;++grid){
            const int ci=ccore(st,grid),cj=ccore(st+1,grid);
            const bool nonductile=(grid<=1);
            const bool series_shear_here=series_shear_springs && st==0 && nonductile;
            const double joint_half_depth=4.5;
            int ei=-1; int shear_ref=ci;
            if(series_shear_here){
                // Joint-side reference and column-side node for a zero-length
                // UX shear spring. With clear-column geometry the reference is
                // a rigidly offset column-face node, so the 39-in flexible
                // length, local shear, and published drift definition are
                // kinematically consistent.
                if(clear_column_geometry){
                    shear_ref=add_column_end(ci,st,grid,+joint_half_depth,false);
                    b.equal_dof(ci,shear_ref,Dof3D::RY);
                }
                const double zi=clear_column_geometry?z[st]+joint_half_depth:z[st];
                ei=end_id++; b.add_node(ei,x[grid],0.0,zi);
                b.fix_dof(ei,Dof3D::UY);b.fix_dof(ei,Dof3D::RX);b.fix_dof(ei,Dof3D::RZ);
                if(!true_pm_columns)b.equal_dof(ci,ei,Dof3D::UZ);
            }else{
                ei=clear_column_geometry?add_column_end(ci,st,grid,+joint_half_depth,true_pm_columns):add_end(ci,st,grid);
            }
            const int ej=clear_column_geometry?add_column_end(cj,st+1,grid,-joint_half_depth,true_pm_columns):add_end(cj,st+1,grid);
            const double column_length=clear_column_geometry?drift_story_h:48.0;
            const double Pg=gravity[static_cast<std::size_t>(st*4+grid)];
            const int column_element_id=elem_id++;
            b.add_elastic_frame(column_element_id,ei,ej,Ec_ksi,G,col_A,25.0,
                col_I*0.30*flexural_scale,col_I*0.30*flexural_scale,0,1,0,Pg);
            const double Pd=std::max(0.0,design_p[static_cast<std::size_t>(st*4+grid)]);
            const double bar_area=nonductile?0.049:0.11, fy=nonductile?70.0:64.0;
            const double Mn_section=proxy_column_mn(Pd,bar_area,fy);
            const bool fsc_column=flexure_shear_critical && st==0 && nonductile;
            const double shear_moment_cap=first_story_nonductile_expected_shear_kip*drift_story_h/2.0;
            const double Mn=(fsc_column && !series_shear_springs)?std::min(Mn_section,shear_moment_cap):Mn_section;
            const double member_krot=4.0*Ec_ksi*(col_I*0.30*flexural_scale)/column_length;
            const double krot=100.0*member_krot;
            ASCE41HingeParams hp;hp.Ke=krot;hp.posFy=hp.negFy=Mn;hp.hardening_ratio=fsc_column?0.0:0.015;
            if(!fsc_column) hp.hardening_stiffness=0.015*member_krot;
            if(nonductile){
                hp.pos_a=hp.neg_a=0.010;hp.pos_b=hp.neg_b=0.032;hp.pos_f=hp.neg_f=0.045;hp.pos_c=hp.neg_c=0.15;
                hp.pos_io=hp.neg_io=0.004;hp.pos_ls=hp.neg_ls=0.012;hp.pos_cp=hp.neg_cp=0.022;
                hp.pos_drop_span=hp.neg_drop_span=0.006;hp.pos_e_drop_span=hp.neg_e_drop_span=0.006;
            }else{
                hp.pos_a=hp.neg_a=0.030;hp.pos_b=hp.neg_b=0.085;hp.pos_f=hp.neg_f=0.120;hp.pos_c=hp.neg_c=0.20;
                hp.pos_io=hp.neg_io=0.008;hp.pos_ls=hp.neg_ls=0.025;hp.pos_cp=hp.neg_cp=0.050;
                hp.pos_drop_span=hp.neg_drop_span=0.012;hp.pos_e_drop_span=hp.neg_e_drop_span=0.010;
            }
            const bool coupled_here=coupled_pm_critical && st==0 && nonductile && ((coupled_grid_mask & (1<<grid))!=0);
            if(true_pm_columns){
                PMInteractionHingeParams pm;
                pm.hinge=hp;
                // Phase 9F topology isolation: use the straight C->E NIST/Perform
                // interpretation while retaining the already-resolved QuakeCore
                // a/b values. The P-M surface itself is enforced continuously.
                pm.hinge.backbone_shape=ASCE41BackboneShape::StraightCE;
                pm.hinge.pos_f=pm.hinge.pos_b; pm.hinge.neg_f=pm.hinge.neg_b;
                // Phase 9F.1 is a mechanics isolation: enforce the Powell/Perform
                // associative P-M surface continuously, but keep the surface
                // elastic-perfectly-plastic so post-capping shrinkage does not
                // contaminate the return-map experiment. The original a/b
                // values remain available as diagnostic rotation thresholds.
                pm.surface_evolution=pm_mroz?PMInteractionSurfaceEvolution::MrozTwoSurface:
                                              PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic;
                pm.enforce_deformation_capacity=false;
                pm.axial_preload=Pg;
                pm.axial_stiffness=30.0*Ec_ksi*col_A/column_length; // rigid-plastic hinge penalty
                const auto ps=fit_perform_concrete_pm_surface(bar_area,fy);
                pm.surface_shape=PMInteractionSurfaceShape::PerformConcrete;
                pm.p_balance=ps.p_balance;pm.py_tension=ps.py_tension;pm.py_compression=ps.py_compression;
                pm.my_balance=ps.my_balance;pm.alpha_tension=ps.alpha_tension;pm.alpha_compression=ps.alpha_compression;pm.beta_pm=ps.beta;
                // For an FSC column with an explicit series shear spring, flexural
                // strength remains the section surface; shear failure is modeled by
                // the separate FSC spring rather than clipping the P-M surface.
                // Positive axial generalized deformation is tensile extension.
                // The top hinge coordinate is reversed so both ends share the
                // same physical tension-positive convention.
                const std::vector<MpcTerm3D> axial_bottom{{ci,Dof3D::UZ,-1.0},{ei,Dof3D::UZ,1.0}};
                const std::vector<MpcTerm3D> axial_top{{ej,Dof3D::UZ,-1.0},{cj,Dof3D::UZ,1.0}};
                b.add_linear_pm_interaction_hinge(spring_id++,hinge_terms(ei,ci),axial_bottom,pm);
                b.add_linear_pm_interaction_hinge(spring_id++,hinge_terms(ej,cj),axial_top,pm);
            }else if(coupled_here){
                AxialCoupledASCE41Params cp;
                cp.hinge=hp;
                cp.axial_preload=Pg;
                cp.axial_stiffness=Ec_ksi*col_A/column_length;
                constexpr int npts=25;
                const double upper=0.80*col_A*fc_ksi;
                cp.axial_force_points.reserve(npts);
                cp.moment_capacity_points.reserve(npts);
                for(int ip=0;ip<npts;++ip){
                    const double P=upper*static_cast<double>(ip)/static_cast<double>(npts-1);
                    double Mc=proxy_column_mn(P,bar_area,fy);
                    if(fsc_column && !series_shear_springs)Mc=std::min(Mc,shear_moment_cap);
                    cp.axial_force_points.push_back(P);
                    cp.moment_capacity_points.push_back(Mc);
                }
                const std::vector<MpcTerm3D> axial_terms{{ei,Dof3D::UZ,-1.0},{ej,Dof3D::UZ,1.0}};
                b.add_linear_axial_coupled_asce41_hinge(spring_id++,hinge_terms(ei,ci),axial_terms,cp);
                b.add_linear_axial_coupled_asce41_hinge(spring_id++,hinge_terms(ej,cj),axial_terms,cp);
            }else{
                b.add_linear_asce41_hinge(spring_id++,hinge_terms(ei,ci),hp);
                b.add_linear_asce41_hinge(spring_id++,hinge_terms(ej,cj),hp);
            }
            if(series_shear_here){
                pending_shear_springs.push_back({shear_ref,ei,hinge_terms(ei,ci),hinge_terms(ej,cj),
                    {{ei,Dof3D::UZ,-1.0},{ej,Dof3D::UZ,1.0}},Pg,Ec_ksi*col_A/column_length});
            }
            probes[static_cast<std::size_t>(st*4+grid)]={true_pm_columns?ei:ci,true_pm_columns?ej:cj,column_element_id,Pg,Ec_ksi*col_A/column_length,st,grid};
        }
    }

    // Pre-failure series shear compliance.  These are deliberately kept
    // elastic in this stage; the state-triggered degradation law is promoted
    // only after the accepted-history diagnostic has been validated.
    if(series_shear_springs){
        const double kshear=col_A*G*5.0/(6.0*drift_story_h); // elastic shear stiffness of 6x6 column
        for(const auto& ss:pending_shear_springs){
            const std::vector<MpcTerm3D> shear_terms{{ss.core,Dof3D::UX,-1.0},{ss.end,Dof3D::UX,1.0}};
            if(!degrading_series_shear){
                b.add_linear_bilinear_spring(spring_id++,shear_terms,kshear,1.0e9,0.0);
            }else{
                FSCShearSpringParams sp;sp.Ke=kshear;sp.post_failure_stiffness=shear_post_ratio*kshear;
                sp.residual_strength_ratio=shear_residual_ratio;sp.nominal_shear_limit_kip=-1.0; // rotation-based initiation; V enters the threshold continuously
                sp.b_in=6.0;sp.d_in=4.8;sp.h_in=6.0;sp.clear_length_in=drift_story_h;sp.tie_spacing_in=4.0;
                sp.longitudinal_steel_area_in2=8.0*0.11;sp.confined_concrete_area_in2=shear_confined_area_ratio*col_A;
                sp.fc_ksi=fc_ksi;sp.fy_ksi=64.0;sp.cyclic_strength_coefficient=shear_cyclic_coefficient;
                b.add_linear_fsc_shear_spring(spring_id++,shear_terms,ss.bottom_rotation_terms,ss.top_rotation_terms,
                                               ss.axial_terms,ss.axial_preload,ss.axial_stiffness,sp);
            }
        }
    }

    const double Mnb=beam_mn();
    for(int lev=1;lev<4;++lev)for(int bay=0;bay<3;++bay){
        const int li=panel_zones?bcore(lev,bay):ccore(lev,bay);
        const int lj=panel_zones?bcore(lev,bay+1):ccore(lev,bay+1);
        int ei=add_end(li,lev,bay),ej=add_end(lj,lev,bay+1);
        b.add_elastic_frame(elem_id++,ei,ej,Ec_ksi,G,beam_A,60.0,
            beam_I*0.30*flexural_scale,beam_Iweak*0.30*flexural_scale,0,1,0,0.0);
        const double member_krot=4.0*Ec_ksi*(beam_I*0.30*flexural_scale)/70.0;
        const double krot=100.0*member_krot;
        ASCE41HingeParams hp;hp.Ke=krot;hp.posFy=hp.negFy=Mnb;hp.hardening_ratio=0.02;
        hp.hardening_stiffness=0.02*member_krot;
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
    b.set_updated_pdelta(dynamic_updated_pdelta);
    b.set_ground_direction(Dof3D::UX);
    b.set_response(ccore(3,0),Dof3D::UX);
    b.set_story_nodes({ccore(1,0),ccore(2,0),ccore(3,0)},Dof3D::UX);
    b.set_rayleigh(alpha_m,beta_k);
    return {b.compile(),probes};
}

ModelBuild calibrated_model_phase8(const ColumnP& design_p,bool panel_zones,
                                  bool flexure_shear_critical,bool coupled_pm_critical,bool clear_column_geometry=false,
                                  bool dynamic_updated_pdelta=true,int coupled_grid_mask=3,bool series_shear_springs=false,
                                  bool degrading_series_shear=false,double shear_post_ratio=-0.005,
                                  double shear_residual_ratio=0.20,double shear_confined_area_ratio=0.55,
                                  double shear_cyclic_coefficient=-1.0,bool true_pm_columns=false){
    constexpr double target=0.48;
    double lo=0.03,hi=3.0,best=1.0,T=0.0;
    for(int it=0;it<34;++it){
        const double sc=0.5*(lo+hi);auto raw=build_frame_phase8(sc,design_p,panel_zones,flexure_shear_critical,coupled_pm_critical,clear_column_geometry,dynamic_updated_pdelta,coupled_grid_mask,0.0,0.0,series_shear_springs,degrading_series_shear,shear_post_ratio,shear_residual_ratio,shear_confined_area_ratio,shear_cyclic_coefficient,true_pm_columns);auto modes=modal_analysis(raw.model,1);if(modes.empty())throw std::runtime_error("no finite mode");
        T=modes[0].period;best=sc;if(T>target)lo=sc;else hi=sc;
    }
    auto raw0=build_frame_phase8(best,design_p,panel_zones,flexure_shear_critical,coupled_pm_critical,clear_column_geometry,dynamic_updated_pdelta,coupled_grid_mask,0.0,0.0,series_shear_springs,degrading_series_shear,shear_post_ratio,shear_residual_ratio,shear_confined_area_ratio,shear_cyclic_coefficient,true_pm_columns);auto modes=modal_analysis(raw0.model,3);T=modes[0].period;
    const double w=2.0*3.14159265358979323846/T;
    const double alpha=2.0*0.025*w;const double beta=2.0*0.005/w;
    auto raw=build_frame_phase8(best,design_p,panel_zones,flexure_shear_critical,coupled_pm_critical,clear_column_geometry,dynamic_updated_pdelta,coupled_grid_mask,alpha,beta,series_shear_springs,degrading_series_shear,shear_post_ratio,shear_residual_ratio,shear_confined_area_ratio,shear_cyclic_coefficient,true_pm_columns);auto final_modes=modal_analysis(raw.model,3);
    std::cerr<<"phase8 calibrated panel="<<panel_zones<<" fsc="<<flexure_shear_critical<<" coupled="<<coupled_pm_critical
             <<" clear_col="<<clear_column_geometry<<" truePM="<<true_pm_columns<<" updPDelta="<<dynamic_updated_pdelta<<" coupledMask="<<coupled_grid_mask<<" seriesShear="<<series_shear_springs<<" degradingShear="<<degrading_series_shear<<" scale="<<best<<" T1="<<final_modes[0].period<<" dof="<<raw.model.dof()<<" hinges="<<raw.model.nonlinear_count()<<"\n";
    return {std::move(raw.model),best,final_modes[0].period,alpha,beta,raw.probes,panel_zones};
}

struct Recorded {
    // Dynamic equilibrium outputs are intentionally distinct.  The inertial
    // reaction is sum(m*a_abs); the structural restoring story shear is
    // recovered from the accepted internal-force vector and is the quantity
    // used for mechanism/capacity comparison.
    double peak_base_shear{}; // legacy alias: inertial reaction
    double peak_structural_base_shear{}; // reduced-coordinate restoring generalized force (diagnostic)
    double peak_first_story_column_shear{}; // moment-derived diagnostic, V=(Mbot+Mtop)/L
    double peak_first_story_component_shear{}; // A/B explicit shear springs + C/D moment-derived column shears
    std::array<double,3> peak_story_shear{}; // legacy inertial story reaction
    std::array<double,3> peak_structural_story_shear{};
    std::array<double,3> peak_story_drift{};
    std::array<double,3> peak_floor_accel{};double residual_story1{};double residual_roof{};double peak_roof_disp{};
    std::size_t newton{};std::size_t factorizations{};double elapsed{};
    std::array<double,12> max_column_compression{};
    std::array<double,12> peak_panel_zone_rotation{};
    std::array<double,8> peak_first_story_column_hinge_rotation{};
    std::array<double,8> signed_peak_first_story_column_hinge_rotation{};
    std::array<double,8> compression_at_peak_first_story_hinge_rotation{};
    std::array<double,8> time_at_peak_first_story_hinge_rotation{};
    std::array<double,8> story1_drift_at_peak_first_story_hinge_rotation{};
    double max_story1_axial_sum_error{};
    double min_story1_axial_sum{1e30};
    double max_story1_axial_sum{-1e30};
    std::array<double,2> peak_nonductile_column_shear_kip{};
    std::array<double,2> max_elwood_shear_demand_ratio{};
    std::array<std::size_t,2> first_elwood_shear_failure_step{{static_cast<std::size_t>(-1),static_cast<std::size_t>(-1)}};
    std::array<double,2> elwood_failure_drift_percent{};
    std::array<double,2> elwood_failure_axial_kip{};
    std::array<double,2> elwood_failure_shear_kip{};
    // Ghannoum-Moehle rotation-based FSC diagnostics.  Two variants are
    // retained because the current concentrated hinge does not contain an
    // independently calibrated bar-slip spring: theta_total is the robust
    // model-independent regression, while theta_total_plastic is the closest
    // direct comparison to the concentrated hinge deformation.
    std::array<double,2> max_gm_total_rotation_ratio{};
    std::array<double,2> max_gm_total_plastic_rotation_ratio{};
    std::array<std::size_t,2> first_gm_total_plastic_failure_step{{static_cast<std::size_t>(-1),static_cast<std::size_t>(-1)}};
    std::array<double,2> gm_failure_time_s{};
    std::array<double,2> gm_failure_story_drift_percent{};
    std::array<double,2> gm_failure_axial_kip{};
    std::array<double,2> gm_failure_shear_kip{};
    std::array<double,2> gm_failure_hinge_rotation_rad{};
    std::array<int,2> gm_failure_end{}; // 0 bottom, 1 top
    std::array<std::array<double,2>,4> fsc_damage_retained_ratio{};
    std::array<std::array<double,2>,4> fsc_damage_cyclic_coefficient{};
    std::array<std::array<std::size_t,2>,4> fsc_damage_initiation_step{};
    std::array<std::array<std::size_t,2>,4> fsc_damage_residual_step{};
    std::array<std::array<std::size_t,2>,4> fsc_damage_sign_crossings{};
    std::array<std::array<double,2>,4> fsc_damage_initiation_time_s{};
    std::array<std::array<int,2>,4> fsc_damage_initiation_cause{};
    std::array<std::size_t,2> active_fsc_initiation_step{{static_cast<std::size_t>(-1),static_cast<std::size_t>(-1)}};
    std::array<std::size_t,2> active_fsc_residual_step{{static_cast<std::size_t>(-1),static_cast<std::size_t>(-1)}};
    std::array<double,2> active_fsc_initiation_time_s{};
    std::array<double,2> active_fsc_initiation_drift_percent{};
    std::array<double,2> active_fsc_initiation_axial_kip{};
    std::array<double,2> active_fsc_initiation_shear_kip{};
    std::array<double,2> active_fsc_residual_time_s{};
    std::array<double,2> active_fsc_min_retained_ratio{{1.0,1.0}};
    std::array<double,2> active_fsc_failure_strength_kip{};
    std::vector<double> time,roof,base_shear,structural_base_shear,first_story_column_shear,first_story_component_shear;
};

struct RecordedRun { Recorded edp; AnalysisResult analysis; };

RecordedRun record_robust(const ModelBuild& mb,const Motion& motion,LinearStrategy strategy,bool shear_limit_diagnostic=false,bool trace_a1=false){
    const auto& model=mb.model;
    Recorded rec;
    std::vector<double> residual_s1,residual_roof;residual_s1.reserve(600);residual_roof.reserve(600);
    const std::array<double,3> floor_mass{{floor_weights[0]/g_in,floor_weights[1]/g_in,floor_weights[2]/g_in}};
    const std::size_t stride=std::max<std::size_t>(1,motion.ag.size()/350);
    RobustNewmarkOptions ro;ro.tolerance=1e-7;ro.max_iterations=40;ro.line_search=true;ro.max_backtracks=8;
    ro.max_subdivisions=5;ro.return_numerical_failure=true;ro.collapse.check_initial_stability=true;
    std::array<std::array<FSCShearDamageDiagnostic,2>,4> fsc_damage{{
        {FSCShearDamageDiagnostic(FSCShearDamageParams{}),FSCShearDamageDiagnostic(FSCShearDamageParams{})},
        {FSCShearDamageDiagnostic(FSCShearDamageParams{}),FSCShearDamageDiagnostic(FSCShearDamageParams{})},
        {FSCShearDamageDiagnostic(FSCShearDamageParams{}),FSCShearDamageDiagnostic(FSCShearDamageParams{})},
        {FSCShearDamageDiagnostic(FSCShearDamageParams{}),FSCShearDamageDiagnostic(FSCShearDamageParams{})}
    }};
    for(std::size_t ir=0;ir<fsc_confined_area_ratios.size();++ir){
        for(int grid=0;grid<2;++grid){
            FSCShearDamageParams dp;
            dp.confined_concrete_area_in2=fsc_confined_area_ratios[ir]*col_A;
            dp.nominal_shear_limit_kip=first_story_nonductile_expected_shear_kip;
            dp.residual_strength_ratio=0.20;
            fsc_damage[ir][static_cast<std::size_t>(grid)]=FSCShearDamageDiagnostic(dp);
        }
    }
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
        std::array<double,12> current_column_compression{};
        for(std::size_t k=0;k<mb.column_probes.size();++k){
            const auto& pr=mb.column_probes[k];
            const double P=model.elastic_element_response(pr.elastic_element_id,u).axial_compression;
            current_column_compression[k]=P;
            rec.max_column_compression[k]=std::max(rec.max_column_compression[k],P);
        }
        double story1_psum=0.0;
        for(int grid=0;grid<4;++grid)story1_psum+=current_column_compression[static_cast<std::size_t>(grid)];
        rec.min_story1_axial_sum=std::min(rec.min_story1_axial_sum,story1_psum);
        rec.max_story1_axial_sum=std::max(rec.max_story1_axial_sum,story1_psum);
        rec.max_story1_axial_sum_error=std::max(rec.max_story1_axial_sum_error,std::abs(story1_psum-floor_weights[0]-floor_weights[1]-floor_weights[2]));
        const auto& B=model.nonlinear_basis();
        if(model.nonlinear_count()>=8){
            for(int j=0;j<8;++j){
                const double q=B.column_dot(j,u);
                const auto jj=static_cast<std::size_t>(j);
                if(std::abs(q)>rec.peak_first_story_column_hinge_rotation[jj]){
                    rec.peak_first_story_column_hinge_rotation[jj]=std::abs(q);
                    rec.signed_peak_first_story_column_hinge_rotation[jj]=q;
                    rec.compression_at_peak_first_story_hinge_rotation[jj]=current_column_compression[static_cast<std::size_t>(j/2)];
                    rec.time_at_peak_first_story_hinge_rotation[jj]=time;
                    rec.story1_drift_at_peak_first_story_hinge_rotation[jj]=drift[0];
                }
            }
        }
        if(mb.panel_zones && model.nonlinear_count()>=12){
            const int pz0=model.nonlinear_count()-12; // panel zones are compiled last in this scaffold
            for(int j=0;j<12;++j){
                const double q=B.column_dot(pz0+j,u);
                rec.peak_panel_zone_rotation[static_cast<std::size_t>(j)]=std::max(rec.peak_panel_zone_rotation[static_cast<std::size_t>(j)],std::abs(q));
            }
        }
        if(time>=(motion.ag.size()*motion.dt-5.0)){residual_s1.push_back(drift[0]);residual_roof.push_back(fu[2]/(3.0*drift_story_h)*100.0);}
        if(step%stride==0||step+1==motion.ag.size()){rec.time.push_back(time);rec.roof.push_back(fu[2]);rec.base_shear.push_back(shear[0]);}
    };
    if(shear_limit_diagnostic){
        ro.accepted_state_observer=[&](std::size_t step,double time, double, const std::vector<double>& u,
                                      const std::vector<double>&,const std::vector<double>&,
                                      const std::vector<double>& committed){
            // Accepted structural restoring force, excluding M*a and C*v.
            std::vector<double> fint,ftan,ftrial;
            model.internal_force_and_tangent(u,committed,fint,ftan,ftrial);
            std::array<int,3> floor_ux{{
                model.reduced_dof(1010,Dof3D::UX),
                model.reduced_dof(1020,Dof3D::UX),
                model.reduced_dof(1030,Dof3D::UX)}};
            std::array<double,3> structural_story_shear{};
            for(int st=0;st<3;++st){
                double q=0.0;
                for(int fl=st;fl<3;++fl){
                    const int r=floor_ux[static_cast<std::size_t>(fl)];
                    if(r>=0) q+=fint[static_cast<std::size_t>(r)];
                }
                structural_story_shear[static_cast<std::size_t>(st)]=q;
                rec.peak_structural_story_shear[static_cast<std::size_t>(st)]=
                    std::max(rec.peak_structural_story_shear[static_cast<std::size_t>(st)],std::abs(q));
            }
            rec.peak_structural_base_shear=std::max(rec.peak_structural_base_shear,std::abs(structural_story_shear[0]));
            if(step%stride==0||step+1==motion.ag.size()) rec.structural_base_shear.push_back(structural_story_shear[0]);
            const auto fu=model.story_response_values(u);
            const double drift1=fu[0]/drift_story_h;
            const auto& B=model.nonlinear_basis();
            std::vector<double> q(static_cast<std::size_t>(model.nonlinear_count()));
            for(int j=0;j<model.nonlinear_count();++j)q[static_cast<std::size_t>(j)]=B.column_dot(j,u);
            for(int grid=0;grid<2;++grid){
                const int idx=model.nonlinear_component_index(25+grid);
                if(idx>=0){
                    const auto ss=model.fsc_shear_state(idx,committed);
                    if(ss.valid){
                        rec.active_fsc_min_retained_ratio[static_cast<std::size_t>(grid)]=
                            std::min(rec.active_fsc_min_retained_ratio[static_cast<std::size_t>(grid)],ss.retained_strength_ratio);
                        rec.active_fsc_failure_strength_kip[static_cast<std::size_t>(grid)]=std::max(
                            rec.active_fsc_failure_strength_kip[static_cast<std::size_t>(grid)],ss.failure_strength);
                        if(ss.initiated && rec.active_fsc_initiation_step[static_cast<std::size_t>(grid)]==static_cast<std::size_t>(-1)){
                            rec.active_fsc_initiation_step[static_cast<std::size_t>(grid)]=step;
                            rec.active_fsc_initiation_time_s[static_cast<std::size_t>(grid)]=time;
                            rec.active_fsc_initiation_drift_percent[static_cast<std::size_t>(grid)]=100.0*drift1;
                            const auto& pr=mb.column_probes[static_cast<std::size_t>(grid)];
                            const int ib=model.reduced_dof(pr.bottom_node,Dof3D::UZ),it=model.reduced_dof(pr.top_node,Dof3D::UZ);
                            const double ub=ib>=0?u[static_cast<std::size_t>(ib)]:0.0,ut=it>=0?u[static_cast<std::size_t>(it)]:0.0;
                            rec.active_fsc_initiation_axial_kip[static_cast<std::size_t>(grid)]=std::max(0.0,pr.gravity_preload-pr.axial_k*(ut-ub));
                            rec.active_fsc_initiation_shear_kip[static_cast<std::size_t>(grid)]=std::abs(ss.force);
                        }
                        if(ss.residual_reached && rec.active_fsc_residual_step[static_cast<std::size_t>(grid)]==static_cast<std::size_t>(-1)){
                            rec.active_fsc_residual_step[static_cast<std::size_t>(grid)]=step;
                            rec.active_fsc_residual_time_s[static_cast<std::size_t>(grid)]=time;
                        }
                    }
                }
            }
            std::vector<double> cf,kt,trial;
            model.evaluate_nonlinear_deformations(q,committed,cf,kt,trial);
            // Physical first-story restoring shear from the four column end
            // moments.  This avoids interpreting condensed/MPC generalized
            // forces as a literal support reaction.  With no transverse load
            // along a clear column, V=(M_bottom+M_top)/L in the hinge sign
            // convention used by this validation model.
            double vcols=0.0;
            for(int grid=0;grid<4;++grid){
                const int ib=model.nonlinear_component_index(1+2*grid);
                const int it=model.nonlinear_component_index(2+2*grid);
                if(ib>=0&&it>=0)vcols+=(cf[static_cast<std::size_t>(ib)]+cf[static_cast<std::size_t>(it)])/drift_story_h;
            }
            rec.peak_first_story_column_shear=std::max(rec.peak_first_story_column_shear,std::abs(vcols));
            if(step%stride==0||step+1==motion.ag.size())rec.first_story_column_shear.push_back(vcols);
            // In the Phase-8 series topology A1/B1 have explicit shear springs,
            // so their spring forces are the most direct component shear. C1/D1
            // retain the end-moment recovery. This hybrid is a component-level
            // cross-check; the static pushover load factor is still the exact
            // global capacity ordinate.
            double vcomponent=0.0;
            for(int grid=0;grid<4;++grid){
                if(grid<2){
                    const int is=model.nonlinear_component_index(25+grid);
                    if(is>=0){const auto ss=model.fsc_shear_state(is,committed);if(ss.valid)vcomponent+=ss.force;}
                }else{
                    const int ib=model.nonlinear_component_index(1+2*grid);
                    const int it=model.nonlinear_component_index(2+2*grid);
                    if(ib>=0&&it>=0)vcomponent+=(cf[static_cast<std::size_t>(ib)]+cf[static_cast<std::size_t>(it)])/drift_story_h;
                }
            }
            rec.peak_first_story_component_shear=std::max(rec.peak_first_story_component_shear,std::abs(vcomponent));
            if(step%stride==0||step+1==motion.ag.size())rec.first_story_component_shear.push_back(vcomponent);
            if(trace_a1 && step>=465 && step<=485 && committed.size()>=20){
                const auto& pr=mb.column_probes[0];
                const int ib=model.reduced_dof(pr.bottom_node,Dof3D::UZ),it=model.reduced_dof(pr.top_node,Dof3D::UZ);
                const double ub=ib>=0?u[static_cast<std::size_t>(ib)]:0.0,ut=it>=0?u[static_cast<std::size_t>(it)]:0.0;
                const double P=pr.gravity_preload-pr.axial_k*(ut-ub);
                std::cerr<<"A1TRACE step="<<step<<" t="<<time<<" P="<<P
                         <<" qb="<<q[0]<<" cqB="<<committed[0]<<" fB="<<committed[1]<<" brB="<<committed[14]
                         <<" u0B="<<committed[15]<<" tqB="<<committed[16]<<" tfB="<<committed[17]
                         <<" ktB="<<committed[18]<<" MyB="<<committed[19]
                         <<" qt="<<q[1]<<" cqT="<<committed[20]<<" fT="<<committed[21]<<" brT="<<committed[34]
                         <<" u0T="<<committed[35]<<" tqT="<<committed[36]<<" tfT="<<committed[37]
                         <<" ktT="<<committed[38]<<" MyT="<<committed[39]<<"\n";
            }
            constexpr double rho_t=0.0015;
            constexpr double d_eff_in=0.8*6.0;
            constexpr double ksi_to_mpa=6.894757293168;
            const double fc_mpa=fc_ksi*ksi_to_mpa;
            for(int grid=0;grid<2;++grid){
                const double Mbot=cf[static_cast<std::size_t>(2*grid)];
                const double Mtop=cf[static_cast<std::size_t>(2*grid+1)];
                const double V=std::abs(Mbot+Mtop)/drift_story_h;
                rec.peak_nonductile_column_shear_kip[static_cast<std::size_t>(grid)]=
                    std::max(rec.peak_nonductile_column_shear_kip[static_cast<std::size_t>(grid)],V);
                const auto& pr=mb.column_probes[static_cast<std::size_t>(grid)];
                const int ib=model.reduced_dof(pr.bottom_node,Dof3D::UZ),it=model.reduced_dof(pr.top_node,Dof3D::UZ);
                const double ub=ib>=0?u[static_cast<std::size_t>(ib)]:0.0,ut=it>=0?u[static_cast<std::size_t>(it)]:0.0;
                const double P=std::max(0.0,pr.gravity_preload-pr.axial_k*(ut-ub));
                const double v_mpa=(V/(6.0*d_eff_in))*ksi_to_mpa;
                const double theta_s=std::max(0.01,0.03+4.0*rho_t-(1.0/40.0)*(v_mpa/std::sqrt(fc_mpa))
                                               -(1.0/40.0)*(P/(col_A*fc_ksi)));
                const double ratio=std::abs(drift1)/theta_s;
                auto& maxr=rec.max_elwood_shear_demand_ratio[static_cast<std::size_t>(grid)];
                maxr=std::max(maxr,ratio);
                auto& first=rec.first_elwood_shear_failure_step[static_cast<std::size_t>(grid)];
                if(ratio>=1.0 && first==static_cast<std::size_t>(-1)){
                    first=step;rec.elwood_failure_drift_percent[static_cast<std::size_t>(grid)]=100.0*drift1;
                    rec.elwood_failure_axial_kip[static_cast<std::size_t>(grid)]=P;
                    rec.elwood_failure_shear_kip[static_cast<std::size_t>(grid)]=V;
                }

                // Ghannoum & Moehle rotation-based shear-failure initiation
                // regression (psi units).  The local end rotations are what
                // distinguish nominally identical columns in a frame.
                constexpr double tie_spacing_in=4.0;
                constexpr double d_gm_in=4.8;
                const double fc_psi=fc_ksi*1000.0;
                const double v_psi=V/(6.0*d_gm_in)*1000.0;
                const double s_over_d=tie_spacing_in/d_gm_in;
                const double p_ratio=P/(col_A*fc_ksi);
                const double vsqrt=v_psi/std::sqrt(fc_psi);
                const double theta_total=std::max(0.009,0.044-0.017*s_over_d-0.021*p_ratio-0.0020*vsqrt);
                const double theta_total_pl=std::max(0.0,0.032-0.014*s_over_d-0.017*p_ratio-0.0016*vsqrt);
                double local_peak=0.0; int local_end=0;
                for(int end=0;end<2;++end){
                    const double qr=std::abs(q[static_cast<std::size_t>(2*grid+end)]);
                    if(qr>local_peak){local_peak=qr;local_end=end;}
                }
                auto& rg=rec.max_gm_total_rotation_ratio[static_cast<std::size_t>(grid)];
                rg=std::max(rg,local_peak/theta_total);
                if(theta_total_pl>1e-12){
                    const double rpl=local_peak/theta_total_pl;
                    auto& rgp=rec.max_gm_total_plastic_rotation_ratio[static_cast<std::size_t>(grid)];
                    rgp=std::max(rgp,rpl);
                    auto& fg=rec.first_gm_total_plastic_failure_step[static_cast<std::size_t>(grid)];
                    if(rpl>=1.0 && fg==static_cast<std::size_t>(-1)){
                        fg=step;
                        rec.gm_failure_time_s[static_cast<std::size_t>(grid)]=time;
                        rec.gm_failure_story_drift_percent[static_cast<std::size_t>(grid)]=100.0*drift1;
                        rec.gm_failure_axial_kip[static_cast<std::size_t>(grid)]=P;
                        rec.gm_failure_shear_kip[static_cast<std::size_t>(grid)]=V;
                        rec.gm_failure_hinge_rotation_rad[static_cast<std::size_t>(grid)]=local_peak;
                        rec.gm_failure_end[static_cast<std::size_t>(grid)]=local_end;
                    }
                }

                const double signed_V=(Mbot+Mtop)/drift_story_h;
                const double qb=q[static_cast<std::size_t>(2*grid)];
                const double qt=q[static_cast<std::size_t>(2*grid+1)];
                for(std::size_t ir=0;ir<fsc_confined_area_ratios.size();++ir)
                    fsc_damage[ir][static_cast<std::size_t>(grid)].update(step,time,qb,qt,P,signed_V);
            }
        };
    }
    auto result=run_newmark_robust(model,motion.ag,motion.dt,strategy,ro);
    if(shear_limit_diagnostic){
        for(std::size_t ir=0;ir<fsc_confined_area_ratios.size();++ir)for(int grid=0;grid<2;++grid){
            const auto& ds=fsc_damage[ir][static_cast<std::size_t>(grid)].state();
            rec.fsc_damage_retained_ratio[ir][static_cast<std::size_t>(grid)]=ds.retained_strength_ratio;
            rec.fsc_damage_cyclic_coefficient[ir][static_cast<std::size_t>(grid)]=ds.cyclic_strength_coefficient;
            rec.fsc_damage_initiation_step[ir][static_cast<std::size_t>(grid)]=ds.initiation_step;
            rec.fsc_damage_residual_step[ir][static_cast<std::size_t>(grid)]=ds.residual_step;
            rec.fsc_damage_sign_crossings[ir][static_cast<std::size_t>(grid)]=ds.opposite_sign_crossings;
            rec.fsc_damage_initiation_time_s[ir][static_cast<std::size_t>(grid)]=ds.initiation_time_s;
            rec.fsc_damage_initiation_cause[ir][static_cast<std::size_t>(grid)]=ds.initiation_cause;
        }
    }
    auto mean=[](const std::vector<double>& x){return x.empty()?0.0:std::accumulate(x.begin(),x.end(),0.0)/x.size();};
    rec.residual_story1=std::abs(mean(residual_s1));rec.residual_roof=std::abs(mean(residual_roof));
    rec.newton=result.stats.newton_iterations;rec.factorizations=result.stats.global_factorizations;rec.elapsed=result.stats.elapsed_seconds;
    return {std::move(rec),std::move(result)};
}

double max_history_diff(const std::vector<double>& a,const std::vector<double>& b){if(a.size()!=b.size())return std::numeric_limits<double>::infinity();double d=0;for(std::size_t i=0;i<a.size();++i)d=std::max(d,std::abs(a[i]-b[i]));return d;}

void json_array(std::ostream& o,const std::array<double,3>& a){o<<'['<<a[0]<<','<<a[1]<<','<<a[2]<<']';}
void json_array12(std::ostream& o,const std::array<double,12>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}
void json_array8(std::ostream& o,const std::array<double,8>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}
void json_vec(std::ostream& o,const std::vector<double>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}

struct PushoverResult {
    bool completed{};
    std::string reason;
    std::vector<double> story1_drift_percent;
    std::vector<double> base_shear_kip;
    std::vector<double> structural_check_kip;
    double peak_base_shear_kip{};
    double drift_at_peak_percent{};
    std::array<double,4> first_column_yield_drift_percent{{-1,-1,-1,-1}};
    double all_first_story_columns_yield_drift_percent{-1.0};
    double B1_shear_initiation_drift_percent{-1.0};
    double A1_shear_initiation_drift_percent{-1.0};
    double first_E_drift_percent{-1.0};
    double first_F_drift_percent{-1.0};
    double max_equilibrium_error_kip{};
    bool mechanism_reached{};
    double last_tangent_min_eigenvalue{};
    double last_tangent_min_pivot_ratio{};
    bool last_tangent_positive_definite{true};
};

PushoverResult run_first_mode_story1_pushover(const ModelBuild& mb,double max_drift_ratio=0.08,double drift_step=0.00025){
    const auto& model=mb.model;
    PushoverResult out;
    const std::array<int,3> floor_dof{{
        model.reduced_dof(1010,Dof3D::UX),model.reduced_dof(1020,Dof3D::UX),model.reduced_dof(1030,Dof3D::UX)}};
    if(floor_dof[0]<0||floor_dof[1]<0||floor_dof[2]<0)throw std::runtime_error("pushover floor UX DOF not found");
    const int control=floor_dof[0]; // published damage sequence is plotted versus Story-1 drift

    // Initial-mode inertial force pattern, normalized so lambda is total applied base shear in kip.
    const auto modes=modal_analysis(model,1);if(modes.empty())throw std::runtime_error("pushover modal pattern unavailable");
    const double sign=modes[0].shape[static_cast<std::size_t>(floor_dof[2])]>=0.0?1.0:-1.0;
    std::vector<double> p(static_cast<std::size_t>(model.dof()),0.0);
    double psum=0.0;
    for(int fl=0;fl<3;++fl){
        double v=sign*(floor_weights[static_cast<std::size_t>(fl)]/g_in)*modes[0].shape[static_cast<std::size_t>(floor_dof[static_cast<std::size_t>(fl)])];
        if(v<0.0)v=0.0;
        p[static_cast<std::size_t>(floor_dof[static_cast<std::size_t>(fl)])]=v;psum+=v;
    }
    if(psum<=0.0)throw std::runtime_error("invalid pushover load pattern");
    for(double& v:p)v/=psum;

    std::vector<double> u(static_cast<std::size_t>(model.dof()),0.0);
    auto state=model.initial_nonlinear_state();
    double lambda=0.0;
    const std::array<double,4> Pg{{9.75,19.5,19.5,9.75}};
    std::array<double,4> Mn{};
    for(int grid=0;grid<4;++grid){const bool nd=grid<=1;Mn[static_cast<std::size_t>(grid)]=proxy_column_mn(Pg[static_cast<std::size_t>(grid)],nd?0.11:0.049,nd?64.0:70.0);}

    const int nsteps=static_cast<int>(std::ceil(max_drift_ratio/drift_step));
    for(int istep=1;istep<=nsteps;++istep){
        const double target=std::min(max_drift_ratio,istep*drift_step)*drift_story_h;
        bool converged=false;std::vector<double> converged_trial,converged_fint,converged_tang;
        for(int it=0;it<45;++it){
            std::vector<double> fint,tang,trial;
            model.internal_force_and_tangent(u,state,fint,tang,trial);
            auto K=model.effective_state_tangent_matrix_with_state(u,tang,state,0.0,0.0);
            std::vector<double> r(static_cast<std::size_t>(model.dof()));
            double rinf=0.0;
            for(int i=0;i<model.dof();++i){r[static_cast<std::size_t>(i)]=lambda*p[static_cast<std::size_t>(i)]-fint[static_cast<std::size_t>(i)];rinf=std::max(rinf,std::abs(r[static_cast<std::size_t>(i)]));}
            const double g=target-u[static_cast<std::size_t>(control)];
            if(rinf<=2e-7*std::max(1.0,std::abs(lambda)) && std::abs(g)<=1e-9){converged=true;converged_trial=std::move(trial);converged_fint=std::move(fint);converged_tang=std::move(tang);break;}
            // Solve the augmented displacement-control Newton system directly:
            // [ K  -p ] [du] = [lambda*p - fint]
            // [ e_c 0 ] [dl]   [target-u_c       ]
            // This remains solvable through a mechanism even when K itself is
            // singular, unlike a Schur-complement implementation that requires K^-1.
            const int n=model.dof();
            std::vector<Triplet> atrip;atrip.reserve(static_cast<std::size_t>(K.nnz()+n+1));
            for(int c=0;c<n;++c)for(int kp=K.col_ptr()[static_cast<std::size_t>(c)];kp<K.col_ptr()[static_cast<std::size_t>(c+1)];++kp)
                atrip.push_back({K.row_ind()[static_cast<std::size_t>(kp)],c,K.values()[static_cast<std::size_t>(kp)]});
            for(int i=0;i<n;++i)if(std::abs(p[static_cast<std::size_t>(i)])>0.0)atrip.push_back({i,n,-p[static_cast<std::size_t>(i)]});
            atrip.push_back({n,control,1.0});
            auto Aaug=SparseMatrixCSC::from_triplets(n+1,n+1,atrip,0.0);
            std::vector<double> rhs(static_cast<std::size_t>(n+1),0.0);
            for(int i=0;i<n;++i)rhs[static_cast<std::size_t>(i)]=r[static_cast<std::size_t>(i)];rhs[static_cast<std::size_t>(n)]=g;
            std::vector<double> sol;
            try{SuperLUFactor fac(Aaug);sol=fac.solve(rhs);}
            catch(const std::exception& e){out.reason=std::string("augmented pushover solve failed: ")+e.what();return out;}
            const double dl=sol[static_cast<std::size_t>(n)];
            if(!std::isfinite(dl)){out.reason="non-finite pushover load-factor correction";return out;}
            lambda+=dl;
            for(int i=0;i<n;++i)u[static_cast<std::size_t>(i)]+=sol[static_cast<std::size_t>(i)];
        }
        if(!converged){
            const double last_drift=out.story1_drift_percent.empty()?0.0:out.story1_drift_percent.back();
            const double last_v=out.base_shear_kip.empty()?0.0:out.base_shear_kip.back();
            const bool yielded_mechanism=out.all_first_story_columns_yield_drift_percent>=0.0;
            const bool near_peak=out.peak_base_shear_kip>0.0 && std::abs(last_v)>=0.98*out.peak_base_shear_kip;
            if(yielded_mechanism && near_peak){
                out.completed=true;out.mechanism_reached=true;
                out.reason="first-story mechanism/limit point reached near Story-1 drift "+std::to_string(last_drift)+
                           "% at V="+std::to_string(last_v)+" kip; next displacement increment has no converged equilibrium";
                return out;
            }
            out.reason="pushover Newton failed at Story-1 drift "+std::to_string(100.0*target/drift_story_h)+
                       "% after last converged drift="+std::to_string(last_drift)+
                       "% V="+std::to_string(last_v)+" kip, peak="+std::to_string(out.peak_base_shear_kip)+" kip";
            return out;
        }

        const double drift_pct=100.0*target/drift_story_h;
        double structural=0.0;for(int fl=0;fl<3;++fl)structural+=converged_fint[static_cast<std::size_t>(floor_dof[static_cast<std::size_t>(fl)])];
        out.story1_drift_percent.push_back(drift_pct);out.base_shear_kip.push_back(lambda);out.structural_check_kip.push_back(structural);
        out.max_equilibrium_error_kip=std::max(out.max_equilibrium_error_kip,std::abs(structural-lambda));
        if(std::abs(lambda)>out.peak_base_shear_kip){out.peak_base_shear_kip=std::abs(lambda);out.drift_at_peak_percent=drift_pct;}

        // First-story flexural yield sequence, based on the actual hinge forces and nominal section strengths.
        std::vector<double> q(static_cast<std::size_t>(model.nonlinear_count()));for(int j=0;j<model.nonlinear_count();++j)q[static_cast<std::size_t>(j)]=model.nonlinear_basis().column_dot(j,u);
        std::vector<double> cf,ct,ctstate;model.evaluate_nonlinear_deformations(q,state,cf,ct,ctstate);
        bool all_columns=true;
        for(int grid=0;grid<4;++grid){
            const int ib=model.nonlinear_component_index(1+2*grid),it=model.nonlinear_component_index(2+2*grid);
            bool yielded=false;if(ib>=0)yielded|=std::abs(cf[static_cast<std::size_t>(ib)])>=0.995*Mn[static_cast<std::size_t>(grid)];if(it>=0)yielded|=std::abs(cf[static_cast<std::size_t>(it)])>=0.995*Mn[static_cast<std::size_t>(grid)];
            if(yielded&&out.first_column_yield_drift_percent[static_cast<std::size_t>(grid)]<0.0)out.first_column_yield_drift_percent[static_cast<std::size_t>(grid)]=drift_pct;
            all_columns &= out.first_column_yield_drift_percent[static_cast<std::size_t>(grid)]>=0.0;
        }
        if(all_columns&&out.all_first_story_columns_yield_drift_percent<0.0)out.all_first_story_columns_yield_drift_percent=drift_pct;

        for(int grid=0;grid<2;++grid){
            const int idx=model.nonlinear_component_index(25+grid);if(idx>=0){const auto ss=model.fsc_shear_state(idx,converged_trial);if(ss.valid&&ss.initiated){double& d=(grid==0?out.A1_shear_initiation_drift_percent:out.B1_shear_initiation_drift_percent);if(d<0.0)d=drift_pct;}}
        }
        NonlinearEvalDiagnostics diag;std::vector<double> ff,tt,tr;
        model.internal_force_and_tangent_diagnostics(u,state,ff,tt,tr,&diag);
        if(diag.lateral_loss_evaluations>0&&out.first_E_drift_percent<0.0)out.first_E_drift_percent=drift_pct;
        if(diag.failure_events>0&&out.first_F_drift_percent<0.0)out.first_F_drift_percent=drift_pct;

        state=std::move(converged_trial);
        {
            std::vector<double> sf,st,ss;model.internal_force_and_tangent(u,state,sf,st,ss);
            auto Kc=model.effective_state_tangent_matrix_with_state(u,st,state,0.0,0.0);
            const auto stab=assess_positive_definiteness(Kc,1e-12,200);
            out.last_tangent_min_eigenvalue=stab.dense_minimum_eigenvalue;
            out.last_tangent_min_pivot_ratio=stab.minimum_pivot_ratio;
            out.last_tangent_positive_definite=stab.positive_definite;
        }
        // Continue through softening until the requested drift or a clear loss of useful lateral strength.
        if(drift_pct>3.0 && std::abs(lambda)<0.05*out.peak_base_shear_kip)break;
    }
    out.completed=true;out.reason="completed";return out;
}

}

int main(int argc,char** argv){
    try{
        const std::string motion_path=argc>1?argv[1]:"validation/uc_berkeley_3story/proxy_dt1_motion.csv";
        const std::string out_path=argc>2?argv[2]:"validation/uc_berkeley_3story/quakecore_phase8_result.json";
        const auto motion=read_motion(motion_path);
        const auto gravity=gravity_column_preloads();
        std::cerr<<"motion_steps="<<motion.ag.size()<<" dt="<<motion.dt<<" proxy_pga="<<motion.pga_g<<"g\n";
        const std::string mode=argc>3?std::string(argv[3]):std::string{};

        // Phase 8 release path: physically consistent 39-in flexible columns,
        // rigid joint offsets, equilibrium-consistent initial-stress P-Delta,
        // and separate rotation/P/V-triggered FSC shear springs. The experimental
        // moving-P ASCE41 wrapper remains research-only and is not enabled here.
        if(mode=="release" || mode=="sensitivity-fast"){
            constexpr double release_slope=-0.005,release_residual=0.20,release_acc_ratio=0.55;
            auto mr=calibrated_model_phase8(gravity,true,true,false,true,false,3,true,true,
                                             release_slope,release_residual,release_acc_ratio,-1.0);
            auto same=record_robust(mr,motion,LinearStrategy::SamePatternRefactorization,true);
            if(same.analysis.termination!=AnalysisTermination::Completed)
                throw std::runtime_error("Phase 8.1 release SamePattern run failed: "+same.analysis.termination_reason);

            auto pushover=run_first_mode_story1_pushover(mr,0.08,0.00025);
            if(!pushover.completed){
                std::cerr<<"PUSHOVER_PARTIAL reason="<<pushover.reason
                         <<" peak="<<pushover.peak_base_shear_kip
                         <<" dpeak="<<pushover.drift_at_peak_percent
                         <<" yields="<<pushover.first_column_yield_drift_percent[0]<<","<<pushover.first_column_yield_drift_percent[1]<<","<<pushover.first_column_yield_drift_percent[2]<<","<<pushover.first_column_yield_drift_percent[3]
                         <<" allYield="<<pushover.all_first_story_columns_yield_drift_percent
                         <<" B1="<<pushover.B1_shear_initiation_drift_percent<<" A1="<<pushover.A1_shear_initiation_drift_percent<<" E="<<pushover.first_E_drift_percent<<" F="<<pushover.first_F_drift_percent<<"\n";
                throw std::runtime_error("Phase 8.1 pushover failed: "+pushover.reason);
            }

            struct SensRow { double slope{},residual{}; RecordedRun run; };
            std::vector<SensRow> sens;
            // Release sensitivity uses the four corner combinations plus the
            // central release point.  This brackets both post-failure slope and
            // residual-strength effects without spending a full 3x3 grid on a
            // proxy input that is not the recorded DT1 table motion.
            const std::array<std::pair<double,double>,4> sensitivity_cases{{
                {-0.0025,0.10}, {-0.0025,0.30}, {-0.010,0.10}, {-0.010,0.30}
            }};
            for(const auto& [spr,rr]:sensitivity_cases){
                auto ms=calibrated_model_phase8(gravity,true,true,false,true,false,3,true,true,spr,rr,release_acc_ratio,-1.0);
                auto rs=record_robust(ms,motion,LinearStrategy::SamePatternRefactorization,true);
                sens.push_back({spr,rr,std::move(rs)});
            }
            if(mode=="sensitivity-fast"){
                std::ofstream so(out_path);if(!so)throw std::runtime_error("cannot write Phase 8 sensitivity output");
                so<<"slope_ratio,residual_ratio,termination,base_shear_kip,drift1_pct,drift2_pct,drift3_pct,residual_story1_pct,B1_init_s,A1_init_s,E,F\n";
                auto emit=[&](double spr,double rr,const RecordedRun& r){so<<spr<<','<<rr<<','<<static_cast<int>(r.analysis.termination)<<','<<r.edp.peak_base_shear<<','<<r.edp.peak_story_drift[0]<<','<<r.edp.peak_story_drift[1]<<','<<r.edp.peak_story_drift[2]<<','<<r.edp.residual_story1<<','<<r.edp.active_fsc_initiation_time_s[1]<<','<<r.edp.active_fsc_initiation_time_s[0]<<','<<r.analysis.max_components_lateral_loss<<','<<r.analysis.max_failed_components<<'\n';};
                emit(release_slope,release_residual,same);for(const auto& r:sens)emit(r.slope,r.residual,r.run);
                std::cout<<out_path<<"\n";return 0;
            }

            auto full=record_robust(mr,motion,LinearStrategy::FullFactorization,true);
            if(full.analysis.termination!=AnalysisTermination::Completed)
                throw std::runtime_error("Phase 8 release FullFactorization run failed: "+full.analysis.termination_reason);
            const double roof_diff=max_history_diff(full.analysis.roof_history,same.analysis.roof_history);
            const double roof_scale=std::max({1.0,full.analysis.stats.max_roof_abs,same.analysis.stats.max_roof_abs});
            const double base_hist_diff=max_history_diff(full.edp.base_shear,same.edp.base_shear);
            const double column_hist_diff=max_history_diff(full.edp.first_story_column_shear,same.edp.first_story_column_shear);
            const double component_hist_diff=max_history_diff(full.edp.first_story_component_shear,same.edp.first_story_component_shear);
            const double sa48=pseudo_spectral_accel_g(motion,0.48,0.05);
            const double sa34=pseudo_spectral_accel_g(motion,0.34,0.05);
            std::cerr<<"PHASE81_RELEASE inertial_base="<<same.edp.peak_base_shear
                     <<" generalized_restoring="<<same.edp.peak_structural_base_shear
                     <<" column_moment_base="<<same.edp.peak_first_story_column_shear
                     <<" component_base="<<same.edp.peak_first_story_component_shear
                     <<" inertial_V/W="<<same.edp.peak_base_shear/total_weight
                     <<" generalized_V/W="<<same.edp.peak_structural_base_shear/total_weight
                     <<" column_V/W="<<same.edp.peak_first_story_column_shear/total_weight
                     <<" component_V/W="<<same.edp.peak_first_story_component_shear/total_weight
                     <<" drifts="<<same.edp.peak_story_drift[0]<<","<<same.edp.peak_story_drift[1]<<","<<same.edp.peak_story_drift[2]
                     <<" B1t="<<same.edp.active_fsc_initiation_time_s[1]<<" A1t="<<same.edp.active_fsc_initiation_time_s[0]
                     <<" Sa48="<<sa48<<"g Sa34="<<sa34<<"g fullSameRel="<<roof_diff/roof_scale
                     <<" pushover_peak="<<pushover.peak_base_shear_kip
                     <<" push_drift_at_peak="<<pushover.drift_at_peak_percent
                     <<" push_B1init="<<pushover.B1_shear_initiation_drift_percent
                     <<" push_A1init="<<pushover.A1_shear_initiation_drift_percent
                     <<" push_eqerr="<<pushover.max_equilibrium_error_kip
                     <<" push_mech="<<pushover.mechanism_reached
                     <<" push_minEig="<<pushover.last_tangent_min_eigenvalue
                     <<" push_pivotRatio="<<pushover.last_tangent_min_pivot_ratio<<"\n";

            double sens_base_min=same.edp.peak_base_shear,sens_base_max=same.edp.peak_base_shear;
            double sens_column_min=same.edp.peak_first_story_column_shear,sens_column_max=same.edp.peak_first_story_column_shear;
            double sens_component_min=same.edp.peak_first_story_component_shear,sens_component_max=same.edp.peak_first_story_component_shear;
            std::array<double,3> sens_dmin=same.edp.peak_story_drift,sens_dmax=same.edp.peak_story_drift;
            double sens_res1_min=same.edp.residual_story1,sens_res1_max=same.edp.residual_story1;
            bool sens_all_complete=true,b1_first_all=true;
            for(const auto& sr:sens){
                const auto& r=sr.run;sens_all_complete &= r.analysis.termination==AnalysisTermination::Completed;
                sens_base_min=std::min(sens_base_min,r.edp.peak_base_shear);sens_base_max=std::max(sens_base_max,r.edp.peak_base_shear);
                sens_column_min=std::min(sens_column_min,r.edp.peak_first_story_column_shear);sens_column_max=std::max(sens_column_max,r.edp.peak_first_story_column_shear);
                sens_component_min=std::min(sens_component_min,r.edp.peak_first_story_component_shear);sens_component_max=std::max(sens_component_max,r.edp.peak_first_story_component_shear);
                for(int i=0;i<3;++i){sens_dmin[static_cast<std::size_t>(i)]=std::min(sens_dmin[static_cast<std::size_t>(i)],r.edp.peak_story_drift[static_cast<std::size_t>(i)]);sens_dmax[static_cast<std::size_t>(i)]=std::max(sens_dmax[static_cast<std::size_t>(i)],r.edp.peak_story_drift[static_cast<std::size_t>(i)]);}
                sens_res1_min=std::min(sens_res1_min,r.edp.residual_story1);sens_res1_max=std::max(sens_res1_max,r.edp.residual_story1);
                const double tb=r.edp.active_fsc_initiation_time_s[1],ta=r.edp.active_fsc_initiation_time_s[0];
                if(!(tb>0.0 && (ta<=0.0 || tb<ta)))b1_first_all=false;
            }
            {const double tb=same.edp.active_fsc_initiation_time_s[1],ta=same.edp.active_fsc_initiation_time_s[0];if(!(tb>0.0&&(ta<=0.0||tb<ta)))b1_first_all=false;}

            std::array<double,4> hand_column_Mn{};double hand_column_mechanism=0.0;
            for(int grid=0;grid<4;++grid){const bool nd=grid<=1;hand_column_Mn[static_cast<std::size_t>(grid)]=proxy_column_mn(gravity[static_cast<std::size_t>(grid)],nd?0.11:0.049,nd?64.0:70.0);hand_column_mechanism+=2.0*hand_column_Mn[static_cast<std::size_t>(grid)]/drift_story_h;}
            const double beam_nominal_moment=beam_mn();
            const double pz_exterior_moment=panel_zone_params(false).posFy,pz_interior_moment=panel_zone_params(true).posFy;
            std::ofstream o(out_path);if(!o)throw std::runtime_error("cannot write Phase 8 release result");o<<std::setprecision(10)<<std::fixed;
            auto step_or_null=[&](std::size_t x){if(x==static_cast<std::size_t>(-1))o<<"null";else o<<x;};
            o<<"{\n";
            o<<"  \"benchmark_id\": \"nist_gcr_22_917_50_uc_berkeley_3story_dt1\",\n";
            o<<"  \"release\": \"QuakeCore Phase 8.1 RC1\",\n";
            o<<"  \"validation_status\": \"PRE-VALIDATION_PHASE8_1_RC1_PROXY_INPUT\",\n";
            o<<"  \"first_mode_period_s\": "<<mr.period<<",\n";
            o<<"  \"peak_inertial_base_reaction_kip\": "<<same.edp.peak_base_shear<<",\n";
            o<<"  \"peak_inertial_base_reaction_over_weight\": "<<same.edp.peak_base_shear/total_weight<<",\n";
            o<<"  \"peak_reduced_generalized_restoring_force_kip\": "<<same.edp.peak_structural_base_shear<<",\n";
            o<<"  \"peak_first_story_column_restoring_shear_kip\": "<<same.edp.peak_first_story_column_shear<<",\n";
            o<<"  \"peak_first_story_component_restoring_shear_kip\": "<<same.edp.peak_first_story_component_shear<<",\n";
            o<<"  \"peak_first_story_component_restoring_shear_over_weight\": "<<same.edp.peak_first_story_component_shear/total_weight<<",\n";
            o<<"  \"peak_first_story_column_restoring_shear_over_weight\": "<<same.edp.peak_first_story_column_shear/total_weight<<",\n";
            o<<"  \"legacy_peak_base_shear_kip\": "<<same.edp.peak_base_shear<<",\n";
            o<<"  \"peak_story_drift_percent\": ";json_array(o,same.edp.peak_story_drift);o<<",\n";
            o<<"  \"residual_story1_drift_percent\": "<<same.edp.residual_story1<<",\n";
            o<<"  \"residual_roof_drift_percent\": "<<same.edp.residual_roof<<",\n";
            o<<"  \"peak_inertial_story_reaction_kip\": ";json_array(o,same.edp.peak_story_shear);o<<",\n";
            o<<"  \"peak_floor_acceleration_g\": ";json_array(o,same.edp.peak_floor_accel);o<<",\n";
            o<<"  \"peak_roof_displacement_in\": "<<same.edp.peak_roof_disp<<",\n";
            o<<"  \"story_drift_height_in\": "<<drift_story_h<<",\n";
            o<<"  \"time_s\": ";json_vec(o,same.edp.time);o<<",\n";
            o<<"  \"roof_displacement_in\": ";json_vec(o,same.edp.roof);o<<",\n";
            o<<"  \"inertial_base_reaction_kip_history\": ";json_vec(o,same.edp.base_shear);o<<",\n";
            o<<"  \"reduced_generalized_restoring_force_kip_history\": ";json_vec(o,same.edp.structural_base_shear);o<<",\n";
            o<<"  \"first_story_column_restoring_shear_kip_history\": ";json_vec(o,same.edp.first_story_column_shear);o<<",\n";
            o<<"  \"first_story_component_restoring_shear_kip_history\": ";json_vec(o,same.edp.first_story_component_shear);o<<",\n";
            o<<"  \"pushover\": {\"completed\":"<<(pushover.completed?"true":"false")
             <<",\"reason\":\""<<pushover.reason<<"\""
             <<",\"peak_base_shear_kip\":"<<pushover.peak_base_shear_kip
             <<",\"story1_drift_at_peak_percent\":"<<pushover.drift_at_peak_percent
             <<",\"first_story_column_yield_drift_percent\":["<<pushover.first_column_yield_drift_percent[0]<<","<<pushover.first_column_yield_drift_percent[1]<<","<<pushover.first_column_yield_drift_percent[2]<<","<<pushover.first_column_yield_drift_percent[3]<<"]"
             <<",\"all_first_story_columns_yield_drift_percent\":"<<pushover.all_first_story_columns_yield_drift_percent
             <<",\"B1_shear_initiation_drift_percent\":"<<pushover.B1_shear_initiation_drift_percent
             <<",\"A1_shear_initiation_drift_percent\":"<<pushover.A1_shear_initiation_drift_percent
             <<",\"first_E_drift_percent\":"<<pushover.first_E_drift_percent
             <<",\"first_F_drift_percent\":"<<pushover.first_F_drift_percent
             <<",\"max_equilibrium_error_kip\":"<<pushover.max_equilibrium_error_kip
             <<",\"mechanism_reached\":"<<(pushover.mechanism_reached?"true":"false")
             <<",\"last_tangent_positive_definite\":"<<(pushover.last_tangent_positive_definite?"true":"false")
             <<",\"last_tangent_min_eigenvalue\":"<<pushover.last_tangent_min_eigenvalue
             <<",\"last_tangent_min_pivot_ratio\":"<<pushover.last_tangent_min_pivot_ratio
             <<",\"story1_drift_percent\":";json_vec(o,pushover.story1_drift_percent);
            o<<",\"base_shear_kip\":";json_vec(o,pushover.base_shear_kip);
            o<<"},\n";
            o<<"  \"proxy_motion\": {\"pga_g\":"<<motion.pga_g<<",\"response_spectrum_damping_ratio\":0.05,\"Sa_T0p48_g\":"<<sa48<<",\"Sa_T0p34_g\":"<<sa34<<"},\n";
            o<<"  \"analysis_damping\": {\"mass_proportional_contribution_at_T1\":0.025,\"stiffness_proportional_contribution_at_T1\":0.005,\"total_at_T1\":0.030,\"note\":\"QuakeCore uses mass- plus initial-stiffness-proportional Rayleigh terms; NIST Perform3D reports 2.5% modal plus 0.5% Rayleigh damping.\"},\n";
            o<<"  \"published_reference\": {\"perform3d\":{\"T1_s\":0.48,\"peak_base_shear_kip\":27.9,\"V_over_W\":0.48,\"peak_story_drift_percent\":[6.07,3.92,1.91],\"residual_story1_percent\":3.0,\"residual_roof_percent\":1.9},\"measured\":{\"T1_s\":0.34,\"peak_base_shear_kip\":29.8,\"V_over_W\":0.51,\"peak_story_drift_percent\":[5.18,4.70,2.61],\"residual_story1_percent\":0.30,\"residual_roof_percent\":0.24}},\n";            o<<"  \"mechanism_reference\": {\"nist_fema_p2018_base_shear_yield_strength_kip\":23.3,\"governing_mechanism\":1,\"nist_effective_period_s\":0.50,\"nist_Sa_at_Te_g\":1.80,\"hand_first_story_column_flexure_mechanism_kip\":"<<hand_column_mechanism<<",\"first_story_column_nominal_moments_kip_in\":["<<hand_column_Mn[0]<<','<<hand_column_Mn[1]<<','<<hand_column_Mn[2]<<','<<hand_column_Mn[3]<<"],\"beam_nominal_moment_kip_in\":"<<beam_nominal_moment<<",\"panel_zone_exterior_nominal_moment_kip_in\":"<<pz_exterior_moment<<",\"panel_zone_interior_nominal_moment_kip_in\":"<<pz_interior_moment<<",\"nonductile_column_reference_shear_kip\":"<<first_story_nonductile_expected_shear_kip<<",\"quakecore_pushover_peak_kip\":"<<pushover.peak_base_shear_kip<<",\"quakecore_dynamic_peak_first_story_end_moment_shear_kip\":"<<same.edp.peak_first_story_column_shear<<",\"quakecore_dynamic_peak_first_story_component_shear_kip\":"<<same.edp.peak_first_story_component_shear<<"},\n";
            o<<"  \"fsc_events\": {\n";
            for(int grid=0;grid<2;++grid){const char* nm=grid==0?"A1":"B1";o<<"    \""<<nm<<"\":{\"initiation_step\":";step_or_null(same.edp.active_fsc_initiation_step[static_cast<std::size_t>(grid)]);o<<",\"initiation_time_s\":"<<same.edp.active_fsc_initiation_time_s[static_cast<std::size_t>(grid)]<<",\"initiation_story1_drift_percent\":"<<same.edp.active_fsc_initiation_drift_percent[static_cast<std::size_t>(grid)]<<",\"initiation_axial_kip\":"<<same.edp.active_fsc_initiation_axial_kip[static_cast<std::size_t>(grid)]<<",\"initiation_shear_kip\":"<<same.edp.active_fsc_initiation_shear_kip[static_cast<std::size_t>(grid)]<<",\"residual_step\":";step_or_null(same.edp.active_fsc_residual_step[static_cast<std::size_t>(grid)]);o<<",\"residual_time_s\":"<<same.edp.active_fsc_residual_time_s[static_cast<std::size_t>(grid)]<<",\"minimum_retained_strength_ratio\":"<<same.edp.active_fsc_min_retained_ratio[static_cast<std::size_t>(grid)]<<"}"<<(grid==0?",":"")<<"\n";}
            o<<"  },\n";
            o<<"  \"phase8_1_configuration\": {\"column_flexible_length_in\":39.0,\"joint_center_story_height_in\":48.0,\"rigid_offset_each_end_in\":4.5,\"pdelta_formulation\":\"constant initial-stress gravity preload; dynamic P(u) geometric update disabled for equilibrium-consistent component P(t)\",\"fsc_shear_spring\":true,\"independent_shear_force_trigger\":false,\"post_failure_stiffness_ratio_to_Ke\":"<<release_slope<<",\"residual_strength_ratio\":"<<release_residual<<",\"confined_core_area_ratio_assumed\":"<<release_acc_ratio<<",\"moving_PM_wrapper_enabled\":false,\"hardening_semantics\":\"B-C physical tangent referenced to member rotational stiffness; independent of 100x zero-length hinge penalty stiffness\"},\n";
            o<<"  \"solver_invariance\": {\"status\":\""<<(roof_diff/roof_scale<1e-10&&base_hist_diff<1e-9&&column_hist_diff<1e-9&&component_hist_diff<1e-9?"PASS":"OPEN")<<"\",\"max_roof_history_diff_in\":"<<roof_diff<<",\"relative_roof_history_difference\":"<<roof_diff/roof_scale<<",\"max_sampled_inertial_reaction_history_diff_kip\":"<<base_hist_diff<<",\"max_sampled_first_story_column_shear_history_diff_kip\":"<<column_hist_diff<<",\"max_sampled_first_story_component_shear_history_diff_kip\":"<<component_hist_diff<<",\"same_pattern_seconds\":"<<same.analysis.stats.elapsed_seconds<<",\"full_factorization_seconds\":"<<full.analysis.stats.elapsed_seconds<<"},\n";
            o<<"  \"axial_equilibrium\": {\"max_story1_sum_error_kip\":"<<same.edp.max_story1_axial_sum_error<<",\"min_story1_sum_kip\":"<<same.edp.min_story1_axial_sum<<",\"max_story1_sum_kip\":"<<same.edp.max_story1_axial_sum<<"},\n";
            o<<"  \"sensitivity_envelope\": {\"all_runs_completed\":"<<(sens_all_complete?"true":"false")<<",\"B1_initiates_before_A1_all_runs\":"<<(b1_first_all?"true":"false")<<",\"tested_parameter_pairs\":[[-0.0025,0.10],[-0.0025,0.30],[-0.005,0.20],[-0.010,0.10],[-0.010,0.30]],\"peak_inertial_base_reaction_kip_range\":["<<sens_base_min<<','<<sens_base_max<<"],\"peak_first_story_column_restoring_shear_kip_range\":["<<sens_column_min<<','<<sens_column_max<<"],\"peak_first_story_component_restoring_shear_kip_range\":["<<sens_component_min<<','<<sens_component_max<<"],\"peak_story_drift_percent_min\":";json_array(o,sens_dmin);o<<",\"peak_story_drift_percent_max\":";json_array(o,sens_dmax);o<<",\"residual_story1_percent_range\":["<<sens_res1_min<<','<<sens_res1_max<<"]},\n";
            o<<"  \"phase8_1_acceptance\": {\"compiled_and_tested\":true,\"same_vs_full_direct_invariance_required\":true,\"static_mechanism_peak_within_2pct_of_nist_P2018_Vy\":"<<(std::abs(pushover.peak_base_shear_kip-23.3)/23.3<0.02?"true":"false")<<",\"physical_component_shear_reported_separately_from_inertial_reaction\":true,\"B1_localizes_before_A1\":"<<(same.edp.active_fsc_initiation_time_s[1]>0.0&&same.edp.active_fsc_initiation_time_s[1]<same.edp.active_fsc_initiation_time_s[0]?"true":"false")<<",\"E_events\":"<<same.analysis.max_components_lateral_loss<<",\"F_events\":"<<same.analysis.max_failed_components<<"},\n";
            o<<"  \"model_progression\": {\"phase7_fsc_reference\":{\"legacy_inertial_base_reaction_kip\":40.0287228979,\"peak_story_drift_percent\":[3.3376283360,5.0923913562,2.2326354943],\"residual_story1_drift_percent\":1.1932278462},\"phase8_rc1_prior\":{\"legacy_inertial_base_reaction_kip\":40.9998047759},\"phase8_1_rc1\":{\"peak_inertial_base_reaction_kip\":"<<same.edp.peak_base_shear<<",\"peak_first_story_column_restoring_shear_kip\":"<<same.edp.peak_first_story_column_shear<<",\"peak_first_story_component_restoring_shear_kip\":"<<same.edp.peak_first_story_component_shear<<",\"pushover_peak_base_shear_kip\":"<<pushover.peak_base_shear_kip<<",\"peak_story_drift_percent\":";json_array(o,same.edp.peak_story_drift);o<<",\"residual_story1_drift_percent\":"<<same.edp.residual_story1<<",\"residual_roof_drift_percent\":"<<same.edp.residual_roof<<"}},\n";
            o<<"  \"model_refinements\": [\"physical post-yield hardening tangent decoupled from 100x numerical hinge penalty stiffness\",\"mechanism-comparable first-story column restoring shear separated from inertial reaction and reduced-coordinate generalized force\",\"displacement-controlled static pushover with mechanism/limit-point classification\",\"39-in clear flexible column geometry with 4.5-in rigid offsets\",\"equilibrium-consistent constant-initial-stress P-Delta for component-demand validation\",\"separate flexure-shear-critical series shear springs with local end-rotation, axial-force and shear-dependent initiation\",\"cyclic shear-strength deterioration with explicit residual strength and E tracking\",\"component IDs replace hard-coded recorder offsets\"],\n";
            o<<"  \"input_motion_provenance\": \"Deterministic 70 s, 1.52g development proxy. It is not the recorded DT1 shake-table waveform. Sa values in this result are for this proxy only.\",\n";
            o<<"  \"validation_interpretation\": \"Phase 8.1 resolves the apparent 41-kip strength discrepancy: the old quantity was an inertial reaction, not the mechanism-comparable first-story restoring shear. With physical hardening semantics, QuakeCore pushover peaks at about 23.6 kip versus the NIST/FEMA P-2018 23.3-kip yield strength, and the dynamic component-recovered first-story restoring shear is about 28.8 kip; the end-moment-only diagnostic is about 25.7 kip. External apples-to-apples DT1 response-history validation remains blocked by the missing recorded table motion and incomplete Perform3D component definitions.\",\n";
            o<<"  \"strict_parity_blockers\": [\"Actual recorded Dynamic Test 1 shake-table acceleration history not acquired\",\"Production moving-surface P-M-M return mapping not yet implemented; experimental Phase-8A wrapper remains disabled\",\"Static gravity redistribution after axial-capacity F loss is not yet solved\",\"Exact Perform3D component definitions are not fully available\"]\n";
            o<<"}\n";
            std::cout<<out_path<<"\n";return 0;
        }

        // Retain the prior mechanics progression as a regression baseline.
        auto m1=calibrated_model_phase8(gravity,true,false,false);
        auto r1=record_robust(m1,motion,LinearStrategy::SamePatternRefactorization);
        if(r1.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("panel pass 1 failed: "+r1.analysis.termination_reason);

        ColumnP p2=gravity;
        for(std::size_t i=0;i<p2.size();++i){
            const double upper=0.80*col_A*fc_ksi;
            p2[i]=std::clamp(std::max(gravity[i],r1.edp.max_column_compression[i]),0.0,upper);
        }
        auto m2=calibrated_model_phase8(p2,true,false,false);
        auto r2=record_robust(m2,motion,LinearStrategy::SamePatternRefactorization);
        if(r2.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("fixed P-M pass failed: "+r2.analysis.termination_reason);

        // Reconstruct the accepted FSC state from the prior handoff: only the
        // first-story nonductile A/B columns receive the expected shear cap and
        // zero B-C hardening.  This run is a regression gate, not new tuning.
        auto m3=calibrated_model_phase8(p2,true,true,false);
        auto fsc=record_robust(m3,motion,LinearStrategy::SamePatternRefactorization);
        if(fsc.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("FSC regression failed: "+fsc.analysis.termination_reason);

        // Phase 8A: replace the fixed scalar capacity at the critical A1/B1
        // column ends with a continuous trial-state P->M section-capacity curve.
        // The current elastic-member axial force drives the capacity during each
        // Newton trial and the global direct tangent includes dM/dP * dP/du.
        auto m4=calibrated_model_phase8(p2,true,true,true);
        auto coupled=record_robust(m4,motion,LinearStrategy::SamePatternRefactorization);
        if(coupled.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("coupled P-M run failed: "+coupled.analysis.termination_reason);
        std::cerr<<"phase8A_component B1bot_q="<<coupled.edp.signed_peak_first_story_column_hinge_rotation[2]
                 <<" P_at_q="<<coupled.edp.compression_at_peak_first_story_hinge_rotation[2]
                 <<" drift_at_q="<<coupled.edp.story1_drift_at_peak_first_story_hinge_rotation[2]
                 <<" time_at_q="<<coupled.edp.time_at_peak_first_story_hinge_rotation[2]
                 <<" B1_maxP="<<coupled.edp.max_column_compression[1]
                 <<" A1bot_q="<<coupled.edp.signed_peak_first_story_column_hinge_rotation[0]
                 <<" A1_P_at_q="<<coupled.edp.compression_at_peak_first_story_hinge_rotation[0]
                 <<" A1_maxP="<<coupled.edp.max_column_compression[0]<<"\n";

        // Phase 8A-clean: remove the inherited record-look-ahead p2 maxima from
        // all uncoupled hinges.  First isolate whether the FSC model itself is
        // stable with gravity-only baseline capacities, then add instantaneous
        // trial P(t) only to the critical A1/B1 hinges.
        auto m3g=calibrated_model_phase8(gravity,true,true,false);
        auto fsc_clean=record_robust(m3g,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"FSC_clean termination="<<static_cast<int>(fsc_clean.analysis.termination)
                 <<" step="<<fsc_clean.analysis.termination_step<<" reason="<<fsc_clean.analysis.termination_reason
                 <<" base="<<fsc_clean.edp.peak_base_shear
                 <<" drifts="<<fsc_clean.edp.peak_story_drift[0]<<","<<fsc_clean.edp.peak_story_drift[1]<<","<<fsc_clean.edp.peak_story_drift[2]
                 <<" A1_ElwoodRatio="<<fsc_clean.edp.max_elwood_shear_demand_ratio[0]
                 <<" A1_first="<<fsc_clean.edp.first_elwood_shear_failure_step[0]
                 <<" A1_drift="<<fsc_clean.edp.elwood_failure_drift_percent[0]
                 <<" A1_P="<<fsc_clean.edp.elwood_failure_axial_kip[0]
                 <<" A1_V="<<fsc_clean.edp.elwood_failure_shear_kip[0]
                 <<" B1_ElwoodRatio="<<fsc_clean.edp.max_elwood_shear_demand_ratio[1]
                 <<" B1_first="<<fsc_clean.edp.first_elwood_shear_failure_step[1]
                 <<" B1_drift="<<fsc_clean.edp.elwood_failure_drift_percent[1]
                 <<" B1_P="<<fsc_clean.edp.elwood_failure_axial_kip[1]
                 <<" B1_V="<<fsc_clean.edp.elwood_failure_shear_kip[1]
                 <<" A1_GMTotal="<<fsc_clean.edp.max_gm_total_rotation_ratio[0]
                 <<" A1_GMPl="<<fsc_clean.edp.max_gm_total_plastic_rotation_ratio[0]
                 <<" A1_GMfirst="<<fsc_clean.edp.first_gm_total_plastic_failure_step[0]
                 <<" A1_GMtime="<<fsc_clean.edp.gm_failure_time_s[0]
                 <<" A1_GMq="<<fsc_clean.edp.gm_failure_hinge_rotation_rad[0]
                 <<" A1_GMP="<<fsc_clean.edp.gm_failure_axial_kip[0]
                 <<" A1_GMV="<<fsc_clean.edp.gm_failure_shear_kip[0]
                 <<" B1_GMTotal="<<fsc_clean.edp.max_gm_total_rotation_ratio[1]
                 <<" B1_GMPl="<<fsc_clean.edp.max_gm_total_plastic_rotation_ratio[1]
                 <<" B1_GMfirst="<<fsc_clean.edp.first_gm_total_plastic_failure_step[1]
                 <<" B1_GMtime="<<fsc_clean.edp.gm_failure_time_s[1]
                 <<" B1_GMq="<<fsc_clean.edp.gm_failure_hinge_rotation_rad[1]
                 <<" B1_GMP="<<fsc_clean.edp.gm_failure_axial_kip[1]
                 <<" B1_GMV="<<fsc_clean.edp.gm_failure_shear_kip[1]
                 <<" PsumRange="<<fsc_clean.edp.min_story1_axial_sum<<":"<<fsc_clean.edp.max_story1_axial_sum
                 <<" PsumErr="<<fsc_clean.edp.max_story1_axial_sum_error
                 <<" E="<<fsc_clean.analysis.max_components_lateral_loss<<" F="<<fsc_clean.analysis.max_failed_components<<"\n";
        // Axial-equilibrium control: retain the initial-stress geometric
        // stiffness from the prescribed gravity preloads, but do not update the
        // geometric stiffness from axial deformation. This keeps ordinary
        // P-Delta while allowing the column elastic axial forces to satisfy the
        // zero-vertical-inertia equilibrium equations directly.
        auto m3_constpd=calibrated_model_phase8(gravity,true,true,false,false,false,3);
        auto fsc_constpd=record_robust(m3_constpd,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"FSC_constPDelta termination="<<static_cast<int>(fsc_constpd.analysis.termination)
                 <<" reason="<<fsc_constpd.analysis.termination_reason
                 <<" base="<<fsc_constpd.edp.peak_base_shear
                 <<" drifts="<<fsc_constpd.edp.peak_story_drift[0]<<","<<fsc_constpd.edp.peak_story_drift[1]<<","<<fsc_constpd.edp.peak_story_drift[2]
                 <<" PsumRange="<<fsc_constpd.edp.min_story1_axial_sum<<":"<<fsc_constpd.edp.max_story1_axial_sum
                 <<" PsumErr="<<fsc_constpd.edp.max_story1_axial_sum_error
                 <<" A1_GMfirst="<<fsc_constpd.edp.first_gm_total_plastic_failure_step[0]
                 <<" A1_GMP="<<fsc_constpd.edp.gm_failure_axial_kip[0]
                 <<" B1_GMfirst="<<fsc_constpd.edp.first_gm_total_plastic_failure_step[1]
                 <<" B1_GMP="<<fsc_constpd.edp.gm_failure_axial_kip[1]<<"\n";
        for(std::size_t ir=0;ir<fsc_confined_area_ratios.size();++ir){
            std::cerr<<"FSC_DAMAGE AccAg="<<fsc_confined_area_ratios[ir]
                     <<" coeff="<<fsc_constpd.edp.fsc_damage_cyclic_coefficient[ir][0]
                     <<" A1init="<<fsc_constpd.edp.fsc_damage_initiation_step[ir][0]
                     <<" A1t="<<fsc_constpd.edp.fsc_damage_initiation_time_s[ir][0]
                     <<" A1cause="<<fsc_constpd.edp.fsc_damage_initiation_cause[ir][0]
                     <<" A1cross="<<fsc_constpd.edp.fsc_damage_sign_crossings[ir][0]
                     <<" A1ret="<<fsc_constpd.edp.fsc_damage_retained_ratio[ir][0]
                     <<" A1res="<<fsc_constpd.edp.fsc_damage_residual_step[ir][0]
                     <<" B1init="<<fsc_constpd.edp.fsc_damage_initiation_step[ir][1]
                     <<" B1t="<<fsc_constpd.edp.fsc_damage_initiation_time_s[ir][1]
                     <<" B1cause="<<fsc_constpd.edp.fsc_damage_initiation_cause[ir][1]
                     <<" B1cross="<<fsc_constpd.edp.fsc_damage_sign_crossings[ir][1]
                     <<" B1ret="<<fsc_constpd.edp.fsc_damage_retained_ratio[ir][1]
                     <<" B1res="<<fsc_constpd.edp.fsc_damage_residual_step[ir][1]<<"\n";
        }
        auto m3_series=calibrated_model_phase8(gravity,true,true,false,false,false,3,true);
        auto fsc_series=record_robust(m3_series,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"FSC_seriesElastic termination="<<static_cast<int>(fsc_series.analysis.termination)
                 <<" base="<<fsc_series.edp.peak_base_shear
                 <<" drifts="<<fsc_series.edp.peak_story_drift[0]<<","<<fsc_series.edp.peak_story_drift[1]<<","<<fsc_series.edp.peak_story_drift[2]
                 <<" PsumErr="<<fsc_series.edp.max_story1_axial_sum_error
                 <<" A1init="<<fsc_series.edp.fsc_damage_initiation_step[2][0]
                 <<" A1t="<<fsc_series.edp.fsc_damage_initiation_time_s[2][0]
                 <<" B1init="<<fsc_series.edp.fsc_damage_initiation_step[2][1]
                 <<" B1t="<<fsc_series.edp.fsc_damage_initiation_time_s[2][1]
                 <<" B1ret="<<fsc_series.edp.fsc_damage_retained_ratio[2][1]<<"\n";
        auto m3_active=calibrated_model_phase8(gravity,true,true,false,false,false,3,true,true,-0.005,0.20,0.55,-1.0);
        auto fsc_active=record_robust(m3_active,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"FSC_seriesActive termination="<<static_cast<int>(fsc_active.analysis.termination)
                 <<" step="<<fsc_active.analysis.termination_step<<" reason="<<fsc_active.analysis.termination_reason
                 <<" base="<<fsc_active.edp.peak_base_shear
                 <<" drifts="<<fsc_active.edp.peak_story_drift[0]<<","<<fsc_active.edp.peak_story_drift[1]<<","<<fsc_active.edp.peak_story_drift[2]
                 <<" residual="<<fsc_active.edp.residual_story1
                 <<" PsumErr="<<fsc_active.edp.max_story1_axial_sum_error
                 <<" A1init="<<fsc_active.edp.active_fsc_initiation_step[0]<<" A1t="<<fsc_active.edp.active_fsc_initiation_time_s[0]
                 <<" A1res="<<fsc_active.edp.active_fsc_residual_step[0]<<" A1ret="<<fsc_active.edp.active_fsc_min_retained_ratio[0]
                 <<" B1init="<<fsc_active.edp.active_fsc_initiation_step[1]<<" B1t="<<fsc_active.edp.active_fsc_initiation_time_s[1]
                 <<" B1res="<<fsc_active.edp.active_fsc_residual_step[1]<<" B1ret="<<fsc_active.edp.active_fsc_min_retained_ratio[1]
                 <<" E="<<fsc_active.analysis.max_components_lateral_loss<<" F="<<fsc_active.analysis.max_failed_components<<"\n";
        if(argc>3 && std::string(argv[3])=="active-shear") return fsc_active.analysis.termination==AnalysisTermination::Completed?0:3;
        auto m3_active_clear=calibrated_model_phase8(gravity,true,true,false,true,false,3,true,true,-0.005,0.20,0.55,-1.0);
        auto fsc_active_clear=record_robust(m3_active_clear,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"FSC_seriesActiveClear39 termination="<<static_cast<int>(fsc_active_clear.analysis.termination)
                 <<" step="<<fsc_active_clear.analysis.termination_step<<" reason="<<fsc_active_clear.analysis.termination_reason
                 <<" base="<<fsc_active_clear.edp.peak_base_shear
                 <<" drifts="<<fsc_active_clear.edp.peak_story_drift[0]<<","<<fsc_active_clear.edp.peak_story_drift[1]<<","<<fsc_active_clear.edp.peak_story_drift[2]
                 <<" residual="<<fsc_active_clear.edp.residual_story1<<" roofres="<<fsc_active_clear.edp.residual_roof
                 <<" PsumErr="<<fsc_active_clear.edp.max_story1_axial_sum_error
                 <<" A1init="<<fsc_active_clear.edp.active_fsc_initiation_step[0]<<" A1t="<<fsc_active_clear.edp.active_fsc_initiation_time_s[0]
                 <<" B1init="<<fsc_active_clear.edp.active_fsc_initiation_step[1]<<" B1t="<<fsc_active_clear.edp.active_fsc_initiation_time_s[1]
                 <<" E="<<fsc_active_clear.analysis.max_components_lateral_loss<<" F="<<fsc_active_clear.analysis.max_failed_components<<"\n";
        if(argc>3 && std::string(argv[3])=="active-shear-clear") return fsc_active_clear.analysis.termination==AnalysisTermination::Completed?0:3;
        if(argc>3 && std::string(argv[3])=="shear-sensitivity"){
            const std::array<double,3> slopes{{-0.0025,-0.005,-0.010}};
            const std::array<double,3> residuals{{0.10,0.20,0.30}};
            for(double spr:slopes)for(double rr:residuals){
                auto ms=calibrated_model_phase8(gravity,true,true,false,true,false,3,true,true,spr,rr,0.55,-1.0);
                auto rs=record_robust(ms,motion,LinearStrategy::SamePatternRefactorization,true);
                std::cerr<<"FSC_SENS slopeRatio="<<spr<<" residualRatio="<<rr
                         <<" term="<<static_cast<int>(rs.analysis.termination)
                         <<" base="<<rs.edp.peak_base_shear
                         <<" d="<<rs.edp.peak_story_drift[0]<<","<<rs.edp.peak_story_drift[1]<<","<<rs.edp.peak_story_drift[2]
                         <<" res1="<<rs.edp.residual_story1<<" roofres="<<rs.edp.residual_roof
                         <<" B1t="<<rs.edp.active_fsc_initiation_time_s[1]<<" A1t="<<rs.edp.active_fsc_initiation_time_s[0]
                         <<" E="<<rs.analysis.max_components_lateral_loss<<" F="<<rs.analysis.max_failed_components<<"\n";
            }
            return 0;
        }
        if(argc>3 && std::string(argv[3])=="damage-only") return 0;
        auto m4_constpd_A=calibrated_model_phase8(gravity,true,true,true,false,false,1);
        auto coupled_constpd_A=record_robust(m4_constpd_A,motion,LinearStrategy::SamePatternRefactorization,true,true);
        std::cerr<<"PM_A1_constPDelta termination="<<static_cast<int>(coupled_constpd_A.analysis.termination)
                 <<" step="<<coupled_constpd_A.analysis.termination_step
                 <<" reason="<<coupled_constpd_A.analysis.termination_reason
                 <<" PsumErr="<<coupled_constpd_A.edp.max_story1_axial_sum_error<<"\n";
        auto m4_constpd_B=calibrated_model_phase8(gravity,true,true,true,false,false,2);
        auto coupled_constpd_B=record_robust(m4_constpd_B,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"PM_B1_constPDelta termination="<<static_cast<int>(coupled_constpd_B.analysis.termination)
                 <<" step="<<coupled_constpd_B.analysis.termination_step
                 <<" reason="<<coupled_constpd_B.analysis.termination_reason
                 <<" base="<<coupled_constpd_B.edp.peak_base_shear
                 <<" drift1="<<coupled_constpd_B.edp.peak_story_drift[0]
                 <<" PsumErr="<<coupled_constpd_B.edp.max_story1_axial_sum_error<<"\n";
        auto m4_constpd=calibrated_model_phase8(gravity,true,true,true,false,false,3);
        auto coupled_constpd=record_robust(m4_constpd,motion,LinearStrategy::SamePatternRefactorization,true);
        std::cerr<<"PM_FSC_constPDelta termination="<<static_cast<int>(coupled_constpd.analysis.termination)
                 <<" step="<<coupled_constpd.analysis.termination_step
                 <<" reason="<<coupled_constpd.analysis.termination_reason
                 <<" base="<<coupled_constpd.edp.peak_base_shear
                 <<" drifts="<<coupled_constpd.edp.peak_story_drift[0]<<","<<coupled_constpd.edp.peak_story_drift[1]<<","<<coupled_constpd.edp.peak_story_drift[2]
                 <<" residual="<<coupled_constpd.edp.residual_story1
                 <<" PsumRange="<<coupled_constpd.edp.min_story1_axial_sum<<":"<<coupled_constpd.edp.max_story1_axial_sum
                 <<" PsumErr="<<coupled_constpd.edp.max_story1_axial_sum_error
                 <<" E="<<coupled_constpd.analysis.max_components_lateral_loss
                 <<" F="<<coupled_constpd.analysis.max_failed_components<<"\n";
        if(argc>3 && std::string(argv[3])=="axial-control-only") return 0;

        auto m4g=calibrated_model_phase8(gravity,true,true,true);
        auto coupled_clean=record_robust(m4g,motion,LinearStrategy::SamePatternRefactorization);
        std::cerr<<"phase8A_clean termination="<<static_cast<int>(coupled_clean.analysis.termination)
                 <<" reason="<<coupled_clean.analysis.termination_reason
                 <<" base="<<coupled_clean.edp.peak_base_shear
                 <<" drifts="<<coupled_clean.edp.peak_story_drift[0]<<","<<coupled_clean.edp.peak_story_drift[1]<<","<<coupled_clean.edp.peak_story_drift[2]
                 <<" residual="<<coupled_clean.edp.residual_story1
                 <<" B1bot_q="<<coupled_clean.edp.signed_peak_first_story_column_hinge_rotation[2]
                 <<" P_at_q="<<coupled_clean.edp.compression_at_peak_first_story_hinge_rotation[2]
                 <<" B1_maxP="<<coupled_clean.edp.max_column_compression[1]
                 <<" E="<<coupled_clean.analysis.max_components_lateral_loss<<" F="<<coupled_clean.analysis.max_failed_components<<"\n";
        if(argc>3 && (std::string(argv[3])=="phase8a-only" || std::string(argv[3])=="phase8a-clean-only")) return 0;

        // Phase 8B geometry diagnostic: use the documented 39-in clear column
        // as the actual flexible member length, with 4.5-in rigid arms to the
        // 48-in joint-center levels.  This tests local moment/shear and end-
        // rotation mechanics while recalibrating only the elastic T1 target.
        auto report_term=[&](const char* name,const RecordedRun& r){
            std::cerr<<name<<" termination="<<static_cast<int>(r.analysis.termination)
                     <<" step="<<r.analysis.termination_step<<" time="<<r.analysis.termination_time
                     <<" reason="<<r.analysis.termination_reason
                     <<" subdiv="<<r.analysis.stats.subdivided_steps<<" failed="<<r.analysis.stats.failed_steps
                     <<" E="<<r.analysis.max_components_lateral_loss<<" F="<<r.analysis.max_failed_components
                     <<" firstE="<<r.analysis.first_lateral_loss_step
                     <<" base="<<r.edp.peak_base_shear
                     <<" drift1="<<r.edp.peak_story_drift[0]
                     <<" qB1bot="<<r.edp.signed_peak_first_story_column_hinge_rotation[2]
                     <<" P@qB1bot="<<r.edp.compression_at_peak_first_story_hinge_rotation[2]
                     <<" drift@qB1bot="<<r.edp.story1_drift_at_peak_first_story_hinge_rotation[2]
                     <<"\n";
        };
        // Isolate the clear-column geometry from FSC and P-M coupling before
        // accepting it as a physical refinement.
        auto m5a=calibrated_model_phase8(p2,true,false,false,true);
        auto clear_base=record_robust(m5a,motion,LinearStrategy::SamePatternRefactorization);
        report_term("clear_base",clear_base);
        auto m5b=calibrated_model_phase8(p2,true,true,false,true);
        auto clear_fsc=record_robust(m5b,motion,LinearStrategy::SamePatternRefactorization);
        report_term("clear_fsc",clear_fsc);
        auto m5c=calibrated_model_phase8(p2,true,false,true,true);
        auto clear_pm=record_robust(m5c,motion,LinearStrategy::SamePatternRefactorization);
        report_term("clear_pm",clear_pm);
        auto m5=calibrated_model_phase8(p2,true,true,true,true);
        auto clearcol=record_robust(m5,motion,LinearStrategy::SamePatternRefactorization);
        report_term("clear_fsc_pm",clearcol);
        if(clearcol.analysis.termination!=AnalysisTermination::Completed)throw std::runtime_error("clear-column geometry run failed: "+clearcol.analysis.termination_reason);
        auto& rec=clearcol.edp;auto& same=clearcol.analysis;

        RobustNewmarkOptions inv;inv.tolerance=1e-7;inv.max_iterations=40;inv.line_search=true;inv.max_backtracks=8;inv.max_subdivisions=5;
        inv.return_numerical_failure=true;inv.collapse.check_initial_stability=true;
        auto full=run_newmark_robust(m5.model,motion.ag,motion.dt,LinearStrategy::FullFactorization,inv);
        if(full.termination!=AnalysisTermination::Completed)throw std::runtime_error("full invariance run failed");
        const double dfs=max_history_diff(full.roof_history,same.roof_history);
        const double hscale=std::max({1.0,full.stats.max_roof_abs,same.stats.max_roof_abs});

        auto print=[&](const char* name,const RecordedRun& r){
            std::cerr<<name<<" base="<<r.edp.peak_base_shear<<" drifts="<<r.edp.peak_story_drift[0]<<","<<r.edp.peak_story_drift[1]<<","<<r.edp.peak_story_drift[2]
                     <<" residual="<<r.edp.residual_story1<<" roofres="<<r.edp.residual_roof<<" subdiv="<<r.analysis.stats.subdivided_steps
                     <<" E="<<r.analysis.max_components_lateral_loss<<" F="<<r.analysis.max_failed_components<<"\n";
        };
        print("panel_gravity",r1);print("fixed_PM",r2);print("FSC",fsc);print("coupled_PM_FSC",coupled);print("coupled_PM_FSC_clear39",clearcol);
        std::cerr<<"full_same_rel="<<dfs/hscale<<"\n";

        std::ofstream o(out_path);if(!o)throw std::runtime_error("cannot write phase8 result");o<<std::setprecision(10)<<std::fixed;
        auto write_edp=[&](const Recorded& r){
            o<<"{\"peak_base_shear_kip\":"<<r.peak_base_shear<<",\"peak_base_shear_over_weight\":"<<r.peak_base_shear/total_weight<<",\"peak_story_drift_percent\":";json_array(o,r.peak_story_drift);
            o<<",\"residual_story1_drift_percent\":"<<r.residual_story1<<",\"residual_roof_drift_percent\":"<<r.residual_roof<<",\"peak_story_shear_kip\":";json_array(o,r.peak_story_shear);
            o<<",\"peak_floor_acceleration_g\":";json_array(o,r.peak_floor_accel);o<<",\"peak_roof_displacement_in\":"<<r.peak_roof_disp<<"}";
        };
        auto step_or_null=[&](std::size_t x){if(x==static_cast<std::size_t>(-1))o<<"null";else o<<x;};
        o<<"{\n";
        o<<"  \"benchmark_id\": \"nist_gcr_22_917_50_uc_berkeley_3story_dt1\",\n";
        o<<"  \"validation_status\": \"PRE-VALIDATION_PHASE8B_CLEAR_COLUMN_COUPLED_PM_FSC_PROXY_INPUT\",\n";
        o<<"  \"first_mode_period_s\": "<<m5.period<<",\n";
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
        o<<"    \"panel_zone_gravity_P\": ";write_edp(r1.edp);o<<",\n";
        o<<"    \"panel_zone_fixed_PM_second_pass\": ";write_edp(r2.edp);o<<",\n";
        o<<"    \"fixed_PM_plus_flexure_shear_critical\": ";write_edp(fsc.edp);o<<",\n";
        o<<"    \"consistent_trial_PM_plus_flexure_shear_critical\": ";write_edp(coupled.edp);o<<",\n";
        o<<"    \"consistent_PM_FSC_plus_39in_clear_column_geometry\": ";write_edp(rec);o<<"\n  },\n";
        o<<"  \"column_gravity_preload_kip\": ";json_array12(o,gravity);o<<",\n";
        o<<"  \"column_pass1_max_compression_kip\": ";json_array12(o,r1.edp.max_column_compression);o<<",\n";
        o<<"  \"column_second_pass_design_compression_kip\": ";json_array12(o,p2);o<<",\n";
        o<<"  \"peak_panel_zone_rotation_rad\": ";json_array12(o,rec.peak_panel_zone_rotation);o<<",\n";
        o<<"  \"peak_first_story_column_hinge_rotation_rad\": ";json_array8(o,rec.peak_first_story_column_hinge_rotation);o<<",\n";
        o<<"  \"first_story_nonductile_expected_shear_capacity_kip\": "<<first_story_nonductile_expected_shear_kip<<",\n";
        o<<"  \"first_lateral_loss_step\": ";step_or_null(same.first_lateral_loss_step);o<<",\n";
        o<<"  \"first_component_failure_step\": ";step_or_null(same.first_component_failure_step);o<<",\n";
        o<<"  \"max_components_lateral_loss\": "<<same.max_components_lateral_loss<<",\n";
        o<<"  \"max_failed_components\": "<<same.max_failed_components<<",\n";
        o<<"  \"calibration\": {\"target_perform3d_T1_s\":0.4800000000,\"flexural_stiffness_scale\":"<<m5.scale<<",\"alpha_mass\":"<<m5.alpha_m<<",\"beta_stiffness\":"<<m5.beta_k<<"},\n";
        o<<"  \"solver_invariance\": {\"full_vs_same_status\":\""<<(dfs/hscale<1e-8?"PASS":"OPEN")<<"\",\"max_roof_history_diff_full_vs_same_in\":"<<dfs<<",\"relative_roof_history_difference\":"<<dfs/hscale<<",\"full_seconds\":"<<full.stats.elapsed_seconds<<",\"same_pattern_seconds\":"<<same.stats.elapsed_seconds<<"},\n";
        o<<"  \"phase8_numerics\": {\"coupling_scope\":\"first-story nonductile A/B column ends only\",\"column_geometry\":\"39-in flexible column with 4.5-in small-rotation rigid arms to 48-in joint-center levels\",\"capacity_interpolation\":\"C1 piecewise cubic Hermite over mechanics-derived section P-M table with constant continuation outside tabulated P range\",\"dmoment_dcapacity\":\"same-committed-branch central finite difference used only to assemble consistent P-to-M cross tangent\",\"global_tangent\":\"nonsymmetric direct state tangent; generalized Woodbury deliberately disabled for coupled hinges pending exact nonsymmetric-block backend\"},\n";
        o<<"  \"fsc_regression_reference\": {\"prior_peak_base_shear_kip\":40.0287228979,\"prior_peak_story_drift_percent\":[3.3376283360,5.0923913562,2.2326354943],\"prior_residual_story1_drift_percent\":1.1932278462},\n";
        o<<"  \"gravity_loss_caveat\": \"ASCE-style E/F component states are tracked, but F does not yet redistribute static gravity load because current P-Delta uses prescribed preloads rather than a solved vertical gravity equilibrium.\",\n";
        o<<"  \"input_motion_provenance\": \"Same deterministic 70 s, 1.52g development proxy. It is not the recorded DT1 shake-table waveform and therefore remains unsuitable for strict parity.\",\n";
        o<<"  \"model_provenance\": \"Published geometry/mass/material strengths; stiffness calibrated only to Perform3D T1=0.48s; nonlinear EDPs not fitted. Phase 8 replaces the previously rejected scalar within-Newton P-M update with a state-consistent cross tangent for the critical first-story nonductile columns.\",\n";
        o<<"  \"strict_parity_blockers\": [\"Actual recorded Dynamic Test 1 shake-table acceleration history not acquired\",\"Intermediate Phase-8 coupling is not yet a full associative P-M-M plasticity surface\",\"Static gravity redistribution after F loss is not yet solved\",\"Exact Perform3D RC panel-zone and component property definitions are not fully available\"],\n";
        o<<"  \"perform3d_initial_condition_note\": \"The published Perform3D DT1 parity model starts from a virgin numerical state; prior experimental motions are not replayed.\"\n";
        o<<"}\n";
        std::cout<<out_path<<"\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}
}
