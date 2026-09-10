#include "quake/frame2d.hpp"
#include "quake/newmark.hpp"
#include "quake/steel2d.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace quake;

namespace {
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double x,double y,double tol,const char* why){
    check(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tol*std::max({1.0,std::abs(x),std::abs(y)}),why);
}
template<class F>void rejects(F f){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,"invalid steel component was accepted");}

SteelMember2DProperties steel(double k=2.0e8,double fy=2.5e5,double b=.02){
    SteelMember2DProperties p;p.E=2.0e8;p.A=.02;p.I=8.0e-4;
    p.hinge_i=NonlinearMaterial(BilinearSpring(k,fy,b));
    p.hinge_j=NonlinearMaterial(BilinearSpring(k,fy,b));return p;
}

void tangent_fd(const SteelMember2D& e,std::array<double,6> u,const std::vector<double>& state,double tol){
    auto r=e.trial(u,state.data());
    for(int j=0;j<6;++j){
        const double h=j%3==2?2e-8:2e-9;auto a=u,b=u;a[j]+=h;b[j]-=h;
        auto ap=e.trial(a,state.data()),bm=e.trial(b,state.data());
        for(int i=0;i<6;++i){
            const double fd=(ap.force[i]-bm.force[i])/(2*h),exact=r.tangent[6*i+j];
            if(std::abs(fd-exact)>tol*std::max({1.0,std::abs(fd),std::abs(exact)})){
                std::cerr<<"steel tangent ("<<i<<','<<j<<") fd="<<fd<<" exact="<<exact<<'\n';
                throw std::runtime_error("steel member consistent tangent mismatch");
            }
        }
    }
}

CompiledFrame2D sdof(bool device,double c){
    Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,2,0,5.0);b.fix(1);b.fix(2,false,true,true);
    b.add_elastic_frame(1,1,2,1000.0,2.0,1.0);
    if(device){b.add_viscous_damper(2,1,2,{c,1.0,0.0});b.set_rayleigh(0.0,0.0);}else b.set_rayleigh(c/5.0,0.0);
    b.set_response_node(2);b.set_story_nodes({2});return b.compile();
}
}

int main(){try{
    rejects([]{auto p=steel();p.E=0;SteelMember2D e(0,0,3,0,p);});
    rejects([]{ViscousDamper2D d({10,.5,0});});
    rejects([]{ViscousDamper2D d({-1,1,0});});

    // Rigid translation and rigid rotation must not create member force.
    SteelMember2D e(1.2,-.5,4.0,2.3,steel());auto s=e.initial_state();
    for(const auto& u:std::vector<std::array<double,6>>{{.4,-.2,0,.4,-.2,0},
         {0,0,.003,-.0084,.0084,.003}}){
        auto r=e.trial(u,s.data());for(double f:r.force)near(f,0,2e-10,"steel rigid-body force");
    }
    auto K=e.initial_tangent();for(int i=0;i<6;++i)for(int j=0;j<6;++j)near(K[6*i+j],K[6*j+i],1e-13,"steel initial tangent symmetry");
    tangent_fd(e,{.001,-.002,.0008,.002,.001,-.0004},s,2e-6);

    // Nonlinear committed branch, rollback, and local condensed Jacobian.
    SteelMember2D vertical(0,0,0,3,steel(1.0e7,120.0,.015));auto sv=vertical.initial_state();
    auto yielded=vertical.trial({0,0,0,.04,0,0},sv.data());check(std::abs(yielded.hinge_rotation[0])>1e-6,"steel hinge did not yield");
    auto committed=yielded.state,copy=committed;
    (void)vertical.trial({0,0,0,-.02,0,0},committed.data());check(committed==copy,"steel trial mutated committed history");
    tangent_fd(vertical,{0,0,0,.032,0,.001},committed,3e-5);

    // A high-stiffness condensed hinge approaches the ordinary elastic frame.
    SteelMember2D stiff(0,0,0,3,steel(1e15,1e14,0));auto ks=stiff.initial_tangent();
    auto ke=frame2d_global_stiffness(0,0,0,3,2e8,.02,8e-4);
    for(int i=0;i<36;++i)near(ks[i],ke[i],2e-7,"condensed rigid hinge did not recover elastic frame");

    // Panel zones and BRBs retain distinct component classification and exact
    // generalized deformations after MPC condensation.
    {
        Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,0,0);b.add_node(3,4,3);
        b.fix(1);b.fix(2,true,true,false);b.fix(3,false,true,true);
        b.add_panel_zone(10,1,2,NonlinearMaterial(BilinearSpring(5000,50,.02)));
        b.add_brb(11,1,3,1000,20,.01);b.set_response_node(3);b.set_story_nodes({3});auto m=b.compile();
        check(m.nonlinear_component_kind(10)==NonlinearComponent2DKind::PanelZone,"panel-zone classification lost");
        check(m.nonlinear_component_kind(11)==NonlinearComponent2DKind::BRB,"BRB classification lost");
        auto st=m.initial_nonlinear_state();std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0);
        u[static_cast<std::size_t>(m.reduced_dof(2,Dof2D::RZ))]=.02;
        u[static_cast<std::size_t>(m.reduced_dof(3,Dof2D::UX))]=.01;
        auto pz=m.nonlinear_component_snapshot(10,u,st),brb=m.nonlinear_component_snapshot(11,u,st);
        near(pz.deformation,.02,1e-13,"panel-zone rotation");near(pz.force,51.0,1e-12,"panel-zone bilinear response");
        near(brb.deformation,.008,1e-13,"BRB projected axial deformation");near(brb.force,8,1e-13,"BRB elastic response");
    }

    // Power-law force and analytic velocity tangent, including the smooth
    // alpha<1 regularization used at velocity reversal.
    {
        ViscousDamper2D linear({125,1,0});auto r=linear.trial(-.4);near(r.force,-50,1e-14,"linear damper force");near(r.tangent,125,1e-14,"linear damper tangent");
        ViscousDamper2D nonlinear({80,.5,.002});for(double v:{-.4,-.01,0.0,.03,.5}){auto x=nonlinear.trial(v);const double h=1e-7;double fd=(nonlinear.trial(v+h).force-nonlinear.trial(v-h).force)/(2*h);near(x.tangent,fd,2e-8,"nonlinear damper tangent");}
    }

    // A linear axial device on a one-DOF oscillator is exactly equivalent to
    // the same damping coefficient supplied as proportional mass damping.
    {
        const double c=3.5;auto device=sdof(true,c),reference=sdof(false,c);
        std::vector<double> ag(500);for(int i=0;i<500;++i)ag[static_cast<std::size_t>(i)]=.8*std::sin(.031*i)+.15*std::sin(.19*i);
        RobustNewmarkOptions o;o.tolerance=1e-11;o.relative_force_tolerance=false;o.kinematic_initial_guess=false;
        auto a=run_newmark_robust(device,ag,.01,LinearStrategy::Woodbury,o),b=run_newmark_robust(reference,ag,.01,LinearStrategy::FullFactorization,o);
        check(a.termination==AnalysisTermination::Completed&&b.termination==AnalysisTermination::Completed,"linear damper NRHA failed");
        near(a.final_displacement[0],b.final_displacement[0],2e-11,"linear damper displacement equivalence");
        near(a.stats.damping_energy,b.stats.damping_energy,2e-11,"linear damper energy equivalence");
    }

    // Nonlinear dashpot exercises the velocity-dependent same-pattern path and
    // must remain identical between requested full and Woodbury strategies.
    {
        Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,2,0,5);b.fix(1);b.fix(2,false,true,true);
        b.add_elastic_frame(1,1,2,1000,2,1);b.add_viscous_damper(2,1,2,{12,.55,.001});b.set_response_node(2);b.set_story_nodes({2});auto m=b.compile();
        check(m.has_velocity_dependent_global_tangent(),"nonlinear damper tangent routing missing");
        std::vector<double> ag(250);for(int i=0;i<250;++i)ag[static_cast<std::size_t>(i)]=.7*std::sin(.09*i);
        RobustNewmarkOptions o;o.tolerance=1e-10;o.relative_force_tolerance=false;
        auto a=run_newmark_robust(m,ag,.01,LinearStrategy::FullFactorization,o),b2=run_newmark_robust(m,ag,.01,LinearStrategy::Woodbury,o);
        check(a.termination==AnalysisTermination::Completed&&b2.termination==AnalysisTermination::Completed,"nonlinear damper NRHA failed");
        near(a.final_displacement[0],b2.final_displacement[0],2e-11,"nonlinear damper solver mismatch");
    }
    std::cout<<"steel members, panel zones, BRBs, viscous dampers, rollback, Jacobians, and NRHA checks passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
