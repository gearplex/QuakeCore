#include "quake/frame2d.hpp"
#include "quake/soil2d.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace quake;
namespace {
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double x,double y,double tol,const char* why){check(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tol*std::max({1.0,std::abs(x),std::abs(y)}),why);}
template<class F>void rejects(F f){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,"invalid soil component was accepted");}
}

int main(){try{
    auto py=py_bilinear_backbone(120.0,1.5,0.03);near(py.nodal_capacity,180.0,1e-14,"p-y tributary force");near(py.initial_stiffness,3000.0,1e-14,"p-y y50 stiffness");
    auto tz=tz_bilinear_backbone(40.0,2.5,1.2,0.02);near(tz.nodal_capacity,120.0,1e-14,"t-z tributary force");near(tz.initial_stiffness,3000.0,1e-14,"t-z z50 stiffness");
    auto qz=qz_bilinear_backbone(500.0,0.8,0.01);near(qz.nodal_capacity,400.0,1e-14,"q-z area force");near(qz.initial_stiffness,20000.0,1e-14,"q-z z50 stiffness");
    auto trib=nodal_tributary_lengths({0,1,3,6});near(trib[0],.5,1e-14,"top tributary");near(trib[1],1.5,1e-14,"interior tributary");near(trib[2],2.5,1e-14,"interior tributary");near(trib[3],1.5,1e-14,"toe tributary");
    rejects([]{(void)py_bilinear_backbone(10,0,1);});rejects([]{(void)nodal_tributary_lengths({0,1,1});});

    // Compression-only return mapping: negative motion carries compression,
    // opening carries no tension and updates only the trial gap state.
    AsymmetricElasticPerfectlyPlasticSpring contact(100.0,0.0,20.0);
    auto c=contact.trial(-.3,0.0);near(c.force,-20,1e-14,"compression cap");near(c.tangent,0,1e-14,"compression cap tangent");near(c.plastic_deformation,-.1,1e-14,"compression plastic set");
    auto unload=contact.trial(-.2,c.plastic_deformation);near(unload.force,-10,1e-14,"compression unloading");
    auto opening=contact.trial(.05,0.0);near(opening.force,0,1e-14,"zero tensile resistance");near(opening.plastic_deformation,.05,1e-14,"opening gap state");

    // Coincident nodes, explicit direction, component classification, and
    // immutable committed history through the compiled low-rank path.
    Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,0,0);b.fix(1);b.fix(2,false,true,true);
    b.add_soil_spring(10,1,2,1,0,NonlinearMaterial(BilinearSpring(100,10,.02)));
    b.add_directional_viscous_damper(11,1,2,1,0,{25,1,0});b.set_response_node(2);b.set_story_nodes({2});auto m=b.compile();
    check(m.nonlinear_component_kind(10)==NonlinearComponent2DKind::SoilSpring,"soil classification lost");
    auto state=m.initial_nonlinear_state(),copy=state;std::vector<double> u{.2},f,t,trial;
    auto x=m.nonlinear_component_snapshot(10,u,state);near(x.deformation,.2,1e-14,"soil projected deformation");near(x.force,10.2,1e-13,"soil bilinear force");
    m.internal_force_and_tangent(u,state,f,t,trial);check(state==copy,"soil trial mutated committed state");
    auto d=m.viscous_damper_response(11,std::vector<double>{-.4});near(d.deformation_rate,-.4,1e-14,"soil dashpot rate");near(d.force,-10,1e-14,"soil dashpot force");

    {Frame2DBuilder r;r.add_node(1,0,0);r.add_node(2,0,0);r.fix(1);r.fix(2,false,true,false);r.add_rotational_soil_spring(20,1,2,NonlinearMaterial(BilinearSpring(500,5,0)));r.add_rotational_soil_dashpot(21,1,2,{30,1,0});r.add_soil_spring(22,1,2,1,0,NonlinearMaterial(BilinearSpring(100,10,0)));r.set_response_node(2);r.set_story_nodes({2});auto rm=r.compile();auto rs=rm.initial_nonlinear_state();auto rx=rm.nonlinear_component_snapshot(20,{0,.02},rs);near(rx.deformation,.02,1e-14,"rotational soil deformation");near(rx.force,5,1e-14,"rotational soil moment");auto rd=rm.viscous_damper_response(21,{0,-.3});near(rd.deformation_rate,-.3,1e-14,"rotational dashpot rate");near(rd.force,-9,1e-14,"rotational dashpot moment");}

    rejects([]{Frame2DBuilder x;x.add_node(1,0,0);x.add_node(2,1,0);x.fix(1);x.fix(2,false,true,true);x.add_soil_spring(1,1,2,1,0,NonlinearMaterial(BilinearSpring(1,1,0)));x.set_response_node(2);x.set_story_nodes({2});(void)x.compile();});

    // Long elastic pile on a uniform Winkler foundation. The finite-element
    // head displacement and rotation converge to the semi-infinite closed form
    // y(0)=P/(2 EI beta^3), theta(0)=-beta*y(0).
    {
        constexpr int ne=200;const double EI=1000.0,ks=400.0,beta=std::pow(ks/(4.0*EI),.25),L=10.0/beta,dz=L/ne,P=10.0;
        Frame2DBuilder x;std::vector<double> depths(static_cast<std::size_t>(ne+1));for(int i=0;i<=ne;++i)depths[static_cast<std::size_t>(i)]=i*dz;
        auto tl=nodal_tributary_lengths(depths);
        for(int i=0;i<=ne;++i){const int pile=1000+i,soil=2000+i;x.add_node(pile,0,20-depths[static_cast<std::size_t>(i)]);x.add_node(soil,0,20-depths[static_cast<std::size_t>(i)]);x.fix(pile,false,true,false);x.fix(soil);x.add_soil_spring(3000+i,soil,pile,1,0,NonlinearMaterial(BilinearSpring(ks*tl[static_cast<std::size_t>(i)],1e12,0)));if(i)x.add_elastic_frame(4000+i,999+i,1000+i,EI,1.0,1.0);}
        x.set_response_node(1000);x.set_story_nodes({1000});auto pile=x.compile();std::vector<double> rhs(static_cast<std::size_t>(pile.dof()),0.0);rhs[static_cast<std::size_t>(pile.reduced_dof(1000,Dof2D::UX))]=P;
        auto sol=superlu_solve_once(pile.K_initial(),rhs);const double exact=P/(2.0*EI*beta*beta*beta);
        const double uh=sol[static_cast<std::size_t>(pile.reduced_dof(1000,Dof2D::UX))],rh=sol[static_cast<std::size_t>(pile.reduced_dof(1000,Dof2D::RZ))];
        near(uh,exact,2e-4,"Winkler head displacement");near(rh,-beta*exact,2e-4,"Winkler head rotation");
        std::cout<<"winkler displacement relative error="<<std::abs(uh-exact)/exact<<", rotation relative error="<<std::abs(rh+beta*exact)/(beta*exact)<<'\n';
    }
    std::cout<<"soil adapters, contact gap, directional spring/dashpot, rollback, and contracts passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
