#include "quake/pm_interaction_hinge.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace quake;
PMInteractionHingeParams params(double L){
 PMInteractionHingeParams p;auto& h=p.hinge;
 h.Ke=5000*L;h.posFy=h.negFy=20*L;h.pos_a=h.neg_a=.02;h.pos_b=h.neg_b=.06;h.pos_f=h.neg_f=.09;h.backbone_shape=ASCE41BackboneShape::StraightCE;h.hardening_stiffness=200*L;
 p.axial_preload=20;p.axial_stiffness=5000/L;p.p_balance=45;p.py_tension=-30;p.py_compression=140;p.my_balance=25*L;
 p.alpha_tension=1.5;p.alpha_compression=1.7;p.beta_pm=1.1;p.surface_shape=PMInteractionSurfaceShape::PerformConcrete;
 p.surface_evolution=PMInteractionSurfaceEvolution::MrozTwoSurface;p.mroz_outer_scale=1.2;p.enforce_deformation_capacity=false;p.max_return_iterations=40;p.return_tolerance=1e-11;return p;
}
int main(){
 std::cout<<std::setprecision(14)<<"length_scale,step,rotation,moment_original_units,axial_plastic_original_units,converged\n";
 for(double L:{1.,25.4}){PMInteractionHinge2D h(params(L));std::vector<double> c(h.kStateSize),s(c.size());h.initialize_state(c.data());int i=0;
 for(double r:{.01,.02,.005,-.01,-.03,.01,.04}){auto t=h.trial(-.0005*L,r,c.data(),s.data());std::cout<<L<<','<<++i<<','<<r<<','<<t.moment/L<<','<<t.plastic_axial_deformation/L<<','<<t.converged<<'\n';c=s;}}
 for(auto evolution:{PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic,PMInteractionSurfaceEvolution::MrozTwoSurface}){
 auto p=params(1);p.surface_evolution=evolution;PMInteractionHinge2D h(p);std::vector<double> c(h.kStateSize),s(c.size());h.initialize_state(c.data());auto t=h.trial(-.05,0,c.data(),s.data());
 std::cerr<<"axial_intercept_test evolution="<<(int)evolution<<" compression="<<t.compression<<" intercept="<<p.py_compression<<" converged="<<t.converged<<'\n';}
 return 0;
}
