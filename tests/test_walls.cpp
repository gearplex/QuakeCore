#include "quake/wall2d.hpp"
#include "quake/frame2d.hpp"
#include "quake/newmark.hpp"
#include "quake/stability.hpp"
#include "quake/reduction.hpp"
#include "quake/superlu_solver.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
using namespace quake;
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double x,double y,double tol,const char* why){check(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tol*std::max({1.,std::abs(x),std::abs(y)}),why);}
template<class F>void rejects(F f){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,"invalid wall definition not rejected");}
MVLEMProperties mv(bool nonlinear=false){MVLEMProperties p;p.c=.4;p.shear=WallUniaxial::elastic(1e5);for(int i=0;i<4;++i)p.fibers.push_back({.5,.2,.02,nonlinear?WallUniaxial::concrete01(-30000,-.002,-6000,-.006):WallUniaxial::elastic(3e7),WallUniaxial::steel(2e8,4e5,.01)});return p;}
SFIMVLEMProperties sf(bool nonlinear=false){SFIMVLEMProperties p;p.c=.4;for(int i=0;i<4;++i){WallPanel panel(3e7,.2);if(nonlinear)panel=WallPanel(1e6,.2,{{.57,1,WallUniaxial::concrete01(-30000,-.002,-6000,-.006)},{.57+std::acos(-1.)/2,1,WallUniaxial::concrete01(-30000,-.002,-6000,-.006)},{0,.01,WallUniaxial::steel(2e8,4e5,.01)},{std::acos(-1.)/2,.02,WallUniaxial::steel(2e8,4e5,.01)}});p.panels.push_back({.5,.2,panel});}return p;}
void tangent(const Wall2D& e,std::array<double,6> u,const std::vector<double>& s,double tol){
 auto r=e.trial(u,s.data());for(int j=0;j<6;++j){double dx=j%3==2?1e-9:3e-9;auto a=u,b=u;a[j]+=dx;b[j]-=dx;auto ap=e.trial(a,s.data()),bm=e.trial(b,s.data());for(int i=0;i<6;++i){double fd=(ap.force[i]-bm.force[i])/(2*dx);if(std::abs(fd-r.tangent[6*i+j])>tol*std::max(1.,std::abs(fd))+1e-9*std::abs(*std::max_element(r.tangent.begin(),r.tangent.end()))) {std::cerr<<"tangent "<<i<<','<<j<<" expected "<<fd<<" got "<<r.tangent[6*i+j]<<'\n';throw std::runtime_error("wall tangent FD mismatch");}}}
}
int main(){try{
 rejects([]{Wall2D w(0,mv());});rejects([]{auto p=mv();p.fibers[0].width=-1;Wall2D w(3,p);});rejects([]{auto p=sf();p.local_max_iterations=0;Wall2D w(3,p);});rejects([]{WallUniaxial::concrete01(30,-.002,-5,-.006);});rejects([]{WallPanel p(-1,.2);});rejects([]{WallUniaxial::elastic(std::numeric_limits<double>::quiet_NaN());});
 {Frame2DBuilder b;rejects([&]{b.add_node(1,std::numeric_limits<double>::quiet_NaN(),0);});rejects([&]{b.add_node(1,0,0,std::numeric_limits<double>::infinity());});}
 for(auto e:{Wall2D(3,mv()),Wall2D(3,sf())}){
  auto s=e.initial_state();auto K=e.initial_tangent();
  for(auto u:{std::array<double,6>{1,0,0,1,0,0},std::array<double,6>{0,1,0,0,1,0},std::array<double,6>{0,0,.01,-.03,0,.01}}){auto r=e.trial(u,s.data());for(double v:r.force)near(v,0,1e-8,"rigid body force");}
  for(int a=0;a<6;++a)for(int b=0;b<6;++b)near(K[6*a+b],K[6*b+a],1e-12,"initial symmetry");
  tangent(e,{.001,-.002,.0001,.003,-.0021,-.0003},s,1e-5);
  auto old=s;auto rejected=e.trial({0,0,0,.025,-.003,-.005},s.data());(void)rejected;check(s==old,"trial mutated history");auto a=e.trial({0,0,0,.002,-.001,-.001},s.data()),b=e.trial({0,0,0,.002,-.001,-.001},old.data());check(a.force==b.force&&a.state==b.state,"rollback mismatch");
 }
 // Closed-form cantilever with discrete-fiber EI and finite shear stiffness.
 auto p=mv();Wall2D e(3,p);double EI=0,EA=0;for(int i=0;i<4;++i){double x=-.75+.5*i,E=.98*3e7+.02*2e8;EI+=.1*x*x*E;EA+=.1*E;}
 auto K=e.initial_tangent();std::vector<Triplet> t;for(int i=0;i<3;++i)for(int j=0;j<3;++j)t.push_back({i,j,K[(i+3)*6+j+3]});auto sol=superlu_solve_once(SparseMatrixCSC::from_triplets(3,3,t),{100,-300,0});near(sol[0],100/1e5+100*.6*.6*27/EI,1e-12,"cantilever shear/flexure compliance");near(sol[1],-300*3/EA,1e-12,"axial compliance");
 // Isotropic panel condensation at free transverse stress gives E_y=E, not E/(1-nu^2).
 auto pe=sf();Wall2D es(3,pe);auto rs=es.trial({0,0,0,0,-.0003,0},es.initial_state().data());for(auto x:rs.panel_strain)near(x[0],.00002,1e-12,"Poisson transverse strain");near(rs.force[4],-.0001*3e7*.4,1e-12,"condensed axial force");
 // Nonlinear tangent on committed unloading branches; the transverse solve is part of FD.
 for(auto wall:{Wall2D(3,mv(true)),Wall2D(3,sf(true))}){auto s=wall.initial_state();auto a=wall.trial({0,0,0,.002,-.009,-.0001},s.data());s=a.state;tangent(wall,{0,0,0,-.0003,-.0087,.00007},s,3e-4);}
 // Coupling: shear response changes with axial strain in the nonlinear RC panel.
 Wall2D coupled(3,sf(true));auto cs=coupled.initial_state();auto r1=coupled.trial({0,0,0,.006,-.002,0},cs.data()),r2=coupled.trial({0,0,0,.006,-.005,0},cs.data());check(std::abs(r1.force[3]-r2.force[3])>1,"missing axial/shear coupling");
 // Wall/frame assembly, lumped self mass, MPC condensation, direct solver routing.
 for(bool sfi:{false,true}){
  Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,0,3,50,20);b.fix(1);auto pm=mv();pm.density=2.5;auto ps=sf();ps.density=2.5;if(sfi)b.add_sfi_mvlem(10,1,2,ps);else b.add_mvlem(10,1,2,pm);b.set_rayleigh(.1,.001);b.set_story_nodes({2});auto model=b.compile();
  near(model.mass()[model.reduced_dof(2,Dof2D::UX)],51.5,1e-12,"wall nodal mass");check(model.wall_count()==1&&model.has_state_dependent_global_tangent(),"wall omitted from coupled assembly");
  std::vector<double> ag(200);for(int i=0;i<200;++i)ag[i]=.02*std::sin(i*.07);RobustNewmarkOptions o;o.kinematic_initial_guess=false;o.tolerance=1e-9;o.relative_force_tolerance=false;
  auto a=run_newmark_robust(model,ag,.005,LinearStrategy::FullFactorization,o),c=run_newmark_robust(model,ag,.005,LinearStrategy::Woodbury,o);check(a.termination==AnalysisTermination::Completed&&c.termination==AnalysisTermination::Completed,"wall NRHA failed");for(int i=0;i<model.dof();++i)near(a.final_displacement[i],c.final_displacement[i],1e-12,"wall direct fallback mismatch");
 }
 // Unit rescaling leaves panel stress and generalized force response invariant.
 auto pmm=mv(true);for(auto& f:pmm.fibers){f.width*=1000;f.thickness*=1000;f.concrete=WallUniaxial::concrete01(-.03,-.002,-.006,-.006);f.steel=WallUniaxial::steel(200, .4,.01);}pmm.shear=WallUniaxial::elastic(100);
 Wall2D wm(3,mv(true)),wmm(3000,pmm);auto sm=wm.initial_state(),smm=wmm.initial_state();std::array<double,6> um{0,0,0,.012,-.009,-.001},umm{0,0,0,12,-9,-.001};auto rm=wm.trial(um,sm.data()),rmm=wmm.trial(umm,smm.data());for(int i=0;i<6;++i)near(rm.force[i],rmm.force[i]/(i%3==2?1000:1),1e-11,"wall unit scaling");
 // Local integration failure must leave committed state unchanged and remain recoverable.
 {auto p=sf(true);p.local_max_iterations=1;Wall2D w(3,p);auto s=w.initial_state(),old=s;try{(void)w.trial({0,0,0,.03,-.007,.001},s.data());}catch(const ConstitutiveIntegrationError&){}check(s==old,"local panel trial mutated committed state");}
 // SFI unit scaling, including the local equilibrium tolerance and condensation.
 {auto p=sf(true);for(auto& f:p.panels){f.width*=1000;f.thickness*=1000;f.material=WallPanel(1.,.2,{{.57,1,WallUniaxial::concrete01(-.03,-.002,-.006,-.006)},{.57+std::acos(-1.)/2,1,WallUniaxial::concrete01(-.03,-.002,-.006,-.006)},{0,.01,WallUniaxial::steel(200,.4,.01)},{std::acos(-1.)/2,.02,WallUniaxial::steel(200,.4,.01)}});}Wall2D w(3,sf(true)),mm(3000,p);auto a=w.trial({0,0,0,.006,-.005,0},w.initial_state().data()),b=mm.trial({0,0,0,6,-5,0},mm.initial_state().data());for(int i=0;i<6;++i)near(a.force[i],b.force[i]/(i%3==2?1000:1),1e-10,"SFI unit scaling");}
 // Stability tracking must use committed panel history rather than virgin material.
 {Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,0,3,10);b.fix(1);b.fix(2,false,true,true);b.add_sfi_mvlem(1,1,2,sf(true));auto model=b.compile();auto s=model.initial_nonlinear_state();std::vector<double>f,t,z;model.internal_force_and_tangent({.02},s,f,t,z);s=z;model.internal_force_and_tangent({.018},s,f,t,z);auto K=model.effective_state_tangent_matrix_with_state({.018},t,s,0,0);auto se=estimate_tangent_stability(model,{.018},s);check(se.factorization_ok,"stateful stability failed");near(se.near_zero_eigenvalue,K.values()[0],1e-12,"stability dropped committed wall history");rejects([&]{ReducedDynamicModel reduced(model,{1.},1);});}
 std::cout<<"wall mechanics, nonlinear Jacobians, rollback, units, and NRHA checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
