#include "quake/frame3d.hpp"
#include "quake/newmark.hpp"
#include <iomanip>
#include <iostream>
using namespace quake;

static CompiledFrame3D cantilever(bool corot){
    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,3.5,1.0,1.0,0.2);b.fix(1);
    b.add_elastic_frame(1,1,2,30000,12000,2.0,.3,.5,.5,1,0,0,25.0);
    b.set_rayleigh(.02,.0005);b.set_ground_direction(Dof3D::UX);b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);b.set_corotational(corot);return b.compile();
}
int main(int argc,char**argv){double amp=argc>1?std::stod(argv[1]):0.25;double dt=.01;auto c=cantilever(true),l=cantilever(false);auto gm=synthetic_ground_motion(120,dt,amp);RobustNewmarkOptions o;o.max_iterations=25;o.max_subdivisions=2;o.return_numerical_failure=true;o.collapse.tangent_check_interval=0;auto rc=run_newmark_robust(c,gm,dt,LinearStrategy::FullFactorization,o);auto rl=run_newmark_robust(l,gm,dt,LinearStrategy::FullFactorization,o);std::cout<<std::fixed<<std::setprecision(8)<<"amp="<<amp<<" corot_roof="<<rc.stats.max_roof_abs<<" linear_roof="<<rl.stats.max_roof_abs<<" ratio="<<(rl.stats.max_roof_abs?rc.stats.max_roof_abs/rl.stats.max_roof_abs:0)<<"\n"<<"corot_newton="<<rc.stats.newton_iterations<<" factorizations="<<rc.stats.global_factorizations<<" elapsed="<<rc.stats.elapsed_seconds<<" termination="<<(int)rc.termination<<"\n";}
