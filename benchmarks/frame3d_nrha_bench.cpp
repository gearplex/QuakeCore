#include "quake/frame3d.hpp"
#include "quake/newmark.hpp"

#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;

static int nid(int s,int ix,int iy){return 10000*s+100*ix+iy;}

static CompiledFrame3D make_space_frame(int stories,int nx,int ny){
    Frame3DBuilder b;
    constexpr double h=3.5,L=6.0,E=25000.0,G=10000.0;
    const int ncols=(nx+1)*(ny+1);
    for(int s=0;s<=stories;++s){
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy){
            const int id=nid(s,ix,iy);const double m=s==0?0.0:1.0/static_cast<double>(ncols);
            b.add_node(id,ix*L,iy*L,s*h,m,m,0.05*m,0,0,0);if(s==0)b.fix(id);
        }
    }
    int eid=1,sid=1,did=500000;
    for(int s=1;s<=stories;++s){
        const double p=6.0*(stories-s+1);
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy)
            b.add_elastic_frame(eid++,nid(s-1,ix,iy),nid(s,ix,iy),E,G,2.0,0.30,0.55,0.55,1,0,0,p);
        for(int ix=0;ix<nx;++ix)for(int iy=0;iy<=ny;++iy){
            const int jl=nid(s,ix,iy),jr=nid(s,ix+1,iy),dl=did++,dr=did++;
            b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,(ix+1)*L,iy*L,s*h);
            for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RX,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}
            b.add_elastic_frame(eid++,dl,dr,E,G,1.5,0.20,0.35,0.35,0,1,0);
            b.add_bilinear_spring(sid++,jl,Dof3D::RY,dl,Dof3D::RY,50000,8.0,0.50);
            b.add_bilinear_spring(sid++,jr,Dof3D::RY,dr,Dof3D::RY,50000,8.0,0.50);
        }
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<ny;++iy){
            const int jl=nid(s,ix,iy),jr=nid(s,ix,iy+1),dl=did++,dr=did++;
            b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,ix*L,(iy+1)*L,s*h);
            for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RY,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}
            b.add_elastic_frame(eid++,dl,dr,E,G,1.5,0.20,0.35,0.35,1,0,0);
            b.add_bilinear_spring(sid++,jl,Dof3D::RX,dl,Dof3D::RX,50000,8.0,0.50);
            b.add_bilinear_spring(sid++,jr,Dof3D::RX,dr,Dof3D::RX,50000,8.0,0.50);
        }
    }
    b.set_rayleigh(0.02,0.0005);b.set_ground_direction(Dof3D::UX);b.set_response(nid(stories,0,0),Dof3D::UX);
    std::vector<int> story_nodes;for(int s=1;s<=stories;++s)story_nodes.push_back(nid(s,0,0));b.set_story_nodes(std::move(story_nodes),Dof3D::UX);
    return b.compile();
}

static double avg_rank(const AnalysisStats& s){return s.reduced_update_solves?static_cast<double>(s.active_rank_sum)/s.reduced_update_solves:0.0;}

int main(int argc,char** argv){
    const int stories=argc>1?std::stoi(argv[1]):10;const int nx=argc>2?std::stoi(argv[2]):1;const int ny=argc>3?std::stoi(argv[3]):1;const int steps=argc>4?std::stoi(argv[4]):400;const double amp=argc>5?std::stod(argv[5]):1.25;const double dt=0.005;
    auto model=make_space_frame(stories,nx,ny);auto gm=synthetic_ground_motion(steps,dt,amp);
    auto full=run_newmark(model,gm,dt,LinearStrategy::FullFactorization,1e-8,35);
    auto same=run_newmark(model,gm,dt,LinearStrategy::SamePatternRefactorization,1e-8,35);
    auto wood=run_newmark(model,gm,dt,LinearStrategy::Woodbury,1e-8,35);
    std::cout<<std::fixed<<std::setprecision(6);
    std::cout<<"stories="<<stories<<" bays_x="<<nx<<" bays_y="<<ny<<" dof="<<model.dof()<<" elastic_elements="<<model.elastic_element_count()<<" nonlinear_springs="<<model.nonlinear_count()<<" basis_nnz="<<model.nonlinear_basis().nnz()<<" steps="<<steps<<"\n";
    std::cout<<"full_seconds="<<full.stats.elapsed_seconds<<" full_factorizations="<<full.stats.global_factorizations<<" full_newton_iters="<<full.stats.newton_iterations<<"\n";
    std::cout<<"same_pattern_seconds="<<same.stats.elapsed_seconds<<" same_pattern_factorizations="<<same.stats.global_factorizations<<" same_pattern_newton_iters="<<same.stats.newton_iterations<<"\n";
    std::cout<<"woodbury_seconds="<<wood.stats.elapsed_seconds<<" woodbury_factorizations="<<wood.stats.global_factorizations<<" woodbury_newton_iters="<<wood.stats.newton_iterations<<"\n";
    std::cout<<"same_pattern_speedup_vs_full="<<full.stats.elapsed_seconds/same.stats.elapsed_seconds<<"x\n";
    std::cout<<"woodbury_speedup_vs_full="<<full.stats.elapsed_seconds/wood.stats.elapsed_seconds<<"x woodbury_speedup_vs_same_pattern="<<same.stats.elapsed_seconds/wood.stats.elapsed_seconds<<"x\n";
    std::cout<<"woodbury_max_active_rank="<<wood.stats.max_active_rank<<" woodbury_mean_active_rank="<<avg_rank(wood.stats)<<" cached_low_rank_columns="<<wood.stats.cached_low_rank_columns<<"\n";
    const double eval_fraction=wood.stats.nonlinear_component_evaluations?static_cast<double>(wood.stats.nonlinear_active_component_evaluations)/wood.stats.nonlinear_component_evaluations:0.0;
    std::cout<<"nonlinear_evaluations="<<wood.stats.nonlinear_component_evaluations<<" active_nonlinear_evaluations="<<wood.stats.nonlinear_active_component_evaluations<<" active_evaluation_fraction="<<eval_fraction<<"\n";
    std::cout<<"full_peak_roof="<<full.stats.max_roof_abs<<" same_pattern_peak_roof="<<same.stats.max_roof_abs<<" woodbury_peak_roof="<<wood.stats.max_roof_abs<<"\n";
}
