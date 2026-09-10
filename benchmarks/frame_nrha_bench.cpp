#include "quake/frame2d.hpp"
#include "quake/newmark.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

using namespace quake;

static CompiledFrame2D make_frame(int stories, int bays) {
    Frame2DBuilder b;
    constexpr double h=3.5, bay=6.0;
    constexpr double E=25000.0, Ac=2.0, Ic=0.55, Ab=1.5, Ib=0.35;
    const int cols=bays+1;
    for(int s=0;s<=stories;++s) {
        for(int c=0;c<cols;++c) {
            const int id=1000*s+c;
            const double mx=s==0?0.0:1.0/static_cast<double>(cols);
            b.add_node(id,c*bay,s*h,mx,0.0,0.0);
            if(s==0) b.fix(id);
        }
        if(s>0) {
            std::vector<int> slaves;
            for(int c=1;c<cols;++c) slaves.push_back(1000*s+c);
            b.rigid_floor_x(1000*s,slaves);
        }
    }
    for(int s=1;s<=stories;++s) {
        const double p=18.0*(stories-s+1);
        for(int c=0;c<cols;++c)
            b.add_elastic_frame(100000+100*s+c,1000*(s-1)+c,1000*s+c,E,Ac,Ic,p);
        for(int bay_i=0;bay_i<bays;++bay_i) {
            const int dl=200000+1000*s+2*bay_i;
            const int dr=dl+1;
            const int jl=1000*s+bay_i, jr=jl+1;
            b.add_node(dl,bay_i*bay,s*h); b.add_node(dr,(bay_i+1)*bay,s*h);
            b.equal_dof(jl,dl,Dof2D::UX); b.equal_dof(jl,dl,Dof2D::UY);
            b.equal_dof(jr,dr,Dof2D::UX); b.equal_dof(jr,dr,Dof2D::UY);
            b.add_elastic_frame(300000+1000*s+bay_i,dl,dr,E,Ab,Ib);
            b.add_rotational_spring(400000+2000*s+2*bay_i,jl,dl,50000.0,8.0,0.50);
            b.add_rotational_spring(400000+2000*s+2*bay_i+1,jr,dr,50000.0,8.0,0.50);
        }
    }
    b.set_rayleigh(0.02,0.0005);
    b.set_response_node(1000*stories);
    std::vector<int> story_nodes;
    for(int s=1;s<=stories;++s) story_nodes.push_back(1000*s);
    b.set_story_nodes(std::move(story_nodes));
    return b.compile();
}

static double avg_rank(const AnalysisStats& s) {
    return s.reduced_update_solves ? static_cast<double>(s.active_rank_sum)/static_cast<double>(s.reduced_update_solves) : 0.0;
}

int main(int argc,char** argv) {
    const int stories=argc>1?std::stoi(argv[1]):20;
    const int bays=argc>2?std::stoi(argv[2]):2;
    const int steps=argc>3?std::stoi(argv[3]):800;
    const double dt=0.005;
    auto model=make_frame(stories,bays);
    auto gm=synthetic_ground_motion(steps,dt,1.5);

    auto full=run_newmark(model,gm,dt,LinearStrategy::FullFactorization,1e-8,35);
    auto modified=run_newmark(model,gm,dt,LinearStrategy::ModifiedNewton,1e-8,50);
    auto wood=run_newmark(model,gm,dt,LinearStrategy::Woodbury,1e-8,35);

    std::cout<<std::fixed<<std::setprecision(6);
    std::cout<<"stories="<<stories<<" bays="<<bays<<" dof="<<model.dof()
             <<" elastic_elements="<<model.elastic_element_count()
             <<" nonlinear_springs="<<model.nonlinear_count()<<" steps="<<steps<<"\n";
    std::cout<<"full_seconds="<<full.stats.elapsed_seconds<<" full_factorizations="<<full.stats.global_factorizations
             <<" full_newton_iters="<<full.stats.newton_iterations<<"\n";
    std::cout<<"modified_seconds="<<modified.stats.elapsed_seconds<<" modified_factorizations="<<modified.stats.global_factorizations
             <<" modified_newton_iters="<<modified.stats.newton_iterations<<"\n";
    std::cout<<"woodbury_seconds="<<wood.stats.elapsed_seconds<<" woodbury_factorizations="<<wood.stats.global_factorizations
             <<" woodbury_newton_iters="<<wood.stats.newton_iterations<<"\n";
    std::cout<<"woodbury_speedup_vs_full="<<full.stats.elapsed_seconds/wood.stats.elapsed_seconds<<"x\n";
    std::cout<<"modified_speedup_vs_full="<<full.stats.elapsed_seconds/modified.stats.elapsed_seconds<<"x\n";
    std::cout<<"woodbury_max_active_rank="<<wood.stats.max_active_rank
             <<" woodbury_mean_active_rank="<<avg_rank(wood.stats)
             <<" cached_low_rank_columns="<<wood.stats.cached_low_rank_columns<<"\n";
    std::cout<<"full_peak_roof="<<full.stats.max_roof_abs<<" woodbury_peak_roof="<<wood.stats.max_roof_abs<<"\n";
    std::cout<<"full_peak_story_drift="<<full.stats.max_story_drift_abs
             <<" woodbury_peak_story_drift="<<wood.stats.max_story_drift_abs<<"\n";
}
