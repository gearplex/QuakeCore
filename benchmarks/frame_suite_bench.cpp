#include "quake/frame2d.hpp"
#include "quake/newmark.hpp"
#include "quake/suite.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread>
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

int main(int argc,char** argv) {
    const int stories=argc>1?std::stoi(argv[1]):20;
    const int bays=argc>2?std::stoi(argv[2]):2;
    const int records=argc>3?std::stoi(argv[3]):12;
    const int steps=argc>4?std::stoi(argv[4]):600;
    constexpr double dt=0.005;
    auto model=make_frame(stories,bays);

    std::vector<std::vector<double>> motions;
    motions.reserve(static_cast<std::size_t>(records));
    for(int r=0;r<records;++r) {
        const double amp=1.10+0.30*static_cast<double>(r)/static_cast<double>(std::max(1,records-1));
        auto gm=synthetic_ground_motion(steps,dt,amp);
        if(r%2) for(std::size_t i=0;i<gm.size();++i) gm[i] *= (i%7<4 ? 1.0 : -0.92);
        motions.push_back(std::move(gm));
    }

    double full_seconds=0.0;
    std::size_t full_factors=0;
    for(const auto& gm:motions) {
        auto r=run_newmark(model,gm,dt,LinearStrategy::FullFactorization,1e-8,35);
        full_seconds += r.stats.elapsed_seconds;
        full_factors += r.stats.global_factorizations;
    }

    PreparedWoodburyNewmark prepared(model,dt);
    double wood_run_seconds=0.0;
    std::size_t wood_iters=0;
    std::size_t max_rank=0;
    for(const auto& gm:motions) {
        auto r=prepared.run(gm,1e-8,35);
        wood_run_seconds += r.stats.elapsed_seconds;
        wood_iters += r.stats.newton_iterations;
        max_rank = std::max(max_rank,r.stats.max_active_rank);
    }
    const double wood_total=prepared.setup_seconds()+wood_run_seconds;

    const int hw=std::max(1,std::min(records,static_cast<int>(std::thread::hardware_concurrency()?std::thread::hardware_concurrency():4)));
    const int parallel_workers=std::min(4,hw);
    auto parallel=run_record_suite_parallel(model,motions,dt,parallel_workers,1e-8,35);

    std::cout<<std::fixed<<std::setprecision(6);
    std::cout<<"stories="<<stories<<" bays="<<bays<<" dof="<<model.dof()
             <<" nonlinear_springs="<<model.nonlinear_count()<<" records="<<records
             <<" steps_per_record="<<steps<<"\n";
    std::cout<<"full_suite_seconds="<<full_seconds<<" full_factorizations="<<full_factors<<"\n";
    std::cout<<"woodbury_setup_seconds="<<prepared.setup_seconds()
             <<" woodbury_run_seconds="<<wood_run_seconds
             <<" woodbury_total_seconds="<<wood_total<<"\n";
    std::cout<<"prepared_woodbury_speedup_vs_full="<<full_seconds/wood_total<<"x\n";
    std::cout<<"woodbury_setup_fraction="<<prepared.setup_seconds()/wood_total
             <<" max_active_rank="<<max_rank<<" woodbury_newton_iters="<<wood_iters<<"\n";
    std::cout<<"parallel_workers="<<parallel.workers<<" parallel_wall_seconds="<<parallel.elapsed_seconds
             <<" parallel_speedup_vs_sequential_woodbury="<<wood_total/parallel.elapsed_seconds<<"x"
             <<" parallel_speedup_vs_full="<<full_seconds/parallel.elapsed_seconds<<"x\n";
}
