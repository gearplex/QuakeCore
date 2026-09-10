#include "quake/newmark.hpp"
#include "quake/shear_building.hpp"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;
int main(int argc,char**argv){
    int stories=argc>1?std::atoi(argv[1]):120; int steps=argc>2?std::atoi(argv[2]):1200; int spacing=argc>3?std::atoi(argv[3]):6;
    std::vector<int> nl; for(int s=0;s<stories;s+=std::max(1,spacing))nl.push_back(s);
    ShearBuilding model(stories,1.0,240.0,nl,120.0,0.22,0.02,0.02,0.001);
    auto gm=synthetic_ground_motion(steps,0.01,5.5);
    auto full=run_newmark(model,gm,0.01,LinearStrategy::FullFactorization,1e-8,30);
    auto modified=run_newmark(model,gm,0.01,LinearStrategy::ModifiedNewton,1e-8,50);
    auto wood=run_newmark(model,gm,0.01,LinearStrategy::Woodbury,1e-8,30);
    std::cout<<std::fixed<<std::setprecision(6)
             <<"stories="<<stories<<" nonlinear_springs="<<nl.size()<<" steps="<<steps<<"\n"
             <<"full_seconds="<<full.stats.elapsed_seconds<<" full_factorizations="<<full.stats.global_factorizations<<" full_newton_iters="<<full.stats.newton_iterations<<"\n"
             <<"modified_seconds="<<modified.stats.elapsed_seconds<<" modified_factorizations="<<modified.stats.global_factorizations<<" modified_newton_iters="<<modified.stats.newton_iterations<<"\n"
             <<"woodbury_seconds="<<wood.stats.elapsed_seconds<<" woodbury_factorizations="<<wood.stats.global_factorizations<<" woodbury_newton_iters="<<wood.stats.newton_iterations<<"\n"
             <<"woodbury_speedup_vs_full="<<(full.stats.elapsed_seconds/wood.stats.elapsed_seconds)<<"x\n"
             <<"modified_speedup_vs_full="<<(full.stats.elapsed_seconds/modified.stats.elapsed_seconds)<<"x\n"
             <<"woodbury_max_active_rank="<<wood.stats.max_active_rank<<" woodbury_mean_active_rank="<<(wood.stats.reduced_update_solves?double(wood.stats.active_rank_sum)/double(wood.stats.reduced_update_solves):0.0)<<"\n"
             <<"full_peak_roof="<<full.stats.max_roof_abs<<" woodbury_peak_roof="<<wood.stats.max_roof_abs<<"\n"
             <<"full_peak_story_drift="<<full.stats.max_story_drift_abs<<" woodbury_peak_story_drift="<<wood.stats.max_story_drift_abs<<"\n";
}
