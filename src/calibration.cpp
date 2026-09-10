#include "quake/calibration.hpp"
#include "quake/low_rank_solver.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace quake {
namespace {
constexpr double beta=0.25, gamma=0.5;

template <class Fn>
double timed(Fn&& fn) {
    const auto t0=std::chrono::steady_clock::now();
    fn();
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
}

std::vector<int> default_ranks(int m) {
    std::vector<int> r;
    for(int x: {0,1,2,5,10,20,40,80,120,200,320,500}) if(x<=m) r.push_back(x);
    if(m>0 && (r.empty() || r.back()!=m)) r.push_back(m);
    std::sort(r.begin(),r.end()); r.erase(std::unique(r.begin(),r.end()),r.end());
    return r;
}
}

SolverCalibrationResult calibrate_solver_crossover(const NonlinearDynamicModel& model,double dt,
                                                     std::vector<int> ranks,int repeats,
                                                     double tangent_ratio) {
    const auto calibration_start=std::chrono::steady_clock::now();
    if(dt<=0.0 || repeats<1 || tangent_ratio<0.0 || tangent_ratio>1.0)
        throw std::invalid_argument("invalid solver calibration options");
    const int n=model.dof(),m=model.nonlinear_count();
    if(ranks.empty()) ranks=default_ranks(m);
    for(int r:ranks) if(r<0||r>m) throw std::out_of_range("solver calibration rank");
    std::sort(ranks.begin(),ranks.end()); ranks.erase(std::unique(ranks.begin(),ranks.end()),ranks.end());

    const double a0=1.0/(beta*dt*dt), a1=gamma/(beta*dt);
    auto baseline=model.effective_initial_matrix(a0,a1);
    LazyLowRankWoodburySolver wood(baseline,model.nonlinear_basis());
    SuperLUSamePatternSolver direct(baseline);
    std::vector<double> rhs(static_cast<std::size_t>(n));
    for(int i=0;i<n;++i) rhs[static_cast<std::size_t>(i)]=std::sin(0.037*(i+1))+0.25*std::cos(0.013*(i+3));
    const auto& init=model.initial_nonlinear_tangents();

    SolverCalibrationResult out; out.dof=n; out.nonlinear_count=m; out.matrix_nnz=static_cast<std::size_t>(baseline.nnz());
    int last_faster=-1;
    for(int r:ranks) {
        std::vector<int> active;
        active.reserve(static_cast<std::size_t>(r));
        if(r>0) {
            // Spread the selected update directions across the whole building
            // rather than biasing calibration toward the first story/bay.
            for(int a=0;a<r;++a) {
                int idx=static_cast<int>((static_cast<long long>(a)*m)/r);
                idx=std::min(m-1,idx); active.push_back(idx);
            }
            std::sort(active.begin(),active.end()); active.erase(std::unique(active.begin(),active.end()),active.end());
            for(int idx=0;static_cast<int>(active.size())<r && idx<m;++idx)
                if(!std::binary_search(active.begin(),active.end(),idx)){active.push_back(idx);std::sort(active.begin(),active.end());}
        }
        std::vector<double> dk(static_cast<std::size_t>(m),0.0), tangents=init;
        for(int j:active){
            tangents[static_cast<std::size_t>(j)]=tangent_ratio*init[static_cast<std::size_t>(j)];
            dk[static_cast<std::size_t>(j)]=tangents[static_cast<std::size_t>(j)]-init[static_cast<std::size_t>(j)];
        }
        // Populate lazy influence cache outside the timed region.
        (void)wood.solve(rhs,dk);
        auto K=model.effective_tangent_matrix(tangents,a0,a1);
        (void)direct.refactor_and_solve(K,rhs);

        volatile double sink=0.0;
        const double tw=timed([&]{for(int k=0;k<repeats;++k){auto x=wood.solve(rhs,dk);sink+=x[static_cast<std::size_t>(k%n)];}});
        const double td=timed([&]{for(int k=0;k<repeats;++k){
            if(r==0){auto x=direct.solve_current(rhs);sink+=x[static_cast<std::size_t>((k+7)%n)];continue;}
            // Slight value perturbation forces a real numeric refactor while
            // preserving the exact sparsity pattern and representative rank.
            const double scale=1.0-1e-6*(k+1);
            auto tk=tangents; for(int j:active) tk[static_cast<std::size_t>(j)]*=scale;
            auto A=model.effective_tangent_matrix(tk,a0,a1); auto x=direct.refactor_and_solve(A,rhs); sink+=x[static_cast<std::size_t>((k+7)%n)];
        }});
        (void)sink;
        const double speed=tw>0.0?td/tw:0.0;
        out.points.push_back({r,tw,td,speed});
    }
    bool crossed=false;
    for(const auto& pt:out.points){
        if(pt.active_rank==0) continue;
        if(!crossed && pt.speedup>1.0) last_faster=pt.active_rank;
        else if(pt.speedup<=1.0) crossed=true;
    }
    // If the model has no nonlinear directions, rank zero is the only valid case.
    if(m==0) last_faster=0;
    out.recommended_woodbury_rank_limit=last_faster;
    out.recommended_rank_fraction=n>0 && last_faster>=0?static_cast<double>(last_faster)/n:0.0;
    out.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-calibration_start).count();
    return out;
}

} // namespace quake
