#include "quake/newmark.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/stability.hpp"
#include "quake/generalized_woodbury.hpp"
#include "quake/analysis_input.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <limits>
#include <stdexcept>
#include <utility>
#include <optional>

namespace quake {
namespace {

constexpr double kBeta=0.25;
constexpr double kGamma=0.5;

double inf_norm(const std::vector<double>& x){double n=0.0;for(double v:x){if(!std::isfinite(v))return std::numeric_limits<double>::infinity();n=std::max(n,std::abs(v));}return n;}
double vec_dot(const std::vector<double>& a,const std::vector<double>& b){double s=0.0;for(std::size_t i=0;i<a.size();++i)s+=a[i]*b[i];return s;}

struct DynamicState {
    std::vector<double> u,v,a;
    std::vector<double> committed;
    std::vector<double> tangents;
};

struct AcceptedSubstepSnapshot {
    int subdivision_depth{};
    double time{};
    double ground_accel{};
    DynamicState state;
};

struct EnergyLedger {
    double input{},internal{},damping{};
    double max_roof{},max_drift{},max_drift_ratio{};
};

struct Eval {
    std::vector<double> residual;
    std::vector<double> tangents;
    std::vector<double> trial_states;
    std::vector<double> a;
    std::vector<double> v;
    double norm{};
    double scale{};
};

class SolverContext {
public:
    SolverContext(const NonlinearDynamicModel& model,double dt,LinearStrategy strategy,
                  AnalysisStats& stats)
        : model_(model),dt_(dt),strategy_(strategy),
          a0_(1.0/(kBeta*dt*dt)),a1_(kGamma/(kBeta*dt)),
          factored_tangents_(static_cast<std::size_t>(model.nonlinear_count())) {
        factored_tangents_=model.initial_nonlinear_tangents();
        if(strategy_==LinearStrategy::Woodbury){
            const auto base=model.effective_initial_matrix(a0_,a1_);
            if(model_.has_generalized_state_update()) generalized_=std::make_unique<GeneralizedWoodburySolver>(base,model_.generalized_state_update_basis());
            else wood_=std::make_unique<LazyLowRankWoodburySolver>(base,model.nonlinear_basis());
            ++stats.global_factorizations;
            if(model_.has_state_dependent_global_tangent()||model_.has_velocity_dependent_global_tangent())
                same_=std::make_unique<SuperLUSamePatternSolver>(base);
        } else if(strategy_==LinearStrategy::SamePatternRefactorization){
            same_=std::make_unique<SuperLUSamePatternSolver>(model.effective_initial_matrix(a0_,a1_));
            ++stats.global_factorizations;
        }
    }

    double a0() const{return a0_;}
    double a1() const{return a1_;}

    std::vector<double> solve(const std::vector<double>& rhs,const std::vector<double>& tangents,const std::vector<double>& u,
                              const std::vector<double>& velocity,
                              const std::vector<double>& committed_state,
                              AnalysisStats& stats,const RobustNewmarkOptions& options){
        ++stats.linear_solves;
        // Some physically dense tangent contributions are nevertheless very
        // low rank. Fixed-mode modal damping is the primary case: its force
        // remains on the residual side, but the exact Newmark velocity
        // Jacobian is a rank-r update a1*sum(2*zeta*w mphi mphi^T). Apply that
        // update by Woodbury around the ordinary sparse structural tangent so
        // Newton convergence does not require a dense global damping matrix.
        if(model_.has_additional_effective_low_rank_update()){
            auto k=(model_.has_state_dependent_global_tangent()||model_.has_velocity_dependent_global_tangent())
                ?model_.effective_state_tangent_matrix_with_state_and_velocity(u,tangents,committed_state,velocity,a0_,a1_)
                :model_.effective_tangent_matrix(tangents,a0_,a1_);
            GeneralizedWoodburySolver low_rank(std::move(k),model_.additional_effective_low_rank_basis());
            auto C=model_.additional_effective_low_rank_coefficients(a0_,a1_);
            ++stats.global_factorizations;++stats.reduced_update_solves;++stats.generalized_update_solves;
            auto x=low_rank.solve(rhs,C);
            stats.minimum_reduced_pivot_ratio=std::min(stats.minimum_reduced_pivot_ratio,low_rank.last_reduced_pivot_ratio());
            stats.max_generalized_update_rank=std::max(stats.max_generalized_update_rank,low_rank.last_active_dimension());
            return x;
        }
        if(model_.has_state_dependent_global_tangent()||model_.has_velocity_dependent_global_tangent()){
            if(!model_.has_velocity_dependent_global_tangent()&&strategy_==LinearStrategy::Woodbury && generalized_){
                auto C=model_.generalized_state_update_coefficients(u,tangents);const int r=generalized_->update_dimension();std::size_t active=0;
                for(int i=0;i<r;++i){double mx=0.0;for(int j=0;j<r;++j)mx=std::max({mx,std::abs(C[static_cast<std::size_t>(i*r+j)]),std::abs(C[static_cast<std::size_t>(j*r+i)])});if(mx>1e-14)++active;}
                stats.active_rank_sum+=active;stats.max_active_rank=std::max(stats.max_active_rank,active);stats.max_generalized_update_rank=std::max(stats.max_generalized_update_rank,active);
                const bool proactive=options.direct_rank_fraction>0.0 && static_cast<double>(active)>options.direct_rank_fraction*model_.dof();
                if(!proactive){try{auto x=generalized_->solve(rhs,C);++stats.reduced_update_solves;++stats.generalized_update_solves;stats.minimum_reduced_pivot_ratio=std::min(stats.minimum_reduced_pivot_ratio,generalized_->last_reduced_pivot_ratio());return x;}catch(const std::exception&){stats.minimum_reduced_pivot_ratio=0.0;}}
            }
            auto k=model_.effective_state_tangent_matrix_with_state_and_velocity(u,tangents,committed_state,velocity,a0_,a1_);++stats.global_factorizations;++stats.direct_fallbacks;if(same_)return same_->refactor_and_solve(k,rhs);return superlu_solve_once(k,rhs);
        }
        if(strategy_==LinearStrategy::FullFactorization){
            auto k=model_.effective_tangent_matrix(tangents,a0_,a1_);++stats.global_factorizations;return superlu_solve_once(k,rhs);
        }
        if(strategy_==LinearStrategy::SamePatternRefactorization){
            bool same=true;for(std::size_t j=0;j<tangents.size();++j)if(tangents[j]!=factored_tangents_[j]){same=false;break;}
            if(same)return same_->solve_current(rhs);
            auto k=model_.effective_tangent_matrix(tangents,a0_,a1_);auto x=same_->refactor_and_solve(k,rhs);factored_tangents_=tangents;++stats.global_factorizations;return x;
        }
        std::vector<double> dk(tangents.size());std::size_t rank=0;
        for(std::size_t j=0;j<tangents.size();++j){dk[j]=tangents[j]-model_.initial_nonlinear_tangents()[j];if(std::abs(dk[j])>1e-14)++rank;}
        stats.active_rank_sum+=rank;stats.max_active_rank=std::max(stats.max_active_rank,rank);++stats.reduced_update_solves;
        const bool proactive=options.direct_rank_fraction>0.0 && static_cast<double>(rank)>options.direct_rank_fraction*model_.dof();
        if(!proactive){
            try{auto x=wood_->solve(rhs,dk);stats.minimum_reduced_pivot_ratio=std::min(stats.minimum_reduced_pivot_ratio,wood_->last_reduced_pivot_ratio());return x;}catch(const std::exception&){stats.minimum_reduced_pivot_ratio=0.0;/* robust direct fallback below */}
        }
        auto k=model_.effective_tangent_matrix(tangents,a0_,a1_);++stats.global_factorizations;++stats.direct_fallbacks;return superlu_solve_once(k,rhs);
    }

    std::size_t cached_columns() const{return generalized_?static_cast<std::size_t>(generalized_->update_dimension()):(wood_?wood_->cached_update_count():0);}
private:
    const NonlinearDynamicModel& model_; double dt_{}; LinearStrategy strategy_{}; double a0_{},a1_{};
    std::vector<double> factored_tangents_;
    std::unique_ptr<LazyLowRankWoodburySolver> wood_;
    std::unique_ptr<GeneralizedWoodburySolver> generalized_;
    std::unique_ptr<SuperLUSamePatternSolver> same_;
};

Eval evaluate(const NonlinearDynamicModel& model,const DynamicState& start,
              const std::vector<double>& u_trial,const std::vector<double>& u_pred,
              const std::vector<double>& v_pred,double ag,const std::vector<double>& constant_load,SolverContext& ctx,
              AnalysisStats& stats,std::string& last_error){
    const int n=model.dof();Eval e;
    e.norm=std::numeric_limits<double>::infinity();e.scale=1.0;
    ++stats.residual_evaluations;
    if(!all_finite(u_trial)){++stats.nonfinite_evaluations;last_error="non-finite trial displacement";return e;}
    std::vector<double> fint; NonlinearEvalDiagnostics nd;
    try{model.internal_force_and_tangent_diagnostics(u_trial,start.committed,fint,e.tangents,e.trial_states,&nd);}
    catch(const ConstitutiveIntegrationError& ex){++stats.constitutive_integration_failures;last_error=ex.what();return e;}
    if(!all_finite(fint)||!all_finite(e.tangents)||!all_finite(e.trial_states)){
        ++stats.nonfinite_evaluations;last_error="non-finite constitutive force, tangent, or state";return e;
    }
    stats.nonlinear_component_evaluations += nd.component_evaluations;
    stats.nonlinear_active_component_evaluations += nd.active_tangent_evaluations;
    stats.nonlinear_fast_path_evaluations += nd.fast_path_evaluations;
    stats.nonlinear_full_state_evaluations += nd.full_state_evaluations;
    stats.nonlinear_transition_events += nd.transition_events;
    stats.nonlinear_reversal_events += nd.reversal_events;
    stats.nonlinear_deterioration_events += nd.deterioration_events;
    stats.nonlinear_failure_events += nd.failure_events;
    stats.nonlinear_io_or_beyond += nd.io_or_beyond_evaluations;
    stats.nonlinear_ls_or_beyond += nd.ls_or_beyond_evaluations;
    stats.nonlinear_cp_or_beyond += nd.cp_or_beyond_evaluations;
    stats.nonlinear_beyond_cp += nd.beyond_cp_evaluations;
    stats.nonlinear_lateral_loss += nd.lateral_loss_evaluations;
    e.a.resize(static_cast<std::size_t>(n));e.v.resize(static_cast<std::size_t>(n));
    for(int i=0;i<n;++i)e.a[static_cast<std::size_t>(i)]=ctx.a0()*(u_trial[static_cast<std::size_t>(i)]-u_pred[static_cast<std::size_t>(i)]);
    // v = v_pred + gamma*dt*a; ctx does not expose dt to keep solve API narrow.
    // Recover gamma*dt = a1/a0 for average-acceleration Newmark.
    const double gamma_dt=ctx.a1()/ctx.a0();
    for(int i=0;i<n;++i)e.v[static_cast<std::size_t>(i)]=v_pred[static_cast<std::size_t>(i)]+gamma_dt*e.a[static_cast<std::size_t>(i)];
    auto cv=model.damping_multiply(e.v),ma=model.mass_multiply(e.a),p=model.base_excitation(ag);
    e.residual.resize(static_cast<std::size_t>(n));
    for(int i=0;i<n;++i)e.residual[static_cast<std::size_t>(i)]=p[static_cast<std::size_t>(i)]+constant_load[static_cast<std::size_t>(i)]-fint[static_cast<std::size_t>(i)]-cv[static_cast<std::size_t>(i)]-ma[static_cast<std::size_t>(i)];
    e.norm=inf_norm(e.residual);e.scale=std::max(1.0,inf_norm(p));
    if(!std::isfinite(e.norm)||!std::isfinite(e.scale)||!all_finite(e.a)||!all_finite(e.v)){
        ++stats.nonfinite_evaluations;last_error="non-finite dynamic residual or kinematics";e.norm=std::numeric_limits<double>::infinity();
    }
    return e;
}

EnergyLedger accepted_energy_increment(const NonlinearDynamicModel& model,const DynamicState& a,
                                             const DynamicState& b,double ag0,double ag1,double dt,
                                             const std::vector<double>& constant_load){
    EnergyLedger e;std::vector<double> du(a.u.size());for(std::size_t i=0;i<du.size();++i)du[i]=b.u[i]-a.u[i];
    std::vector<double> f0,t0,s0,f1,t1,s1;
    model.internal_force_and_tangent(a.u,a.committed,f0,t0,s0);
    model.internal_force_and_tangent(b.u,b.committed,f1,t1,s1);
    auto p0=model.base_excitation(ag0),p1=model.base_excitation(ag1);
    std::vector<double> favg(f0.size()),pavg(p0.size());for(std::size_t i=0;i<favg.size();++i){favg[i]=0.5*(f0[i]+f1[i]);pavg[i]=0.5*(p0[i]+p1[i])+constant_load[i];}
    e.internal=vec_dot(favg,du);e.input=vec_dot(pavg,du);
    auto cv0=model.damping_multiply(a.v),cv1=model.damping_multiply(b.v);
    e.damping=0.5*dt*(vec_dot(a.v,cv0)+vec_dot(b.v,cv1));
    e.max_roof=std::abs(model.response_value(b.u));
    e.max_drift=model.max_drift_measure(b.u);
    const double ratio=model.max_drift_ratio(b.u);
    if(std::isfinite(ratio))e.max_drift_ratio=ratio;
    return e;
}

// We need dt in predictor; keep it outside SolverContext private implementation.
bool attempt_step(const NonlinearDynamicModel& model,const DynamicState& start,double ag,double dt,
                  SolverContext& ctx,const RobustNewmarkOptions& options,
                  const std::vector<double>& constant_load,AnalysisStats& stats,DynamicState& out,std::string& last_error){
    const int n=model.dof();
    std::vector<double> u_pred(static_cast<std::size_t>(n)),v_pred(static_cast<std::size_t>(n));
    for(int i=0;i<n;++i){u_pred[static_cast<std::size_t>(i)]=start.u[static_cast<std::size_t>(i)]+dt*start.v[static_cast<std::size_t>(i)]+dt*dt*(0.5-kBeta)*start.a[static_cast<std::size_t>(i)];v_pred[static_cast<std::size_t>(i)]=start.v[static_cast<std::size_t>(i)]+dt*(1.0-kGamma)*start.a[static_cast<std::size_t>(i)];}
    std::vector<double> u_trial=options.kinematic_initial_guess?u_pred:start.u;
    std::optional<Eval> cached;
    for(int iter=0;iter<options.max_iterations;++iter){
        Eval cur=cached?std::move(*cached):evaluate(model,start,u_trial,u_pred,v_pred,ag,constant_load,ctx,stats,last_error);cached.reset();++stats.newton_iterations;
        stats.last_residual_norm=cur.norm;
        stats.last_residual_tolerance=options.tolerance*(options.relative_force_tolerance?cur.scale:1.0);
        if(!std::isfinite(cur.norm)||!std::isfinite(cur.scale))return false;
        if(cur.norm<=options.tolerance*(options.relative_force_tolerance?cur.scale:1.0)){out.u=std::move(u_trial);out.v=std::move(cur.v);out.a=std::move(cur.a);out.committed=std::move(cur.trial_states);out.tangents=std::move(cur.tangents);return true;}
        std::vector<double> du;
        try{du=ctx.solve(cur.residual,cur.tangents,u_trial,cur.v,start.committed,stats,options);}
        catch(const ConstitutiveIntegrationError& ex){++stats.constitutive_integration_failures;last_error=ex.what();return false;}
        catch(const std::runtime_error& ex){++stats.linear_solve_failures;last_error=ex.what();return false;}
        if(!all_finite(du)){++stats.nonfinite_evaluations;last_error="non-finite Newton correction";return false;}
        if(!options.line_search){for(int i=0;i<n;++i)u_trial[static_cast<std::size_t>(i)]+=du[static_cast<std::size_t>(i)];continue;}
        double alpha=1.0;bool accepted=false;double best=cur.norm;std::vector<double> best_u;
        double best_any=std::numeric_limits<double>::infinity();std::vector<double> best_any_u;int best_any_bt=0;
        std::optional<Eval> best_any_eval;
        for(int bt=0;bt<=options.max_backtracks;++bt){
            std::vector<double> cand=u_trial;for(int i=0;i<n;++i)cand[static_cast<std::size_t>(i)]+=alpha*du[static_cast<std::size_t>(i)];
            Eval ce=evaluate(model,start,cand,u_pred,v_pred,ag,constant_load,ctx,stats,last_error);
            const bool best_candidate=std::isfinite(ce.norm)&&ce.norm<best_any;
            if(best_candidate){best_any=ce.norm;best_any_u=cand;best_any_bt=bt;}
            if(std::isfinite(ce.norm) && ce.norm<best){best=ce.norm;best_u=std::move(cand);cached=std::move(ce);accepted=true;if(bt>0)stats.line_search_backtracks+=static_cast<std::size_t>(bt);break;}
            if(best_candidate)best_any_eval=std::move(ce);
            alpha*=options.backtrack_ratio;
        }
        if(!accepted && options.nonmonotone_line_search_factor>1.0 && std::isfinite(best_any) && best_any<=options.nonmonotone_line_search_factor*cur.norm){best_u=std::move(best_any_u);cached=std::move(best_any_eval);accepted=true;stats.line_search_backtracks+=static_cast<std::size_t>(best_any_bt);}
        if(!accepted){++stats.line_search_failures;return false;}
        u_trial=std::move(best_u);
    }
    ++stats.max_iteration_failures;return false;
}

} // namespace

namespace {

using RobustContextGetter = std::function<SolverContext&(int)>;
using RobustCachedColumnGetter = std::function<std::size_t()>;

AnalysisResult run_newmark_robust_core(const NonlinearDynamicModel& model,const std::vector<double>& ground_accel,
                                       double dt,LinearStrategy strategy,const RobustNewmarkOptions& options,
                                       const RobustContextGetter& get_ctx,
                                       const RobustCachedColumnGetter& cached_columns){
    validate_transient_input(ground_accel,dt,options.tolerance,options.max_iterations);
    if(!std::isfinite(options.backtrack_ratio)||!std::isfinite(options.nonmonotone_line_search_factor)||!std::isfinite(options.direct_rank_fraction)||options.max_subdivisions>30)
        throw std::invalid_argument("invalid finite robust Newmark controls");
    if(options.max_iterations<1||options.max_backtracks<0||options.max_subdivisions<0||options.backtrack_ratio<=0.0||options.backtrack_ratio>=1.0||options.nonmonotone_line_search_factor<1.0)throw std::invalid_argument("invalid robust Newmark options");
    if(strategy==LinearStrategy::ModifiedNewton || strategy==LinearStrategy::Adaptive) throw std::invalid_argument("robust driver does not yet support ModifiedNewton/Adaptive");
    const int n=model.dof();
    auto sized=[&](const std::vector<double>& v,std::size_t count,const char* name){if(!v.empty()&&v.size()!=count)throw std::invalid_argument(std::string(name)+" size mismatch");return v.empty()?std::vector<double>(count,0.0):v;};
    DynamicState state{
        sized(options.initial_displacement,static_cast<std::size_t>(n),"initial displacement"),
        sized(options.initial_velocity,static_cast<std::size_t>(n),"initial velocity"),
        sized(options.initial_acceleration,static_cast<std::size_t>(n),"initial acceleration"),
        options.initial_committed_state.empty()?model.initial_nonlinear_state():options.initial_committed_state,
        model.initial_nonlinear_tangents()};
    if(state.committed.size()!=static_cast<std::size_t>(model.nonlinear_state_size()))throw std::invalid_argument("initial committed state size mismatch");
    {std::vector<double> initial_force,trial_state;model.internal_force_and_tangent(state.u,state.committed,initial_force,state.tangents,trial_state);}
    const std::vector<double> zero_load(static_cast<std::size_t>(n),0.0);
    const auto& constant_load=options.constant_load.empty()?zero_load:options.constant_load;
    if(constant_load.size()!=static_cast<std::size_t>(n))throw std::invalid_argument("constant load size mismatch");
    AnalysisResult result;result.roof_history.reserve(ground_accel.size());result.stats.minimum_dt=dt;
    const auto start_time=std::chrono::steady_clock::now();
    double initial_stability=0.0;
    const bool tangent_tracking=options.collapse.tangent_check_interval>0 && options.collapse.min_tangent_ratio>0.0;
    if(options.collapse.check_initial_stability || tangent_tracking){
        // Initial admissibility and subsequent mechanism tracking are distinct
        // questions. First certify positive-definiteness with sparse symmetric
        // elimination; only then (if requested) use inverse iteration to
        // establish the near-zero eigenvalue reference for the evolving tangent.
        auto ia=assess_initial_stability_from_tangents(model,state.u,state.tangents,1e-12,0);
        result.stats.initial_spd_minimum_pivot_ratio=ia.minimum_pivot_ratio;
        result.stats.initial_symmetry_relative_error=ia.symmetry_relative_error;
        if(ia.status!=InitialStabilityStatus::PositiveDefinite){
            result.termination=AnalysisTermination::InitialInstability;result.collapse_mechanism=CollapseMechanism::InitialInstability;
            result.termination_step=0;result.termination_time=0.0;
            result.termination_reason=ia.status==InitialStabilityStatus::NotPositiveDefinite ?
                "initial structural tangent failed sparse positive-definiteness certification" :
                "initial structural tangent stability certification failed numerically";
            result.stats.minimum_tangent_ratio=1.0;
            result.stats.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start_time).count();
            result.final_displacement=state.u;result.final_state=state.committed;return result;
        }
    }
    if(tangent_tracking){
        auto se=estimate_tangent_stability(model,state.u,state.committed);
        if(!se.factorization_ok || se.near_zero_eigenvalue<=0.0){
            result.termination=AnalysisTermination::InitialInstability;result.collapse_mechanism=CollapseMechanism::InitialInstability;
            result.termination_step=0;result.termination_time=0.0;
            result.termination_reason="initial tangent passed SPD certification but near-zero reference solve failed";
            result.stats.minimum_tangent_ratio=1.0;
            result.stats.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start_time).count();
            result.final_displacement=state.u;result.final_state=state.committed;return result;
        }
        initial_stability=se.absolute_eigenvalue;
        result.stats.minimum_tangent_eigenvalue=se.near_zero_eigenvalue;
    }

    bool current_step_lateral_seen=false;
    using SubstepBuffer=std::vector<AcceptedSubstepSnapshot>;
    std::function<bool(const DynamicState&,double,double,double,int,DynamicState&,EnergyLedger&,SubstepBuffer*)> advance;
    advance=[&](const DynamicState& s,double ag0,double ag1,double t0,int depth,DynamicState& out,EnergyLedger& energy,SubstepBuffer* accepted)->bool{
        const double d=dt/std::pow(2.0,depth);result.stats.minimum_dt=std::min(result.stats.minimum_dt,d);auto& ctx=get_ctx(depth);const std::size_t lat0=result.stats.nonlinear_lateral_loss;
        if(attempt_step(model,s,ag1,d,ctx,options,constant_load,result.stats,out,result.last_integration_error)){
            if(result.stats.nonlinear_lateral_loss>lat0)current_step_lateral_seen=true;
            energy=accepted_energy_increment(model,s,out,ag0,ag1,d,constant_load);++result.stats.internal_substeps;
            if(accepted) accepted->push_back({depth,t0+d,ag1,out});
            return true;
        }
        if(result.stats.nonlinear_lateral_loss>lat0)current_step_lateral_seen=true;
        if(depth>=options.max_subdivisions)return false;
        ++result.stats.subdivided_steps;const double mid=0.5*(ag0+ag1);DynamicState half;EnergyLedger e1,e2;
        SubstepBuffer first,second;
        SubstepBuffer* b1=accepted?&first:nullptr; SubstepBuffer* b2=accepted?&second:nullptr;
        if(!advance(s,ag0,mid,t0,depth+1,half,e1,b1))return false;
        if(!advance(half,mid,ag1,t0+0.5*d,depth+1,out,e2,b2))return false;
        if(accepted){accepted->insert(accepted->end(),std::make_move_iterator(first.begin()),std::make_move_iterator(first.end()));accepted->insert(accepted->end(),std::make_move_iterator(second.begin()),std::make_move_iterator(second.end()));}
        energy.input=e1.input+e2.input;energy.internal=e1.internal+e2.internal;energy.damping=e1.damping+e2.damping;
        energy.max_roof=std::max(e1.max_roof,e2.max_roof);energy.max_drift=std::max(e1.max_drift,e2.max_drift);
        energy.max_drift_ratio=std::max(e1.max_drift_ratio,e2.max_drift_ratio);return true;
    };

    double ag_prev=0.0;
    for(std::size_t step=0;step<ground_accel.size();++step){
        DynamicState next;EnergyLedger step_energy;current_step_lateral_seen=false;
        SubstepBuffer accepted_substeps; SubstepBuffer* accepted_ptr=options.accepted_substep_state_observer?&accepted_substeps:nullptr;
        if(!advance(state,ag_prev,ground_accel[step],step*dt,0,next,step_energy,accepted_ptr)){
            ++result.stats.failed_steps;
            if(options.collapse.classify_unrecoverable_lateral_loss_as_collapse && current_step_lateral_seen){
                result.termination=AnalysisTermination::PhysicalCollapse;result.collapse_mechanism=CollapseMechanism::UnrecoverableLateralLoss;
                result.termination_reason="unrecoverable loss of equilibrium coincident with component E-point lateral resistance loss";
            }else{
                result.termination=AnalysisTermination::NumericalFailure;result.collapse_mechanism=CollapseMechanism::NumericalFailure;
                result.termination_reason="robust Newmark failed after maximum subdivision";
            }
            result.termination_step=step;result.termination_time=(step+1)*dt;
            if(result.termination==AnalysisTermination::PhysicalCollapse || options.return_numerical_failure) break;
            throw std::runtime_error("robust Newmark step failed after subdivision at step "+std::to_string(step));
        }
        state=std::move(next);ag_prev=ground_accel[step];
        if(options.accepted_substep_state_observer){
            for(std::size_t ss=0;ss<accepted_substeps.size();++ss){
                const auto& snap=accepted_substeps[ss];
                options.accepted_substep_state_observer(step,ss,snap.subdivision_depth,snap.time,snap.ground_accel,
                    snap.state.u,snap.state.v,snap.state.a,snap.state.committed);
            }
        }
        if(options.accepted_step_observer)
            options.accepted_step_observer(step,(step+1)*dt,ground_accel[step],state.u,state.v,state.a);
        if(options.accepted_state_observer)
            options.accepted_state_observer(step,(step+1)*dt,ground_accel[step],state.u,state.v,state.a,state.committed);
        result.stats.input_energy+=step_energy.input;result.stats.internal_work+=step_energy.internal;result.stats.damping_energy+=step_energy.damping;
        auto mv_energy=model.mass_multiply(state.v);result.stats.kinetic_energy=0.5*vec_dot(state.v,mv_energy);
        result.stats.energy_balance_error=result.stats.input_energy-result.stats.internal_work-result.stats.damping_energy-result.stats.kinetic_energy;
        const double escale=std::max({1e-12,std::abs(result.stats.input_energy),std::abs(result.stats.internal_work)+std::abs(result.stats.damping_energy)+std::abs(result.stats.kinetic_energy)});
        result.stats.energy_balance_relative_error=std::abs(result.stats.energy_balance_error)/escale;
        const double r=model.response_value(state.u);result.roof_history.push_back(r);
        result.stats.max_roof_abs=std::max(result.stats.max_roof_abs,step_energy.max_roof);
        result.stats.max_story_drift_abs=std::max(result.stats.max_story_drift_abs,step_energy.max_drift);
        const double drift_ratio=step_energy.max_drift_ratio;
        if(std::isfinite(drift_ratio)) result.max_story_drift_ratio=std::max(result.max_story_drift_ratio,drift_ratio);
        ++result.stats.steps;

        // Limit-state assessment is intentionally separate from convergence.
        // CP exceedance is reported but is not physical collapse unless the
        // caller explicitly makes it a termination criterion.
        std::vector<double> f,t,trial; NonlinearEvalDiagnostics ld;
        model.internal_force_and_tangent_diagnostics(state.u,state.committed,f,t,trial,&ld);
        const auto none=static_cast<std::size_t>(-1);
        if(ld.io_or_beyond_evaluations>0 && result.first_io_step==none) result.first_io_step=step;
        if(ld.ls_or_beyond_evaluations>0 && result.first_ls_step==none) result.first_ls_step=step;
        if(ld.cp_or_beyond_evaluations>0 && result.first_cp_step==none) result.first_cp_step=step;
        if(ld.lateral_loss_evaluations>0 && result.first_lateral_loss_step==none) result.first_lateral_loss_step=step;
        if(ld.failure_events>0 && result.first_component_failure_step==none) result.first_component_failure_step=step;
        result.max_components_io_or_beyond=std::max(result.max_components_io_or_beyond,ld.io_or_beyond_evaluations);
        result.max_components_ls_or_beyond=std::max(result.max_components_ls_or_beyond,ld.ls_or_beyond_evaluations);
        result.max_components_cp_or_beyond=std::max(result.max_components_cp_or_beyond,ld.cp_or_beyond_evaluations);
        result.max_components_beyond_cp=std::max(result.max_components_beyond_cp,ld.beyond_cp_evaluations);
        result.max_failed_components=std::max(result.max_failed_components,ld.failure_events);
        result.max_components_lateral_loss=std::max(result.max_components_lateral_loss,ld.lateral_loss_evaluations);

        bool collapsed=false;
        if(options.collapse.tangent_check_interval>0 && options.collapse.min_tangent_ratio>0.0 && initial_stability>0.0 && ((step+1)%static_cast<std::size_t>(options.collapse.tangent_check_interval)==0)){
            auto se=estimate_tangent_stability(model,state.u,state.committed);
            if(!se.factorization_ok){collapsed=true;result.collapse_mechanism=CollapseMechanism::TangentInstability;result.termination_reason="structural tangent became singular during stability check";result.stats.minimum_tangent_ratio=0.0;}
            else{
                const double ratio=se.absolute_eigenvalue/initial_stability;
                result.stats.minimum_tangent_ratio=std::min(result.stats.minimum_tangent_ratio,ratio);
                result.stats.minimum_tangent_eigenvalue=std::min(result.stats.minimum_tangent_eigenvalue,se.near_zero_eigenvalue);
                if(se.near_zero_eigenvalue<=0.0 || ratio<=options.collapse.min_tangent_ratio){collapsed=true;result.collapse_mechanism=CollapseMechanism::TangentInstability;result.termination_reason="near-zero or negative structural tangent criterion exceeded";}
            }
        }
        if(!collapsed && options.collapse.max_story_drift_ratio>0.0 && std::isfinite(drift_ratio) && drift_ratio>=options.collapse.max_story_drift_ratio){collapsed=true;result.collapse_mechanism=CollapseMechanism::DriftLimit;result.termination_reason="story drift ratio collapse criterion exceeded";}
        if(!collapsed && options.collapse.max_lateral_loss_components>0 && ld.lateral_loss_evaluations>=options.collapse.max_lateral_loss_components){collapsed=true;result.collapse_mechanism=CollapseMechanism::LateralResistanceLoss;result.termination_reason="component lateral-resistance-loss criterion exceeded";}
        if(!collapsed && options.collapse.max_failed_components>0 && ld.failure_events>=options.collapse.max_failed_components){collapsed=true;result.collapse_mechanism=CollapseMechanism::GravityResistanceLoss;result.termination_reason="component gravity/effective-resistance-loss criterion exceeded";}
        if(!collapsed && options.collapse.max_beyond_cp_components>0 && ld.beyond_cp_evaluations>=options.collapse.max_beyond_cp_components){collapsed=true;result.collapse_mechanism=CollapseMechanism::BeyondCPConfigured;result.termination_reason="configured beyond-CP component criterion exceeded";}
        if(collapsed){result.termination=AnalysisTermination::PhysicalCollapse;result.termination_step=step;result.termination_time=(step+1)*dt;break;}
    }
    if(result.termination==AnalysisTermination::Completed){result.termination_step=result.stats.steps;result.termination_time=result.stats.steps*dt;result.termination_reason="completed";}
    result.stats.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start_time).count();
    if(cached_columns) result.stats.cached_low_rank_columns=std::max(result.stats.cached_low_rank_columns,cached_columns());
    result.final_displacement=state.u;result.final_state=state.committed;return result;
}


} // namespace

AnalysisResult run_newmark_robust(const NonlinearDynamicModel& model,const std::vector<double>& ground_accel,
                                  double dt,LinearStrategy strategy,const RobustNewmarkOptions& options){
    std::map<int,std::unique_ptr<SolverContext>> contexts;
    AnalysisStats setup_stats;
    auto getter=[&](int depth)->SolverContext&{
        auto it=contexts.find(depth);if(it!=contexts.end())return *it->second;
        const double d=dt/std::pow(2.0,depth);
        auto p=std::make_unique<SolverContext>(model,d,strategy,setup_stats);auto* raw=p.get();
        contexts.emplace(depth,std::move(p));return *raw;
    };
    auto cached=[&](){std::size_t n=0;for(const auto& [depth,ctx]:contexts){(void)depth;n=std::max(n,ctx->cached_columns());}return n;};
    auto result=run_newmark_robust_core(model,ground_accel,dt,strategy,options,getter,cached);
    result.stats.global_factorizations+=setup_stats.global_factorizations;
    return result;
}

namespace {

bool same_sparse_numeric_matrix(const SparseMatrixCSC& a,const SparseMatrixCSC& b){
    return a.rows()==b.rows()&&a.cols()==b.cols()&&a.col_ptr()==b.col_ptr()&&
           a.row_ind()==b.row_ind()&&a.values()==b.values();
}

} // namespace

struct PreparedRobustNewmark::Impl {
    std::map<int,std::unique_ptr<SolverContext>> contexts;
    std::map<int,SparseMatrixCSC> prepared_effective_initial;
    std::size_t preparation_count{};
    std::size_t setup_factorizations{};
    bool last_reused{};
};

PreparedRobustNewmark::PreparedRobustNewmark(const NonlinearDynamicModel& model,double dt,
                                             LinearStrategy strategy,int prepared_subdivision_depth)
    :model_(&model),dt_(dt),strategy_(strategy),prepared_subdivision_depth_(prepared_subdivision_depth),
     impl_(std::make_unique<Impl>()){
    if(!std::isfinite(dt_)||dt_<=0.0||prepared_subdivision_depth_<0||prepared_subdivision_depth_>30)throw std::invalid_argument("invalid prepared robust Newmark configuration");
    if(strategy_==LinearStrategy::ModifiedNewton||strategy_==LinearStrategy::Adaptive)
        throw std::invalid_argument("prepared robust driver does not support ModifiedNewton/Adaptive");
}

PreparedRobustNewmark::~PreparedRobustNewmark()=default;
PreparedRobustNewmark::PreparedRobustNewmark(PreparedRobustNewmark&&) noexcept=default;
PreparedRobustNewmark& PreparedRobustNewmark::operator=(PreparedRobustNewmark&&) noexcept=default;

AnalysisResult PreparedRobustNewmark::run(const std::vector<double>& ground_accel,const RobustNewmarkOptions& options){
    if(options.max_subdivisions>prepared_subdivision_depth_)
        throw std::invalid_argument("prepared robust driver max_subdivisions exceeds prepared depth");
    // Validate only contexts that have actually been used. Subdivision depths
    // are prepared lazily, avoiding unnecessary factorization of half/quarter
    // steps on records that converge at the parent dt.
    bool valid=true;
    for(const auto& [depth,base]:impl_->prepared_effective_initial){
        const double d=dt_/std::pow(2.0,depth);const double a0=1.0/(kBeta*d*d),a1=kGamma/(kBeta*d);
        if(!same_sparse_numeric_matrix(model_->effective_initial_matrix(a0,a1),base)){valid=false;break;}
    }
    const bool had_preparation=!impl_->contexts.empty();
    if(!valid){impl_->contexts.clear();impl_->prepared_effective_initial.clear();++impl_->preparation_count;}
    else if(!had_preparation)++impl_->preparation_count;
    impl_->last_reused=valid&&had_preparation;

    auto getter=[&](int depth)->SolverContext&{
        auto it=impl_->contexts.find(depth);if(it!=impl_->contexts.end())return *it->second;
        const double d=dt_/std::pow(2.0,depth);const double a0=1.0/(kBeta*d*d),a1=kGamma/(kBeta*d);
        AnalysisStats setup;auto p=std::make_unique<SolverContext>(*model_,d,strategy_,setup);auto* raw=p.get();
        impl_->setup_factorizations+=setup.global_factorizations;
        impl_->prepared_effective_initial.emplace(depth,model_->effective_initial_matrix(a0,a1));
        impl_->contexts.emplace(depth,std::move(p));return *raw;
    };
    auto cached=[&](){std::size_t n=0;for(const auto& [depth,ctx]:impl_->contexts){(void)depth;n=std::max(n,ctx->cached_columns());}return n;};
    return run_newmark_robust_core(*model_,ground_accel,dt_,strategy_,options,getter,cached);
}

std::size_t PreparedRobustNewmark::preparation_count() const{return impl_->preparation_count;}
std::size_t PreparedRobustNewmark::setup_factorizations() const{return impl_->setup_factorizations;}
bool PreparedRobustNewmark::last_run_reused_preparation() const{return impl_->last_reused;}

} // namespace quake
