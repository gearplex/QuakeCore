#include "quake/newmark.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/analysis_input.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace quake {

static double norm_inf(const std::vector<double>& x) {
    double v = 0.0;
    for (double a : x) {if(!std::isfinite(a))return std::numeric_limits<double>::infinity();v = std::max(v, std::abs(a));}
    return v;
}

static AnalysisResult run_newmark_core(const NonlinearDynamicModel& model,
                                       const std::vector<double>& ground_accel,
                                       double dt, LinearStrategy strategy,
                                       LowRankSolverBase* prepared_woodbury,
                                       SuperLUSamePatternSolver* same_pattern_solver,
                                       double tolerance, int max_iterations,
                                       int adaptive_rank_limit = -1) {
    validate_transient_input(ground_accel,dt,tolerance,max_iterations);
    if(model.has_state_dependent_global_tangent() || model.has_additional_effective_low_rank_update())
        throw std::invalid_argument("state-dependent or additional low-rank tangents require run_newmark_robust / PreparedRobustNewmark");
    if (strategy == LinearStrategy::Woodbury && prepared_woodbury == nullptr)
        throw std::invalid_argument("Woodbury strategy requires prepared solver");
    if (strategy == LinearStrategy::SamePatternRefactorization && same_pattern_solver == nullptr)
        throw std::invalid_argument("same-pattern strategy requires prepared direct solver");
    if (strategy == LinearStrategy::Adaptive && (prepared_woodbury == nullptr || same_pattern_solver == nullptr || adaptive_rank_limit < 0))
        throw std::invalid_argument("adaptive strategy requires Woodbury/direct solvers and calibrated rank limit");
    const int n = model.dof();
    const int m = model.nonlinear_count();
    constexpr double beta = 0.25;
    constexpr double gamma = 0.5;
    const double a0 = 1.0 / (beta * dt * dt);
    const double a1 = gamma / (beta * dt);

    std::vector<double> u(static_cast<std::size_t>(n), 0.0);
    std::vector<double> v(static_cast<std::size_t>(n), 0.0);
    std::vector<double> acc(static_cast<std::size_t>(n), 0.0);
    std::vector<double> committed = model.initial_nonlinear_state();
    std::vector<double> factored_tangents(static_cast<std::size_t>(m));
    factored_tangents = model.initial_nonlinear_tangents();
    // A prepared adaptive solver may retain the final tangent from its last run.
    // Treat it as unknown until this run explicitly refactors it.
    if(strategy==LinearStrategy::Adaptive)std::fill(factored_tangents.begin(),factored_tangents.end(),std::numeric_limits<double>::quiet_NaN());

    AnalysisResult result;
    result.roof_history.reserve(ground_accel.size());
    const auto start = std::chrono::steady_clock::now();

    for (std::size_t step = 0; step < ground_accel.size(); ++step) {
        std::vector<double> u_pred(static_cast<std::size_t>(n));
        std::vector<double> v_pred(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            u_pred[static_cast<std::size_t>(i)] = u[static_cast<std::size_t>(i)] + dt*v[static_cast<std::size_t>(i)] +
                dt*dt*(0.5-beta)*acc[static_cast<std::size_t>(i)];
            v_pred[static_cast<std::size_t>(i)] = v[static_cast<std::size_t>(i)] + dt*(1.0-gamma)*acc[static_cast<std::size_t>(i)];
        }
        std::vector<double> u_trial = u_pred;
        std::vector<double> trial_states = committed;
        std::vector<double> tangents(static_cast<std::size_t>(m));
        bool converged = false;
        std::unique_ptr<SuperLUFactor> modified_factor;

        for (int iter = 0; iter < max_iterations; ++iter) {
            std::vector<double> fint;
            NonlinearEvalDiagnostics nd;
            model.internal_force_and_tangent_diagnostics(u_trial, committed, fint, tangents, trial_states, &nd);
            if(!all_finite(fint)||!all_finite(tangents)||!all_finite(trial_states))
                throw ConstitutiveIntegrationError("non-finite constitutive response in Newmark");
            result.stats.nonlinear_component_evaluations += nd.component_evaluations;
            result.stats.nonlinear_active_component_evaluations += nd.active_tangent_evaluations;
            result.stats.nonlinear_fast_path_evaluations += nd.fast_path_evaluations;
            result.stats.nonlinear_full_state_evaluations += nd.full_state_evaluations;
            result.stats.nonlinear_transition_events += nd.transition_events;
            result.stats.nonlinear_reversal_events += nd.reversal_events;
            result.stats.nonlinear_deterioration_events += nd.deterioration_events;
            result.stats.nonlinear_failure_events += nd.failure_events;
            result.stats.nonlinear_io_or_beyond += nd.io_or_beyond_evaluations;
            result.stats.nonlinear_ls_or_beyond += nd.ls_or_beyond_evaluations;
            result.stats.nonlinear_cp_or_beyond += nd.cp_or_beyond_evaluations;
            result.stats.nonlinear_beyond_cp += nd.beyond_cp_evaluations;
            result.stats.nonlinear_lateral_loss += nd.lateral_loss_evaluations;
            std::vector<double> a_trial(static_cast<std::size_t>(n));
            std::vector<double> v_trial(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                a_trial[static_cast<std::size_t>(i)] = a0 * (u_trial[static_cast<std::size_t>(i)] - u_pred[static_cast<std::size_t>(i)]);
                v_trial[static_cast<std::size_t>(i)] = v_pred[static_cast<std::size_t>(i)] + gamma*dt*a_trial[static_cast<std::size_t>(i)];
            }
            auto cv = model.damping_multiply(v_trial);
            auto ma = model.mass_multiply(a_trial);
            auto p = model.base_excitation(ground_accel[step]);
            std::vector<double> residual(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                residual[static_cast<std::size_t>(i)] = p[static_cast<std::size_t>(i)] - fint[static_cast<std::size_t>(i)] -
                    cv[static_cast<std::size_t>(i)] - ma[static_cast<std::size_t>(i)];
            }
            ++result.stats.newton_iterations;
            const double scale = std::max(1.0, norm_inf(p));
            if (!std::isfinite(norm_inf(residual))||!std::isfinite(scale))
                throw ConstitutiveIntegrationError("non-finite dynamic residual in Newmark");
            if (norm_inf(residual) <= tolerance * scale) {
                acc = std::move(a_trial);
                v = std::move(v_trial);
                converged = true;
                break;
            }

            std::vector<double> du;
            if (strategy == LinearStrategy::FullFactorization) {
                auto keff = model.effective_tangent_matrix(tangents, a0, a1);
                du = superlu_solve_once(keff, residual);
                ++result.stats.global_factorizations;
            } else if (strategy == LinearStrategy::SamePatternRefactorization) {
                bool same=true;
                for(int j=0;j<m;++j) if(tangents[static_cast<std::size_t>(j)]!=factored_tangents[static_cast<std::size_t>(j)]) { same=false; break; }
                if(same) {
                    du=same_pattern_solver->solve_current(residual);
                } else {
                    auto keff=model.effective_tangent_matrix(tangents,a0,a1);
                    du=same_pattern_solver->refactor_and_solve(keff,residual);
                    factored_tangents=tangents;
                    ++result.stats.global_factorizations;
                }
            } else if (strategy == LinearStrategy::ModifiedNewton) {
                if (!modified_factor) {
                    auto keff = model.effective_tangent_matrix(tangents, a0, a1);
                    modified_factor = std::make_unique<SuperLUFactor>(keff);
                    ++result.stats.global_factorizations;
                }
                du = modified_factor->solve(residual);
            } else {
                std::vector<double> delta_k(static_cast<std::size_t>(m));
                std::size_t active_rank = 0;
                for (int j = 0; j < m; ++j) {
                    delta_k[static_cast<std::size_t>(j)] = tangents[static_cast<std::size_t>(j)] -
                        model.initial_nonlinear_tangents()[static_cast<std::size_t>(j)];
                    if (std::abs(delta_k[static_cast<std::size_t>(j)]) > 1e-14) ++active_rank;
                }
                result.stats.active_rank_sum += active_rank;
                result.stats.max_active_rank = std::max(result.stats.max_active_rank, active_rank);
                bool direct_matches=true;
                if(strategy==LinearStrategy::Adaptive) for(int j=0;j<m;++j) if(tangents[static_cast<std::size_t>(j)]!=factored_tangents[static_cast<std::size_t>(j)]){direct_matches=false;break;}
                // At rank zero, an already-baseline direct factor avoids the
                // Woodbury active-set scan/reduced-system overhead. At high rank
                // use the calibrated direct path; intermediate ranks use Woodbury.
                const bool use_direct = strategy==LinearStrategy::Adaptive &&
                    ((active_rank==0 && direct_matches) || static_cast<int>(active_rank)>adaptive_rank_limit);
                if(use_direct){
                    if(direct_matches) du=same_pattern_solver->solve_current(residual);
                    else {auto keff=model.effective_tangent_matrix(tangents,a0,a1);du=same_pattern_solver->refactor_and_solve(keff,residual);factored_tangents=tangents;++result.stats.global_factorizations;}
                    ++result.stats.adaptive_direct_solves;
                } else {
                    ++result.stats.reduced_update_solves;
                    du = prepared_woodbury->solve(residual, delta_k);
                    if(strategy==LinearStrategy::Adaptive) ++result.stats.adaptive_woodbury_solves;
                }
            }
            ++result.stats.linear_solves;
            for (int i = 0; i < n; ++i) u_trial[static_cast<std::size_t>(i)] += du[static_cast<std::size_t>(i)];
        }

        if (!converged) {
            ++result.stats.failed_steps;
            throw std::runtime_error("Newmark step failed to converge at step " + std::to_string(step));
        }
        u = std::move(u_trial);
        committed = std::move(trial_states);
        const double response = model.response_value(u);
        result.roof_history.push_back(response);
        result.stats.max_roof_abs = std::max(result.stats.max_roof_abs, std::abs(response));
        result.stats.max_story_drift_abs = std::max(result.stats.max_story_drift_abs,
                                                    model.max_drift_measure(u));
        ++result.stats.steps;
    }
    result.stats.elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if (prepared_woodbury) result.stats.cached_low_rank_columns = prepared_woodbury->cached_update_count();
    result.final_displacement = u;
    result.final_state = committed;
    return result;
}

AnalysisResult run_newmark(const NonlinearDynamicModel& model,
                           const std::vector<double>& ground_accel,
                           double dt, LinearStrategy strategy,
                           double tolerance, int max_iterations) {
    if (strategy == LinearStrategy::SamePatternRefactorization) {
        validate_transient_input({0.0},dt,1e-8,1);
        constexpr double beta=0.25, gamma=0.5;
        const double a0=1.0/(beta*dt*dt), a1=gamma/(beta*dt);
        const auto setup_start=std::chrono::steady_clock::now();
        SuperLUSamePatternSolver solver(model.effective_initial_matrix(a0,a1));
        const double setup=std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_start).count();
        auto result=run_newmark_core(model,ground_accel,dt,strategy,nullptr,&solver,tolerance,max_iterations);
        result.stats.elapsed_seconds += setup;
        result.stats.global_factorizations += 1;
        return result;
    }
    if (strategy == LinearStrategy::Adaptive)
        throw std::invalid_argument("use PreparedAdaptiveNewmark with a calibrated rank limit");
    if (strategy != LinearStrategy::Woodbury)
        return run_newmark_core(model,ground_accel,dt,strategy,nullptr,nullptr,tolerance,max_iterations);

    validate_transient_input({0.0},dt,1e-8,1);
    constexpr double beta = 0.25;
    constexpr double gamma = 0.5;
    const double a0 = 1.0 / (beta * dt * dt);
    const double a1 = gamma / (beta * dt);
    const auto setup_start = std::chrono::steady_clock::now();
    LazyLowRankWoodburySolver solver(model.effective_initial_matrix(a0,a1),model.nonlinear_basis());
    const double setup = std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_start).count();
    auto result = run_newmark_core(model,ground_accel,dt,strategy,&solver,nullptr,tolerance,max_iterations);
    result.stats.elapsed_seconds += setup;
    result.stats.global_factorizations += 1;
    return result;
}

PreparedWoodburyNewmark::PreparedWoodburyNewmark(const NonlinearDynamicModel& model, double dt)
    : model_(&model), dt_(dt) {
    validate_transient_input({0.0},dt,1e-8,1);
    constexpr double beta = 0.25;
    constexpr double gamma = 0.5;
    const double a0 = 1.0 / (beta * dt * dt);
    const double a1 = gamma / (beta * dt);
    const auto start = std::chrono::steady_clock::now();
    solver_ = std::make_unique<LazyLowRankWoodburySolver>(model.effective_initial_matrix(a0,a1),
                                                           model.nonlinear_basis());
    setup_seconds_ = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}

AnalysisResult PreparedWoodburyNewmark::run(const std::vector<double>& ground_accel,
                                             double tolerance,
                                             int max_iterations) const {
    return run_newmark_core(*model_,ground_accel,dt_,LinearStrategy::Woodbury,solver_.get(),nullptr,tolerance,max_iterations);
}

PreparedAdaptiveNewmark::PreparedAdaptiveNewmark(const NonlinearDynamicModel& model,double dt,int rank_limit)
    :model_(&model),dt_(dt),rank_limit_(rank_limit){
    validate_transient_input({0.0},dt,1e-8,1);
    if(rank_limit<0)throw std::invalid_argument("invalid adaptive solver setup");
    constexpr double beta=.25,gamma=.5;const double a0=1.0/(beta*dt*dt),a1=gamma/(beta*dt);
    auto A=model.effective_initial_matrix(a0,a1);const auto t0=std::chrono::steady_clock::now();
    wood_=std::make_unique<LazyLowRankWoodburySolver>(A,model.nonlinear_basis());
    direct_=std::make_unique<SuperLUSamePatternSolver>(A);
    setup_seconds_=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
}
AnalysisResult PreparedAdaptiveNewmark::run(const std::vector<double>& gm,double tol,int maxit) const{
    return run_newmark_core(*model_,gm,dt_,LinearStrategy::Adaptive,wood_.get(),direct_.get(),tol,maxit,rank_limit_);
}

std::vector<double> synthetic_ground_motion(int steps, double dt, double amplitude) {
    std::vector<double> a(static_cast<std::size_t>(steps));
    constexpr double pi = 3.14159265358979323846;
    for (int i = 0; i < steps; ++i) {
        const double t = i * dt;
        const double duration = steps * dt;
        const double env = std::sin(pi * std::min(1.0, t / std::max(1e-9, duration)));
        const double carrier = 0.62*std::sin(2*pi*1.05*t) + 0.28*std::sin(2*pi*2.4*t+0.7) +
                               0.10*std::sin(2*pi*4.7*t+1.2);
        a[static_cast<std::size_t>(i)] = amplitude * env * carrier;
    }
    return a;
}

} // namespace quake
