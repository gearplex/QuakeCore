#include "quake/static_analysis.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace quake {
namespace {
double inf_norm(const std::vector<double>& x) {
    double out = 0.0;
    for (double v : x) out = std::max(out, std::abs(v));
    return out;
}
}

StaticAnalysisResult solve_static_load(const NonlinearDynamicModel& model,
                                        const std::vector<double>& load,
                                        int load_steps,
                                        double tolerance,
                                        int max_iterations) {
    if (load.size() != static_cast<std::size_t>(model.dof()))
        throw std::invalid_argument("static load size mismatch");
    if (load_steps < 1 || max_iterations < 1 || !std::isfinite(tolerance) || tolerance <= 0.0)
        throw std::invalid_argument("invalid static-analysis controls");
    StaticAnalysisResult result;
    result.displacement.assign(static_cast<std::size_t>(model.dof()), 0.0);
    result.committed_state = model.initial_nonlinear_state();
    result.tangents = model.initial_nonlinear_tangents();
    for (int step = 1; step <= load_steps; ++step) {
        const double factor_load = static_cast<double>(step) / load_steps;
        bool converged = false;
        std::vector<double> trial_state;
        for (int iter = 0; iter < max_iterations; ++iter) {
            std::vector<double> internal;
            model.internal_force_and_tangent(result.displacement, result.committed_state,
                                             internal, result.tangents, trial_state);
            std::vector<double> residual(load.size());
            for (std::size_t i = 0; i < load.size(); ++i)
                residual[i] = factor_load * load[i] - internal[i];
            result.residual_inf_norm = inf_norm(residual);
            ++result.newton_iterations;
            if (result.residual_inf_norm <= tolerance) {
                result.committed_state = std::move(trial_state);
                converged = true;
                break;
            }
            auto tangent = model.effective_state_tangent_matrix_with_state(
                result.displacement, result.tangents, result.committed_state, 0.0, 0.0);
            const auto du = superlu_solve_once(tangent, residual);
            for (std::size_t i = 0; i < du.size(); ++i) result.displacement[i] += du[i];
        }
        if (!converged) return result;
        result.load_steps_completed = step;
    }
    result.converged = true;
    return result;
}

}  // namespace quake
