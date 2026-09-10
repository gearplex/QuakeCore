#pragma once

#include "quake/dynamic_model.hpp"
#include <vector>
#include <cstddef>

namespace quake {

struct SolverCalibrationPoint {
    int active_rank{};
    double woodbury_seconds{};
    double same_pattern_refactor_seconds{};
    double speedup{}; // direct / Woodbury
};

struct SolverCalibrationResult {
    int dof{};
    int nonlinear_count{};
    std::size_t matrix_nnz{};
    std::vector<SolverCalibrationPoint> points;
    // Largest sampled rank for which Woodbury is faster. -1 means none;
    // nonlinear_count means no crossover was observed in the sampled range.
    int recommended_woodbury_rank_limit{-1};
    double recommended_rank_fraction{};
    double elapsed_seconds{};

    bool prefer_woodbury(int active_rank) const {
        return recommended_woodbury_rank_limit >= 0 && active_rank <= recommended_woodbury_rank_limit;
    }
};

// Calibrate tangent-refresh cost on the actual compiled model. The benchmark
// reuses Woodbury influence columns but forces a same-pattern numeric
// refactorization on every repeat, matching the expensive event when a tangent
// state changes during NRHA rather than the cheap unchanged-tangent solve.
SolverCalibrationResult calibrate_solver_crossover(
    const NonlinearDynamicModel& model, double dt,
    std::vector<int> ranks = {}, int repeats = 4,
    double tangent_ratio = 0.05);

} // namespace quake
