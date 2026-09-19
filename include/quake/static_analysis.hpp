#pragma once

#include "quake/dynamic_model.hpp"
#include <cstddef>
#include <vector>

namespace quake {

struct StaticAnalysisResult {
    bool converged{false};
    int load_steps_completed{0};
    std::size_t newton_iterations{0};
    double residual_inf_norm{0.0};
    std::vector<double> displacement;
    std::vector<double> committed_state;
    std::vector<double> tangents;
};

StaticAnalysisResult solve_static_load(const NonlinearDynamicModel& model,
                                        const std::vector<double>& load,
                                        int load_steps = 100,
                                        double tolerance = 1e-8,
                                        int max_iterations = 100);

}  // namespace quake
