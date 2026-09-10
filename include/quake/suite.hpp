#pragma once

#include "quake/newmark.hpp"

#include <vector>

namespace quake {

struct SuiteRunResult {
    std::vector<AnalysisResult> records;
    double elapsed_seconds{};
    int workers{};
    std::size_t setup_factorizations{};
};

// Shared-memory ground-motion parallelism. Each worker owns a prepared
// Woodbury solver, so lazy influence caches are never shared/mutated across
// threads. Model data are read-only.
SuiteRunResult run_record_suite_parallel(const NonlinearDynamicModel& model,
                                         const std::vector<std::vector<double>>& ground_motions,
                                         double dt,
                                         int workers,
                                         double tolerance = 1e-8,
                                         int max_iterations = 20);

} // namespace quake
