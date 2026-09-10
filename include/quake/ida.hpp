#pragma once
#include "quake/newmark.hpp"
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace quake {

struct GroundMotionRecord {
    std::string name;
    std::vector<double> acceleration;
    double dt{};
};

struct IDARunResult {
    std::size_t record_index{};
    std::string record_name;
    double scale_factor{};
    double pga{};
    AnalysisTermination termination{AnalysisTermination::Completed};
    CollapseMechanism collapse_mechanism{CollapseMechanism::None};
    std::size_t termination_step{};
    double termination_time{};
    double max_story_drift_ratio{};
    double max_roof_abs{};
    std::size_t newton_iterations{};
    std::size_t global_factorizations{};
    std::size_t direct_fallbacks{};
    double elapsed_seconds{};
    std::string termination_reason;
};

struct IDARecordSummary {
    std::size_t record_index{};
    std::string record_name;
    double base_pga{};
    bool collapsed{false};
    bool right_censored{false};
    bool numerical_failure{false};
    bool bracket_has_numerical_gap{false};
    bool nonmonotonic_response{false};
    double last_noncollapse_scale{0.0};
    double first_collapse_scale{std::numeric_limits<double>::quiet_NaN()};
    double collapse_pga{std::numeric_limits<double>::quiet_NaN()};
    CollapseMechanism collapse_mechanism{CollapseMechanism::None};
};

struct IDAOptions {
    std::vector<double> scale_factors{0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0};
    LinearStrategy strategy{LinearStrategy::Woodbury};
    RobustNewmarkOptions analysis{};
    int workers{1};
    bool stop_after_first_collapse{true};
    // Refine the first noncollapse/collapse bracket in log scale.
    int collapse_refinement_steps{3};
    // One independently owned context per record; reuse it across scale runs.
    bool reuse_preparation{true};
};

struct IDASuiteResult {
    std::vector<IDARunResult> runs;
    std::vector<IDARecordSummary> records;
    double elapsed_seconds{};
    int workers{};
    std::size_t physical_collapses{};
    std::size_t numerical_failures{};
    std::size_t right_censored_records{};
    // Simple lognormal capacity fit over uncensored physical-collapse PGAs.
    // These are descriptive sample estimates, not a censored-data MLE.
    double collapse_pga_median{std::numeric_limits<double>::quiet_NaN()};
    double collapse_pga_log_stddev{std::numeric_limits<double>::quiet_NaN()};
};

IDASuiteResult run_ida_suite(const NonlinearDynamicModel& model,
                             const std::vector<GroundMotionRecord>& records,
                             const IDAOptions& options={});

} // namespace quake
