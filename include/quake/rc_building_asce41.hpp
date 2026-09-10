#pragma once

#include "quake/rc_column_asce41.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace quake {

// One physical RC column participating in a building-level code-model update.
// component_id must be unique within a coordinated run.
struct RCBuildingColumnDefinition {
    RCColumnSectionInput section;
    ASCE41BackboneShape topology{ASCE41BackboneShape::StraightCE};
};

struct RCBuildingColumnModel {
    std::string component_id;
    RCColumnModelSpec model;
};

struct RCBuildingColumnDemandObservation {
    std::string component_id;
    RCColumnDemandState demand;
};

// Returned by one complete building analysis. The analysis driver is responsible
// for harvesting the P/V envelopes of every coordinated RC column from the same
// global response history. The coordinator verifies exact component coverage.
struct RCBuildingAnalysisObservation {
    bool analysis_succeeded{true};
    // A physical collapse is a valid structural outcome, but it terminates the
    // code-parameter fixed-point loop because the full record demand envelope
    // was not completed. Keep it distinct from numerical analysis failure.
    bool physical_collapse{false};
    std::string status;
    std::vector<RCBuildingColumnDemandObservation> column_demands;
};

using RCBuildingResponseRunner = std::function<RCBuildingAnalysisObservation(
    const std::vector<RCBuildingColumnModel>& models,
    std::size_t iteration)>;

struct RCBuildingColumnIterationMetric {
    std::string component_id;
    RCColumnDemandState demand_input;
    RCColumnDemandState demand_observed;
    double demand_relative_change{};
    double parameter_relative_change{};
};

struct RCBuildingASCE41IterationStep {
    std::size_t iteration{};
    std::vector<RCBuildingColumnModel> models;
    std::vector<RCBuildingColumnIterationMetric> component_metrics;
    double global_demand_relative_change{};
    double global_parameter_relative_change{};
    std::string analysis_status;
};

enum class RCBuildingASCE41Termination : int {
    Converged = 0,
    IterationLimit = 1,
    AnalysisFailure = 2,
    PhysicalCollapse = 3
};

struct RCBuildingASCE41IterationResult {
    bool converged{false};
    RCBuildingASCE41Termination termination{RCBuildingASCE41Termination::IterationLimit};
    std::string termination_detail;
    std::vector<RCBuildingASCE41IterationStep> steps;
    std::vector<RCBuildingColumnModel> final_models;
    std::vector<RCBuildingColumnDemandObservation> final_observed_demands;
};

// Whole-building synchronous code-model coordinator. All columns are resolved,
// analyzed and regenerated as one parameter field. No column is updated using
// another column's demand from the same iteration, making the update invariant
// to input ordering and suitable for deterministic validation/regression work.
class RCBuildingASCE41Iteration {
public:
    static RCBuildingASCE41IterationResult run_two_pass_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        const RCBuildingResponseRunner& response_runner);

    static RCBuildingASCE41IterationResult run_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        const RCBuildingResponseRunner& response_runner,
        const RCColumnAxialIterationOptions& options = {});
};

} // namespace quake
