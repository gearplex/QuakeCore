#pragma once

#include "quake/frame2d.hpp"
#include "quake/frame3d.hpp"
#include "quake/modal.hpp"
#include "quake/newmark.hpp"
#include "quake/rc_building_asce41.hpp"
#include "quake/rc_frame_asce41.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace quake {

// One constitutive alternative applied to an otherwise identical compiled
// building. replacements is intentionally a delta rather than a rebuilt model;
// the comparison harness restores the affected components to their captured
// baseline before applying each variant.
struct ModelComparisonVariant {
    std::string name;
    std::string provenance;
    std::vector<std::pair<int,NonlinearMaterial>> replacements;
};

// Convenience adapter from a whole-building RC code-model field to a comparison
// variant. Each physical column model is assigned to all hinge IDs in its native
// frame binding. Unbound beams/joints retain the captured baseline material.
ModelComparisonVariant make_rc_column_model_comparison_variant(
    std::string name,
    std::string provenance,
    const std::vector<RCColumnElementBinding>& bindings,
    const std::vector<RCBuildingColumnModel>& models);

struct ModelComparisonComponent {
    std::string label;
    int component_id{};
};

struct ModelComparisonOptions {
    // Output-step histories are aligned across variants by the common input dt.
    bool record_story_history{true};
    bool record_component_history{true};
    // Exact envelopes additionally observe rollback-safe accepted subdivision
    // states so local peaks between record samples are not missed.
    bool include_accepted_substeps_in_envelopes{true};
    // Tangent modes are evaluated after the analysis at the committed states
    // that produced each story's peak |drift|. Duplicate peak steps are merged.
    bool tangent_modes_at_story_peaks{true};
    int tangent_mode_count{3};
    std::vector<std::size_t> additional_tangent_mode_steps;
    // Optional native column-force histories. Bindings use the same physical
    // column-to-element map as the ASCE demand coordinator.
    std::vector<RCColumnElementBinding> column_force_bindings;
};

struct ModelComparisonColumnSample {
    std::size_t step{};
    double time{};
    double max_compression_kip{};
    double max_tension_kip{};
    double max_abs_shear_kip{};
};

struct ModelComparisonColumnResult {
    std::string component_id;
    std::vector<ModelComparisonColumnSample> history;
    double max_compression_kip{};
    double time_at_max_compression{};
    double max_tension_kip{};
    double time_at_max_tension{};
    double max_abs_shear_kip{};
    double time_at_max_abs_shear{};
};

struct ModelComparisonComponentSample {
    std::size_t step{};
    double time{};
    double deformation{};
    double force{};
    double tangent{};
    MaterialEvalDiagnostics diagnostics{};
};

struct ModelComparisonComponentResult {
    std::string label;
    int component_id{};
    std::vector<ModelComparisonComponentSample> history;
    double max_abs_deformation{};
    double time_at_max_abs_deformation{};
    double max_abs_force{};
    double time_at_max_abs_force{};
    double minimum_tangent{};
    double time_at_minimum_tangent{};
    double first_io_time{-1.0};
    double first_ls_time{-1.0};
    double first_cp_time{-1.0};
    double first_lateral_loss_time{-1.0};
    double first_failure_time{-1.0};
};

struct ModelComparisonStorySample {
    std::size_t step{};
    double time{};
    std::vector<double> floor_displacement;
    std::vector<double> story_drift_ratio;
};

struct ModelComparisonTangentMode {
    double period{};
    double eigenvalue{};
    std::vector<double> reduced_shape;
    std::vector<double> story_shape;
};

struct ModelComparisonTangentSnapshot {
    std::size_t step{};
    double time{};
    std::string trigger;
    std::vector<ModelComparisonTangentMode> modes;
    // Static condensation of the same committed structural tangent onto the
    // configured floor-response coordinates, then transformed to interstory
    // displacement coordinates. Matrices are row-major.
    bool story_tangent_available{false};
    std::vector<double> floor_tangent_matrix;
    std::vector<double> interstory_tangent_matrix;
    std::vector<double> interstory_tangent_diagonal;
};

struct ModelComparisonVariantResult {
    std::string name;
    std::string provenance;
    AnalysisResult analysis;
    bool reused_prepared_solver{};
    std::size_t preparation_count_after_run{};

    std::vector<ModelComparisonStorySample> story_history;
    std::vector<double> peak_abs_story_drift_ratio;
    std::vector<double> time_at_peak_abs_story_drift;
    std::vector<double> residual_story_drift_ratio;
    double peak_abs_roof_response{};
    double time_at_peak_abs_roof_response{};
    double residual_roof_response{};

    std::vector<ModelComparisonComponentResult> components;
    std::vector<ModelComparisonColumnResult> columns;
    std::vector<ModelComparisonTangentSnapshot> tangent_snapshots;
};

struct ModelComparisonResult {
    std::vector<ModelComparisonVariantResult> variants;
    std::size_t prepared_solver_preparations{};
    std::size_t prepared_solver_setup_factorizations{};
};

class RCFrameModelComparison {
public:
    static ModelComparisonResult run_frame2d(
        CompiledFrame2D& frame,
        const std::vector<ModelComparisonVariant>& variants,
        const std::vector<ModelComparisonComponent>& monitored_components,
        const RCNativeFrameAnalysisSettings& settings,
        const ModelComparisonOptions& options = {});

    static ModelComparisonResult run_frame3d(
        CompiledFrame3D& frame,
        const std::vector<ModelComparisonVariant>& variants,
        const std::vector<ModelComparisonComponent>& monitored_components,
        const RCNativeFrameAnalysisSettings& settings,
        const ModelComparisonOptions& options = {});
};

} // namespace quake
