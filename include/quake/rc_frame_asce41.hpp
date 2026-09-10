#pragma once

#include "quake/frame2d.hpp"
#include "quake/frame3d.hpp"
#include "quake/newmark.hpp"
#include "quake/rc_building_asce41.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace quake {

// Maps one physical RC column to the elastic member segment(s) that carry its
// section forces in a compiled frame. Multiple elastic elements are permitted
// for segmented physical columns; demand envelopes are taken over all segments.
// Hinge IDs are optional but, when supplied, are validated at model creation so
// a stale application-level binding fails before the NLRH starts.
struct RCColumnElementBinding {
    std::string component_id;
    std::vector<int> elastic_element_ids;
    std::vector<int> hinge_component_ids;
};

// Native transient settings shared by the 2D and 3D adapters. Robust Newmark is
// used deliberately because the committed-substep observer is rollback-safe.
struct RCNativeFrameAnalysisSettings {
    std::vector<double> ground_accel;
    double dt{};
    LinearStrategy strategy{LinearStrategy::SamePatternRefactorization};
    RobustNewmarkOptions newmark_options{};
};

struct RCNativeFrame2DBuild {
    CompiledFrame2D model;
    std::vector<RCColumnElementBinding> bindings;
};

struct RCNativeFrame3DBuild {
    CompiledFrame3D model;
    std::vector<RCColumnElementBinding> bindings;
};

// Phase 9D.5 production path: topology is compiled once and the ASCE material
// field is replaced in-place between global analyses. The model state itself is
// still reinitialized for each NRHA, but nodes/elements/MPCs/sparse patterns and
// nonlinear basis directions are not rebuilt.
struct RCNativeCompiledFrame2D {
    CompiledFrame2D model;
    std::vector<RCColumnElementBinding> bindings;
};

struct RCNativeCompiledFrame3D {
    CompiledFrame3D model;
    std::vector<RCColumnElementBinding> bindings;
};

using RCNativeFrame2DFactory = std::function<RCNativeFrame2DBuild(
    const std::vector<RCBuildingColumnModel>& models,std::size_t iteration)>;
using RCNativeFrame3DFactory = std::function<RCNativeFrame3DBuild(
    const std::vector<RCBuildingColumnModel>& models,std::size_t iteration)>;

// Exact envelope recorders for committed states. They are intentionally small
// state objects so applications can also use them outside the ASCE coordinator.
class RCFrame2DColumnDemandRecorder {
public:
    RCFrame2DColumnDemandRecorder(const CompiledFrame2D& model,
                                  std::vector<RCColumnElementBinding> bindings);
    void observe(const std::vector<double>& displacement,
                 const std::vector<double>& committed_state);
    std::vector<RCBuildingColumnDemandObservation> observations() const;
    std::size_t committed_samples() const { return committed_samples_; }

private:
    const CompiledFrame2D* model_{};
    std::vector<RCColumnElementBinding> bindings_;
    std::vector<RCColumnDemandState> demands_;
    std::size_t committed_samples_{};
};

class RCFrame3DColumnDemandRecorder {
public:
    RCFrame3DColumnDemandRecorder(const CompiledFrame3D& model,
                                  std::vector<RCColumnElementBinding> bindings);
    void observe(const std::vector<double>& displacement,
                 const std::vector<double>& committed_state);
    std::vector<RCBuildingColumnDemandObservation> observations() const;
    std::size_t committed_samples() const { return committed_samples_; }

private:
    const CompiledFrame3D* model_{};
    std::vector<RCColumnElementBinding> bindings_;
    std::vector<RCColumnDemandState> demands_;
    std::size_t committed_samples_{};
};

// Native building coordinator adapters. The application still owns geometry
// construction through the factory, but it no longer owns response probing,
// accepted-step bookkeeping, rollback filtering, or P/V envelope assembly.
class RCBuildingASCE41NativeFrame {
public:
    // Compiled-topology overloads. The supplied model is mutated only in its
    // nonlinear material field and is left carrying the final generated field.
    static RCBuildingASCE41IterationResult run_compiled_frame2d_two_pass_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        RCNativeCompiledFrame2D& frame,
        const RCNativeFrameAnalysisSettings& settings);

    static RCBuildingASCE41IterationResult run_compiled_frame2d_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        RCNativeCompiledFrame2D& frame,
        const RCNativeFrameAnalysisSettings& settings,
        const RCColumnAxialIterationOptions& iteration_options = {});

    static RCBuildingASCE41IterationResult run_compiled_frame3d_two_pass_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        RCNativeCompiledFrame3D& frame,
        const RCNativeFrameAnalysisSettings& settings);

    static RCBuildingASCE41IterationResult run_compiled_frame3d_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        RCNativeCompiledFrame3D& frame,
        const RCNativeFrameAnalysisSettings& settings,
        const RCColumnAxialIterationOptions& iteration_options = {});

    static RCBuildingASCE41IterationResult run_frame2d_two_pass_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        const RCNativeFrame2DFactory& factory,
        const RCNativeFrameAnalysisSettings& settings);

    static RCBuildingASCE41IterationResult run_frame2d_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        const RCNativeFrame2DFactory& factory,
        const RCNativeFrameAnalysisSettings& settings,
        const RCColumnAxialIterationOptions& iteration_options = {});

    static RCBuildingASCE41IterationResult run_frame3d_two_pass_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        const RCNativeFrame3DFactory& factory,
        const RCNativeFrameAnalysisSettings& settings);

    static RCBuildingASCE41IterationResult run_frame3d_asce41_23_aci369_1_22(
        const std::vector<RCBuildingColumnDefinition>& columns,
        const RCColumnRulesResolver& resolver,
        const RCNativeFrame3DFactory& factory,
        const RCNativeFrameAnalysisSettings& settings,
        const RCColumnAxialIterationOptions& iteration_options = {});
};

} // namespace quake
