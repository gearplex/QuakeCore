#pragma once

#include "quake/asce41_hinge.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace quake {

// Edition is metadata/provenance, not a hidden switch for proprietary table
// values. QuakeCore intentionally does not embed copyrighted ASCE/ACI tables.
enum class RCColumnCodeEdition : int {
    NIST_ASCE41_17_Benchmark = 0,
    ASCE41_23_ACI369_1_22 = 1
};

// Physical/member information supplied to an edition-specific engineering
// rules resolver. Compression is positive. These fields deliberately describe
// the column rather than pre-resolving ASCE/ACI table parameters in the solver.
struct RCColumnSectionInput {
    std::string component_id;

    double width_in{};
    double depth_in{};
    double clear_length_in{};
    double expected_fc_ksi{};
    double expected_fy_long_ksi{};
    double expected_fy_trans_ksi{};
    double steel_modulus_ksi{29000.0};

    int longitudinal_bar_count{};
    double longitudinal_bar_area_in2{};
    double longitudinal_steel_area_in2{}; // optional; 0 => count * bar area

    double cover_to_tie_center_in{};
    double transverse_bar_area_in2{};
    double transverse_spacing_in{};
    int transverse_legs_x{};
    int transverse_legs_y{};

    // Analysis/design state supplied before the first nonlinear pass.
    double gravity_axial_compression_kip{};
    double initial_shear_demand_kip{};
    double expected_shear_capacity_kip{};

    // Detailing/condition flags are intentionally factual inputs. Their effect
    // is owned by the authorized edition-specific rules resolver.
    bool conforming_transverse_reinforcement{false};
    bool lap_splice_in_hinge_region{false};
    bool development_or_anchorage_critical{false};
};

// Analysis demand used to resolve axial/shear-sensitive modeling parameters.
// Compression and tension are separately retained so a future P-M surface can
// use the full range without changing this interface.
struct RCColumnDemandState {
    double max_compression_kip{};
    double max_tension_kip{};
    double max_abs_shear_kip{};
};

// Values resolved from an authorized code/table workflow or a benchmark source.
// All rotations are plastic hinge rotations. Strengths are positive magnitudes.
struct RCColumnResolvedParameters {
    double numerical_hinge_Ke{};
    double posMy{}, negMy{};
    double posMc{}, negMc{};
    double pos_a{}, neg_a{};
    double pos_b{}, neg_b{};
    double pos_c{}, neg_c{};
    double pos_io{}, pos_ls{}, pos_cp{};
    double neg_io{}, neg_ls{}, neg_cp{};

    // Optional separate effective/gravity failure deformation. Set to zero to
    // make F coincide with E, which is the NIST benchmark column convention.
    double pos_f{}, neg_f{};

    // Only used by ResearchExtendedCDE.
    double pos_drop_span{}, neg_drop_span{};
    double pos_e_drop_span{}, neg_e_drop_span{};

    // Cyclic policy remains separate from code-backbone generation.
    double lambda_strength{1e30};
    double lambda_unloading{1e30};
    double cyclic_exponent{1.0};
};

// Fine-grained provenance makes a generated hinge reviewable without requiring
// the nonlinear material itself to understand the governing standard.
struct RCColumnParameterAuditEntry {
    std::string parameter;
    double value{};
    std::string source_reference;
    std::string controlling_condition;
};

struct RCColumnRuleResolution {
    RCColumnResolvedParameters resolved;
    RCColumnDemandState demand_used;
    std::string provenance;
    std::vector<RCColumnParameterAuditEntry> audit;
};

struct RCColumnModelSpec {
    RCColumnCodeEdition edition{RCColumnCodeEdition::NIST_ASCE41_17_Benchmark};
    ASCE41HingeParams hinge;
    RCColumnDemandState demand_used;
    std::string provenance;
    std::vector<RCColumnParameterAuditEntry> audit;
    bool code_coefficients_embedded{false};
};

// Boundary between QuakeCore and an authorized implementation of ASCE 41 / ACI
// 369 rules. This can be backed by a licensed table service, a firm-controlled
// rules library, or a benchmark resolver. The solver never guesses coefficients.
class RCColumnRulesResolver {
public:
    virtual ~RCColumnRulesResolver() = default;
    virtual RCColumnRuleResolution resolve(const RCColumnSectionInput& section,
                                           const RCColumnDemandState& demand) const = 0;
};

// Lightweight adapter useful for applications that already have a rules engine
// and simply want to provide it as a callback.
class CallbackRCColumnRulesResolver final : public RCColumnRulesResolver {
public:
    using Function = std::function<RCColumnRuleResolution(
        const RCColumnSectionInput&, const RCColumnDemandState&)>;

    explicit CallbackRCColumnRulesResolver(Function fn);
    RCColumnRuleResolution resolve(const RCColumnSectionInput& section,
                                   const RCColumnDemandState& demand) const override;

private:
    Function fn_;
};

class RCColumnASCE41Provider {
public:
    // Reconstruct the NIST ASCE 41-17 benchmark topology from already-resolved
    // component parameters: explicit B/C strengths, straight C->E, and F=E.
    static RCColumnModelSpec nist_asce41_17_benchmark(
        const RCColumnResolvedParameters& resolved,
        std::string provenance = {});
    static RCColumnModelSpec nist_asce41_17_benchmark(
        const RCColumnRuleResolution& resolution);

    // Production ASCE 41-23 / ACI 369.1-22 adapter. The caller must supply the
    // resolved modeling and acceptance parameters from an authorized source.
    // No table coefficients are guessed or silently embedded here.
    static RCColumnModelSpec asce41_23_aci369_1_22(
        const RCColumnResolvedParameters& resolved,
        ASCE41BackboneShape topology,
        std::string provenance);
    static RCColumnModelSpec asce41_23_aci369_1_22(
        const RCColumnRuleResolution& resolution,
        ASCE41BackboneShape topology);

private:
    static ASCE41HingeParams common_params(const RCColumnResolvedParameters& r);
};

struct RCColumnAxialIterationOptions {
    std::size_t max_iterations{5};
    std::size_t min_iterations{2};
    double demand_relative_tolerance{0.01};
    double parameter_relative_tolerance{0.01};
    // 1.0 uses the newly observed demand directly; lower values under-relax
    // difficult demand/parameter fixed-point loops.
    double relaxation{1.0};
};

struct RCColumnAxialIterationStep {
    std::size_t iteration{};
    RCColumnDemandState demand_input;
    RCColumnDemandState demand_observed;
    RCColumnModelSpec model;
    double demand_relative_change{};
    double parameter_relative_change{};
};

struct RCColumnAxialIterationResult {
    bool converged{false};
    std::vector<RCColumnAxialIterationStep> steps;
    RCColumnModelSpec final_model;
    RCColumnDemandState final_observed_demand;
};

// Generic response callback. The application runs its NLRH/NDP analysis using
// the supplied model and returns the column's observed demand envelope.
using RCColumnResponseRunner = std::function<RCColumnDemandState(
    const RCColumnModelSpec&, std::size_t iteration)>;

// Practical axial-demand-sensitive workflow used before a full moving-surface
// P-M implementation is available. It is a fixed-point iteration:
//   demand -> code parameters -> NLRH -> observed demand -> regenerated params.
// The workflow is solver-agnostic and therefore testable without embedding a
// particular frame-analysis driver here.
class RCColumnAxialIteration {
public:
    // Convenience wrapper for the practical two-analysis workflow: pass 1 uses
    // gravity/preanalysis demand, pass 2 regenerates the hinge at the maximum
    // compression observed in pass 1. The result still reports whether pass 2
    // is self-consistent, but it always stops after the second analysis.
    static RCColumnAxialIterationResult run_two_pass_asce41_23_aci369_1_22(
        const RCColumnSectionInput& section,
        const RCColumnRulesResolver& resolver,
        const RCColumnResponseRunner& response_runner,
        ASCE41BackboneShape topology);

    static RCColumnAxialIterationResult run_asce41_23_aci369_1_22(
        const RCColumnSectionInput& section,
        const RCColumnRulesResolver& resolver,
        const RCColumnResponseRunner& response_runner,
        ASCE41BackboneShape topology,
        const RCColumnAxialIterationOptions& options = {});
};

} // namespace quake
