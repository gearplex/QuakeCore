#pragma once

#include "quake/asce41_hinge.hpp"

#include <vector>

namespace quake {

// Phase 9F two-force P-M hinge. The hinge owns a zero-length axial deformation
// and a rotational deformation. Compression is positive in P-space, while the
// generalized axial hinge force N is positive in tension, so P=P0-N.
//
// The yield surface is |M| = Y(P,kappa), where the axial-force-dependent base
// moment capacity is supplied as a monotone table and the ASCE41 hinge
// parameters govern directional post-yield hardening/capping/strength loss as
// a function of accumulated plastic rotation on each bending side.
//
// Plastic flow is associative in the two generalized force coordinates [N,M].
// Therefore flexural yielding can generate axial plastic deformation as well as
// rotation, unlike the legacy one-way P->M capacity wrapper.
enum class PMInteractionSurfaceShape {
    // Legacy Phase-9F development surface supplied as My(P) points.
    TabulatedMomentCapacity,
    // PERFORM concrete-type surface in a P-M plane:
    // |(P-PB)/(PY0-PB)|^alpha + |M/MYB|^beta = 1.
    PerformConcrete
};

enum class PMInteractionSurfaceEvolution {
    // Capacity follows the supplied ASCE41 directional backbone as plastic
    // rotation accumulates. This is the Phase-9F research/degrading mode.
    BackboneScaled,
    // The axial-force-dependent P-M surface is fixed after first yield. This
    // isolates Powell/Perform-style associative P-M mechanics without mixing
    // in post-capping surface shrinkage.
    ElasticPerfectlyPlastic,
    // Phase 9F.2: two-surface kinematic hardening patterned after PERFORM's
    // documented Mroz formulation.  A same-shape inner Y surface translates
    // toward a fixed, homothetic U surface; when the two become tangent the
    // response proceeds EPP on U.  Translation follows the Mroz conjugate-
    // point direction.  The prescribed U-point plastic rotation calibrates
    // the translation rate under proportional bending.
    MrozTwoSurface
};

struct PMInteractionHingeParams {
    ASCE41HingeParams hinge;
    double axial_preload{};          // compression-positive total P at zero incremental deformation
    double axial_stiffness{};        // zero-length hinge elastic Ka (force/length)
    std::vector<double> axial_force_points;      // compression-positive P
    std::vector<double> moment_capacity_points;  // positive yield moment My(P)
    double return_tolerance{1e-10};
    int max_return_iterations{30};
    PMInteractionSurfaceEvolution surface_evolution{PMInteractionSurfaceEvolution::BackboneScaled};
    PMInteractionSurfaceShape surface_shape{PMInteractionSurfaceShape::TabulatedMomentCapacity};
    // Parameters for the PERFORM concrete-type P-M surface. Compression is
    // positive.  py_tension should be below p_balance and py_compression above.
    double p_balance{};
    double py_tension{};
    double py_compression{};
    double my_balance{};
    double alpha_tension{1.5};
    double alpha_compression{1.5};
    double beta_pm{1.1};
    // When false, a/b/f remain available as diagnostic rotation thresholds
    // but do not remove lateral resistance. Used by the Phase-9F.1 mechanics
    // isolation so the experiment tests P-M coupling rather than collapse
    // semantics.
    bool enforce_deformation_capacity{true};
    // U/Y scale for the Mroz outer surface. <=1 derives the scale from the
    // directional point-C strength (or B-C physical tangent) and a-value.
    double mroz_outer_scale{0.0};
    // Reproduction only: Phase 9F.2 measured distances in mixed P/M units.
    // New analyses use dimensionless coordinates based on the surface scales.
    bool legacy_mroz_mixed_unit_metric{false};
};

enum PMInteractionEvent : unsigned {
    PM_EVENT_NONE          = 0u,
    PM_EVENT_ELASTIC       = 1u << 0,
    PM_EVENT_YIELD         = 1u << 1,
    PM_EVENT_CAP           = 1u << 2,
    PM_EVENT_STRENGTH_DROP = 1u << 3,
    PM_EVENT_LATERAL_LOSS  = 1u << 4,
    PM_EVENT_FAILURE       = 1u << 5,
    PM_EVENT_RETURN_FAIL   = 1u << 6
};

struct PMInteractionTrialResult {
    double axial_force{};       // generalized N, positive tension (incremental about preload)
    double compression{};       // total compression-positive P = P0-N
    double moment{};
    // Row-major generalized tangent d[N,M]/d[axial_deformation,rotation].
    double tangent[4]{};
    double plastic_axial_deformation{};
    double plastic_rotation{};
    double positive_plastic_rotation{};
    double negative_plastic_rotation{};
    double active_capacity{};
    unsigned events{PM_EVENT_NONE};
    bool converged{true};
};

class PMInteractionHinge2D {
public:
    // ep_axial, ep_rotation, kappa_pos, kappa_neg, last_N, last_M,
    // lateral_loss, failure, last_capacity, Mroz center-P, center-M,
    // Mroz outer-surface flag.
    static constexpr int kStateSize = 12;

    explicit PMInteractionHinge2D(PMInteractionHingeParams params);

    int state_size() const { return kStateSize; }
    double axial_initial_stiffness() const { return p_.axial_stiffness; }
    double rotational_initial_stiffness() const { return p_.hinge.Ke; }
    const PMInteractionHingeParams& params() const { return p_; }

    void initialize_state(double* state) const;
    PMInteractionTrialResult trial(double axial_deformation,double rotation,
                                   const double* committed,double* trial_state) const;

private:
    struct CapacityEval { double value{},d1{},d2{}; };
    struct BackboneEval {
        double value{},dP{},dPP{},dKappa{},dPdKappa{};
        unsigned events{};
        bool lateral_loss{};
    };

    CapacityEval base_capacity(double compression) const;
    CapacityEval base_capacity_scaled(double compression,double scale) const;
    BackboneEval directional_capacity(double compression,double kappa,bool positive) const;
    PMInteractionTrialResult trial_mroz(double axial_deformation,double rotation,
                                        const double* committed,double* trial_state,
                                        bool compute_tangent) const;

    PMInteractionHingeParams p_;
    std::vector<double> capacity_slopes_;
};

} // namespace quake
