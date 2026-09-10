#pragma once

namespace quake {

// Envelope topology is deliberately separated from component parameter
// generation.  The legacy/research formulation retains an explicit C-D
// degradation ramp, residual branch, and final drop to E.  StraightCE is the
// NIST benchmark interpretation used for its ASCE 41-17 Perform3D/OpenSees
// comparisons: once point C is reached, strength degrades linearly to zero at E.
enum class ASCE41BackboneShape : int {
    ResearchExtendedCDE = 0,
    StraightCE = 1
};

// Research implementation of the generalized ASCE 41 / FEMA component
// force-deformation envelope.  The user supplies the component-specific
// modeling/acceptance parameters; this class does not embed copyrighted
// tabulated values from the standard.
//
// A-B: elastic to yield
// B-C: strain hardening to the capping point
// C-D: strength degradation (finite width for numerical robustness)
// D-E: residual-strength branch
// E-F: zero lateral/seismic strength while gravity/effective resistance may remain
// beyond F: effective/gravity component failure
//
// a and b are plastic deformations measured from effective yield. c is the
// residual-strength ratio.  The finite CD transition is controlled by
// drop_span and must fit between a and b.
struct ASCE41HingeParams {
    double Ke{};
    double posFy{}, negFy{};            // positive magnitudes at B (effective yield)
    // Optional explicit point-C strengths. <=0 derives C from the B-C slope.
    // This lets a code/rules provider define the backbone independently of the
    // numerical zero-length hinge penalty stiffness.
    double posMc{}, negMc{};
    double hardening_ratio{0.02};       // legacy B-C slope / Ke when hardening_stiffness <= 0
    // Optional physical B-C tangent in force/deformation units.  This is
    // deliberately independent of Ke so penalty-stiff zero-length hinges do
    // not amplify physical post-yield hardening.  <=0 preserves the legacy
    // hardening_ratio*Ke behavior for backward compatibility.
    double hardening_stiffness{-1.0};
    double pos_a{}, neg_a{};            // plastic deformation at C
    double pos_b{}, neg_b{};            // plastic deformation at E (loss of seismic/lateral resistance)
    double pos_f{}, neg_f{};            // optional plastic deformation at F (gravity/effective component loss); 0 => E
    double pos_c{}, neg_c{};            // residual strength / Fy
    ASCE41BackboneShape backbone_shape{ASCE41BackboneShape::ResearchExtendedCDE};
    double pos_drop_span{}, neg_drop_span{}; // C-D deformation width; ignored by StraightCE
    double pos_e_drop_span{}, neg_e_drop_span{}; // final residual-to-zero ramp immediately before E; 0 => automatic

    // Optional plastic-deformation acceptance limits. <=0 disables a limit.
    double pos_io{}, pos_ls{}, pos_cp{};
    double neg_io{}, neg_ls{}, neg_cp{};

    // Optional energy-based cyclic degradation.  Very large lambda values
    // effectively disable cyclic deterioration while preserving envelope
    // degradation.  Strength and unloading degradation are deliberately
    // simple research controls until external parity tests are available.
    double lambda_strength{1e30};
    double lambda_unloading{1e30};
    double cyclic_exponent{1.0};
};

enum ASCE41Event : unsigned {
    ASCE41_EVENT_NONE          = 0u,
    ASCE41_EVENT_FAST_ELASTIC  = 1u << 0,
    ASCE41_EVENT_YIELD         = 1u << 1,
    ASCE41_EVENT_CAP           = 1u << 2,
    ASCE41_EVENT_STRENGTH_DROP = 1u << 3,
    ASCE41_EVENT_RESIDUAL      = 1u << 4,
    ASCE41_EVENT_REVERSAL      = 1u << 5,
    ASCE41_EVENT_DETERIORATION = 1u << 6,
    ASCE41_EVENT_FAILURE       = 1u << 7,
    ASCE41_EVENT_RELOAD        = 1u << 8,
    ASCE41_EVENT_LATERAL_LOSS  = 1u << 9
};

enum class ASCE41PerformanceLevel : int {
    Elastic = 0,
    InelasticBelowIO = 1,
    IO = 2,
    LS = 3,
    CP = 4,
    BeyondCP = 5,
    Failed = 6
};

struct ASCE41TrialResult {
    double force{};
    double tangent{};
    unsigned events{ASCE41_EVENT_NONE};
};

class ASCE41HingeMaterial {
public:
    // State is intentionally explicit so nearby Newton trial points cannot
    // infer different hysteretic branches from the same committed history.
    // q, f, cumulative/excursion energy, unloading K, last travel direction,
    // +/- envelope peaks, +/- strength scale, F-failure flag, persistent
    // E-lateral-loss flag, branch flag, zero-force intercept, reload target
    // q/f, and last converged material tangent.
    static constexpr int kStateSize = 19;

    explicit ASCE41HingeMaterial(ASCE41HingeParams params);

    int state_size() const { return kStateSize; }
    double initial_stiffness() const { return p_.Ke; }
    const ASCE41HingeParams& params() const { return p_; }

    void initialize_state(double* state) const;
    ASCE41TrialResult trial(double deformation, const double* committed,
                            double* trial_state) const;
    ASCE41PerformanceLevel performance_level(double deformation) const;

private:
    struct EnvelopePoint { double force{}, tangent{}; unsigned events{}; bool failed{}; };
    EnvelopePoint envelope(double deformation, double strength_scale) const;
    double degradation_beta(double excursion_energy, double cumulative_energy,
                            double lambda) const;

    ASCE41HingeParams p_;
};

} // namespace quake
