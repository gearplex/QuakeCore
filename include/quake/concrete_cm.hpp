#pragma once

namespace quake {

struct ConcreteCMParameters {
    double fc{};
    double epsc{};
    double Ec{};
    double rc{};
    double xcrn{};
    double ft{};
    double et{};
    double rt{};
    double xcrp{};
    bool gap_close{false};
};

struct ConcreteCMResponse {
    double stress{};
    double tangent{};
};

enum class ConcreteCMRule : int {
    Initial = 0,
    CompressionEnvelope = 1,
    TensionEnvelope = 2,
    CompressionUnloading = 3,
    TensionUnloading = 4,
    TensionCutoff = 6,
    CompressionRejoining = 7,
    TensionRejoining = 8,
    CompressionToTension = 9,
    TensionToCompression = 10,
    NestedNegativeTarget = 11,
    NestedPositiveTarget = 12,
    CompressionReversalTransition = 77,
};

struct ConcreteCMState {
    double strain{};
    double stress{};
    double tangent{};
    double increment{};
    ConcreteCMRule rule{ConcreteCMRule::Initial};

    // Negative-to-positive reversal landmarks. These are refreshed when a
    // new compression-envelope extreme is admitted and then drive rules
    // 3, 9, 8, and the shifted positive envelope.
    double unloading_strain{};
    double unloading_stress{};
    double zero_stress_strain{};
    double zero_stress_tangent{};
    double tension_zero_strain{};
    double tension_peak_strain{};
    double tension_peak_stress{};
    double tension_new_stress{};
    double tension_new_tangent{};
    double tension_rejoin_strain{};
    double tension_rejoin_stress{};
    double tension_rejoin_tangent{};

    // First positive-to-negative reversal history. These landmarks drive
    // rules 4, 10, 7 and the return to the negative envelope. For OpenSees
    // rule 77, the same storage is reused for the reversal point and the
    // two transition endpoints while has_positive_to_negative_reversal
    // remains false, so it is not confused with the rule-4/10/7 path.
    bool has_positive_to_negative_reversal{false};
    double positive_reversal_strain{};
    double positive_reversal_stress{};
    double positive_zero_stress_strain{};
    double positive_zero_stress_tangent{};
    double compression_new_stress{};
    double compression_new_tangent{};
    double compression_rejoin_strain{};
    double compression_rejoin_stress{};
    double compression_rejoin_tangent{};

    // OpenSees rule-12 transition created when rule 77 reverses before the
    // prior compression unloading point is exceeded. The origin (Teb) and
    // positive target (Tea) remain fixed while the material advances.
    double nested_positive_origin_strain{};
    double nested_positive_origin_stress{};
    double nested_positive_target_strain{};

    // OpenSees rule-11 transition created when that rule-12 path reverses
    // negative. For the admitted rule77->12->11 provenance, rule 11 targets
    // the stored Teb point on the prior rule-77 path.
    double nested_negative_origin_strain{};
    double nested_negative_origin_stress{};

    // Second negative-to-positive reversal. Gate 4 currently admits only
    // the retained-prior-tension-history branch exercised by the frozen
    // YORi story-1 ConcreteCM protocol. Deeper nested reversals remain out
    // of scope and are rejected.
    bool has_second_negative_to_positive_reversal{false};
};

struct ConcreteCMTrial {
    ConcreteCMResponse response{};
    ConcreteCMState state{};
};

class ConcreteCMEnvelope {
public:
    explicit ConcreteCMEnvelope(ConcreteCMParameters parameters);

    const ConcreteCMParameters& parameters() const noexcept { return parameters_; }
    ConcreteCMResponse compression(double strain) const;
    ConcreteCMResponse tension(double strain, double zero_strain = 0.0) const;

private:
    ConcreteCMParameters parameters_;
};

// Incremental Gate 4 ConcreteCM state machine. The admitted cyclic scope
// covers the complete frozen YORi story-1 protocol plus the OpenSees rule-77
// reversal paths needed by the reconstructed 10001 system path. Deeper
// nested reversals remain rejected until independently admitted against an
// OpenSees 3.8.0 oracle.
class ConcreteCM {
public:
    explicit ConcreteCM(ConcreteCMParameters parameters);

    const ConcreteCMParameters& parameters() const noexcept { return envelope_.parameters(); }
    ConcreteCMState initial_state() const noexcept;
    ConcreteCMTrial trial(double strain, const ConcreteCMState& committed) const;
    ConcreteCMTrial trial_legacy(double strain, const ConcreteCMState& committed) const;

private:
    ConcreteCMEnvelope envelope_;
};

} // namespace quake

#ifdef QUAKE_CONCRETE_CM_RENAME_TRIAL
#define trial trial_legacy
#endif
