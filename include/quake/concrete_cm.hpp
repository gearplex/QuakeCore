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
};

struct ConcreteCMState {
    double strain{};
    double stress{};
    double tangent{};
    double increment{};
    ConcreteCMRule rule{ConcreteCMRule::Initial};

    // First negative-to-positive reversal history. These landmarks are frozen
    // at the reversal and carried through the first positive excursion.
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
    // rules 4, 10, 7 and the return to the negative envelope.
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

// Incremental Gate 4 ConcreteCM state machine. The admitted cyclic scope now
// covers the first complete compression -> tension -> compression excursion:
// rules 3/9/8/2 followed by 4/10/7/1. A subsequent reversal from the second
// compression extreme is rejected until nested-cycle rules are admitted.
class ConcreteCM {
public:
    explicit ConcreteCM(ConcreteCMParameters parameters);

    const ConcreteCMParameters& parameters() const noexcept { return envelope_.parameters(); }
    ConcreteCMState initial_state() const noexcept;
    ConcreteCMTrial trial(double strain, const ConcreteCMState& committed) const;

private:
    ConcreteCMEnvelope envelope_;
};

} // namespace quake
