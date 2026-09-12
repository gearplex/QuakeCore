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
    TensionCutoff = 6,
    TensionRejoining = 8,
    CompressionToTension = 9,
};

struct ConcreteCMState {
    double strain{};
    double stress{};
    double tangent{};
    double increment{};
    ConcreteCMRule rule{ConcreteCMRule::Initial};

    // First negative-to-positive reversal history. These landmarks are frozen
    // at the reversal and carried explicitly through rules 3, 9, 8, and the
    // shifted positive envelope.
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

// Incremental Gate 4 ConcreteCM state machine. The admitted cyclic scope is
// the first compression-to-tension excursion: rules 3, 9, 8, then the shifted
// positive envelope (rules 2/6). Reversal back toward compression is rejected
// until the corresponding OpenSees rules are independently admitted.
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
