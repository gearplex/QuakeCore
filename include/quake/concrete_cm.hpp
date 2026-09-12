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
    CompressionUnloading = 3,
};

struct ConcreteCMState {
    double strain{};
    double stress{};
    double tangent{};
    double increment{};
    ConcreteCMRule rule{ConcreteCMRule::Initial};

    // First negative-to-positive reversal landmarks. These remain fixed while
    // traversing ConcreteCM rule 3 and are carried explicitly with the state.
    double unloading_strain{};
    double unloading_stress{};
    double zero_stress_strain{};
    double zero_stress_tangent{};
};

struct ConcreteCMTrial {
    ConcreteCMResponse response{};
    ConcreteCMState state{};
};

// Monotonic Chang-Mander envelope kernel used by the Gate 4 ConcreteCM
// implementation. Cyclic history rules are kept separate so their explicit
// committed/trial state can be validated independently.
class ConcreteCMEnvelope {
public:
    explicit ConcreteCMEnvelope(ConcreteCMParameters parameters);

    const ConcreteCMParameters& parameters() const noexcept { return parameters_; }
    ConcreteCMResponse compression(double strain) const;
    ConcreteCMResponse tension(double strain, double zero_strain = 0.0) const;

private:
    ConcreteCMParameters parameters_;
};

// Incremental Gate 4 ConcreteCM state machine. This milestone intentionally
// admits only monotonic compression and the first negative-to-positive
// unloading branch (ConcreteCM rule 3). Later cyclic rules are rejected until
// independently admitted against the frozen OpenSees 3.8.0 oracle.
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
