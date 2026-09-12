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

// Monotonic Chang-Mander envelope kernel used by the Gate 4 ConcreteCM
// implementation. Cyclic history rules are intentionally kept out of this
// class so their commit/revert state can be validated independently.
class ConcreteCMEnvelope {
public:
    explicit ConcreteCMEnvelope(ConcreteCMParameters parameters);

    const ConcreteCMParameters& parameters() const noexcept { return parameters_; }
    ConcreteCMResponse compression(double strain) const;
    ConcreteCMResponse tension(double strain, double zero_strain = 0.0) const;

private:
    ConcreteCMParameters parameters_;
};

} // namespace quake
