#include "quake/concrete_cm.hpp"

#include <cmath>
#include <stdexcept>

namespace quake {
namespace {

struct Shape {
    double y{};
    double z{};
};

struct CompressionUnloadingLandmarks {
    double zero_stress_strain{};
    double zero_stress_tangent{};
};

Shape chang_mander_shape(double x, double n, double r) {
    double D = 0.0;
    if (r != 1.0) {
        D = 1.0 + (n - r / (r - 1.0)) * x + std::pow(x, r) / (r - 1.0);
    } else {
        if (x <= 0.0) return {0.0, 1.0};
        D = 1.0 + (n - 1.0 + std::log10(x)) * x;
    }
    return {n * x / D, (1.0 - std::pow(x, r)) / (D * D)};
}

ConcreteCMResponse envelope(double strain,
                            double zero_strain,
                            double peak_strain,
                            double peak_stress,
                            double Ec,
                            double r,
                            double xcr) {
    const double scale = std::abs(peak_strain);
    const double x = std::abs((strain - zero_strain) / scale);
    const double n = std::abs(Ec * peak_strain / peak_stress);
    const auto critical = chang_mander_shape(xcr, n, r);
    const double cutoff = std::abs(xcr - critical.y / (n * critical.z));

    if (x > cutoff) return {0.0, 0.0};
    if (x < xcr) {
        const auto point = chang_mander_shape(x, n, r);
        return {peak_stress * point.y, Ec * point.z};
    }

    return {
        peak_stress * (critical.y + n * critical.z * (x - xcr)),
        Ec * critical.z,
    };
}

CompressionUnloadingLandmarks compression_unloading_landmarks(
    const ConcreteCMParameters& p, double eunn, double funn) {
    const double Esecn = p.Ec * ((std::abs(funn / (p.Ec * p.epsc)) + 0.57) /
                                 (std::abs(eunn / p.epsc) + 0.57));
    return {
        eunn - funn / Esecn,
        0.1 * p.Ec * std::exp(-2.0 * std::abs(eunn / p.epsc)),
    };
}

ConcreteCMResponse smooth_transition(double strain,
                                     double start_strain,
                                     double start_stress,
                                     double start_tangent,
                                     double end_strain,
                                     double end_stress,
                                     double end_tangent) {
    constexpr double huge = 1.797e308;
    const double Esec = (end_stress - start_stress) / (end_strain - start_strain);
    const double R = (end_tangent - Esec) / (Esec - start_tangent);
    const double span = std::abs(end_strain - start_strain);
    const double check = std::pow(span, R);

    double A = 0.0;
    if (check == 0.0 || check > huge || check < -huge || Esec == start_tangent) {
        A = 1.0e-300;
    } else {
        A = (Esec - start_tangent) / std::pow(span, R);
        if (A > huge || A < -huge) A = 1.0e300;
    }

    const double delta = std::abs(strain - start_strain);
    const double inverse_power = std::pow(delta, -R);
    const bool use_secant =
        A == 1.0e300 || A == 0.0 ||
        inverse_power == 0.0 || inverse_power > huge || inverse_power < -huge ||
        (start_tangent >= Esec && end_tangent >= Esec) ||
        (start_tangent <= Esec && end_tangent <= Esec);

    if (use_secant) {
        return {
            start_stress + Esec * (strain - start_strain),
            Esec,
        };
    }

    const double power = std::pow(delta, R);
    const double stress = start_stress +
        (strain - start_strain) * (start_tangent + A * power);
    const double tangent = start_tangent + A * (R + 1.0) * power;
    if (tangent > huge || tangent < -huge) {
        return {
            start_stress + Esec * (strain - start_strain),
            Esec,
        };
    }
    return {stress, tangent};
}

ConcreteCMResponse compression_rule3(const ConcreteCMParameters& p,
                                     const ConcreteCMState& state,
                                     double strain) {
    return smooth_transition(
        strain,
        state.unloading_strain,
        state.unloading_stress,
        p.Ec,
        state.zero_stress_strain,
        0.0,
        state.zero_stress_tangent);
}

ConcreteCMTrial make_trial(const ConcreteCMState& base,
                           double strain,
                           ConcreteCMResponse response,
                           double increment,
                           ConcreteCMRule rule) {
    ConcreteCMState next = base;
    next.strain = strain;
    next.stress = response.stress;
    next.tangent = response.tangent;
    next.increment = increment;
    next.rule = rule;
    return {response, next};
}

} // namespace

ConcreteCMEnvelope::ConcreteCMEnvelope(ConcreteCMParameters parameters)
    : parameters_(parameters) {
    if (!(parameters_.fc < 0.0) || !(parameters_.epsc < 0.0) ||
        !(parameters_.Ec > 0.0) || !(parameters_.rc > 0.0) ||
        !(parameters_.xcrn > 0.0) || !(parameters_.ft > 0.0) ||
        !(parameters_.et > 0.0) || !(parameters_.rt > 0.0) ||
        !(parameters_.xcrp > 0.0)) {
        throw std::invalid_argument("ConcreteCM parameters have invalid sign or zero value");
    }
}

ConcreteCMResponse ConcreteCMEnvelope::compression(double strain) const {
    if (strain > 0.0) throw std::invalid_argument("ConcreteCM compression envelope requires strain <= 0");
    if (strain == 0.0) return {0.0, parameters_.Ec};
    return envelope(strain, 0.0, parameters_.epsc, parameters_.fc,
                    parameters_.Ec, parameters_.rc, parameters_.xcrn);
}

ConcreteCMResponse ConcreteCMEnvelope::tension(double strain, double zero_strain) const {
    if (strain < zero_strain) {
        throw std::invalid_argument("ConcreteCM tension envelope requires strain >= zero strain");
    }
    if (strain == zero_strain) return {0.0, parameters_.Ec};
    return envelope(strain, zero_strain, parameters_.et, parameters_.ft,
                    parameters_.Ec, parameters_.rt, parameters_.xcrp);
}

ConcreteCM::ConcreteCM(ConcreteCMParameters parameters)
    : envelope_(parameters) {}

ConcreteCMState ConcreteCM::initial_state() const noexcept {
    ConcreteCMState state;
    state.tangent = parameters().Ec;
    return state;
}

ConcreteCMTrial ConcreteCM::trial(double strain, const ConcreteCMState& committed) const {
    if (strain > 0.0) {
        throw std::logic_error("ConcreteCM positive-strain cyclic rules are not yet admitted in Gate 4");
    }

    if (committed.rule == ConcreteCMRule::Initial) {
        if (strain == 0.0) {
            return make_trial(committed, 0.0, {0.0, parameters().Ec}, 0.0,
                              ConcreteCMRule::Initial);
        }
        const auto response = envelope_.compression(strain);
        return make_trial(committed, strain, response, -1.0,
                          ConcreteCMRule::CompressionEnvelope);
    }

    if (committed.rule == ConcreteCMRule::CompressionEnvelope) {
        if (strain <= committed.strain) {
            const auto response = envelope_.compression(strain);
            return make_trial(committed, strain, response, -1.0,
                              ConcreteCMRule::CompressionEnvelope);
        }

        if (!(committed.strain < 0.0) || !(committed.stress < 0.0)) {
            throw std::logic_error("ConcreteCM rule 3 requires a committed compression reversal point");
        }

        ConcreteCMState reversal = committed;
        reversal.unloading_strain = committed.strain;
        reversal.unloading_stress = committed.stress;
        const auto landmarks = compression_unloading_landmarks(
            parameters(), reversal.unloading_strain, reversal.unloading_stress);
        reversal.zero_stress_strain = landmarks.zero_stress_strain;
        reversal.zero_stress_tangent = landmarks.zero_stress_tangent;

        if (strain > reversal.zero_stress_strain) {
            throw std::logic_error("ConcreteCM rule 9 is not yet admitted in Gate 4");
        }
        const auto response = compression_rule3(parameters(), reversal, strain);
        return make_trial(reversal, strain, response, 1.0,
                          ConcreteCMRule::CompressionUnloading);
    }

    if (committed.rule == ConcreteCMRule::CompressionUnloading) {
        if (strain < committed.strain) {
            throw std::logic_error("ConcreteCM reversal from rule 3 is not yet admitted in Gate 4");
        }
        if (strain > committed.zero_stress_strain) {
            throw std::logic_error("ConcreteCM rule 9 is not yet admitted in Gate 4");
        }
        const auto response = compression_rule3(parameters(), committed, strain);
        return make_trial(committed, strain, response, 1.0,
                          ConcreteCMRule::CompressionUnloading);
    }

    throw std::logic_error("ConcreteCM committed rule is not implemented in Gate 4");
}

} // namespace quake
