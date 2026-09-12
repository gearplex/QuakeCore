#include "quake/concrete_cm.hpp"

#include <algorithm>
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

double tension_secant(const ConcreteCMParameters& p,
                      double e0,
                      double eunp,
                      double funp,
                      double espln) {
    double Esecp = p.Ec * ((std::abs(funp / (p.Ec * p.et)) + 0.67) /
                           (std::abs((eunp - e0) / p.et) + 0.67));
    const double lower_bound = std::abs(funp / std::abs(eunp - espln));
    if (Esecp < lower_bound) Esecp = lower_bound;
    return Esecp;
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
        return {start_stress + Esec * (strain - start_strain), Esec};
    }

    const double power = std::pow(delta, R);
    const double stress = start_stress +
        (strain - start_strain) * (start_tangent + A * power);
    const double tangent = start_tangent + A * (R + 1.0) * power;
    if (tangent > huge || tangent < -huge) {
        return {start_stress + Esec * (strain - start_strain), Esec};
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

ConcreteCMResponse rule9(const ConcreteCMState& state, double strain) {
    return smooth_transition(
        strain,
        state.zero_stress_strain,
        0.0,
        state.zero_stress_tangent,
        state.tension_peak_strain,
        state.tension_new_stress,
        state.tension_new_tangent);
}

ConcreteCMResponse rule8(const ConcreteCMState& state, double strain) {
    return smooth_transition(
        strain,
        state.tension_peak_strain,
        state.tension_new_stress,
        state.tension_new_tangent,
        state.tension_rejoin_strain,
        state.tension_rejoin_stress,
        state.tension_rejoin_tangent);
}

ConcreteCMResponse rule4(const ConcreteCMParameters& p,
                         const ConcreteCMState& state,
                         double strain) {
    return smooth_transition(
        strain,
        state.positive_reversal_strain,
        state.positive_reversal_stress,
        p.Ec,
        state.positive_zero_stress_strain,
        0.0,
        state.positive_zero_stress_tangent);
}

ConcreteCMResponse rule10(const ConcreteCMState& state, double strain) {
    return smooth_transition(
        strain,
        state.positive_zero_stress_strain,
        0.0,
        state.positive_zero_stress_tangent,
        state.unloading_strain,
        state.compression_new_stress,
        state.compression_new_tangent);
}

ConcreteCMResponse rule7(const ConcreteCMState& state, double strain) {
    return smooth_transition(
        strain,
        state.unloading_strain,
        state.compression_new_stress,
        state.compression_new_tangent,
        state.compression_rejoin_strain,
        state.compression_rejoin_stress,
        state.compression_rejoin_tangent);
}

void populate_positive_rejoin_landmarks(const ConcreteCMEnvelope& envelope_kernel,
                                        ConcreteCMState& reversal) {
    const auto& p = envelope_kernel.parameters();
    const double Esecp = tension_secant(
        p, reversal.tension_zero_strain, reversal.tension_peak_strain,
        reversal.tension_peak_stress, reversal.zero_stress_strain);
    const double esplp = reversal.tension_peak_strain -
                         reversal.tension_peak_stress / Esecp;
    const double delfp =
        reversal.tension_peak_strain >= reversal.tension_zero_strain + p.et / 2.0
            ? 0.15 * reversal.tension_peak_stress
            : 0.0;
    reversal.tension_new_stress = reversal.tension_peak_stress - delfp;
    reversal.tension_new_tangent =
        reversal.tension_peak_strain == esplp
            ? p.Ec
            : std::min(p.Ec, reversal.tension_new_stress /
                                 (reversal.tension_peak_strain - esplp));

    reversal.tension_rejoin_strain =
        reversal.tension_peak_strain +
        0.22 * std::abs(reversal.tension_peak_strain - reversal.tension_zero_strain);
    const auto rejoin = envelope_kernel.tension(
        reversal.tension_rejoin_strain, reversal.tension_zero_strain);
    reversal.tension_rejoin_stress = rejoin.stress;
    reversal.tension_rejoin_tangent = rejoin.tangent;
}

ConcreteCMState first_reversal_state(const ConcreteCMEnvelope& envelope_kernel,
                                     const ConcreteCMState& committed) {
    const auto& p = envelope_kernel.parameters();
    ConcreteCMState reversal = committed;
    reversal.unloading_strain = committed.strain;
    reversal.unloading_stress = committed.stress;

    const auto compression = compression_unloading_landmarks(
        p, reversal.unloading_strain, reversal.unloading_stress);
    reversal.zero_stress_strain = compression.zero_stress_strain;
    reversal.zero_stress_tangent = compression.zero_stress_tangent;

    const double xup = std::abs(reversal.unloading_strain / p.epsc);
    const double reference_peak_strain = xup * p.et;
    const double reference_peak_stress =
        envelope_kernel.tension(reference_peak_strain, 0.0).stress;
    const double reference_secant = tension_secant(
        p, 0.0, reference_peak_strain, reference_peak_stress,
        reversal.zero_stress_strain);
    const double dele0 = 2.0 * reference_peak_stress /
                         (reference_secant + reversal.zero_stress_tangent);

    reversal.tension_zero_strain =
        reversal.zero_stress_strain + dele0 - xup * p.et;
    reversal.tension_peak_strain =
        xup * p.et + reversal.tension_zero_strain;
    reversal.tension_peak_stress = envelope_kernel.tension(
        reversal.tension_peak_strain, reversal.tension_zero_strain).stress;
    populate_positive_rejoin_landmarks(envelope_kernel, reversal);
    return reversal;
}

ConcreteCMState second_reversal_state(const ConcreteCMEnvelope& envelope_kernel,
                                      const ConcreteCMState& committed) {
    const auto& p = envelope_kernel.parameters();
    ConcreteCMState reversal = committed;
    reversal.has_positive_to_negative_reversal = true;
    reversal.positive_reversal_strain = committed.strain;
    reversal.positive_reversal_stress = envelope_kernel.tension(
        committed.strain, committed.tension_zero_strain).stress;

    const double Esecp = tension_secant(
        p,
        committed.tension_zero_strain,
        reversal.positive_reversal_strain,
        reversal.positive_reversal_stress,
        committed.zero_stress_strain);
    reversal.positive_zero_stress_strain =
        reversal.positive_reversal_strain - reversal.positive_reversal_stress / Esecp;
    reversal.positive_zero_stress_tangent = p.gap_close
        ? p.Ec / (std::pow(std::abs(
              (reversal.positive_reversal_strain - committed.tension_zero_strain) / p.et), 1.1) + 1.0)
        : 0.0;

    const double delfn = committed.unloading_strain <= p.epsc / 10.0
        ? 0.09 * committed.unloading_stress *
              std::pow(std::abs(committed.unloading_strain / p.epsc), 0.5)
        : 0.0;
    reversal.compression_new_stress = committed.unloading_stress - delfn;
    reversal.compression_new_tangent =
        committed.unloading_strain == committed.zero_stress_strain
            ? p.Ec
            : std::min(p.Ec,
                       reversal.compression_new_stress /
                           (committed.unloading_strain - committed.zero_stress_strain));

    const double delen = committed.unloading_strain /
        (1.15 + 2.75 * std::abs(committed.unloading_strain / p.epsc));
    reversal.compression_rejoin_strain = committed.unloading_strain + delen;
    const auto rejoin = envelope_kernel.compression(reversal.compression_rejoin_strain);
    reversal.compression_rejoin_stress = rejoin.stress;
    reversal.compression_rejoin_tangent = rejoin.tangent;
    return reversal;
}

ConcreteCMState second_negative_to_positive_reversal_state(
    const ConcreteCMEnvelope& envelope_kernel,
    const ConcreteCMState& committed) {
    const auto& p = envelope_kernel.parameters();
    ConcreteCMState reversal = committed;
    reversal.has_second_negative_to_positive_reversal = true;
    reversal.unloading_strain = committed.strain;
    reversal.unloading_stress = committed.stress;

    const auto compression = compression_unloading_landmarks(
        p, reversal.unloading_strain, reversal.unloading_stress);
    reversal.zero_stress_strain = compression.zero_stress_strain;
    reversal.zero_stress_tangent = compression.zero_stress_tangent;

    // OpenSees e0eunpfunpf: the frozen YORi protocol reaches a much larger
    // prior positive excursion than the normalized second compression
    // extreme, so xup >= xun and the prior tension history is retained.
    const double xun = std::abs(reversal.unloading_strain / p.epsc);
    const double xup = std::abs(
        (committed.positive_reversal_strain - committed.tension_zero_strain) / p.et);
    if (xup < xun) {
        throw std::logic_error(
            "ConcreteCM second rebound compression-dominant branch is not yet admitted in Gate 4");
    }

    const double reference_peak_strain = committed.positive_reversal_strain;
    const double reference_peak_stress = committed.positive_reversal_stress;
    const double reference_secant = tension_secant(
        p, committed.tension_zero_strain, reference_peak_strain,
        reference_peak_stress, reversal.zero_stress_strain);
    const double dele0 = 2.0 * reference_peak_stress /
                         (reference_secant + reversal.zero_stress_tangent);

    reversal.tension_zero_strain =
        reversal.zero_stress_strain + dele0 - xup * p.et;
    reversal.tension_peak_strain = xup * p.et + reversal.tension_zero_strain;
    reversal.tension_peak_stress = envelope_kernel.tension(
        reversal.tension_peak_strain, reversal.tension_zero_strain).stress;
    populate_positive_rejoin_landmarks(envelope_kernel, reversal);
    return reversal;
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

ConcreteCMTrial positive_path_trial(const ConcreteCMEnvelope& envelope_kernel,
                                    const ConcreteCMState& base,
                                    double strain) {
    const auto& p = envelope_kernel.parameters();
    if (strain <= base.zero_stress_strain) {
        return make_trial(base, strain, compression_rule3(p, base, strain), 1.0,
                          ConcreteCMRule::CompressionUnloading);
    }
    if (strain <= base.tension_peak_strain) {
        return make_trial(base, strain, rule9(base, strain), 1.0,
                          ConcreteCMRule::CompressionToTension);
    }
    if (strain <= base.tension_rejoin_strain) {
        return make_trial(base, strain, rule8(base, strain), 1.0,
                          ConcreteCMRule::TensionRejoining);
    }

    const auto response = envelope_kernel.tension(strain, base.tension_zero_strain);
    const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
        ? ConcreteCMRule::TensionCutoff
        : ConcreteCMRule::TensionEnvelope;
    return make_trial(base, strain, response, 1.0, rule);
}

ConcreteCMTrial first_negative_return_trial(const ConcreteCMEnvelope& envelope_kernel,
                                            const ConcreteCMState& base,
                                            double strain) {
    const auto& p = envelope_kernel.parameters();
    if (strain >= base.positive_zero_stress_strain) {
        return make_trial(base, strain, rule4(p, base, strain), -1.0,
                          ConcreteCMRule::TensionUnloading);
    }
    if (strain >= base.unloading_strain) {
        return make_trial(base, strain, rule10(base, strain), -1.0,
                          ConcreteCMRule::TensionToCompression);
    }
    if (strain >= base.compression_rejoin_strain) {
        return make_trial(base, strain, rule7(base, strain), -1.0,
                          ConcreteCMRule::CompressionRejoining);
    }
    return make_trial(base, strain, envelope_kernel.compression(strain), -1.0,
                      ConcreteCMRule::CompressionEnvelope);
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
    if (strain > 0.0) {
        throw std::invalid_argument("ConcreteCM compression envelope requires strain <= 0");
    }
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
    if (committed.rule == ConcreteCMRule::Initial) {
        if (strain == 0.0) {
            return make_trial(committed, 0.0, {0.0, parameters().Ec}, 0.0,
                              ConcreteCMRule::Initial);
        }
        if (strain < 0.0) {
            return make_trial(committed, strain, envelope_.compression(strain), -1.0,
                              ConcreteCMRule::CompressionEnvelope);
        }
        ConcreteCMState positive = committed;
        positive.tension_zero_strain = 0.0;
        const auto response = envelope_.tension(strain, 0.0);
        return make_trial(positive, strain, response, 1.0,
                          ConcreteCMRule::TensionEnvelope);
    }

    if (committed.rule == ConcreteCMRule::CompressionEnvelope) {
        if (strain <= committed.strain) {
            return make_trial(committed, strain, envelope_.compression(strain), -1.0,
                              ConcreteCMRule::CompressionEnvelope);
        }
        if (!(committed.strain < 0.0) || !(committed.stress < 0.0)) {
            throw std::logic_error(
                "ConcreteCM negative-to-positive reversal requires a committed compression point");
        }
        if (committed.has_second_negative_to_positive_reversal) {
            throw std::logic_error(
                "ConcreteCM deeper nested reversal is not yet admitted in Gate 4");
        }
        if (committed.has_positive_to_negative_reversal) {
            const auto reversal = second_negative_to_positive_reversal_state(
                envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
        const auto reversal = first_reversal_state(envelope_, committed);
        return positive_path_trial(envelope_, reversal, strain);
    }

    const bool first_rebound = committed.unloading_strain < 0.0;
    if (first_rebound && !committed.has_positive_to_negative_reversal &&
        !committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::CompressionUnloading ||
         committed.rule == ConcreteCMRule::CompressionToTension ||
         committed.rule == ConcreteCMRule::TensionRejoining ||
         committed.rule == ConcreteCMRule::TensionEnvelope ||
         committed.rule == ConcreteCMRule::TensionCutoff)) {
        if (strain < committed.strain) {
            if (committed.rule != ConcreteCMRule::TensionEnvelope &&
                committed.rule != ConcreteCMRule::TensionRejoining) {
                throw std::logic_error(
                    "ConcreteCM reversal from this positive-going rule is not yet admitted in Gate 4");
            }
            const auto reversal = second_reversal_state(envelope_, committed);
            return first_negative_return_trial(envelope_, reversal, strain);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    if (committed.has_positive_to_negative_reversal &&
        !committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::TensionUnloading ||
         committed.rule == ConcreteCMRule::TensionToCompression ||
         committed.rule == ConcreteCMRule::CompressionRejoining)) {
        if (strain > committed.strain) {
            throw std::logic_error(
                "ConcreteCM reversal from the first negative return is not yet admitted in Gate 4");
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }

    if (committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::CompressionUnloading ||
         committed.rule == ConcreteCMRule::CompressionToTension ||
         committed.rule == ConcreteCMRule::TensionRejoining ||
         committed.rule == ConcreteCMRule::TensionEnvelope ||
         committed.rule == ConcreteCMRule::TensionCutoff)) {
        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM reversal from the second rebound is not yet admitted in Gate 4");
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::TensionEnvelope ||
        committed.rule == ConcreteCMRule::TensionCutoff) {
        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM virgin tension reversal rules are not yet admitted in Gate 4");
        }
        const auto response = envelope_.tension(strain, committed.tension_zero_strain);
        const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
            ? ConcreteCMRule::TensionCutoff
            : ConcreteCMRule::TensionEnvelope;
        return make_trial(committed, strain, response, 1.0, rule);
    }

    throw std::logic_error("ConcreteCM committed rule is not implemented in Gate 4");
}

} // namespace quake
