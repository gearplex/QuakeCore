#include "quake/concrete_cm.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {

struct CompressionUnloadingLandmarks {
    double zero_stress_strain{};
    double zero_stress_tangent{};
};

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

ConcreteCMState post_extreme_reversal_state(
    const ConcreteCMEnvelope& envelope_kernel,
    const ConcreteCMState& committed) {
    const auto& p = envelope_kernel.parameters();
    ConcreteCMState reversal = committed;
    reversal.unloading_strain = committed.strain;
    reversal.unloading_stress = committed.stress;

    const auto compression = compression_unloading_landmarks(
        p, reversal.unloading_strain, reversal.unloading_stress);
    reversal.zero_stress_strain = compression.zero_stress_strain;
    reversal.zero_stress_tangent = compression.zero_stress_tangent;

    const double xun = std::abs(reversal.unloading_strain / p.epsc);
    double xup = std::abs(
        (committed.tension_peak_strain - committed.tension_zero_strain) / p.et);

    double e0ref = committed.tension_zero_strain;
    double eunpref = committed.tension_peak_strain;
    double funpref = committed.tension_peak_stress;
    if (xup < xun) {
        xup = xun;
        e0ref = 0.0;
        eunpref = xup * p.et;
        funpref = envelope_kernel.tension(eunpref, e0ref).stress;
    }

    const double Esecp = tension_secant(
        p, e0ref, eunpref, funpref, reversal.zero_stress_strain);
    const double dele0 = 2.0 * funpref /
                         (Esecp + reversal.zero_stress_tangent);
    reversal.tension_zero_strain =
        reversal.zero_stress_strain + dele0 - xup * p.et;
    reversal.tension_peak_strain =
        xup * p.et + reversal.tension_zero_strain;
    reversal.tension_peak_stress = envelope_kernel.tension(
        reversal.tension_peak_strain, reversal.tension_zero_strain).stress;
    populate_positive_rejoin_landmarks(envelope_kernel, reversal);
    return reversal;
}

ConcreteCMTrial make_trial(const ConcreteCMState& base,
                           double strain,
                           ConcreteCMResponse response,
                           ConcreteCMRule rule) {
    ConcreteCMState next = base;
    next.strain = strain;
    next.stress = response.stress;
    next.tangent = response.tangent;
    next.increment = 1.0;
    next.rule = rule;
    return {response, next};
}

ConcreteCMTrial positive_path_trial(const ConcreteCMEnvelope& envelope_kernel,
                                    const ConcreteCMState& base,
                                    double strain) {
    const auto& p = envelope_kernel.parameters();
    if (strain <= base.zero_stress_strain) {
        return make_trial(base, strain, compression_rule3(p, base, strain),
                          ConcreteCMRule::CompressionUnloading);
    }
    if (strain <= base.tension_peak_strain) {
        return make_trial(base, strain, rule9(base, strain),
                          ConcreteCMRule::CompressionToTension);
    }
    if (strain <= base.tension_rejoin_strain) {
        return make_trial(base, strain, rule8(base, strain),
                          ConcreteCMRule::TensionRejoining);
    }

    const auto response = envelope_kernel.tension(strain, base.tension_zero_strain);
    const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
        ? ConcreteCMRule::TensionCutoff
        : ConcreteCMRule::TensionEnvelope;
    return make_trial(base, strain, response, rule);
}

} // namespace

ConcreteCMTrial ConcreteCM::trial(double strain, const ConcreteCMState& committed) const {
    if (committed.rule == ConcreteCMRule::CompressionReversalTransition &&
        strain > committed.strain) {
        if (committed.strain < committed.unloading_strain) {
            const auto reversal = post_extreme_reversal_state(envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
        throw std::logic_error(
            "ConcreteCM rule77 positive reversal nested-targeting branch is not yet admitted in Gate 4");
    }
    return trial_legacy(strain, committed);
}

} // namespace quake
