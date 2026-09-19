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
    auto response = smooth_transition(
        strain,
        state.positive_zero_stress_strain,
        0.0,
        state.positive_zero_stress_tangent,
        state.unloading_strain,
        state.compression_new_stress,
        state.compression_new_tangent);

    // OpenSees 3.8.0 ConcreteCM "Fix 2": when fcEturf collapses onto
    // the endpoint secant for rule 10, replace that degenerate branch by
    // either the zero-stress gap or the negative-side Enewn line through
    // espln. This is an explicit post-processing step in ConcreteCM.cpp.
    const double endpoint_secant =
        state.compression_new_stress /
        (state.unloading_strain - state.positive_zero_stress_strain);
    if (response.tangent == endpoint_secant) {
        if (strain >= state.zero_stress_strain) {
            return {0.0, 0.0};
        }
        return {
            state.compression_new_tangent * (strain - state.zero_stress_strain),
            state.compression_new_tangent,
        };
    }
    return response;
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

    // OpenSees 3.8.0 Crule 2/8 negative reversal promotes Cstrain/Cstress
    // to the primary positive history point Teunp/Tfunp.  The earlier Gate-4
    // implementation only saved this point in positive_reversal_* and left
    // tension_peak_* stale, which later corrupts ea1112f on Crule 10 -> 12.
    reversal.tension_peak_strain = committed.strain;
    reversal.tension_peak_stress = envelope_kernel.tension(
        committed.strain, committed.tension_zero_strain).stress;
    reversal.positive_reversal_strain = reversal.tension_peak_strain;
    reversal.positive_reversal_stress = reversal.tension_peak_stress;

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

    // OpenSees reconstructs fnewp/Enewp/esrep/frep from the newly committed
    // Teunp/Tfunp at the start of every subsequent cyclic trial.
    populate_positive_rejoin_landmarks(envelope_kernel, reversal);

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

    // OpenSees e0eunpfunpf is cycle-count agnostic. It compares the new
    // compression extreme xun with the current primary Teunp/Tfunp history,
    // not with an auxiliary earlier reversal point. If compression dominates,
    // it resets the tensile reference to the virgin envelope at xup=xun.
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

    const double reference_secant = tension_secant(
        p, e0ref, eunpref, funpref, reversal.zero_stress_strain);
    const double dele0 = 2.0 * funpref /
                         (reference_secant + reversal.zero_stress_tangent);

    reversal.tension_zero_strain =
        reversal.zero_stress_strain + dele0 - xup * p.et;
    reversal.tension_peak_strain = xup * p.et + reversal.tension_zero_strain;
    reversal.tension_peak_stress = envelope_kernel.tension(
        reversal.tension_peak_strain, reversal.tension_zero_strain).stress;
    populate_positive_rejoin_landmarks(envelope_kernel, reversal);
    return reversal;
}

ConcreteCMState rule77_reversal_state(const ConcreteCMEnvelope& envelope_kernel,
                                      const ConcreteCMState& committed) {
    const auto& p = envelope_kernel.parameters();
    ConcreteCMState reversal = committed;

    const double eunn = committed.unloading_strain;
    const double funn = committed.unloading_stress;
    const double er0n = committed.strain;
    const double fr0n = committed.stress;
    const double espln = committed.zero_stress_strain;

    if (!(eunn < er0n) || !(eunn < espln)) {
        throw std::logic_error("ConcreteCM rule77 landmarks are not ordered");
    }

    const double delfn = eunn <= p.epsc / 10.0
        ? 0.09 * funn * std::pow(std::abs(eunn / p.epsc), 0.5)
        : 0.0;
    const double fnewstn =
        funn - delfn * ((eunn - er0n) / (eunn - espln));
    const double Enewstn = (fnewstn - fr0n) / (eunn - er0n);
    const double delen = eunn /
        (1.15 + 2.75 * std::abs(eunn / p.epsc));
    const double esrestn =
        eunn + delen * (eunn - er0n) / (eunn - espln);
    const auto rest = envelope_kernel.compression(esrestn);

    reversal.positive_reversal_strain = er0n;
    reversal.positive_reversal_stress = fr0n;
    reversal.compression_new_stress = fnewstn;
    reversal.compression_new_tangent = Enewstn;
    reversal.compression_rejoin_strain = esrestn;
    reversal.compression_rejoin_stress = rest.stress;
    reversal.compression_rejoin_tangent = rest.tangent;
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

struct PositiveReturnLandmarks {
    double shortened_peak_stress{};
    double shortened_peak_tangent{};
    double shortened_rejoin_strain{};
    double shortened_rejoin_stress{};
    double shortened_rejoin_tangent{};
};

PositiveReturnLandmarks positive_return_landmarks(
    const ConcreteCMEnvelope& envelope_kernel,
    const ConcreteCMState& state) {
    const auto& p = envelope_kernel.parameters();
    const double Esecp = tension_secant(
        p, state.tension_zero_strain, state.tension_peak_strain,
        state.tension_peak_stress, state.zero_stress_strain);
    const double esplp = state.tension_peak_strain - state.tension_peak_stress / Esecp;
    const double delfp =
        state.tension_peak_strain >= state.tension_zero_strain + p.et / 2.0
            ? 0.15 * state.tension_peak_stress
            : 0.0;
    const double ratio =
        (state.tension_peak_strain - state.positive_return_reversal_strain) /
        (state.tension_peak_strain - esplp);
    const double fnewstp = state.tension_peak_stress - delfp * ratio;
    const double Enewstp =
        (fnewstp - state.positive_return_reversal_stress) /
        (state.tension_peak_strain - state.positive_return_reversal_strain);
    const double delep = 0.22 * std::abs(
        state.tension_peak_strain - state.tension_zero_strain);
    const double esrestp = state.tension_peak_strain + delep * ratio;
    const auto rest = envelope_kernel.tension(esrestp, state.tension_zero_strain);
    return {fnewstp, Enewstp, esrestp, rest.stress, rest.tangent};
}

ConcreteCMTrial rule88_trial(const ConcreteCMEnvelope& envelope_kernel,
                             const ConcreteCMState& base,
                             double strain) {
    const auto& p = envelope_kernel.parameters();
    const auto lm = positive_return_landmarks(envelope_kernel, base);
    if (strain <= base.tension_peak_strain) {
        return make_trial(
            base, strain,
            smooth_transition(
                strain,
                base.positive_return_reversal_strain,
                base.positive_return_reversal_stress,
                p.Ec,
                base.tension_peak_strain,
                lm.shortened_peak_stress,
                lm.shortened_peak_tangent),
            1.0, ConcreteCMRule::TensionReversalTransition);
    }
    if (strain < lm.shortened_rejoin_strain) {
        return make_trial(
            base, strain,
            smooth_transition(
                strain,
                base.tension_peak_strain,
                lm.shortened_peak_stress,
                lm.shortened_peak_tangent,
                lm.shortened_rejoin_strain,
                lm.shortened_rejoin_stress,
                lm.shortened_rejoin_tangent),
            1.0, ConcreteCMRule::TensionReversalTransition);
    }
    const auto response = envelope_kernel.tension(strain, base.tension_zero_strain);
    const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
        ? ConcreteCMRule::TensionCutoff
        : ConcreteCMRule::TensionEnvelope;
    return make_trial(base, strain, response, 1.0, rule);
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

ConcreteCMTrial rule77_trial(const ConcreteCMEnvelope& envelope_kernel,
                             const ConcreteCMState& base,
                             double strain) {
    const auto& p = envelope_kernel.parameters();
    if (strain >= base.unloading_strain) {
        const auto response = smooth_transition(
            strain,
            base.positive_reversal_strain,
            base.positive_reversal_stress,
            p.Ec,
            base.unloading_strain,
            base.compression_new_stress,
            base.compression_new_tangent);
        return make_trial(base, strain, response, -1.0,
                          ConcreteCMRule::CompressionReversalTransition);
    }
    if (strain > base.compression_rejoin_strain) {
        const auto response = smooth_transition(
            strain,
            base.unloading_strain,
            base.compression_new_stress,
            base.compression_new_tangent,
            base.compression_rejoin_strain,
            base.compression_rejoin_stress,
            base.compression_rejoin_tangent);
        return make_trial(base, strain, response, -1.0,
                          ConcreteCMRule::CompressionReversalTransition);
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
        if (committed.has_positive_to_negative_reversal ||
            committed.has_second_negative_to_positive_reversal) {
            const auto reversal = second_negative_to_positive_reversal_state(
                envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
        const auto reversal = first_reversal_state(envelope_, committed);
        return positive_path_trial(envelope_, reversal, strain);
    }

    if (committed.rule == ConcreteCMRule::CompressionReversalTransition) {
        if (strain > committed.strain) {
            if (committed.strain < committed.unloading_strain) {
                const auto& p = parameters();
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
                    funpref = envelope_.tension(eunpref, e0ref).stress;
                }

                const double Esecp = tension_secant(
                    p, e0ref, eunpref, funpref, reversal.zero_stress_strain);
                const double dele0 = 2.0 * funpref /
                                     (Esecp + reversal.zero_stress_tangent);
                reversal.tension_zero_strain =
                    reversal.zero_stress_strain + dele0 - xup * p.et;
                reversal.tension_peak_strain =
                    xup * p.et + reversal.tension_zero_strain;
                reversal.tension_peak_stress = envelope_.tension(
                    reversal.tension_peak_strain, reversal.tension_zero_strain).stress;
                populate_positive_rejoin_landmarks(envelope_, reversal);
                return positive_path_trial(envelope_, reversal, strain);
            }

            // OpenSees 3.8.0 Crule=77, Cstrain>=Teunn. Tea is the original
            // rule-77 reversal strain (Ter0n), while Teb is this reversal
            // point. Rule 12 transitions from Ter/Tfr to Tea with Ec as the
            // initial tangent; Tea and Teb remain fixed across nested cycles.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain;
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_positive_target_strain = committed.positive_reversal_strain; // Tea=Ter0n
            reversal.nested_negative_target_strain = committed.strain; // Teb (rule77)
            reversal.nested_negative_target_uses_rule77 = true;
            const auto target = positive_path_trial(
                envelope_, committed, reversal.nested_positive_target_strain);
            if (strain <= reversal.nested_positive_target_strain) {
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    reversal.nested_positive_target_strain,
                    target.response.stress,
                    target.response.tangent);
                return make_trial(reversal, strain, response, 1.0,
                                  ConcreteCMRule::NestedPositiveTarget);
            }
            return positive_path_trial(envelope_, reversal, strain);
        }
        return rule77_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::NestedPositiveTarget) {
        const auto positive_target = positive_path_trial(
            envelope_, committed, committed.nested_positive_target_strain);
        if (strain < committed.strain) {
            // OpenSees 3.8.0 Crule=12 reversal for the rule77->12 provenance:
            // Tea == Ter0n and Teb is the saved rule-77 reversal point. The
            // negative reversal creates rule 11 from current Ter/Tfr to Teb.
            ConcreteCMState reversal = committed;
            reversal.nested_negative_origin_strain = committed.strain;
            reversal.nested_negative_origin_stress = committed.stress;
            const double target_strain = committed.nested_negative_target_strain; // Teb
            const bool from_rule77 =
                committed.nested_negative_target_uses_rule77;
            const auto target = from_rule77
                ? rule77_trial(envelope_, committed, target_strain)
                : first_negative_return_trial(envelope_, committed, target_strain);
            if (strain >= target_strain) {
                const auto normalized_target = smooth_transition(
                    target_strain,
                    reversal.nested_negative_origin_strain,
                    reversal.nested_negative_origin_stress,
                    parameters().Ec,
                    target_strain,
                    target.response.stress,
                    target.response.tangent);
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_negative_origin_strain,
                    reversal.nested_negative_origin_stress,
                    parameters().Ec,
                    target_strain,
                    normalized_target.stress,
                    normalized_target.tangent);
                return make_trial(reversal, strain, response, -1.0,
                                  ConcreteCMRule::NestedNegativeTarget);
            }
            return from_rule77
                ? rule77_trial(envelope_, reversal, strain)
                : first_negative_return_trial(envelope_, reversal, strain);
        }
        if (strain <= committed.nested_positive_target_strain) {
            const auto response = smooth_transition(
                strain,
                committed.nested_positive_origin_strain,
                committed.nested_positive_origin_stress,
                parameters().Ec,
                committed.nested_positive_target_strain,
                positive_target.response.stress,
                positive_target.response.tangent);
            return make_trial(committed, strain, response, 1.0,
                              ConcreteCMRule::NestedPositiveTarget);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::NestedNegativeTarget) {
        const double target_strain = committed.nested_negative_target_strain; // Teb
        const bool from_rule77 =
            committed.nested_negative_target_uses_rule77;
        const auto negative_target = from_rule77
            ? rule77_trial(envelope_, committed, target_strain)
            : first_negative_return_trial(envelope_, committed, target_strain);
        if (strain > committed.strain) {
            // OpenSees 3.8.0 Crule=11 positive reversal for this provenance:
            // keep Tea/Teb fixed, move Ter/Tfr to the current point, and
            // create rule 12 targeting Tea on the established positive path.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain;
            reversal.nested_positive_origin_stress = committed.stress;
            const auto positive_target = positive_path_trial(
                envelope_, committed, committed.nested_positive_target_strain);
            if (strain <= committed.nested_positive_target_strain) {
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    committed.nested_positive_target_strain,
                    positive_target.response.stress,
                    positive_target.response.tangent);
                return make_trial(reversal, strain, response, 1.0,
                                  ConcreteCMRule::NestedPositiveTarget);
            }
            return positive_path_trial(envelope_, reversal, strain);
        }
        if (strain >= target_strain) {
            const auto response = smooth_transition(
                strain,
                committed.nested_negative_origin_strain,
                committed.nested_negative_origin_stress,
                parameters().Ec,
                target_strain,
                negative_target.response.stress,
                negative_target.response.tangent);
            return make_trial(committed, strain, response, -1.0,
                              ConcreteCMRule::NestedNegativeTarget);
        }
        return from_rule77
            ? rule77_trial(envelope_, committed, strain)
            : first_negative_return_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::CompressionToTension) {
        if (strain < committed.strain) {
            // OpenSees 3.8.0 Crule=9 negative reversal. This branch is driven
            // by Crule, not by QuakeCore's cycle-history flags. Reconstruct
            // the derived esplp/Eplp and negative-return landmarks from the
            // primary Teunp/Tfunp and Teunn/Tfunn history exactly as OpenSees
            // does at the start of each cyclic trial, then apply eb1112f.
            const auto& p = parameters();
            ConcreteCMState negative_base = committed;

            const double Esecp = tension_secant(
                p,
                committed.tension_zero_strain,
                committed.tension_peak_strain,
                committed.tension_peak_stress,
                committed.zero_stress_strain);
            const double esplp =
                committed.tension_peak_strain - committed.tension_peak_stress / Esecp;
            const double Eplp = p.gap_close
                ? p.Ec / (std::pow(std::abs(
                      (committed.tension_peak_strain - committed.tension_zero_strain) /
                      p.et), 1.1) + 1.0)
                : 0.0;
            negative_base.positive_reversal_strain = committed.tension_peak_strain;
            negative_base.positive_reversal_stress = committed.tension_peak_stress;
            negative_base.positive_zero_stress_strain = esplp;
            negative_base.positive_zero_stress_tangent = Eplp;

            const double delfn = committed.unloading_strain <= p.epsc / 10.0
                ? 0.09 * committed.unloading_stress *
                      std::pow(std::abs(committed.unloading_strain / p.epsc), 0.5)
                : 0.0;
            negative_base.compression_new_stress = committed.unloading_stress - delfn;
            negative_base.compression_new_tangent =
                committed.unloading_strain == committed.zero_stress_strain
                    ? p.Ec
                    : std::min(
                          p.Ec,
                          negative_base.compression_new_stress /
                              (committed.unloading_strain - committed.zero_stress_strain));
            const double delen = committed.unloading_strain /
                (1.15 + 2.75 * std::abs(committed.unloading_strain / p.epsc));
            negative_base.compression_rejoin_strain = committed.unloading_strain + delen;
            const auto compression_rejoin = envelope_.compression(
                negative_base.compression_rejoin_strain);
            negative_base.compression_rejoin_stress = compression_rejoin.stress;
            negative_base.compression_rejoin_tangent = compression_rejoin.tangent;

            ConcreteCMState reversal = negative_base;
            reversal.nested_negative_origin_strain = committed.strain; // Ter
            reversal.nested_negative_origin_stress = committed.stress; // Tfr
            reversal.nested_positive_target_strain = committed.strain; // Tea
            const double denom =
                committed.tension_peak_strain - committed.zero_stress_strain;
            reversal.nested_negative_target_strain =
                committed.unloading_strain -
                ((committed.strain - committed.zero_stress_strain) / denom) *
                    (committed.unloading_strain - esplp); // Teb = eb1112f(...)
            reversal.nested_negative_target_uses_rule77 = false;

            const double target_strain = reversal.nested_negative_target_strain;
            const auto target = first_negative_return_trial(
                envelope_, negative_base, target_strain);
            if (strain >= target_strain) {
                // OpenSees r11f normalizes the endpoint at Teb before the
                // caller evaluates the actual trial strain.
                const auto normalized_target = smooth_transition(
                    target_strain,
                    reversal.nested_negative_origin_strain,
                    reversal.nested_negative_origin_stress,
                    p.Ec,
                    target_strain,
                    target.response.stress,
                    target.response.tangent);
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_negative_origin_strain,
                    reversal.nested_negative_origin_stress,
                    p.Ec,
                    target_strain,
                    normalized_target.stress,
                    normalized_target.tangent);
                return make_trial(reversal, strain, response, -1.0,
                                  ConcreteCMRule::NestedNegativeTarget);
            }
            return first_negative_return_trial(envelope_, reversal, strain);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::TensionReversalTransition) {
        if (strain < committed.strain) {
            // OpenSees 3.8.0 Crule=88 negative reversal.
            // If the committed point is still at/below Teunp, preserve the
            // rule-4 reversal point as Teb=Ter0p and create rule 11 from the
            // current Ter/Tfr back to that point.  Ter0p was created while
            // rule 4 was active, so this demanded provenance has Teb>=esplp.
            if (committed.strain <= committed.tension_peak_strain) {
                ConcreteCMState reversal = committed;
                reversal.nested_negative_origin_strain = committed.strain; // Ter
                reversal.nested_negative_origin_stress = committed.stress; // Tfr
                reversal.nested_positive_target_strain = committed.strain; // Tea=Ter
                reversal.nested_negative_target_strain =
                    committed.positive_return_reversal_strain; // Teb=Ter0p
                reversal.nested_negative_target_uses_rule77 = false;

                const double target_strain = reversal.nested_negative_target_strain;
                const auto target = first_negative_return_trial(
                    envelope_, committed, target_strain); // rule 4 at Teb
                if (strain >= target_strain) {
                    // OpenSees r11f evaluates the raw target once at Teb,
                    // promotes that evaluated stress/tangent to ff/Ef, then
                    // the caller evaluates the actual trial point.
                    const auto normalized_target = smooth_transition(
                        target_strain,
                        reversal.nested_negative_origin_strain,
                        reversal.nested_negative_origin_stress,
                        parameters().Ec,
                        target_strain,
                        target.response.stress,
                        target.response.tangent);
                    const auto response = smooth_transition(
                        strain,
                        reversal.nested_negative_origin_strain,
                        reversal.nested_negative_origin_stress,
                        parameters().Ec,
                        target_strain,
                        normalized_target.stress,
                        normalized_target.tangent);
                    return make_trial(reversal, strain, response, -1.0,
                                      ConcreteCMRule::NestedNegativeTarget);
                }
                return first_negative_return_trial(envelope_, reversal, strain);
            }

            // OpenSees Crule=88 with Cstrain>Teunp promotes the current point
            // to the new Teunp/Tfunp, rebuilds esplp/Eplp, and then follows
            // rules 4,10,7,1/5.  second_reversal_state is the QuakeCore
            // representation of that primary-history refresh.
            const auto reversal = second_reversal_state(envelope_, committed);
            return first_negative_return_trial(envelope_, reversal, strain);
        }
        return rule88_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::CompressionRejoining) {
        if (strain > committed.strain) {
            const auto reversal = second_negative_to_positive_reversal_state(
                envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::TensionToCompression) {
        // OpenSees 3.8.0 dispatches Crule=10 by rule identity, independent of
        // how the material arrived there. Reconstruct the derived cyclic
        // landmarks from the primary Teunp/Tfunp and Teunn/Tfunn history
        // before either continuing negative or reversing positive.
        const auto& p = parameters();
        ConcreteCMState refreshed = committed;

        const double Esecp = tension_secant(
            p,
            committed.tension_zero_strain,
            committed.tension_peak_strain,
            committed.tension_peak_stress,
            committed.zero_stress_strain);
        const double esplp =
            committed.tension_peak_strain - committed.tension_peak_stress / Esecp;
        const double Eplp = p.gap_close
            ? p.Ec / (std::pow(std::abs(
                  (committed.tension_peak_strain - committed.tension_zero_strain) /
                  p.et), 1.1) + 1.0)
            : 0.0;
        refreshed.positive_zero_stress_strain = esplp;
        refreshed.positive_zero_stress_tangent = Eplp;
        refreshed.positive_reversal_strain = committed.tension_peak_strain;
        refreshed.positive_reversal_stress = committed.tension_peak_stress;

        const double delfn = committed.unloading_strain <= p.epsc / 10.0
            ? 0.09 * committed.unloading_stress *
                  std::pow(std::abs(committed.unloading_strain / p.epsc), 0.5)
            : 0.0;
        refreshed.compression_new_stress = committed.unloading_stress - delfn;
        refreshed.compression_new_tangent =
            committed.unloading_strain == committed.zero_stress_strain
                ? p.Ec
                : std::min(
                      p.Ec,
                      refreshed.compression_new_stress /
                          (committed.unloading_strain - committed.zero_stress_strain));
        const double delen = committed.unloading_strain /
            (1.15 + 2.75 * std::abs(committed.unloading_strain / p.epsc));
        refreshed.compression_rejoin_strain = committed.unloading_strain + delen;
        const auto compression_rejoin = envelope_.compression(
            refreshed.compression_rejoin_strain);
        refreshed.compression_rejoin_stress = compression_rejoin.stress;
        refreshed.compression_rejoin_tangent = compression_rejoin.tangent;
        populate_positive_rejoin_landmarks(envelope_, refreshed);

        if (strain > committed.strain) {
            // OpenSees Crule=10 positive reversal: Teb=current, Tea=ea1112f.
            ConcreteCMState reversal = refreshed;
            reversal.nested_positive_origin_strain = committed.strain; // Teb/Ter
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_negative_target_strain = committed.strain; // Teb
            reversal.nested_negative_target_uses_rule77 = false;
            const double denom = committed.unloading_strain - esplp;
            reversal.nested_positive_target_strain =
                committed.zero_stress_strain +
                ((committed.unloading_strain - committed.strain) / denom) *
                    (committed.tension_peak_strain - committed.zero_stress_strain); // Tea

            const auto target = positive_path_trial(
                envelope_, refreshed, reversal.nested_positive_target_strain);
            if (strain <= reversal.nested_positive_target_strain) {
                // r12f first normalizes its endpoint at Tea, then the caller
                // evaluates the actual trial strain with a second RAf/fcEturf.
                const auto normalized_target = smooth_transition(
                    reversal.nested_positive_target_strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    p.Ec,
                    reversal.nested_positive_target_strain,
                    target.response.stress,
                    target.response.tangent);
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    p.Ec,
                    reversal.nested_positive_target_strain,
                    normalized_target.stress,
                    normalized_target.tangent);
                return make_trial(reversal, strain, response, 1.0,
                                  ConcreteCMRule::NestedPositiveTarget);
            }
            return positive_path_trial(envelope_, reversal, strain);
        }

        return first_negative_return_trial(envelope_, refreshed, strain);
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
            if (committed.rule == ConcreteCMRule::CompressionUnloading) {
                const auto reversal = rule77_reversal_state(envelope_, committed);
                return rule77_trial(envelope_, reversal, strain);
            }
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
        (committed.rule == ConcreteCMRule::TensionEnvelope ||
         committed.rule == ConcreteCMRule::TensionRejoining)) {
        if (strain < committed.strain) {
            const auto reversal = second_reversal_state(envelope_, committed);
            return first_negative_return_trial(envelope_, reversal, strain);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    // OpenSees dispatch remains rule-driven after later rebounds.  A committed
    // rule 4 still continues on the ordinary negative-return path; if it
    // reverses positive, it creates a fresh Ter0p/Tfr0p and rule 88 exactly as
    // it does on the first cycle.
    if (committed.has_positive_to_negative_reversal &&
        committed.has_second_negative_to_positive_reversal &&
        committed.rule == ConcreteCMRule::TensionUnloading) {
        if (strain > committed.strain) {
            ConcreteCMState reversal = committed;
            reversal.positive_return_reversal_strain = committed.strain;
            reversal.positive_return_reversal_stress = committed.stress;
            reversal.nested_negative_target_strain = committed.strain;
            return rule88_trial(envelope_, reversal, strain);
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }

    if (committed.has_positive_to_negative_reversal &&
        !committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::TensionUnloading ||
         committed.rule == ConcreteCMRule::TensionToCompression ||
         committed.rule == ConcreteCMRule::CompressionRejoining)) {
        if (strain > committed.strain) {
            // OpenSees 3.8.0: a reversal from rule 7 follows the same
            // negative-to-positive refresh used by rules 1/5: update eunn/funn,
            // rebuild the shifted positive path with e0eunpfunpf, then select
            // rules 3/9/8/2/6 at the requested strain. QuakeCore's existing
            // second-negative-to-positive helper is the direct mapping.
            if (committed.rule == ConcreteCMRule::CompressionRejoining) {
                if (committed.has_second_negative_to_positive_reversal) {
                    throw std::logic_error(
                        "ConcreteCM deeper nested reversal is not yet admitted in Gate 4");
                }
                const auto reversal = second_negative_to_positive_reversal_state(
                    envelope_, committed);
                return positive_path_trial(envelope_, reversal, strain);
            }
            if (committed.rule == ConcreteCMRule::TensionUnloading) {
                ConcreteCMState reversal = committed;
                reversal.positive_return_reversal_strain = committed.strain; // Ter0p
                reversal.positive_return_reversal_stress = committed.stress; // Tfr0p
                reversal.nested_negative_target_strain = committed.strain; // Teb=Ter0p
                return rule88_trial(envelope_, reversal, strain);
            }

            // OpenSees 3.8.0 Crule=10 positive reversal. Save Teb at the
            // current point, derive Tea with ea1112f, and create rule 12 from
            // Ter/Tfr to the established positive path at Tea.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain; // Teb
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_negative_target_strain = committed.strain; // Teb
            reversal.nested_negative_target_uses_rule77 = false;
            const double denom =
                committed.unloading_strain - committed.positive_zero_stress_strain;
            reversal.nested_positive_target_strain =
                committed.zero_stress_strain +
                ((committed.unloading_strain - committed.strain) / denom) *
                    (committed.tension_peak_strain - committed.zero_stress_strain); // Tea
            const auto target = positive_path_trial(
                envelope_, committed, reversal.nested_positive_target_strain);
            if (strain <= reversal.nested_positive_target_strain) {
                // OpenSees r12f first evaluates the raw transition at Tea and
                // replaces ff/Ef with that evaluated endpoint. The caller then
                // runs RAf again before evaluating the actual trial strain.
                const auto normalized_target = smooth_transition(
                    reversal.nested_positive_target_strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    reversal.nested_positive_target_strain,
                    target.response.stress,
                    target.response.tangent);
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    reversal.nested_positive_target_strain,
                    normalized_target.stress,
                    normalized_target.tangent);
                return make_trial(reversal, strain, response, 1.0,
                                  ConcreteCMRule::NestedPositiveTarget);
            }
            return positive_path_trial(envelope_, reversal, strain);
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }

    if (committed.has_second_negative_to_positive_reversal &&
        committed.rule == ConcreteCMRule::TensionToCompression) {
        if (strain > committed.strain) {
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain; // Teb
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_negative_target_strain = committed.strain; // Teb
            reversal.nested_negative_target_uses_rule77 = false;
            const double denom =
                committed.unloading_strain - committed.positive_zero_stress_strain;
            reversal.nested_positive_target_strain =
                committed.zero_stress_strain +
                ((committed.unloading_strain - committed.strain) / denom) *
                    (committed.tension_peak_strain - committed.zero_stress_strain); // Tea
            const auto target = positive_path_trial(
                envelope_, committed, reversal.nested_positive_target_strain);
            if (strain <= reversal.nested_positive_target_strain) {
                const auto normalized_target = smooth_transition(
                    reversal.nested_positive_target_strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    reversal.nested_positive_target_strain,
                    target.response.stress,
                    target.response.tangent);
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    reversal.nested_positive_target_strain,
                    normalized_target.stress,
                    normalized_target.tangent);
                return make_trial(reversal, strain, response, 1.0,
                                  ConcreteCMRule::NestedPositiveTarget);
            }
            return positive_path_trial(envelope_, reversal, strain);
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
            // OpenSees 3.8.0 negative reversals after the second rebound are
            // governed by the committed rule, not by a separate history law.
            if (committed.rule == ConcreteCMRule::CompressionUnloading) {
                const auto reversal = rule77_reversal_state(envelope_, committed);
                return rule77_trial(envelope_, reversal, strain);
            }
            if (committed.rule == ConcreteCMRule::CompressionToTension) {
                throw std::logic_error("ConcreteCM second rebound negative reversal from rule9");
            }
            if (committed.rule == ConcreteCMRule::TensionRejoining ||
                committed.rule == ConcreteCMRule::TensionEnvelope) {
                // OpenSees Crule 2/8 always promotes the current point to
                // Teunp/Tfunp and restarts the ordinary 4->10->7->1/5 return.
                const auto reversal = second_reversal_state(envelope_, committed);
                return first_negative_return_trial(envelope_, reversal, strain);
            }
            throw std::logic_error("ConcreteCM second rebound negative reversal from rule6");
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::TensionEnvelope ||
        committed.rule == ConcreteCMRule::TensionCutoff) {
        if (strain < committed.strain) {
            if (committed.rule == ConcreteCMRule::TensionEnvelope) {
                // OpenSees 3.8.0 Crule=2 reversal: save the current positive
                // extreme as Teunp/Tfunp, construct the positive unloading
                // landmarks, then follow rules 4, 10, 7, and the compression
                // envelope. second_reversal_state + first_negative_return_trial
                // is QuakeCore's direct representation of that sequence.
                const auto reversal = second_reversal_state(envelope_, committed);
                return first_negative_return_trial(envelope_, reversal, strain);
            }
            throw std::logic_error(
                "ConcreteCM virgin tension reversal from rule6");
        }
        const auto response = envelope_.tension(strain, committed.tension_zero_strain);
        const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
            ? ConcreteCMRule::TensionCutoff
            : ConcreteCMRule::TensionEnvelope;
        return make_trial(committed, strain, response, 1.0, rule);
    }

    if (committed.rule == ConcreteCMRule::Initial) {
        throw std::logic_error("ConcreteCM fallback rule0 initial");
    }
    if (committed.rule == ConcreteCMRule::CompressionEnvelope) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule1 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule1 p2n history"
                : "ConcreteCM fallback rule1"));
    }
    if (committed.rule == ConcreteCMRule::TensionEnvelope) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule2 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule2 p2n history"
                : "ConcreteCM fallback rule2"));
    }
    if (committed.rule == ConcreteCMRule::CompressionUnloading) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule3 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule3 p2n history"
                : "ConcreteCM fallback rule3"));
    }
    if (committed.rule == ConcreteCMRule::TensionUnloading) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule4 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule4 p2n history"
                : "ConcreteCM fallback rule4"));
    }
    if (committed.rule == ConcreteCMRule::TensionCutoff) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule6 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule6 p2n history"
                : "ConcreteCM fallback rule6"));
    }
    if (committed.rule == ConcreteCMRule::CompressionRejoining) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule7 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule7 p2n history"
                : "ConcreteCM fallback rule7"));
    }
    if (committed.rule == ConcreteCMRule::TensionRejoining) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule8 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule8 p2n history"
                : "ConcreteCM fallback rule8"));
    }
    if (committed.rule == ConcreteCMRule::CompressionToTension) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule9 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule9 p2n history"
                : "ConcreteCM fallback rule9"));
    }
    if (committed.rule == ConcreteCMRule::TensionToCompression) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule10 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule10 p2n history"
                : "ConcreteCM fallback rule10"));
    }
    if (committed.rule == ConcreteCMRule::NestedNegativeTarget) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule11 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule11 p2n history"
                : "ConcreteCM fallback rule11"));
    }
    if (committed.rule == ConcreteCMRule::NestedPositiveTarget) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule12 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule12 p2n history"
                : "ConcreteCM fallback rule12"));
    }
    if (committed.rule == ConcreteCMRule::CompressionReversalTransition) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule77 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule77 p2n history"
                : "ConcreteCM fallback rule77"));
    }
    throw std::logic_error("ConcreteCM fallback enum value outside declared rules");
}

} // namespace quake