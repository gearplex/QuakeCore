from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()

old_helper = '''ConcreteCMState second_negative_to_positive_reversal_state(
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
'''
new_helper = '''ConcreteCMState second_negative_to_positive_reversal_state(
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
'''
if text.count(old_helper) != 1:
    raise SystemExit("expected second-negative-to-positive helper not found exactly once")
text = text.replace(old_helper, new_helper, 1)

old_envelope_rebound = '''        if (committed.has_second_negative_to_positive_reversal) {
            throw std::logic_error(
                "ConcreteCM deeper nested reversal is not yet admitted in Gate 4");
        }
        if (committed.has_positive_to_negative_reversal) {
            const auto reversal = second_negative_to_positive_reversal_state(
                envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
'''
new_envelope_rebound = '''        if (committed.has_positive_to_negative_reversal ||
            committed.has_second_negative_to_positive_reversal) {
            const auto reversal = second_negative_to_positive_reversal_state(
                envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
'''
if text.count(old_envelope_rebound) != 1:
    raise SystemExit("expected compression-envelope deeper-rebound guard not found exactly once")
text = text.replace(old_envelope_rebound, new_envelope_rebound, 1)

# Rule 7 is also cycle-count agnostic in OpenSees: continued negative motion
# remains on the negative-return path; any positive reversal refreshes
# Teunn/Tfunn and the shifted positive path with e0eunpfunpf.
anchor = '''    const bool first_rebound = committed.unloading_strain < 0.0;
'''
rule7_dispatch = '''    if (committed.rule == ConcreteCMRule::CompressionRejoining) {
        if (strain > committed.strain) {
            const auto reversal = second_negative_to_positive_reversal_state(
                envelope_, committed);
            return positive_path_trial(envelope_, reversal, strain);
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }

'''
if text.count(anchor) != 1:
    raise SystemExit("expected first_rebound anchor not found exactly once")
text = text.replace(anchor, rule7_dispatch + anchor, 1)

path.write_text(text)
