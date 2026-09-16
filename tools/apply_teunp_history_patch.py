from pathlib import Path

header_path = Path("include/quake/concrete_cm.hpp")
header = header_path.read_text()
old_enum = '''    NestedPositiveTarget = 12,
    CompressionReversalTransition = 77,
'''
new_enum = '''    NestedPositiveTarget = 12,
    CompressionReversalTransition = 77,
    TensionReversalTransition = 88,
'''
if header.count(old_enum) != 1:
    raise SystemExit("expected ConcreteCM enum anchor not found exactly once")
header = header.replace(old_enum, new_enum, 1)
old_state_anchor = '''    double compression_rejoin_tangent{};

    // OpenSees nested rule-12/rule-11 history for the rule77 provenance.
'''
new_state_anchor = '''    double compression_rejoin_tangent{};

    // OpenSees primary Ter0p/Tfr0p history introduced by a positive reversal
    // from rule 4. Derived rule-88 landmarks are reconstructed from these
    // values and Teunp/Tfunp rather than persisted.
    double positive_return_reversal_strain{};
    double positive_return_reversal_stress{};

    // OpenSees nested rule-12/rule-11 history for the rule77 provenance.
'''
if header.count(old_state_anchor) != 1:
    raise SystemExit("expected ConcreteCM state anchor not found exactly once")
header_path.write_text(header.replace(old_state_anchor, new_state_anchor, 1))

path = Path("src/concrete_cm.cpp")
text = path.read_text()

old = '''    reversal.has_positive_to_negative_reversal = true;
    reversal.positive_reversal_strain = committed.strain;
    reversal.positive_reversal_stress = envelope_kernel.tension(
        committed.strain, committed.tension_zero_strain).stress;

    const double Esecp = tension_secant(
'''
new = '''    reversal.has_positive_to_negative_reversal = true;

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
'''
if text.count(old) != 1:
    raise SystemExit("expected second_reversal_state positive-history block not found exactly once")
text = text.replace(old, new, 1)

old_landmarks = '''    reversal.positive_zero_stress_tangent = p.gap_close
        ? p.Ec / (std::pow(std::abs(
              (reversal.positive_reversal_strain - committed.tension_zero_strain) / p.et), 1.1) + 1.0)
        : 0.0;

    const double delfn = committed.unloading_strain <= p.epsc / 10.0
'''
new_landmarks = '''    reversal.positive_zero_stress_tangent = p.gap_close
        ? p.Ec / (std::pow(std::abs(
              (reversal.positive_reversal_strain - committed.tension_zero_strain) / p.et), 1.1) + 1.0)
        : 0.0;

    // OpenSees reconstructs fnewp/Enewp/esrep/frep from the newly committed
    // Teunp/Tfunp at the start of every subsequent cyclic trial.
    populate_positive_rejoin_landmarks(envelope_kernel, reversal);

    const double delfn = committed.unloading_strain <= p.epsc / 10.0
'''
if text.count(old_landmarks) != 1:
    raise SystemExit("expected second_reversal_state landmark block not found exactly once")
text = text.replace(old_landmarks, new_landmarks, 1)

# Rule 88 is OpenSees' shortened positive return created by a positive reversal
# from rule 4. Persist only Ter0p/Tfr0p; reconstruct fnewstp/Enewstp/esrestp
# and the rejoin response from the primary history on every trial.
rule88_anchor = '''ConcreteCMTrial positive_path_trial(const ConcreteCMEnvelope& envelope_kernel,
'''
rule88_helpers = '''struct PositiveReturnLandmarks {
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

'''
if text.count(rule88_anchor) != 1:
    raise SystemExit("expected positive_path_trial anchor not found exactly once")
text = text.replace(rule88_anchor, rule88_helpers + rule88_anchor, 1)

# OpenSees uses the same Crule 2/8 negative-reversal block regardless of the
# earlier cyclic provenance: Teunp/Tfunp are promoted to the current point,
# then the return proceeds through rules 4, 10, 7, and the compression
# envelope.  After a nested rule-12 excursion QuakeCore can therefore be on
# rule 8 while has_positive_to_negative_reversal is already true.
anchor = '''    if (committed.has_positive_to_negative_reversal &&
        !committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::TensionUnloading ||
'''
rule28_rereversal = '''    if (committed.has_positive_to_negative_reversal &&
        !committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::TensionEnvelope ||
         committed.rule == ConcreteCMRule::TensionRejoining)) {
        if (strain < committed.strain) {
            const auto reversal = second_reversal_state(envelope_, committed);
            return first_negative_return_trial(envelope_, reversal, strain);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

'''
if text.count(anchor) != 1:
    raise SystemExit("expected first-negative-return dispatch anchor not found exactly once")
text = text.replace(anchor, rule28_rereversal + anchor, 1)

# OpenSees Crule=4 positive reversal saves Ter0p/Tfr0p and constructs rule 88.
# The rule-10 patch intentionally left this branch as a diagnostic until the
# primary Teunp/Tfunp history was corrected.
rule4_guard = '''            if (committed.rule == ConcreteCMRule::TensionUnloading) {
                throw std::logic_error(
                    "ConcreteCM positive reversal from rule4 is not yet admitted in Gate 4");
            }
'''
rule4_port = '''            if (committed.rule == ConcreteCMRule::TensionUnloading) {
                ConcreteCMState reversal = committed;
                reversal.positive_return_reversal_strain = committed.strain; // Ter0p
                reversal.positive_return_reversal_stress = committed.stress; // Tfr0p
                reversal.nested_negative_target_strain = committed.strain; // Teb=Ter0p
                return rule88_trial(envelope_, reversal, strain);
            }
'''
if text.count(rule4_guard) != 1:
    raise SystemExit("expected rule4 positive-reversal diagnostic not found exactly once")
text = text.replace(rule4_guard, rule4_port, 1)

# Continue rule 88 in the positive direction. Its negative reversal has a
# separate OpenSees branch and remains explicitly diagnosed until demanded.
trial_anchor = '''    const bool first_rebound = committed.unloading_strain < 0.0;
'''
rule88_dispatch = '''    if (committed.rule == ConcreteCMRule::TensionReversalTransition) {
        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM negative reversal from rule88 is not yet admitted in Gate 4");
        }
        return rule88_trial(envelope_, committed, strain);
    }

'''
if text.count(trial_anchor) != 1:
    raise SystemExit("expected first_rebound trial anchor not found exactly once")
text = text.replace(trial_anchor, rule88_dispatch + trial_anchor, 1)

path.write_text(text)
