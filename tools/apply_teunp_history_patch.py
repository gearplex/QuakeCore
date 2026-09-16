from pathlib import Path

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

path.write_text(text)
