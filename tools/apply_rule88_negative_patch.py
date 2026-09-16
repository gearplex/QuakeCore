from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()

old = '''    if (committed.rule == ConcreteCMRule::TensionReversalTransition) {
        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM negative reversal from rule88 is not yet admitted in Gate 4");
        }
        return rule88_trial(envelope_, committed, strain);
    }
'''
new = '''    if (committed.rule == ConcreteCMRule::TensionReversalTransition) {
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
'''
if text.count(old) != 1:
    raise SystemExit("expected rule88 negative-reversal diagnostic not found exactly once")
path.write_text(text.replace(old, new, 1))
