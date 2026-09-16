from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()

anchor = '''    if (committed.rule == ConcreteCMRule::TensionReversalTransition) {
'''
insert = '''    if (committed.has_positive_to_negative_reversal &&
        !committed.has_second_negative_to_positive_reversal &&
        committed.rule == ConcreteCMRule::CompressionToTension) {
        if (strain < committed.strain) {
            // OpenSees 3.8.0 Crule=9 negative reversal. Tea is the current
            // reversal point; Teb is obtained from eb1112f and becomes the
            // rule-11 target on the established negative-return path.
            ConcreteCMState reversal = committed;
            reversal.nested_negative_origin_strain = committed.strain; // Ter
            reversal.nested_negative_origin_stress = committed.stress; // Tfr
            reversal.nested_positive_target_strain = committed.strain; // Tea
            const double denom =
                committed.tension_peak_strain - committed.zero_stress_strain;
            reversal.nested_negative_target_strain =
                committed.unloading_strain -
                ((committed.strain - committed.zero_stress_strain) / denom) *
                    (committed.unloading_strain - committed.positive_zero_stress_strain); // Teb

            const double target_strain = reversal.nested_negative_target_strain;
            const auto target = first_negative_return_trial(
                envelope_, committed, target_strain);
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
            return first_negative_return_trial(envelope_, reversal, strain);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

'''
if text.count(anchor) != 1:
    raise SystemExit("expected rule88 dispatch anchor not found exactly once")
path.write_text(text.replace(anchor, insert + anchor, 1))
