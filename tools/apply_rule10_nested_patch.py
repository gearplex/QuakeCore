from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()

old_rule10 = '''            if (committed.rule == ConcreteCMRule::TensionUnloading) {
                throw std::logic_error(
                    "ConcreteCM positive reversal from rule4 is not yet admitted in Gate 4");
            }
            throw std::logic_error(
                "ConcreteCM positive reversal from rule10 is not yet admitted in Gate 4");
'''
new_rule10 = '''            if (committed.rule == ConcreteCMRule::TensionUnloading) {
                throw std::logic_error(
                    "ConcreteCM positive reversal from rule4 is not yet admitted in Gate 4");
            }

            // OpenSees 3.8.0 Crule=10 positive reversal. Save Teb at the
            // current point, derive Tea with ea1112f, and create rule 12 from
            // Ter/Tfr to the established positive path at Tea.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain; // Teb
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_negative_target_strain = committed.strain; // Teb
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
'''
if text.count(old_rule10) != 1:
    raise SystemExit("expected rule10 diagnostic not found exactly once")
text = text.replace(old_rule10, new_rule10, 1)

old_rule12_negative = '''            const double target_strain = committed.nested_negative_target_strain; // Teb
            const auto target = rule77_trial(envelope_, committed, target_strain);
            if (strain >= target_strain) {
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_negative_origin_strain,
                    reversal.nested_negative_origin_stress,
                    parameters().Ec,
                    target_strain,
                    target.response.stress,
                    target.response.tangent);
                return make_trial(reversal, strain, response, -1.0,
                                  ConcreteCMRule::NestedNegativeTarget);
            }
            return rule77_trial(envelope_, reversal, strain);
'''
new_rule12_negative = '''            const double target_strain = committed.nested_negative_target_strain; // Teb
            const bool from_first_negative_return =
                committed.has_positive_to_negative_reversal;
            const auto target = from_first_negative_return
                ? first_negative_return_trial(envelope_, committed, target_strain)
                : rule77_trial(envelope_, committed, target_strain);
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
            return from_first_negative_return
                ? first_negative_return_trial(envelope_, reversal, strain)
                : rule77_trial(envelope_, reversal, strain);
'''
if text.count(old_rule12_negative) != 1:
    raise SystemExit("expected rule12 negative target block not found exactly once")
text = text.replace(old_rule12_negative, new_rule12_negative, 1)

old_rule11_target = '''        const auto negative_target = rule77_trial(
            envelope_, committed, target_strain);
'''
new_rule11_target = '''        const bool from_first_negative_return =
            committed.has_positive_to_negative_reversal;
        const auto negative_target = from_first_negative_return
            ? first_negative_return_trial(envelope_, committed, target_strain)
            : rule77_trial(envelope_, committed, target_strain);
'''
if text.count(old_rule11_target) != 1:
    raise SystemExit("expected rule11 negative target block not found exactly once")
text = text.replace(old_rule11_target, new_rule11_target, 1)

old_rule11_continue = '''        return rule77_trial(envelope_, committed, strain);
    }

    const bool first_rebound = committed.unloading_strain < 0.0;
'''
new_rule11_continue = '''        return from_first_negative_return
            ? first_negative_return_trial(envelope_, committed, strain)
            : rule77_trial(envelope_, committed, strain);
    }

    const bool first_rebound = committed.unloading_strain < 0.0;
'''
if text.count(old_rule11_continue) != 1:
    raise SystemExit("expected rule11 continuation block not found exactly once")
text = text.replace(old_rule11_continue, new_rule11_continue, 1)

# Rule identity, rather than QuakeCore's history flag, governs OpenSees Crule=10.
# A second-rebound history can re-enter rule 10 after a nested rule-11 target.
# Continue down the established first-negative-return path for negative motion;
# on positive reversal, apply the same ea1112f/rule-12 construction used above.
second_history_anchor = '''    if (committed.has_second_negative_to_positive_reversal &&
        (committed.rule == ConcreteCMRule::CompressionUnloading ||
'''
second_history_rule10 = '''    if (committed.has_second_negative_to_positive_reversal &&
        committed.rule == ConcreteCMRule::TensionToCompression) {
        if (strain > committed.strain) {
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain; // Teb
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_negative_target_strain = committed.strain; // Teb
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

'''
if text.count(second_history_anchor) != 1:
    raise SystemExit("expected second-rebound block anchor not found exactly once")
text = text.replace(second_history_anchor, second_history_rule10 + second_history_anchor, 1)

path.write_text(text)
