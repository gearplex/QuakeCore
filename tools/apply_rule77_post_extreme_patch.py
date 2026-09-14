from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()
old = '''        if (strain > committed.strain) {
            if (committed.strain < committed.unloading_strain) {
                throw std::logic_error(
                    "ConcreteCM rule77 positive reversal: Cstrain < Teunn");
            }
            throw std::logic_error(
                "ConcreteCM rule77 positive reversal: Cstrain >= Teunn");
        }
'''
new = '''        if (strain > committed.strain) {
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

            // OpenSees 3.8.0 Crule=77, Cstrain>=Teunn.  Tea is the original
            // rule-77 reversal strain (Ter0n); rule 12 transitions from the
            // current reversal point (Ter/Tfr) to the response of the existing
            // positive path at Tea, with Ec as the initial tangent.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain;
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_positive_target_strain = committed.positive_reversal_strain;
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
'''
if text.count(old) != 1:
    raise SystemExit("expected rule77 diagnostic block not found exactly once")
text = text.replace(old, new)

anchor = '''    const bool first_rebound = committed.unloading_strain < 0.0;
'''
insert = '''    if (committed.rule == ConcreteCMRule::NestedPositiveTarget) {
        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM reversal from rule12 is not yet admitted in Gate 4");
        }
        const auto target = positive_path_trial(
            envelope_, committed, committed.nested_positive_target_strain);
        if (strain <= committed.nested_positive_target_strain) {
            const auto response = smooth_transition(
                strain,
                committed.nested_positive_origin_strain,
                committed.nested_positive_origin_stress,
                parameters().Ec,
                committed.nested_positive_target_strain,
                target.response.stress,
                target.response.tangent);
            return make_trial(committed, strain, response, 1.0,
                              ConcreteCMRule::NestedPositiveTarget);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

'''
if text.count(anchor) != 1:
    raise SystemExit("expected first_rebound anchor not found exactly once")
text = text.replace(anchor, insert + anchor)
path.write_text(text)
