from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()

anchor = '''    const bool first_rebound = committed.unloading_strain < 0.0;
'''
insert = '''    if (committed.rule == ConcreteCMRule::TensionToCompression) {
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

'''
if text.count(anchor) != 1:
    raise SystemExit("expected first_rebound anchor not found exactly once")
path.write_text(text.replace(anchor, insert + anchor, 1))
