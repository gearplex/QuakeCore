from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()

anchor = '''    if (committed.rule == ConcreteCMRule::TensionReversalTransition) {
'''
insert = '''    if (committed.rule == ConcreteCMRule::CompressionToTension) {
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

'''
if text.count(anchor) != 1:
    raise SystemExit("expected rule88 dispatch anchor not found exactly once")
path.write_text(text.replace(anchor, insert + anchor, 1))
