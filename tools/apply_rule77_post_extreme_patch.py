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
            throw std::logic_error(
                "ConcreteCM rule77 positive reversal nested-targeting branch is not yet admitted in Gate 4");
        }
'''
if text.count(old) != 1:
    raise SystemExit("expected rule77 diagnostic block not found exactly once")
path.write_text(text.replace(old, new))
