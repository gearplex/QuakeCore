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

            // OpenSees 3.8.0 Crule=77, Cstrain>=Teunn. Tea is the original
            // rule-77 reversal strain (Ter0n), while Teb is this reversal
            // point. Rule 12 transitions from Ter/Tfr to Tea with Ec as the
            // initial tangent; Tea and Teb remain fixed across nested cycles.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain;
            reversal.nested_positive_origin_stress = committed.stress;
            reversal.nested_positive_target_strain = committed.positive_reversal_strain; // Tea=Ter0n
            reversal.nested_negative_target_strain = committed.strain; // Teb
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
        const auto positive_target = positive_path_trial(
            envelope_, committed, committed.nested_positive_target_strain);
        if (strain < committed.strain) {
            // OpenSees 3.8.0 Crule=12 reversal for the rule77->12 provenance:
            // Tea == Ter0n and Teb is the saved rule-77 reversal point. The
            // negative reversal creates rule 11 from current Ter/Tfr to Teb.
            ConcreteCMState reversal = committed;
            reversal.nested_negative_origin_strain = committed.strain;
            reversal.nested_negative_origin_stress = committed.stress;
            const double target_strain = committed.nested_negative_target_strain; // Teb
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
        }
        if (strain <= committed.nested_positive_target_strain) {
            const auto response = smooth_transition(
                strain,
                committed.nested_positive_origin_strain,
                committed.nested_positive_origin_stress,
                parameters().Ec,
                committed.nested_positive_target_strain,
                positive_target.response.stress,
                positive_target.response.tangent);
            return make_trial(committed, strain, response, 1.0,
                              ConcreteCMRule::NestedPositiveTarget);
        }
        return positive_path_trial(envelope_, committed, strain);
    }

    if (committed.rule == ConcreteCMRule::NestedNegativeTarget) {
        const double target_strain = committed.nested_negative_target_strain; // Teb
        const auto negative_target = rule77_trial(
            envelope_, committed, target_strain);
        if (strain > committed.strain) {
            // OpenSees 3.8.0 Crule=11 positive reversal for this provenance:
            // keep Tea/Teb fixed, move Ter/Tfr to the current point, and
            // create rule 12 targeting Tea on the established positive path.
            ConcreteCMState reversal = committed;
            reversal.nested_positive_origin_strain = committed.strain;
            reversal.nested_positive_origin_stress = committed.stress;
            const auto positive_target = positive_path_trial(
                envelope_, committed, committed.nested_positive_target_strain);
            if (strain <= committed.nested_positive_target_strain) {
                const auto response = smooth_transition(
                    strain,
                    reversal.nested_positive_origin_strain,
                    reversal.nested_positive_origin_stress,
                    parameters().Ec,
                    committed.nested_positive_target_strain,
                    positive_target.response.stress,
                    positive_target.response.tangent);
                return make_trial(reversal, strain, response, 1.0,
                                  ConcreteCMRule::NestedPositiveTarget);
            }
            return positive_path_trial(envelope_, reversal, strain);
        }
        if (strain >= target_strain) {
            const auto response = smooth_transition(
                strain,
                committed.nested_negative_origin_strain,
                committed.nested_negative_origin_stress,
                parameters().Ec,
                target_strain,
                negative_target.response.stress,
                negative_target.response.tangent);
            return make_trial(committed, strain, response, -1.0,
                              ConcreteCMRule::NestedNegativeTarget);
        }
        return rule77_trial(envelope_, committed, strain);
    }

'''
if text.count(anchor) != 1:
    raise SystemExit("expected first_rebound anchor not found exactly once")
text = text.replace(anchor, insert + anchor)

first_negative_return = '''    if (committed.has_positive_to_negative_reversal &&
        (committed.rule == ConcreteCMRule::TensionUnloading ||
         committed.rule == ConcreteCMRule::TensionToCompression ||
         committed.rule == ConcreteCMRule::CompressionRejoining)) {
        if (strain > committed.strain) {
            throw std::logic_error(
                "ConcreteCM reversal from the first negative return is not yet admitted in Gate 4");
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }
'''
first_negative_return_replacement = '''    if (committed.has_positive_to_negative_reversal &&
        (committed.rule == ConcreteCMRule::TensionUnloading ||
         committed.rule == ConcreteCMRule::TensionToCompression ||
         committed.rule == ConcreteCMRule::CompressionRejoining)) {
        if (strain > committed.strain) {
            // OpenSees 3.8.0: a reversal from rule 7 follows the same
            // negative-to-positive refresh used by rules 1/5: update eunn/funn,
            // rebuild the shifted positive path with e0eunpfunpf, then select
            // rules 3/9/8/2/6 at the requested strain. QuakeCore's existing
            // second-negative-to-positive helper is the direct mapping.
            if (committed.rule == ConcreteCMRule::CompressionRejoining) {
                if (committed.has_second_negative_to_positive_reversal) {
                    throw std::logic_error(
                        "ConcreteCM deeper nested reversal is not yet admitted in Gate 4");
                }
                const auto reversal = second_negative_to_positive_reversal_state(
                    envelope_, committed);
                return positive_path_trial(envelope_, reversal, strain);
            }
            if (committed.rule == ConcreteCMRule::TensionUnloading) {
                throw std::logic_error(
                    "ConcreteCM positive reversal from rule4 is not yet admitted in Gate 4");
            }
            throw std::logic_error(
                "ConcreteCM positive reversal from rule10 is not yet admitted in Gate 4");
        }
        return first_negative_return_trial(envelope_, committed, strain);
    }
'''
if text.count(first_negative_return) != 1:
    raise SystemExit("expected first-negative-return guard not found exactly once")
text = text.replace(first_negative_return, first_negative_return_replacement)
path.write_text(text)

wall_path = Path("src/wall_material.cpp")
wall = wall_path.read_text()
old_size = "constexpr int concrete_cm_state_size=28;"
if wall.count(old_size) != 1:
    raise SystemExit("expected ConcreteCM state-size anchor not found exactly once")
wall = wall.replace(old_size, "constexpr int concrete_cm_state_size=34;", 1)

old_encode = "    v[27]=s.has_second_negative_to_positive_reversal?1.0:0.0;\n"
new_encode = (
    "    v[27]=s.has_second_negative_to_positive_reversal?1.0:0.0;\n"
    "    v[28]=s.nested_positive_origin_strain;v[29]=s.nested_positive_origin_stress;"
    "v[30]=s.nested_positive_target_strain;\n"
    "    v[31]=s.nested_negative_origin_strain;v[32]=s.nested_negative_origin_stress;"
    "v[33]=s.nested_negative_target_strain;\n"
)
if wall.count(old_encode) != 1:
    raise SystemExit("expected ConcreteCM encode anchor not found exactly once")
wall = wall.replace(old_encode, new_encode, 1)

old_decode = "    s.has_second_negative_to_positive_reversal=v[27]!=0.0;return s;\n"
new_decode = (
    "    s.has_second_negative_to_positive_reversal=v[27]!=0.0;"
    "s.nested_positive_origin_strain=v[28];s.nested_positive_origin_stress=v[29];"
    "s.nested_positive_target_strain=v[30];"
    "s.nested_negative_origin_strain=v[31];s.nested_negative_origin_stress=v[32];"
    "s.nested_negative_target_strain=v[33];return s;\n"
)
if wall.count(old_decode) != 1:
    raise SystemExit("expected ConcreteCM decode anchor not found exactly once")
wall = wall.replace(old_decode, new_decode, 1)
wall_path.write_text(wall)
