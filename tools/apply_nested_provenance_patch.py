from pathlib import Path

# Runtime-only Gate-4 patch: distinguish rule77-origin nested targets from
# ordinary first-negative-return targets. OpenSees decides the rule-11/12
# target path from Tea/Teb provenance, not from whether any earlier p2n event
# has occurred. The latter can remain true across later rule77 cycles.

header_path = Path("include/quake/concrete_cm.hpp")
header = header_path.read_text()
old_header = '''    double nested_negative_origin_strain{};
    double nested_negative_origin_stress{};

    // Second negative-to-positive reversal. Gate 4 currently admits only
'''
new_header = '''    double nested_negative_origin_strain{};
    double nested_negative_origin_stress{};
    bool nested_negative_target_uses_rule77{false};

    // Second negative-to-positive reversal. Gate 4 currently admits only
'''
if header.count(old_header) != 1:
    raise SystemExit("expected nested provenance header anchor not found exactly once")
header_path.write_text(header.replace(old_header, new_header, 1))

path = Path("src/concrete_cm.cpp")
text = path.read_text()

# Mark the rule77 -> rule12 path explicitly. Change the comment as well so the
# generic rule10 replacement below cannot overwrite this provenance.
old_rule77 = '''            reversal.nested_positive_target_strain = committed.positive_reversal_strain; // Tea=Ter0n
            reversal.nested_negative_target_strain = committed.strain; // Teb
'''
new_rule77 = '''            reversal.nested_positive_target_strain = committed.positive_reversal_strain; // Tea=Ter0n
            reversal.nested_negative_target_strain = committed.strain; // Teb (rule77)
            reversal.nested_negative_target_uses_rule77 = true;
'''
if text.count(old_rule77) != 1:
    raise SystemExit("expected rule77 nested-target provenance anchor not found exactly once")
text = text.replace(old_rule77, new_rule77, 1)

# Rule10 -> rule12 targets the ordinary first-negative-return path. There are
# two demanded runtime dispatches (ordinary and later-history rule10).
old_rule10_target = '''            reversal.nested_negative_target_strain = committed.strain; // Teb
'''
rule10_count = text.count(old_rule10_target)
if rule10_count < 1:
    raise SystemExit("expected at least one rule10 nested target")
text = text.replace(
    old_rule10_target,
    '''            reversal.nested_negative_target_strain = committed.strain; // Teb
            reversal.nested_negative_target_uses_rule77 = false;
''',
)

# Rule88 negative reversal -> rule11 targets the ordinary first-negative-return
# path at Ter0p (rule 4 for the demanded provenance).
old_rule88 = '''                reversal.nested_negative_target_strain =
                    committed.positive_return_reversal_strain; // Teb=Ter0p
'''
new_rule88 = '''                reversal.nested_negative_target_strain =
                    committed.positive_return_reversal_strain; // Teb=Ter0p
                reversal.nested_negative_target_uses_rule77 = false;
'''
if text.count(old_rule88) != 1:
    raise SystemExit("expected rule88 target provenance anchor not found exactly once")
text = text.replace(old_rule88, new_rule88, 1)

# Rule9 negative reversal -> rule11 also targets the ordinary negative-return
# path at eb1112f's Teb. Explicitly clear any rule77 provenance inherited in
# the copied state.
old_rule9 = '''            reversal.nested_negative_target_strain =
                committed.unloading_strain -
                ((committed.strain - committed.zero_stress_strain) / denom) *
                    (committed.unloading_strain - esplp); // Teb = eb1112f(...)
'''
new_rule9 = '''            reversal.nested_negative_target_strain =
                committed.unloading_strain -
                ((committed.strain - committed.zero_stress_strain) / denom) *
                    (committed.unloading_strain - esplp); // Teb = eb1112f(...)
            reversal.nested_negative_target_uses_rule77 = false;
'''
if text.count(old_rule9) != 1:
    raise SystemExit("expected rule9 target provenance anchor not found exactly once")
text = text.replace(old_rule9, new_rule9, 1)

# Replace the heuristic target-path selection in rule12 negative reversal.
old_rule12_select = '''            const bool from_first_negative_return =
                committed.has_positive_to_negative_reversal;
            const auto target = from_first_negative_return
                ? first_negative_return_trial(envelope_, committed, target_strain)
                : rule77_trial(envelope_, committed, target_strain);
'''
new_rule12_select = '''            const bool from_rule77 =
                committed.nested_negative_target_uses_rule77;
            const auto target = from_rule77
                ? rule77_trial(envelope_, committed, target_strain)
                : first_negative_return_trial(envelope_, committed, target_strain);
'''
if text.count(old_rule12_select) != 1:
    raise SystemExit("expected rule12 provenance selector not found exactly once")
text = text.replace(old_rule12_select, new_rule12_select, 1)
old_rule12_continue = '''            return from_first_negative_return
                ? first_negative_return_trial(envelope_, reversal, strain)
                : rule77_trial(envelope_, reversal, strain);
'''
new_rule12_continue = '''            return from_rule77
                ? rule77_trial(envelope_, reversal, strain)
                : first_negative_return_trial(envelope_, reversal, strain);
'''
if text.count(old_rule12_continue) != 1:
    raise SystemExit("expected rule12 provenance continuation not found exactly once")
text = text.replace(old_rule12_continue, new_rule12_continue, 1)

# Replace the same heuristic in rule11 continuation/reversal.
old_rule11_select = '''        const bool from_first_negative_return =
            committed.has_positive_to_negative_reversal;
        const auto negative_target = from_first_negative_return
            ? first_negative_return_trial(envelope_, committed, target_strain)
            : rule77_trial(envelope_, committed, target_strain);
'''
new_rule11_select = '''        const bool from_rule77 =
            committed.nested_negative_target_uses_rule77;
        const auto negative_target = from_rule77
            ? rule77_trial(envelope_, committed, target_strain)
            : first_negative_return_trial(envelope_, committed, target_strain);
'''
if text.count(old_rule11_select) != 1:
    raise SystemExit("expected rule11 provenance selector not found exactly once")
text = text.replace(old_rule11_select, new_rule11_select, 1)
old_rule11_continue = '''        return from_first_negative_return
            ? first_negative_return_trial(envelope_, committed, strain)
            : rule77_trial(envelope_, committed, strain);
'''
new_rule11_continue = '''        return from_rule77
            ? rule77_trial(envelope_, committed, strain)
            : first_negative_return_trial(envelope_, committed, strain);
'''
if text.count(old_rule11_continue) != 1:
    raise SystemExit("expected rule11 provenance continuation not found exactly once")
text = text.replace(old_rule11_continue, new_rule11_continue, 1)

path.write_text(text)

# Extend runtime material-state persistence by one scalar. The preceding
# rule77 patch expands state to 36 values and the rule88 patch uses 34/35.
wall_path = Path("src/wall_material.cpp")
wall = wall_path.read_text()
old_size = "constexpr int concrete_cm_state_size=36;"
if wall.count(old_size) != 1:
    raise SystemExit("expected 36-value ConcreteCM state size not found exactly once")
wall = wall.replace(old_size, "constexpr int concrete_cm_state_size=37;", 1)
old_encode = '''    v[34]=s.positive_return_reversal_strain;v[35]=s.positive_return_reversal_stress;\n'''
new_encode = '''    v[34]=s.positive_return_reversal_strain;v[35]=s.positive_return_reversal_stress;\n    v[36]=s.nested_negative_target_uses_rule77?1.0:0.0;\n'''
if wall.count(old_encode) != 1:
    raise SystemExit("expected rule88 encode tail not found exactly once")
wall = wall.replace(old_encode, new_encode, 1)
old_decode = '''s.nested_negative_target_strain=v[33];s.positive_return_reversal_strain=v[34];s.positive_return_reversal_stress=v[35];return s;\n'''
new_decode = '''s.nested_negative_target_strain=v[33];s.positive_return_reversal_strain=v[34];s.positive_return_reversal_stress=v[35];s.nested_negative_target_uses_rule77=v[36]!=0.0;return s;\n'''
if wall.count(old_decode) != 1:
    raise SystemExit("expected rule88 decode tail not found exactly once")
wall_path.write_text(wall.replace(old_decode, new_decode, 1))
