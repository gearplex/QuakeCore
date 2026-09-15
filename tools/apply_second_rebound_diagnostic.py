from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()
old = '''        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM reversal from the second rebound is not yet admitted in Gate 4");
        }
'''
new = '''        if (strain < committed.strain) {
            // OpenSees 3.8.0 Crule=3 negative reversal: capture the current
            // rule-3 point as Ter0n/Tfr0n, recompute the shortened negative
            // rejoin landmarks, then follow rule 77 or the compression envelope.
            if (committed.rule == ConcreteCMRule::CompressionUnloading) {
                const auto reversal = rule77_reversal_state(envelope_, committed);
                return rule77_trial(envelope_, reversal, strain);
            }
            if (committed.rule == ConcreteCMRule::CompressionToTension) {
                throw std::logic_error("ConcreteCM second rebound negative reversal from rule9");
            }
            if (committed.rule == ConcreteCMRule::TensionRejoining) {
                throw std::logic_error("ConcreteCM second rebound negative reversal from rule8");
            }
            if (committed.rule == ConcreteCMRule::TensionEnvelope) {
                throw std::logic_error("ConcreteCM second rebound negative reversal from rule2");
            }
            throw std::logic_error("ConcreteCM second rebound negative reversal from rule6");
        }
'''
if text.count(old) != 1:
    raise SystemExit("expected second-rebound guard not found exactly once")
text = text.replace(old, new, 1)

old_virgin = '''    if (committed.rule == ConcreteCMRule::TensionEnvelope ||
        committed.rule == ConcreteCMRule::TensionCutoff) {
        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM virgin tension reversal rules are not yet admitted in Gate 4");
        }
        const auto response = envelope_.tension(strain, committed.tension_zero_strain);
        const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
            ? ConcreteCMRule::TensionCutoff
            : ConcreteCMRule::TensionEnvelope;
        return make_trial(committed, strain, response, 1.0, rule);
    }
'''
new_virgin = '''    if (committed.rule == ConcreteCMRule::TensionEnvelope ||
        committed.rule == ConcreteCMRule::TensionCutoff) {
        if (strain < committed.strain) {
            if (committed.rule == ConcreteCMRule::TensionEnvelope) {
                // OpenSees 3.8.0 Crule=2 reversal: save the current positive
                // extreme as Teunp/Tfunp, construct the positive unloading
                // landmarks, then follow rules 4, 10, 7, and the compression
                // envelope. second_reversal_state + first_negative_return_trial
                // is QuakeCore's direct representation of that sequence.
                const auto reversal = second_reversal_state(envelope_, committed);
                return first_negative_return_trial(envelope_, reversal, strain);
            }
            throw std::logic_error(
                "ConcreteCM virgin tension reversal from rule6");
        }
        const auto response = envelope_.tension(strain, committed.tension_zero_strain);
        const auto rule = (response.stress == 0.0 && response.tangent == 0.0)
            ? ConcreteCMRule::TensionCutoff
            : ConcreteCMRule::TensionEnvelope;
        return make_trial(committed, strain, response, 1.0, rule);
    }
'''
if text.count(old_virgin) != 1:
    raise SystemExit("expected virgin tension reversal block not found exactly once")
text = text.replace(old_virgin, new_virgin, 1)
path.write_text(text)
