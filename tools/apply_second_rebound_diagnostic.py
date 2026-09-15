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
path.write_text(text.replace(old, new, 1))
