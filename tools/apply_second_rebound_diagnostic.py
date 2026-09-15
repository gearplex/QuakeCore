from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()
old = '''        if (strain < committed.strain) {
            throw std::logic_error(
                "ConcreteCM reversal from the second rebound is not yet admitted in Gate 4");
        }
'''
new = '''        if (strain < committed.strain) {
            if (committed.rule == ConcreteCMRule::CompressionUnloading) {
                throw std::logic_error("ConcreteCM second rebound negative reversal from rule3");
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
