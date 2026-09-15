from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()
old = '''    throw std::logic_error("ConcreteCM committed rule is not implemented in Gate 4");
'''
new = '''    if (committed.rule == ConcreteCMRule::TensionUnloading) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule4 second-rebound history"
            : "ConcreteCM fallback rule4");
    }
    if (committed.rule == ConcreteCMRule::TensionToCompression) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule10 second-rebound history"
            : "ConcreteCM fallback rule10");
    }
    if (committed.rule == ConcreteCMRule::CompressionRejoining) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule7 second-rebound history"
            : "ConcreteCM fallback rule7");
    }
    if (committed.rule == ConcreteCMRule::NestedNegativeTarget) {
        throw std::logic_error("ConcreteCM fallback rule11");
    }
    if (committed.rule == ConcreteCMRule::NestedPositiveTarget) {
        throw std::logic_error("ConcreteCM fallback rule12");
    }
    if (committed.rule == ConcreteCMRule::CompressionReversalTransition) {
        throw std::logic_error("ConcreteCM fallback rule77");
    }
    throw std::logic_error("ConcreteCM fallback unknown rule");
'''
if text.count(old) != 1:
    raise SystemExit("expected ConcreteCM fallback not found exactly once")
path.write_text(text.replace(old, new, 1))
