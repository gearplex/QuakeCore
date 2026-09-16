from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()
old = '''    throw std::logic_error("ConcreteCM committed rule is not implemented in Gate 4");
'''
new = '''    if (committed.rule == ConcreteCMRule::Initial) {
        throw std::logic_error("ConcreteCM fallback rule0 initial");
    }
    if (committed.rule == ConcreteCMRule::CompressionEnvelope) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule1 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule1 p2n history"
                : "ConcreteCM fallback rule1"));
    }
    if (committed.rule == ConcreteCMRule::TensionEnvelope) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule2 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule2 p2n history"
                : "ConcreteCM fallback rule2"));
    }
    if (committed.rule == ConcreteCMRule::CompressionUnloading) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule3 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule3 p2n history"
                : "ConcreteCM fallback rule3"));
    }
    if (committed.rule == ConcreteCMRule::TensionUnloading) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule4 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule4 p2n history"
                : "ConcreteCM fallback rule4"));
    }
    if (committed.rule == ConcreteCMRule::TensionCutoff) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule6 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule6 p2n history"
                : "ConcreteCM fallback rule6"));
    }
    if (committed.rule == ConcreteCMRule::CompressionRejoining) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule7 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule7 p2n history"
                : "ConcreteCM fallback rule7"));
    }
    if (committed.rule == ConcreteCMRule::TensionRejoining) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule8 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule8 p2n history"
                : "ConcreteCM fallback rule8"));
    }
    if (committed.rule == ConcreteCMRule::CompressionToTension) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule9 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule9 p2n history"
                : "ConcreteCM fallback rule9"));
    }
    if (committed.rule == ConcreteCMRule::TensionToCompression) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule10 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule10 p2n history"
                : "ConcreteCM fallback rule10"));
    }
    if (committed.rule == ConcreteCMRule::NestedNegativeTarget) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule11 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule11 p2n history"
                : "ConcreteCM fallback rule11"));
    }
    if (committed.rule == ConcreteCMRule::NestedPositiveTarget) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule12 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule12 p2n history"
                : "ConcreteCM fallback rule12"));
    }
    if (committed.rule == ConcreteCMRule::CompressionReversalTransition) {
        throw std::logic_error(committed.has_second_negative_to_positive_reversal
            ? "ConcreteCM fallback rule77 second-rebound history"
            : (committed.has_positive_to_negative_reversal
                ? "ConcreteCM fallback rule77 p2n history"
                : "ConcreteCM fallback rule77"));
    }
    throw std::logic_error("ConcreteCM fallback enum value outside declared rules");
'''
if text.count(old) != 1:
    raise SystemExit("expected ConcreteCM fallback not found exactly once")
path.write_text(text.replace(old, new, 1))
