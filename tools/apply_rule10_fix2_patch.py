from pathlib import Path

path = Path("src/concrete_cm.cpp")
text = path.read_text()
old = '''ConcreteCMResponse rule10(const ConcreteCMState& state, double strain) {
    return smooth_transition(
        strain,
        state.positive_zero_stress_strain,
        0.0,
        state.positive_zero_stress_tangent,
        state.unloading_strain,
        state.compression_new_stress,
        state.compression_new_tangent);
}
'''
new = '''ConcreteCMResponse rule10(const ConcreteCMState& state, double strain) {
    auto response = smooth_transition(
        strain,
        state.positive_zero_stress_strain,
        0.0,
        state.positive_zero_stress_tangent,
        state.unloading_strain,
        state.compression_new_stress,
        state.compression_new_tangent);

    // OpenSees 3.8.0 ConcreteCM "Fix 2": when fcEturf collapses onto
    // the endpoint secant for rule 10, replace that degenerate branch by
    // either the zero-stress gap or the negative-side Enewn line through
    // espln. This is an explicit post-processing step in ConcreteCM.cpp.
    const double endpoint_secant =
        state.compression_new_stress /
        (state.unloading_strain - state.positive_zero_stress_strain);
    if (response.tangent == endpoint_secant) {
        if (strain >= state.zero_stress_strain) {
            return {0.0, 0.0};
        }
        return {
            state.compression_new_tangent * (strain - state.zero_stress_strain),
            state.compression_new_tangent,
        };
    }
    return response;
}
'''
if text.count(old) != 1:
    raise SystemExit("expected rule10 helper not found exactly once")
path.write_text(text.replace(old, new, 1))
