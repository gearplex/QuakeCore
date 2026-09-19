#include "quake/wall_material.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// These are analytic and internal lifecycle regressions, not an OpenSees oracle.
// No assert(): all checks remain active in Release builds with NDEBUG defined.
namespace {
using Material = quake::WallUniaxial;
using State = std::vector<double>;
constexpr double guard = 9.87654321e123;

void require(bool condition, const std::string& context) {
    if (!condition) throw std::runtime_error(context);
}
void near(double actual, double expected, const std::string& context,
          double atol = 1e-12, double rtol = 1e-11) {
    require(std::isfinite(actual) && std::isfinite(expected), context + ": nonfinite value");
    if (std::abs(actual - expected) > atol + rtol * std::max(std::abs(actual), std::abs(expected))) {
        throw std::runtime_error(context + ": actual=" + std::to_string(actual) +
                                 ", expected=" + std::to_string(expected));
    }
}
State initial(const Material& material) {
    const auto n = static_cast<std::size_t>(material.state_size());
    State buffer(n + 2, guard);
    material.initialize(buffer.data() + 1);
    require(buffer.front() == guard && buffer.back() == guard, "initialize overwrote state guards");
    State state(buffer.begin() + 1, buffer.end() - 1);
    for (double value : state) require(std::isfinite(value) && value != guard, "uninitialized state entry");
    return state;
}
struct Step { Material::Result response; State state; };
Step evaluate(const Material& material, double strain, const State& committed) {
    require(committed.size() == static_cast<std::size_t>(material.state_size()), "wrong input state size");
    const State snapshot = committed;
    State buffer(committed.size() + 2, guard);
    auto response = material.trial(strain, committed.data(), buffer.data() + 1);
    require(committed == snapshot, "trial evaluation mutated committed history");
    require(buffer.front() == guard && buffer.back() == guard, "trial overwrote state guards");
    State state(buffer.begin() + 1, buffer.end() - 1);
    for (double value : state) require(std::isfinite(value) && value != guard, "unwritten trial state entry");
    require(std::isfinite(response.stress) && std::isfinite(response.tangent), "nonfinite response");
    return {response, std::move(state)};
}
void equal_step(const Step& a, const Step& b, const std::string& context) {
    near(a.response.stress, b.response.stress, context + " stress");
    near(a.response.tangent, b.response.tangent, context + " tangent");
    require(a.state == b.state, context + " internal state");
}
template <typename Function>
void invalid(Function&& function, const std::string& context) {
    try { function(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error(context + ": expected invalid_argument");
}

void initialization_and_state_layout() {
    const auto steel = Material::steel(1000, 10, .1);
    const auto concrete = Material::concrete01(-4, -.002, -.8, -.006);
    const auto nested = Material::parallel({Material::minmax(steel, -.05, .05), concrete});
    require(nested.state_size() == 13, "heterogeneous state size");
    const auto s = initial(nested);
    require(s.size() == 13 && s[6] == 0, "MinMax flag offset");
    near(s[9], 4000, "Concrete01 unloading modulus offset");
    near(nested.initial_tangent(), 5000, "nested initial tangent");
    (void)evaluate(nested, -.001, s);
}
void exact_failure_boundaries() {
    const auto m = Material::minmax(Material::elastic(1024), -.03125, .015625);
    const auto s = initial(m);
    for (double strain : {-.03125, .015625}) {
        const auto t = evaluate(m, strain, s);
        near(t.response.stress, 0, "failure boundary stress");
        near(t.response.tangent, 1e-8 * 1024, "failure boundary tangent", 1e-16);
        require(t.state.back() == 1, "inclusive boundary must fail");
    }
}
void adjacent_failure_boundaries() {
    const double lo = -.03125, hi = .015625;
    const auto m = Material::minmax(Material::elastic(1024), lo, hi);
    const auto s = initial(m);
    for (double strain : {std::nextafter(lo, 0.), std::nextafter(hi, 0.)}) {
        const auto t = evaluate(m, strain, s);
        require(t.state.back() == 0, "one ULP inside must remain intact");
        near(t.response.stress, 1024 * strain, "inside boundary stress");
        near(t.response.tangent, 1024, "inside boundary tangent");
    }
    for (double strain : {std::nextafter(lo, -1.), std::nextafter(hi, 1.)}) {
        require(evaluate(m, strain, s).state.back() == 1, "one ULP outside must fail");
    }
}
void failure_uses_initial_not_current_tangent() {
    const auto m = Material::minmax(Material::steel(1000, 10, .1), -.05, .05);
    const auto yielded = evaluate(m, .03, initial(m));
    near(yielded.response.tangent, 100, "yielded child tangent");
    const auto failed = evaluate(m, .05, yielded.state);
    near(failed.response.tangent, 1e-5, "failure tangent uses initial child E", 1e-16);
    require(std::equal(yielded.state.begin(), yielded.state.end() - 1, failed.state.begin()),
            "failure must preserve child history");
}
void rejected_failure_after_yield() {
    const auto m = Material::minmax(Material::steel(1000, 10, .1), -.05, .05);
    const auto committed = evaluate(m, .03, initial(m)).state;
    const auto reference = evaluate(m, .025, committed);
    const auto rejected = evaluate(m, .06, committed);
    require(rejected.state.back() == 1, "rejected trial did not fail");
    equal_step(evaluate(m, .025, committed), reference, "discard failed trial");
    near(reference.response.stress, 7, "nonlinear-history unload stress");
    near(reference.response.tangent, 1000, "nonlinear-history unload tangent");
}
void committed_failure_is_irreversible() {
    const auto m = Material::minmax(Material::steel(1000, 10, .1), -.05, .05);
    const auto yielded = evaluate(m, .03, initial(m));
    State committed = evaluate(m, .05, yielded.state).state;
    const State failed_snapshot = committed;
    for (double strain : {.02, -.02, 0., -.1, .1}) {
        const auto t = evaluate(m, strain, committed);
        near(t.response.stress, 0, "committed failure stress");
        near(t.response.tangent, 1e-5, "committed failure tangent", 1e-16);
        require(t.state == failed_snapshot, "failed child history changed");
        committed = t.state;
    }
}
void reinitialization_clears_nested_failure() {
    const auto m = Material::parallel({Material::minmax(Material::steel(1000, 10, .1), -.05, .05),
                                     Material::elastic(25)});
    const auto virgin = initial(m);
    State failed = evaluate(m, .06, virgin).state;
    require(failed != virgin, "test must begin with changed state");
    m.initialize(failed.data());
    require(failed == virgin, "initialize did not clear failure/history");
    equal_step(evaluate(m, .025, failed), evaluate(m, .025, virgin), "reset replay");
}
void bilinear_analytic_trace() {
    const auto m = Material::steel(1000, 10, .1);
    State committed = initial(m);
    struct Point { double strain, stress, tangent; };
    const std::vector<Point> points{{0,0,1000},{.005,5,1000},{.03,12,100},
                                    {.025,7,1000},{0,-9,100},{-.03,-12,100},{0,9,100}};
    for (const auto& point : points) {
        const auto t = evaluate(m, point.strain, committed);
        near(t.response.stress, point.stress, "analytic kinematic-hardening stress");
        near(t.response.tangent, point.tangent, "analytic kinematic-hardening tangent");
        committed = t.state;
    }
}
void exact_yield_tangent_convention() {
    const auto m = Material::steel(1024, 8, .125);
    const auto virgin = initial(m);
    for (double strain : {8./1024, -8./1024}) {
        const auto t = evaluate(m, strain, virgin);
        near(t.response.stress, 1024 * strain, "yield force");
        near(t.response.tangent, 128, "exact-yield post-yield tangent");
    }
}
void nested_minmax_failure_scope() {
    const auto inner = Material::minmax(Material::elastic(1000), -.02, .02);
    const auto m = Material::minmax(inner, -.05, .05);
    auto t = evaluate(m, .03, initial(m));
    require(t.state[6] == 1 && t.state[7] == 0, "inner-only failure flags");
    near(t.response.tangent, 1e-5, "nested residual stiffness is not multiplied twice", 1e-16);
    t = evaluate(m, .01, t.state);
    require(t.state[6] == 1 && t.state[7] == 0, "inner failure must remain committed");
    near(t.response.stress, 0, "nested irreversible failure stress");
}
void parallel_analytic_elastic_response() {
    const auto m = Material::parallel({Material::elastic(1024), Material::elastic(256), Material::elastic(64)});
    auto committed = initial(m);
    for (double strain : {0., .01, -.02, .04, 0.}) {
        const auto t = evaluate(m, strain, committed);
        near(t.response.stress, 1344 * strain, "parallel elastic stress");
        near(t.response.tangent, 1344, "parallel elastic tangent");
        committed = t.state;
    }
}
void heterogeneous_parallel_history_and_work() {
    const std::vector<Material> children{
        Material::minmax(Material::steel(1000, 10, .1), -.045, .035),
        Material::concrete01(-4, -.002, -.8, -.006),
        Material::elastic(17)};
    const auto m = Material::parallel(children);
    State committed = initial(m);
    std::vector<State> states;
    for (const auto& child : children) states.push_back(initial(child));
    double previous_strain = 0, previous_stress = 0, total_work = 0;
    std::vector<double> last_stress(children.size(), 0), work(children.size(), 0);
    const std::vector<double> path{0,.001,.015,-.001,-.003,-.0005,.03,-.025,.04,.015,-.02,0};
    for (double strain : path) {
        const auto t = evaluate(m, strain, committed);
        double stress = 0, tangent = 0;
        State packed;
        for (std::size_t i = 0; i < children.size(); ++i) {
            const auto c = evaluate(children[i], strain, states[i]);
            stress += c.response.stress;
            tangent += c.response.tangent;
            packed.insert(packed.end(), c.state.begin(), c.state.end());
            work[i] += .5 * (last_stress[i] + c.response.stress) * (strain - previous_strain);
            last_stress[i] = c.response.stress;
            states[i] = c.state;
        }
        near(t.response.stress, stress, "heterogeneous stress additivity");
        near(t.response.tangent, tangent, "heterogeneous tangent additivity");
        require(t.state == packed, "heterogeneous state packing");
        total_work += .5 * (previous_stress + t.response.stress) * (strain - previous_strain);
        double child_work = 0;
        for (double value : work) child_work += value;
        near(total_work, child_work, "trapezoidal work additivity");
        previous_strain = strain;
        previous_stress = t.response.stress;
        committed = t.state;
    }
    // Work is signed trapezoidal stress-strain work, not inferred dissipated energy.
}
void wrapper_placement_changes_failure_scope() {
    const auto fail_all = Material::minmax(Material::parallel({Material::elastic(100), Material::elastic(20)}), -.02, .02);
    const auto fail_one = Material::parallel({Material::minmax(Material::elastic(100), -.02, .02), Material::elastic(20)});
    const auto a = evaluate(fail_all, .02, initial(fail_all));
    const auto b = evaluate(fail_one, .02, initial(fail_one));
    near(a.response.stress, 0, "outer wrapper fails whole composite");
    near(a.response.tangent, 120e-8, "outer wrapper residual tangent", 1e-16);
    near(b.response.stress, .4, "inner wrapper retains companion load path");
    near(b.response.tangent, 20 + 100e-8, "inner wrapper companion tangent");
}
void discarded_trials_leave_accepted_trace_unchanged() {
    const auto m = Material::parallel({Material::minmax(Material::steel(1000, 10, .1), -.05, .05),
                                     Material::concrete01(-4, -.002, 0, -.006)});
    State clean = initial(m), noisy = clean;
    const std::vector<double> path{0,.003,.02,-.001,-.004,.01,-.03,.04,0,.051,.01};
    double clean_work = 0, noisy_work = 0, previous_strain = 0, a_last = 0, b_last = 0;
    for (double strain : path) {
        for (double rejected : {.2,-.2,.007,-.003,0.,strain}) (void)evaluate(m, rejected, noisy);
        const auto a = evaluate(m, strain, clean);
        const auto b = evaluate(m, strain, noisy);
        equal_step(a, b, "accepted trace after rejected trials");
        clean_work += .5 * (a_last + a.response.stress) * (strain - previous_strain);
        noisy_work += .5 * (b_last + b.response.stress) * (strain - previous_strain);
        near(noisy_work, clean_work, "rejected trials must not contribute accepted work");
        clean = a.state; noisy = b.state;
        previous_strain = strain; a_last = a.response.stress; b_last = b.response.stress;
    }
}
void tangent_finite_difference_away_from_corners() {
    const auto m = Material::parallel({Material::minmax(Material::steel(1000, 10, .1), -.05, .05),
                                     Material::elastic(20)});
    const auto committed = evaluate(m, .03, initial(m)).state;
    constexpr double h = 1e-7;
    for (double strain : {.025, .04, -.02}) {
        const auto center = evaluate(m, strain, committed);
        const auto left = evaluate(m, strain - h, committed);
        const auto right = evaluate(m, strain + h, committed);
        near(center.response.tangent, (right.response.stress - left.response.stress) / (2*h),
             "algorithmic tangent finite difference", 1e-6, 1e-8);
    }
}
void copied_material_has_independent_explicit_histories() {
    const auto original = Material::parallel({Material::minmax(Material::steel(1000, 10, .1), -.05, .05),
                                             Material::elastic(20)});
    const auto copied = original; // Shared configuration is immutable; history is external.
    State a = initial(original), b = initial(copied);
    a = evaluate(original, .06, a).state;
    const auto intact = evaluate(copied, .005, b);
    near(intact.response.stress, 5.1, "copy must not inherit another instance's failure");
    equal_step(evaluate(copied, -.02, a), evaluate(original, -.02, a), "copy same-history replay");
    require(b == initial(original), "copied instance's virgin history changed");
}
void invalid_input_and_nonfinite_state() {
    invalid([] { (void)Material::parallel({}); }, "empty Parallel");
    invalid([] { (void)Material::minmax(Material::elastic(1), 0, 0); }, "equal limits");
    invalid([] { (void)Material::minmax(Material::elastic(1), 1, -1); }, "reversed limits");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double bad : {nan, inf, -inf}) {
        invalid([&] { (void)Material::minmax(Material::elastic(1), bad, 1); }, "nonfinite minimum");
        invalid([&] { (void)Material::minmax(Material::elastic(1), -1, bad); }, "nonfinite maximum");
        const auto m = Material::parallel({Material::minmax(Material::steel(1000, 10, .1), -.05, .05),
                                         Material::concrete01(-4, -.002, 0, -.006)});
        const auto s = initial(m);
        invalid([&] { (void)evaluate(m, bad, s); }, "nonfinite trial strain");
        for (std::size_t i = 0; i < s.size(); ++i) {
            auto corrupt = s;
            corrupt[i] = bad;
            invalid([&] { (void)evaluate(m, .001, corrupt); }, "nonfinite committed state entry");
        }
    }
}
} // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests{
        {"initialization_and_state_layout", initialization_and_state_layout},
        {"exact_failure_boundaries", exact_failure_boundaries},
        {"adjacent_failure_boundaries", adjacent_failure_boundaries},
        {"failure_uses_initial_not_current_tangent", failure_uses_initial_not_current_tangent},
        {"rejected_failure_after_yield", rejected_failure_after_yield},
        {"committed_failure_is_irreversible", committed_failure_is_irreversible},
        {"reinitialization_clears_nested_failure", reinitialization_clears_nested_failure},
        {"bilinear_analytic_trace", bilinear_analytic_trace},
        {"exact_yield_tangent_convention", exact_yield_tangent_convention},
        {"nested_minmax_failure_scope", nested_minmax_failure_scope},
        {"parallel_analytic_elastic_response", parallel_analytic_elastic_response},
        {"heterogeneous_parallel_history_and_work", heterogeneous_parallel_history_and_work},
        {"wrapper_placement_changes_failure_scope", wrapper_placement_changes_failure_scope},
        {"discarded_trials_leave_accepted_trace_unchanged", discarded_trials_leave_accepted_trace_unchanged},
        {"tangent_finite_difference_away_from_corners", tangent_finite_difference_away_from_corners},
        {"copied_material_has_independent_explicit_histories", copied_material_has_independent_explicit_histories},
        {"invalid_input_and_nonfinite_state", invalid_input_and_nonfinite_state}};
    std::size_t passed = 0;
    for (const auto& [name, test] : tests) {
        try { test(); ++passed; std::cout << "PASS " << name << '\n'; }
        catch (const std::exception& error) {
            std::cerr << "FAIL " << name << ": " << error.what() << '\n';
            return 1;
        }
    }
    std::cout << "Passed " << passed << '/' << tests.size() << " lifecycle/analytic groups.\n";
    return 0;
}
