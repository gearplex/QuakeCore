#include "quake/concrete_cm.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quake::ConcreteCMEnvelope;
using quake::ConcreteCMParameters;

namespace {

bool near(double a, double b, double rel = 1.0e-11, double abs = 1.0e-13) {
    return std::abs(a - b) <= abs + rel * std::max(std::abs(a), std::abs(b));
}

void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

std::vector<double> compression_protocol_prefix() {
    std::vector<double> strains{0.0};
    const double targets[] = {-0.0005, -0.0020, -0.0040};
    for (double target : targets) {
        const double start = strains.back();
        for (int i = 1; i <= 20; ++i) {
            strains.push_back(start + (target - start) * static_cast<double>(i) / 20.0);
        }
    }
    return strains;
}

void verify_oracle_prefix(const ConcreteCMEnvelope& material,
                          const std::vector<double>& expected_stress,
                          const std::vector<double>& expected_tangent,
                          double expected_work) {
    const auto strains = compression_protocol_prefix();
    const int steps[] = {0, 1, 20, 40, 41, 60};
    double work = 0.0;
    auto previous = material.compression(strains.front());
    for (std::size_t i = 1; i < strains.size(); ++i) {
        const auto current = material.compression(strains[i]);
        work += 0.5 * (previous.stress + current.stress) * (strains[i] - strains[i - 1]);
        previous = current;
    }
    for (std::size_t i = 0; i < std::size(steps); ++i) {
        const auto response = material.compression(strains[static_cast<std::size_t>(steps[i])]);
        check(near(response.stress, expected_stress[i]), "ConcreteCM compression stress disagrees with frozen OpenSees oracle");
        check(near(response.tangent, expected_tangent[i]), "ConcreteCM compression tangent disagrees with frozen OpenSees oracle");
    }
    check(near(work, expected_work, 1.0e-10, 1.0e-13), "ConcreteCM compression work disagrees with frozen OpenSees oracle");
}

ConcreteCMParameters compression_parameters(double fc, double epsc, double Ec, double xcrn) {
    // Tension fields do not participate in this Gate 4 compression-envelope
    // regression. They are populated with valid neutral values until the
    // tension/cyclic oracle is admitted in the next component step.
    return {fc, epsc, Ec, 7.0, xcrn, 1.0, 1.0 / Ec, 1.2, 10000.0, true};
}

} // namespace

int main() try {
    const ConcreteCMEnvelope unconfined(compression_parameters(
        -6.5, -0.002, 4595.486916530173, 1.030));
    verify_oracle_prefix(
        unconfined,
        {0.0, -0.11453308056925031, -2.163921164843064, -6.5,
         -6.465420005089847, -5.52346369744675},
        {4595.486916530173, 4567.2031827561395, 4075.5367970685274, 0.0,
         -495.76647770689294, -495.76647770689294},
        0.019733877567269662);

    const ConcreteCMEnvelope confined(compression_parameters(
        -8.01435, -0.00432976, 5102.805419570689, 1.015));
    verify_oracle_prefix(
        confined,
        {0.0, -0.12640952032925176, -2.155578462184863, -5.881318217785812,
         -6.046356343523506, -7.958377184018189},
        {5102.805419570689, 5010.378571312529, 3642.3236647637405,
         1687.047305323382, 1614.318243147867, 330.2002092395013},
        0.021221688703467297);

    std::cout << "ConcreteCM monotonic compression envelope matches frozen OpenSees 3.8.0 oracle.\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}
