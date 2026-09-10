#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace quake {
inline bool all_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double x){ return std::isfinite(x); });
}
inline void validate_transient_input(const std::vector<double>& acceleration,
                                     double dt, double tolerance, int iterations) {
    if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(1.0/(dt*dt)) ||
        acceleration.empty() || !all_finite(acceleration) ||
        !std::isfinite(tolerance) || tolerance <= 0.0 || iterations < 1)
        throw std::invalid_argument("transient input requires finite acceleration, positive finite dt/tolerance, and iterations >= 1");
}
}
