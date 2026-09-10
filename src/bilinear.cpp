#include "quake/bilinear.hpp"
#include <cmath>
#include <stdexcept>

namespace quake {

BilinearSpring::BilinearSpring(double k0, double yield_force, double post_yield_ratio)
    : k0_(k0), fy_(yield_force) {
    if (k0 <= 0.0 || yield_force <= 0.0 || post_yield_ratio < 0.0 || post_yield_ratio >= 1.0)
        throw std::invalid_argument("Invalid bilinear spring parameters");
    H_ = post_yield_ratio == 0.0 ? 0.0 : post_yield_ratio * k0_ / (1.0 - post_yield_ratio);
}

BilinearTrial BilinearSpring::trial(double deformation, const BilinearState& c) const {
    const double sigma_trial = k0_ * (deformation - c.plastic);
    const double xi = sigma_trial - c.backstress;
    const double f = std::abs(xi) - fy_;
    // Steel01 selects the post-yield branch at the exact yield surface.  The
    // force is continuous either way, but matching the tangent convention
    // avoids a one-iteration discrepancy when a protocol lands exactly on Fy.
    if (f < 0.0) return {sigma_trial, k0_, c};
    const double sign = xi >= 0.0 ? 1.0 : -1.0;
    const double dgamma = f / (k0_ + H_);
    BilinearState s = c;
    s.plastic += dgamma * sign;
    s.backstress += H_ * dgamma * sign;
    const double force = sigma_trial - k0_ * dgamma * sign;
    const double tangent = k0_ * H_ / (k0_ + H_);
    return {force, tangent, s};
}

} // namespace quake
