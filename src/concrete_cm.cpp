#include "quake/concrete_cm.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {

struct Shape {
    double y{};
    double z{};
};

Shape chang_mander_shape(double x, double n, double r) {
    double D = 0.0;
    if (r != 1.0) {
        D = 1.0 + (n - r / (r - 1.0)) * x + std::pow(x, r) / (r - 1.0);
    } else {
        if (x <= 0.0) return {0.0, 1.0};
        D = 1.0 + (n - 1.0 + std::log10(x)) * x;
    }
    return {n * x / D, (1.0 - std::pow(x, r)) / (D * D)};
}

ConcreteCMResponse envelope(double strain,
                            double zero_strain,
                            double peak_strain,
                            double peak_stress,
                            double Ec,
                            double r,
                            double xcr) {
    const double scale = std::abs(peak_strain);
    const double x = std::abs((strain - zero_strain) / scale);
    const double n = std::abs(Ec * peak_strain / peak_stress);
    const auto critical = chang_mander_shape(xcr, n, r);
    const double cutoff = std::abs(xcr - critical.y / (n * critical.z));

    if (x > cutoff) return {0.0, 0.0};
    if (x < xcr) {
        const auto point = chang_mander_shape(x, n, r);
        return {peak_stress * point.y, Ec * point.z};
    }

    return {
        peak_stress * (critical.y + n * critical.z * (x - xcr)),
        Ec * critical.z,
    };
}

} // namespace

ConcreteCMEnvelope::ConcreteCMEnvelope(ConcreteCMParameters parameters)
    : parameters_(parameters) {
    if (!(parameters_.fc < 0.0) || !(parameters_.epsc < 0.0) ||
        !(parameters_.Ec > 0.0) || !(parameters_.rc > 0.0) ||
        !(parameters_.xcrn > 0.0) || !(parameters_.ft > 0.0) ||
        !(parameters_.et > 0.0) || !(parameters_.rt > 0.0) ||
        !(parameters_.xcrp > 0.0)) {
        throw std::invalid_argument("ConcreteCM parameters have invalid sign or zero value");
    }
}

ConcreteCMResponse ConcreteCMEnvelope::compression(double strain) const {
    if (strain > 0.0) throw std::invalid_argument("ConcreteCM compression envelope requires strain <= 0");
    if (strain == 0.0) return {0.0, parameters_.Ec};
    return envelope(strain, 0.0, parameters_.epsc, parameters_.fc,
                    parameters_.Ec, parameters_.rc, parameters_.xcrn);
}

ConcreteCMResponse ConcreteCMEnvelope::tension(double strain, double zero_strain) const {
    if (strain < zero_strain) {
        throw std::invalid_argument("ConcreteCM tension envelope requires strain >= zero strain");
    }
    if (strain == zero_strain) return {0.0, parameters_.Ec};
    return envelope(strain, zero_strain, parameters_.et, parameters_.ft,
                    parameters_.Ec, parameters_.rt, parameters_.xcrp);
}

} // namespace quake
