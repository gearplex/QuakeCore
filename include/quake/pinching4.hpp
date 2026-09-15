#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>

namespace quake {

struct Pinching4Point {
    double deformation{};
    double force{};
};

struct Pinching4EnvelopeParameters {
    std::array<Pinching4Point, 4> positive{};
    std::array<Pinching4Point, 4> negative{};
};

struct Pinching4Response {
    double force{};
    double tangent{};
};

// OpenSees-compatible Pinching4 virgin/damaged-envelope geometry only.
// Cyclic transition states, pinching paths, and damage evolution are admitted
// separately under Gate 4 after this envelope kernel is verified.
class Pinching4Envelope {
public:
    explicit Pinching4Envelope(Pinching4EnvelopeParameters parameters)
        : parameters_(parameters) {
        validate_positive(parameters_.positive);
        validate_negative(parameters_.negative);

        const double k_pos = parameters_.positive[0].force / parameters_.positive[0].deformation;
        const double k_neg = parameters_.negative[0].force / parameters_.negative[0].deformation;
        const double k = std::max(k_pos, k_neg);
        const double u = std::max(parameters_.positive[0].deformation,
                                  -parameters_.negative[0].deformation) * 1.0e-4;

        positive_envelope_[0] = {u, u * k};
        negative_envelope_[0] = {-u, -u * k};
        for (std::size_t i = 0; i < 4; ++i) {
            positive_envelope_[i + 1] = parameters_.positive[i];
            negative_envelope_[i + 1] = parameters_.negative[i];
        }

        const double k_last_pos = slope(parameters_.positive[2], parameters_.positive[3]);
        const double k_last_neg = slope(parameters_.negative[2], parameters_.negative[3]);
        positive_envelope_[5].deformation = 1.0e6 * parameters_.positive[3].deformation;
        positive_envelope_[5].force = k_last_pos > 0.0
            ? parameters_.positive[3].force +
                  k_last_pos * (positive_envelope_[5].deformation - parameters_.positive[3].deformation)
            : 1.1 * parameters_.positive[3].force;
        negative_envelope_[5].deformation = 1.0e6 * parameters_.negative[3].deformation;
        negative_envelope_[5].force = k_last_neg > 0.0
            ? parameters_.negative[3].force +
                  k_last_neg * (negative_envelope_[5].deformation - parameters_.negative[3].deformation)
            : 1.1 * parameters_.negative[3].force;

        positive_elastic_tangent_ =
            positive_envelope_[1].force / positive_envelope_[1].deformation;
        negative_elastic_tangent_ =
            negative_envelope_[1].force / negative_envelope_[1].deformation;
    }

    const Pinching4EnvelopeParameters& parameters() const noexcept { return parameters_; }
    Pinching4Response response(double deformation) const {
        return deformation >= 0.0 ? positive(deformation) : negative(deformation);
    }
    Pinching4Response positive(double deformation) const {
        if (deformation < 0.0) {
            throw std::invalid_argument("Pinching4 positive envelope requires deformation >= 0");
        }
        if (deformation == 0.0) return {0.0, positive_elastic_tangent_};
        if (deformation <= positive_envelope_[0].deformation) {
            return {positive_elastic_tangent_ * deformation, positive_elastic_tangent_};
        }
        for (std::size_t i = 0; i < 5; ++i) {
            if (deformation <= positive_envelope_[i + 1].deformation) {
                return interpolate(deformation, positive_envelope_[i], positive_envelope_[i + 1]);
            }
        }
        return interpolate(deformation, positive_envelope_[4], positive_envelope_[5]);
    }
    Pinching4Response negative(double deformation) const {
        if (deformation > 0.0) {
            throw std::invalid_argument("Pinching4 negative envelope requires deformation <= 0");
        }
        if (deformation == 0.0) return {0.0, negative_elastic_tangent_};
        if (deformation >= negative_envelope_[0].deformation) {
            return {negative_elastic_tangent_ * deformation, negative_elastic_tangent_};
        }
        for (std::size_t i = 0; i < 5; ++i) {
            if (deformation >= negative_envelope_[i + 1].deformation) {
                return interpolate(deformation, negative_envelope_[i], negative_envelope_[i + 1]);
            }
        }
        return interpolate(deformation, negative_envelope_[4], negative_envelope_[5]);
    }
    double positive_elastic_tangent() const noexcept { return positive_elastic_tangent_; }
    double negative_elastic_tangent() const noexcept { return negative_elastic_tangent_; }

private:
    static double slope(const Pinching4Point& a, const Pinching4Point& b) {
        return (b.force - a.force) / (b.deformation - a.deformation);
    }
    static Pinching4Response interpolate(double deformation,
                                         const Pinching4Point& a,
                                         const Pinching4Point& b) {
        const double tangent = slope(a, b);
        return {a.force + (deformation - a.deformation) * tangent, tangent};
    }
    static void validate_positive(const std::array<Pinching4Point, 4>& points) {
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (!(points[i].deformation > 0.0) || !(points[i].force > 0.0)) {
                throw std::invalid_argument("Pinching4 positive envelope points must be positive");
            }
            if (i > 0 && !(points[i].deformation > points[i - 1].deformation)) {
                throw std::invalid_argument("Pinching4 positive envelope deformation must increase");
            }
        }
    }
    static void validate_negative(const std::array<Pinching4Point, 4>& points) {
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (!(points[i].deformation < 0.0) || !(points[i].force < 0.0)) {
                throw std::invalid_argument("Pinching4 negative envelope points must be negative");
            }
            if (i > 0 && !(points[i].deformation < points[i - 1].deformation)) {
                throw std::invalid_argument("Pinching4 negative envelope deformation must decrease");
            }
        }
    }

    std::array<Pinching4Point, 6> positive_envelope_{};
    std::array<Pinching4Point, 6> negative_envelope_{};
    Pinching4EnvelopeParameters parameters_{};
    double positive_elastic_tangent_{};
    double negative_elastic_tangent_{};
};

} // namespace quake
