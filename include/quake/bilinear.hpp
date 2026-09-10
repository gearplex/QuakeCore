#pragma once

namespace quake {

struct BilinearState {
    double plastic{0.0};
    double backstress{0.0};
};

struct BilinearTrial {
    double force{0.0};
    double tangent{0.0};
    BilinearState state{};
};

class BilinearSpring {
public:
    BilinearSpring(double k0, double yield_force, double post_yield_ratio);
    BilinearTrial trial(double deformation, const BilinearState& committed) const;
    double initial_stiffness() const { return k0_; }
    double yield_force() const { return fy_; }
    double post_yield_modulus() const { return H_; }

private:
    double k0_{};
    double fy_{};
    double H_{};
};

} // namespace quake
