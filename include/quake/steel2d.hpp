#pragma once

#include "quake/nonlinear_material.hpp"

#include <array>
#include <vector>

namespace quake {

// Concentrated-plasticity, small-displacement steel frame member.  The two
// end-hinge rotations are local element variables and are statically
// condensed, so the external topology remains the ordinary six-DOF frame
// topology [UXi, UYi, RZi, UXj, UYj, RZj].
struct SteelMember2DProperties {
    double E{};
    double A{};
    double I{};
    double axial_compression{}; // positive compression; tangent-only preload
    NonlinearMaterial hinge_i{BilinearSpring(1.0, 1.0, 0.0)};
    NonlinearMaterial hinge_j{BilinearSpring(1.0, 1.0, 0.0)};
    int local_max_iterations{50};
    double local_relative_tolerance{1e-11};
};

struct SteelMember2DResponse {
    std::array<double,6> force{};
    std::array<double,36> tangent{};
    std::vector<double> state;
    double axial_deformation{};
    double axial_force{};
    double chord_rotation{};
    std::array<double,2> end_rotation{};
    std::array<double,2> hinge_rotation{};
    std::array<double,2> end_moment{};
    std::array<MaterialEvalDiagnostics,2> hinge_diagnostics{};
    int local_iterations{};
};

class SteelMember2D {
public:
    SteelMember2D(double xi, double yi, double xj, double yj,
                  SteelMember2DProperties properties);

    int state_size() const { return state_size_; }
    std::vector<double> initial_state() const;
    std::array<double,36> initial_tangent() const;
    SteelMember2DResponse trial(const std::array<double,6>& u,
                                const double* committed) const;

private:
    SteelMember2DProperties properties_;
    double length_{}, c_{}, s_{};
    int state_size_{}, hinge_i_offset_{2}, hinge_j_offset_{};
    std::array<double,12> bending_B_{}; // 2 x 6, row-major
    std::array<double,6> axial_B_{};
    std::array<double,36> geometric_K_{};

    std::array<double,4> elastic_bending_stiffness() const;
    std::array<double,4> condensed_bending_stiffness(double ki,
                                                      double kj) const;
};

// Memoryless axial power-law dashpot.  For regularization_velocity > 0,
// F = C*v*(v^2+v_r^2)^((alpha-1)/2), which is smooth and supplies a finite,
// consistent tangent at zero velocity.  alpha=1 is the exact linear dashpot.
struct ViscousDamper2DProperties {
    double coefficient{};
    double alpha{1.0};
    double regularization_velocity{};
};

struct ViscousDamper2DTrial {
    double deformation_rate{};
    double force{};
    double tangent{};
};

class ViscousDamper2D {
public:
    explicit ViscousDamper2D(ViscousDamper2DProperties properties);
    ViscousDamper2DTrial trial(double deformation_rate) const;
    double initial_tangent() const { return trial(0.0).tangent; }
    const ViscousDamper2DProperties& properties() const { return properties_; }
private:
    ViscousDamper2DProperties properties_;
};

} // namespace quake
