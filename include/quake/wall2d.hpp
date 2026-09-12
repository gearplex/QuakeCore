#pragma once
#include "quake/wall_material.hpp"
#include <array>
#include <variant>
#include <vector>

namespace quake {
struct MVLEMFiber {
    double width{}, thickness{}, reinforcement_ratio{};
    WallUniaxial concrete{WallUniaxial::elastic(1)};
    WallUniaxial steel{WallUniaxial::elastic(1)};
};
struct SFIMVLEMPanel { double width{}, thickness{}; WallPanel material{1,0}; };
struct MVLEMProperties {
    double c{.4}; double density{}; // mass / volume; half to each node, UX and UY
    std::vector<MVLEMFiber> fibers;
    // Force versus shear DISPLACEMENT, not stress versus strain.
    WallUniaxial shear{WallUniaxial::elastic(1)};
};
struct SFIMVLEMProperties {
    double c{.4}; double density{};
    std::vector<SFIMVLEMPanel> panels;
    int local_max_iterations{60};
    double local_relative_tolerance{1e-10};
};
struct Wall2DResponse {
    std::array<double,6> force{};
    std::array<double,36> tangent{};
    std::vector<double> state;
    double curvature{}, shear_deformation{};
    std::vector<double> fiber_strain, concrete_stress, steel_stress;
    std::vector<std::array<double,3>> panel_strain, panel_stress;
    double maximum_transverse_stress_residual{};
    int local_iterations{};
};
// Vertical two-node, six external DOFs [UXi,UYi,RZi,UXj,UYj,RZj].
// Small displacements; no geometric stiffness or implicit gravity preload.
// SFI internal horizontal panel extensions are solved at sigma_x=0, then
// statically condensed with the FULL panel constitutive Jacobian.
class Wall2D {
public:
    Wall2D(double height,MVLEMProperties properties);
    Wall2D(double height,SFIMVLEMProperties properties);
    int state_size() const { return state_size_; }
    double total_mass() const { return total_mass_; }
    bool is_sfi() const { return std::holds_alternative<SFIMVLEMProperties>(properties_); }
    std::vector<double> initial_state() const;
    std::array<double,36> initial_tangent() const;
    Wall2DResponse trial(const std::array<double,6>& u,const double* committed) const;
private:
    double h_{},c_{},total_mass_{};
    int state_size_{};
    int shear_offset_{-1};
    std::variant<MVLEMProperties,SFIMVLEMProperties> properties_;
    std::vector<double> x_;
    std::vector<int> offsets_;
    std::array<double,6> shear_B() const;
    std::array<double,6> axial_B(int fiber) const;
};
} // namespace quake
