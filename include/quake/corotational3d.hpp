#pragma once
#include <array>

namespace quake {

// Research/reference 3D corotational Euler-Bernoulli frame element.
// The 12 generalized element coordinates are global translations/rotation
// vectors [ux,uy,uz,rx,ry,rz] at end i followed by end j.  The formulation
// is objective under finite rigid-body translation/rotation and uses a
// co-rotated six-deformation basic system.  It is intentionally implemented
// from a scalar strain-energy potential with numerical differentiation so it
// can serve as a high-confidence correctness oracle before an optimized
// analytic tangent is introduced.
struct CorotationalFrame3DProperties {
    double xi{}, yi{}, zi{};
    double xj{}, yj{}, zj{};
    double E{}, G{}, A{}, J{}, Iy{}, Iz{};
    std::array<double,3> reference{1.0,0.0,0.0};
    // Positive compression. The initial prestress contribution is subtracted
    // from the residual at the undeformed configuration while its geometric
    // stiffness remains in the tangent.
    double axial_compression{0.0};
};

// Objective local section-force result recovered from the co-rotated basic
// system. The values are invariant to superposed rigid-body translation and
// rotation. Compression is positive; end shear signs follow the current local
// member axes.
struct CorotationalFrame3DSectionResponse {
    double axial_compression{};
    double shear_y_i{}, shear_y_j{};
    double shear_z_i{}, shear_z_j{};
    double torsion_i{}, torsion_j{};
    double moment_y_i{}, moment_y_j{};
    double moment_z_i{}, moment_z_j{};
    double current_length{};
};

struct CorotationalFrame3DResponse {
    std::array<double,12> force{};
    std::array<double,144> tangent{}; // row-major
    std::array<double,6> basic_deformation{};
    double strain_energy{0.0};
    double current_length{0.0};
};

std::array<double,6> corotational3d_basic_deformation(
    const CorotationalFrame3DProperties& p,
    const std::array<double,12>& dof);

double corotational3d_strain_energy(
    const CorotationalFrame3DProperties& p,
    const std::array<double,12>& dof);

CorotationalFrame3DSectionResponse corotational3d_section_response(
    const CorotationalFrame3DProperties& p,
    const std::array<double,12>& dof);

std::array<double,12> corotational3d_internal_force(
    const CorotationalFrame3DProperties& p,
    const std::array<double,12>& dof,
    double relative_step = 2e-7);

CorotationalFrame3DResponse corotational3d_response(
    const CorotationalFrame3DProperties& p,
    const std::array<double,12>& dof,
    double relative_step = 2e-6);

} // namespace quake
