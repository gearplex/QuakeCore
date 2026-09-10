#pragma once
#include "quake/dynamic_model.hpp"
#include <vector>

namespace quake {

struct ModeShape {
    double eigenvalue{};      // omega^2
    double omega{};           // rad / time
    double frequency_hz{};    // cycles / time
    double period{};          // time / cycle
    std::vector<double> shape;
};

// Dense validation-oriented generalized eigen solution of
//     K phi = lambda M phi
// using LAPACK DGGEV. Unlike Cholesky-based symmetric drivers, this tolerates
// the singular generalized mass matrices common in frame models with massless
// rotational/internal DOFs. Finite positive real modes are returned in
// ascending eigenvalue order. This is intentionally a reference utility, not
// the production sparse eigensolver.
std::vector<ModeShape> modal_analysis(const NonlinearDynamicModel& model,
                                      int requested_modes,
                                      double beta_tolerance = 1e-12,
                                      double imaginary_tolerance = 1e-8);

// Validation/diagnostic companion for an instantaneous structural tangent. The
// mass operator is still taken from model, but K may be any assembled tangent
// with the same reduced DOF ordering (for example a committed-state tangent at
// a drift peak).
std::vector<ModeShape> modal_analysis_from_stiffness(const NonlinearDynamicModel& model,
                                                     const SparseMatrixCSC& stiffness,
                                                     int requested_modes,
                                                     double beta_tolerance = 1e-12,
                                                     double imaginary_tolerance = 1e-8);

} // namespace quake
