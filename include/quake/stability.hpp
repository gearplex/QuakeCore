#pragma once
#include "quake/dynamic_model.hpp"
#include "quake/sparse.hpp"
#include <cstddef>
#include <limits>
#include <vector>

namespace quake {

struct TangentStabilityEstimate {
    bool factorization_ok{false};
    // Signed Rayleigh estimate of the tangent eigenvalue nearest zero.
    double near_zero_eigenvalue{0.0};
    double absolute_eigenvalue{0.0};
    // Backward-compatible alias.
    double rayleigh_eigenvalue{0.0};
};

enum class InitialStabilityStatus {
    PositiveDefinite,
    NotPositiveDefinite,
    NumericalFailure
};

struct InitialStabilityAssessment {
    InitialStabilityStatus status{InitialStabilityStatus::NumericalFailure};
    bool positive_definite{false};
    int dof{0};
    int failing_pivot{-1};
    double minimum_pivot{0.0};
    double minimum_pivot_ratio{0.0};
    double symmetry_relative_error{0.0};
    std::size_t factor_nonzeros{0};
    // Dense cross-check is populated only when requested and the matrix is
    // below dense_check_limit. It is validation-oriented, not the scalable path.
    bool dense_check_performed{false};
    int dense_negative_eigenvalues{0};
    int dense_near_zero_eigenvalues{0};
    double dense_minimum_eigenvalue{std::numeric_limits<double>::quiet_NaN()};
};

// Scalable initial-stability certification. A symmetric matrix is positive
// definite iff symmetric elimination/Cholesky can proceed with strictly
// positive pivots. The implementation applies a reverse-Cuthill-McKee
// symmetric permutation and performs sparse Schur updates; it never forms an
// n-by-n dense matrix. For small validation models an optional LAPACK dense
// eigen cross-check can also be requested.
InitialStabilityAssessment assess_positive_definiteness(const SparseMatrixCSC& matrix,
                                                        double relative_pivot_tolerance=1e-12,
                                                        int dense_check_limit=0);
InitialStabilityAssessment assess_initial_stability_from_tangents(const NonlinearDynamicModel& model,
                                                                  const std::vector<double>& u,
                                                                  const std::vector<double>& tangents,
                                                                  double relative_pivot_tolerance=1e-12,
                                                                  int dense_check_limit=0);

// Inverse-iteration estimate of the tangent eigenvalue nearest zero. Intended
// to track loss of stiffness from a certified stable initial configuration,
// not to compute the complete inertia of an arbitrarily indefinite matrix.
TangentStabilityEstimate estimate_tangent_stability(const NonlinearDynamicModel& model,
                                                     const std::vector<double>& u,
                                                     const std::vector<double>& committed_state,
                                                     int inverse_iterations=8);
TangentStabilityEstimate estimate_tangent_stability_from_tangents(const NonlinearDynamicModel& model,
                                                                  const std::vector<double>& u,
                                                                  const std::vector<double>& tangents,
                                                                  int inverse_iterations=8);

} // namespace quake
