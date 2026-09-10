#pragma once
#include "quake/sparse.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/update_basis.hpp"
#include "quake/low_rank_solver.hpp"
#include <vector>

namespace quake {

// Exact solver for A + U C U^T where U is a fixed sparse basis and C is a
// state-dependent dense symmetric coefficient matrix. Unlike the diagonal-D
// hinge Woodbury solver, this supports coupled low-rank blocks such as
// localized geometric-stiffness updates.
class GeneralizedWoodburySolver {
public:
    GeneralizedWoodburySolver(const SparseMatrixCSC& baseline,
                              SparseUpdateBasis basis);

    std::vector<double> solve(const std::vector<double>& rhs,
                              const std::vector<double>& C_row_major,
                              double active_tol=1e-14);

    int size() const { return n_; }
    int update_dimension() const { return r_; }
    double last_reduced_pivot_ratio() const { return last_pivot_ratio_; }
    std::size_t last_active_dimension() const { return last_active_dimension_; }
    const SparseUpdateBasis& basis() const { return basis_; }

private:
    int n_{};
    int r_{};
    SuperLUFactor factor_;
    SparseUpdateBasis basis_;
    std::vector<double> W_; // A^-1 U, n x r column-major
    std::vector<double> H_; // U^T A^-1 U, r x r row-major
    double last_pivot_ratio_{1.0};
    std::size_t last_active_dimension_{0};
};

} // namespace quake
