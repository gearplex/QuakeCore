#pragma once
#include "quake/sparse.hpp"
#include <memory>
#include <vector>

namespace quake {

class SuperLUFactor {
public:
    explicit SuperLUFactor(const SparseMatrixCSC& matrix);
    ~SuperLUFactor();
    SuperLUFactor(const SuperLUFactor&) = delete;
    SuperLUFactor& operator=(const SuperLUFactor&) = delete;
    SuperLUFactor(SuperLUFactor&&) noexcept;
    SuperLUFactor& operator=(SuperLUFactor&&) noexcept;

    std::vector<double> solve(const std::vector<double>& rhs) const;
    // Column-major RHS/result: n rows x nrhs columns.
    std::vector<double> solve_multiple(const std::vector<double>& rhs, int nrhs) const;
    int size() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};


// Repeated sparse direct solver for matrices with identical CSC sparsity.
// The first call performs ordering/symbolic analysis; subsequent numerical
// factorizations reuse SuperLU's column ordering and elimination tree via
// SamePattern while allowing new row pivots as the tangent degrades. Intended as a fair optimized full-Newton baseline.
class SuperLUSamePatternSolver {
public:
    struct Impl;
    explicit SuperLUSamePatternSolver(const SparseMatrixCSC& initial_matrix);
    ~SuperLUSamePatternSolver();
    SuperLUSamePatternSolver(const SuperLUSamePatternSolver&) = delete;
    SuperLUSamePatternSolver& operator=(const SuperLUSamePatternSolver&) = delete;
    SuperLUSamePatternSolver(SuperLUSamePatternSolver&&) noexcept;
    SuperLUSamePatternSolver& operator=(SuperLUSamePatternSolver&&) noexcept;

    std::vector<double> solve_current(const std::vector<double>& rhs) const;
    std::vector<double> refactor_and_solve(const SparseMatrixCSC& matrix,
                                           const std::vector<double>& rhs);
    int size() const;
    std::size_t factorizations() const;
    std::size_t symbolic_analyses() const { return 1; }

private:
    std::unique_ptr<Impl> impl_;
};

std::vector<double> superlu_solve_once(const SparseMatrixCSC& matrix,
                                       const std::vector<double>& rhs);

} // namespace quake
