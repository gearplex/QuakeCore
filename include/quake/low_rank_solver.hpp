#pragma once
#include "quake/sparse.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/update_basis.hpp"
#include <cstddef>
#include <vector>

namespace quake {

struct DenseSolveDiagnostics { double min_abs_pivot{0.0}; double max_abs_pivot{0.0}; double pivot_ratio{1.0}; };

class LowRankSolverBase {
public:
    virtual ~LowRankSolverBase() = default;
    virtual std::vector<double> solve(const std::vector<double>& rhs,
                                      const std::vector<double>& delta_k,
                                      double active_tol = 1e-14) = 0;
    virtual int size() const = 0;
    virtual int update_count() const = 0;
    virtual std::size_t cached_update_count() const = 0;
    virtual double last_reduced_pivot_ratio() const { return 1.0; }
};

// Eager exact Woodbury solver. Precomputes A^-1 B and the full Gram matrix
// G=B^T A^-1 B for all potential update directions.
class LowRankWoodburySolver final : public LowRankSolverBase {
public:
    LowRankWoodburySolver(const SparseMatrixCSC& baseline,
                          SparseUpdateBasis basis);
    // Backward-compatible research/benchmark constructor.
    LowRankWoodburySolver(const SparseMatrixCSC& baseline,
                          std::vector<double> B_column_major,
                          int update_count);

    std::vector<double> solve(const std::vector<double>& rhs,
                              const std::vector<double>& delta_k,
                              double active_tol = 1e-14) override;

    int size() const override { return n_; }
    int update_count() const override { return m_; }
    std::size_t cached_update_count() const override { return static_cast<std::size_t>(m_); }
    std::size_t baseline_factorizations() const { return 1; }
    const std::vector<double>& influence_matrix() const { return W_; }
    double last_reduced_pivot_ratio() const override { return last_pivot_ratio_; }

private:
    int n_{};
    int m_{};
    SuperLUFactor factor_;
    SparseUpdateBasis basis_;
    std::vector<double> W_; // A^-1 B, n x m column-major
    std::vector<double> G_; // B^T A^-1 B, m x m row-major
    double last_pivot_ratio_{1.0};
};

// Lazy exact Woodbury solver for building NRHA. The baseline A is factorized
// once, but A^-1 b_j columns are generated and cached only when update j first
// becomes active. Candidate B columns themselves remain sparse.
class LazyLowRankWoodburySolver final : public LowRankSolverBase {
public:
    LazyLowRankWoodburySolver(const SparseMatrixCSC& baseline,
                              SparseUpdateBasis basis);
    // Backward-compatible research/benchmark constructor.
    LazyLowRankWoodburySolver(const SparseMatrixCSC& baseline,
                              std::vector<double> B_column_major,
                              int update_count);

    std::vector<double> solve(const std::vector<double>& rhs,
                              const std::vector<double>& delta_k,
                              double active_tol = 1e-14) override;

    int size() const override { return n_; }
    int update_count() const override { return m_; }
    std::size_t cached_update_count() const override { return cached_indices_.size(); }
    std::size_t baseline_factorizations() const { return 1; }
    std::size_t basis_nnz() const { return static_cast<std::size_t>(basis_.nnz()); }
    double last_reduced_pivot_ratio() const override { return last_pivot_ratio_; }

private:
    void cache_update(int update_index);

    int n_{};
    int m_{};
    SuperLUFactor factor_;
    SparseUpdateBasis basis_;
    std::vector<int> index_to_cache_; // candidate index -> cache column, -1 if unseen
    std::vector<int> cached_indices_; // cache column -> candidate index
    std::vector<double> W_cache_; // n x k column-major
    std::vector<double> G_cache_; // k x k row-major
    double last_pivot_ratio_{1.0};
};

std::vector<double> dense_solve(std::vector<double> A_row_major,
                                std::vector<double> b, int n, DenseSolveDiagnostics* diagnostics=nullptr);

} // namespace quake
