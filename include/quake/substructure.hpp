#pragma once

#include "quake/frame3d.hpp"
#include "quake/low_rank_solver.hpp"
#include "quake/sparse.hpp"
#include "quake/superlu_solver.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace quake {

struct BlockPartition {
    std::vector<int> interface_dofs;
    std::vector<std::vector<int>> interior_blocks;
};

// Build a floor-banded exact-elimination partition. Every Nth active elevation
// is retained so interior DOFs can be eliminated in local story blocks without
// a global dense transformation. Nonlinear-support DOFs need not remain on the
// interface because Woodbury only requires exact baseline inverse actions A^-1 b_j;
// retain_nonlinear_support=true is available as a diagnostic/conservative option.
BlockPartition make_story_block_partition(const CompiledFrame3D& model,
                                          int stories_per_block,
                                          bool retain_nonlinear_support=false);

struct BlockSchurStats {
    int full_dof{};
    int interface_dof{};
    int interior_dof{};
    int blocks{};
    int schur_nnz{};
    double setup_seconds{};
};

// Exact factorization A x=b using independent interior blocks and a global
// sparse Schur complement on retained interface DOFs. No model reduction or
// dynamic approximation is introduced.
class BlockSchurFactor {
public:
    BlockSchurFactor(const SparseMatrixCSC& matrix, BlockPartition partition);
    ~BlockSchurFactor();
    BlockSchurFactor(BlockSchurFactor&&) noexcept;
    BlockSchurFactor& operator=(BlockSchurFactor&&) noexcept;
    BlockSchurFactor(const BlockSchurFactor&) = delete;
    BlockSchurFactor& operator=(const BlockSchurFactor&) = delete;

    int size() const;
    std::vector<double> solve(const std::vector<double>& rhs) const;
    std::vector<double> solve_multiple(const std::vector<double>& rhs, int nrhs) const;
    const BlockSchurStats& stats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Exact lazy Woodbury solver whose baseline A^-1 operations use the local
// story-block Schur factor instead of one monolithic sparse triangular solve.
class SubstructuredLazyWoodburySolver final : public LowRankSolverBase {
public:
    SubstructuredLazyWoodburySolver(const SparseMatrixCSC& baseline,
                                    SparseUpdateBasis basis,
                                    BlockPartition partition);
    std::vector<double> solve(const std::vector<double>& rhs,
                              const std::vector<double>& delta_k,
                              double active_tol=1e-14) override;
    int size() const override { return n_; }
    int update_count() const override { return m_; }
    std::size_t cached_update_count() const override { return cached_indices_.size(); }
    const BlockSchurStats& block_stats() const { return factor_.stats(); }
private:
    void cache_update(int update_index);
    int n_{},m_{};
    BlockSchurFactor factor_;
    SparseUpdateBasis basis_;
    std::vector<int> index_to_cache_;
    std::vector<int> cached_indices_;
    std::vector<double> W_cache_;
    std::vector<double> G_cache_;
};

} // namespace quake
