#include "quake/superlu_solver.hpp"

#include <stdexcept>
#include <utility>

extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda,
             int* ipiv, int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs,
             const double* a, const int* lda, const int* ipiv,
             double* b, const int* ldb, int* info);
}

namespace quake {
namespace {

std::vector<double> dense_column_major(const SparseMatrixCSC& matrix) {
    if (matrix.rows() != matrix.cols())
        throw std::invalid_argument("Factor requires square matrix");
    const int n = matrix.rows();
    std::vector<double> a(static_cast<std::size_t>(n * n), 0.0);
    for (int col = 0; col < n; ++col) {
        for (int k = matrix.col_ptr()[static_cast<std::size_t>(col)];
             k < matrix.col_ptr()[static_cast<std::size_t>(col + 1)]; ++k) {
            const int row = matrix.row_ind()[static_cast<std::size_t>(k)];
            a[static_cast<std::size_t>(row + col * n)] += matrix.values()[static_cast<std::size_t>(k)];
        }
    }
    return a;
}

void factor(std::vector<double>& a, std::vector<int>& piv, int n) {
    int info = 0;
    dgetrf_(&n, &n, a.data(), &n, piv.data(), &info);
    if (info != 0)
        throw std::runtime_error("LAPACK dense factorization failed, info=" + std::to_string(info));
}

std::vector<double> solve_factored(const std::vector<double>& a,
                                   const std::vector<int>& piv, int n,
                                   const std::vector<double>& rhs, int nrhs) {
    if (nrhs <= 0 || static_cast<int>(rhs.size()) != n * nrhs)
        throw std::invalid_argument("Dense fallback RHS dimension mismatch");
    std::vector<double> x = rhs;
    const char trans = 'N';
    int info = 0;
    dgetrs_(&trans, &n, &nrhs, a.data(), &n, piv.data(), x.data(), &n, &info);
    if (info != 0)
        throw std::runtime_error("LAPACK dense solve failed, info=" + std::to_string(info));
    return x;
}

void validate_pattern(const SparseMatrixCSC& matrix, int n,
                      const std::vector<int>& cols, const std::vector<int>& rows) {
    if (matrix.rows() != n || matrix.cols() != n || matrix.col_ptr() != cols ||
        matrix.row_ind() != rows)
        throw std::invalid_argument("Same-pattern solver received different sparsity pattern");
}
}  // namespace

struct SuperLUFactor::Impl {
    int n{};
    std::vector<double> lu;
    std::vector<int> piv;
};

SuperLUFactor::SuperLUFactor(const SparseMatrixCSC& matrix)
    : impl_(std::make_unique<Impl>()) {
    impl_->n = matrix.rows();
    impl_->lu = dense_column_major(matrix);
    impl_->piv.resize(static_cast<std::size_t>(impl_->n));
    factor(impl_->lu, impl_->piv, impl_->n);
}
SuperLUFactor::~SuperLUFactor() = default;
SuperLUFactor::SuperLUFactor(SuperLUFactor&&) noexcept = default;
SuperLUFactor& SuperLUFactor::operator=(SuperLUFactor&&) noexcept = default;
int SuperLUFactor::size() const { return impl_->n; }
std::vector<double> SuperLUFactor::solve(const std::vector<double>& rhs) const {
    return solve_multiple(rhs, 1);
}
std::vector<double> SuperLUFactor::solve_multiple(const std::vector<double>& rhs, int nrhs) const {
    return solve_factored(impl_->lu, impl_->piv, impl_->n, rhs, nrhs);
}

struct SuperLUSamePatternSolver::Impl {
    int n{};
    std::vector<int> cols;
    std::vector<int> rows;
    std::vector<double> lu;
    std::vector<int> piv;
    bool factored{false};
    std::size_t factorization_count{};
};

SuperLUSamePatternSolver::SuperLUSamePatternSolver(const SparseMatrixCSC& matrix)
    : impl_(std::make_unique<Impl>()) {
    impl_->n = matrix.rows();
    impl_->cols = matrix.col_ptr();
    impl_->rows = matrix.row_ind();
    impl_->lu = dense_column_major(matrix);
    impl_->piv.resize(static_cast<std::size_t>(impl_->n));
    factor(impl_->lu, impl_->piv, impl_->n);
    impl_->factored = true;
    impl_->factorization_count = 1;
}
SuperLUSamePatternSolver::~SuperLUSamePatternSolver() = default;
SuperLUSamePatternSolver::SuperLUSamePatternSolver(SuperLUSamePatternSolver&&) noexcept = default;
SuperLUSamePatternSolver& SuperLUSamePatternSolver::operator=(SuperLUSamePatternSolver&&) noexcept = default;
int SuperLUSamePatternSolver::size() const { return impl_->n; }
std::size_t SuperLUSamePatternSolver::factorizations() const { return impl_->factorization_count; }
std::vector<double> SuperLUSamePatternSolver::solve_current(const std::vector<double>& rhs) const {
    if (!impl_->factored)
        throw std::runtime_error("Dense fallback factors invalid after failed refactorization");
    return solve_factored(impl_->lu, impl_->piv, impl_->n, rhs, 1);
}
std::vector<double> SuperLUSamePatternSolver::refactor_and_solve(
    const SparseMatrixCSC& matrix, const std::vector<double>& rhs) {
    validate_pattern(matrix, impl_->n, impl_->cols, impl_->rows);
    auto next_lu = dense_column_major(matrix);
    auto next_piv = impl_->piv;
    try {
        factor(next_lu, next_piv, impl_->n);
    } catch (...) {
        impl_->factored = false;
        throw;
    }
    impl_->lu = std::move(next_lu);
    impl_->piv = std::move(next_piv);
    impl_->factored = true;
    ++impl_->factorization_count;
    return solve_current(rhs);
}

std::vector<double> superlu_solve_once(const SparseMatrixCSC& matrix,
                                       const std::vector<double>& rhs) {
    SuperLUFactor factorization(matrix);
    return factorization.solve(rhs);
}

}  // namespace quake
