#pragma once
#include <cstddef>
#include <tuple>
#include <vector>

namespace quake {

struct Triplet {
    int row{};
    int col{};
    double value{};
};

class SparseMatrixCSC {
public:
    SparseMatrixCSC() = default;
    SparseMatrixCSC(int rows, int cols, std::vector<int> col_ptr,
                    std::vector<int> row_ind, std::vector<double> values);

    static SparseMatrixCSC from_triplets(int rows, int cols,
                                         const std::vector<Triplet>& triplets,
                                         double drop_tol = 0.0);

    int rows() const { return rows_; }
    int cols() const { return cols_; }
    int nnz() const { return static_cast<int>(values_.size()); }

    const std::vector<int>& col_ptr() const { return col_ptr_; }
    const std::vector<int>& row_ind() const { return row_ind_; }
    const std::vector<double>& values() const { return values_; }
    std::vector<double>& values() { return values_; }

    std::vector<double> multiply(const std::vector<double>& x) const;
    std::vector<double> diagonal() const;
    double quadratic_form(const std::vector<double>& x) const;

private:
    int rows_{0};
    int cols_{0};
    std::vector<int> col_ptr_;
    std::vector<int> row_ind_;
    std::vector<double> values_;
};

SparseMatrixCSC add(const SparseMatrixCSC& a, const SparseMatrixCSC& b,
                    double scale_a = 1.0, double scale_b = 1.0);
SparseMatrixCSC add_diagonal(const SparseMatrixCSC& a,
                             const std::vector<double>& diagonal);

} // namespace quake
