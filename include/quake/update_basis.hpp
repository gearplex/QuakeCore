#pragma once
#include <utility>
#include <vector>

namespace quake {

// Sparse column-oriented collection of generalized nonlinear deformation
// operators. Each column b_j maps global generalized displacement u to one
// nonlinear deformation q_j=b_j^T u. Concentrated-plasticity columns are
// normally extremely sparse, so storing n*m dense values is intentionally
// avoided.
class SparseUpdateBasis {
public:
    SparseUpdateBasis() = default;
    SparseUpdateBasis(int rows, int cols, std::vector<int> col_ptr,
                      std::vector<int> row_ind, std::vector<double> values);

    static SparseUpdateBasis from_dense(int rows, int cols,
                                        const std::vector<double>& column_major,
                                        double drop_tol = 0.0);
    static SparseUpdateBasis from_columns(
        int rows,
        const std::vector<std::vector<std::pair<int,double>>>& columns,
        double drop_tol = 0.0);

    int rows() const { return rows_; }
    int cols() const { return cols_; }
    int nnz() const { return static_cast<int>(values_.size()); }
    const std::vector<int>& col_ptr() const { return col_ptr_; }
    const std::vector<int>& row_ind() const { return row_ind_; }
    const std::vector<double>& values() const { return values_; }

    double column_dot(int col, const std::vector<double>& x) const;
    std::vector<double> dense_column(int col) const;
    std::vector<double> dense_column_major() const;
    void axpy_column(int col, double scale, std::vector<double>& y) const;
    double value(int row, int col) const;

private:
    int rows_{0};
    int cols_{0};
    std::vector<int> col_ptr_;
    std::vector<int> row_ind_;
    std::vector<double> values_;
};

} // namespace quake
