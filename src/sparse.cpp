#include "quake/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace quake {

SparseMatrixCSC::SparseMatrixCSC(int rows, int cols, std::vector<int> col_ptr,
                                 std::vector<int> row_ind, std::vector<double> values)
    : rows_(rows), cols_(cols), col_ptr_(std::move(col_ptr)),
      row_ind_(std::move(row_ind)), values_(std::move(values)) {
    if (rows_ < 0 || cols_ < 0 || static_cast<int>(col_ptr_.size()) != cols_ + 1 ||
        row_ind_.size() != values_.size()) {
        throw std::invalid_argument("Invalid CSC matrix dimensions/storage");
    }
}

SparseMatrixCSC SparseMatrixCSC::from_triplets(int rows, int cols,
                                                const std::vector<Triplet>& triplets,
                                                double drop_tol) {
    std::vector<std::map<int, double>> columns(static_cast<std::size_t>(cols));
    for (const auto& t : triplets) {
        if (t.row < 0 || t.row >= rows || t.col < 0 || t.col >= cols) {
            throw std::out_of_range("Triplet index out of range");
        }
        columns[static_cast<std::size_t>(t.col)][t.row] += t.value;
    }

    std::vector<int> col_ptr(static_cast<std::size_t>(cols + 1), 0);
    std::vector<int> row_ind;
    std::vector<double> values;
    for (int c = 0; c < cols; ++c) {
        for (const auto& [r, v] : columns[static_cast<std::size_t>(c)]) {
            if (std::abs(v) > drop_tol) {
                row_ind.push_back(r);
                values.push_back(v);
            }
        }
        col_ptr[static_cast<std::size_t>(c + 1)] = static_cast<int>(values.size());
    }
    return SparseMatrixCSC(rows, cols, std::move(col_ptr), std::move(row_ind), std::move(values));
}

std::vector<double> SparseMatrixCSC::multiply(const std::vector<double>& x) const {
    if (static_cast<int>(x.size()) != cols_) {
        throw std::invalid_argument("Sparse multiply dimension mismatch");
    }
    std::vector<double> y(static_cast<std::size_t>(rows_), 0.0);
    for (int c = 0; c < cols_; ++c) {
        for (int p = col_ptr_[static_cast<std::size_t>(c)];
             p < col_ptr_[static_cast<std::size_t>(c + 1)]; ++p) {
            y[static_cast<std::size_t>(row_ind_[static_cast<std::size_t>(p)])] +=
                values_[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(c)];
        }
    }
    return y;
}

std::vector<double> SparseMatrixCSC::diagonal() const {
    const int n = std::min(rows_, cols_);
    std::vector<double> d(static_cast<std::size_t>(n), 0.0);
    for (int c = 0; c < n; ++c) {
        for (int p = col_ptr_[static_cast<std::size_t>(c)];
             p < col_ptr_[static_cast<std::size_t>(c + 1)]; ++p) {
            if (row_ind_[static_cast<std::size_t>(p)] == c) {
                d[static_cast<std::size_t>(c)] = values_[static_cast<std::size_t>(p)];
                break;
            }
        }
    }
    return d;
}

double SparseMatrixCSC::quadratic_form(const std::vector<double>& x) const {
    auto y = multiply(x);
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += x[i] * y[i];
    return s;
}

SparseMatrixCSC add(const SparseMatrixCSC& a, const SparseMatrixCSC& b,
                    double scale_a, double scale_b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("Sparse add dimension mismatch");
    }
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(a.nnz() + b.nnz()));
    for (int c = 0; c < a.cols(); ++c) {
        for (int p = a.col_ptr()[static_cast<std::size_t>(c)];
             p < a.col_ptr()[static_cast<std::size_t>(c + 1)]; ++p) {
            triplets.push_back({a.row_ind()[static_cast<std::size_t>(p)], c,
                                scale_a * a.values()[static_cast<std::size_t>(p)]});
        }
        for (int p = b.col_ptr()[static_cast<std::size_t>(c)];
             p < b.col_ptr()[static_cast<std::size_t>(c + 1)]; ++p) {
            triplets.push_back({b.row_ind()[static_cast<std::size_t>(p)], c,
                                scale_b * b.values()[static_cast<std::size_t>(p)]});
        }
    }
    return SparseMatrixCSC::from_triplets(a.rows(), a.cols(), triplets, 1e-18);
}

SparseMatrixCSC add_diagonal(const SparseMatrixCSC& a,
                             const std::vector<double>& diagonal) {
    if (static_cast<int>(diagonal.size()) != a.rows() || a.rows() != a.cols()) {
        throw std::invalid_argument("add_diagonal dimension mismatch");
    }
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(a.nnz() + a.rows()));
    for (int c = 0; c < a.cols(); ++c) {
        for (int p = a.col_ptr()[static_cast<std::size_t>(c)];
             p < a.col_ptr()[static_cast<std::size_t>(c + 1)]; ++p) {
            triplets.push_back({a.row_ind()[static_cast<std::size_t>(p)], c,
                                a.values()[static_cast<std::size_t>(p)]});
        }
        triplets.push_back({c, c, diagonal[static_cast<std::size_t>(c)]});
    }
    return SparseMatrixCSC::from_triplets(a.rows(), a.cols(), triplets, 1e-18);
}

} // namespace quake
