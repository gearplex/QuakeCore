#include "quake/low_rank_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace quake {

std::vector<double> dense_solve(std::vector<double> A, std::vector<double> b, int n, DenseSolveDiagnostics* diagnostics) {
    if (n == 0) { if(diagnostics)*diagnostics={0.0,0.0,1.0}; return {}; }
    if (static_cast<int>(A.size()) != n * n || static_cast<int>(b.size()) != n)
        throw std::invalid_argument("dense_solve dimension mismatch");
    double minp=std::numeric_limits<double>::infinity(),maxp=0.0;
    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double best = std::abs(A[static_cast<std::size_t>(k * n + k)]);
        for (int i = k + 1; i < n; ++i) {
            const double v = std::abs(A[static_cast<std::size_t>(i * n + k)]);
            if (v > best) { best = v; pivot = i; }
        }
        if (best < 1e-18) { if(diagnostics)*diagnostics={0.0,maxp,0.0}; throw std::runtime_error("dense_solve singular matrix"); }
        minp=std::min(minp,best);maxp=std::max(maxp,best);
        if (pivot != k) {
            for (int j = k; j < n; ++j)
                std::swap(A[static_cast<std::size_t>(k*n+j)], A[static_cast<std::size_t>(pivot*n+j)]);
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(pivot)]);
        }
        const double akk = A[static_cast<std::size_t>(k*n+k)];
        for (int i = k + 1; i < n; ++i) {
            const double f = A[static_cast<std::size_t>(i*n+k)] / akk;
            A[static_cast<std::size_t>(i*n+k)] = 0.0;
            for (int j = k + 1; j < n; ++j)
                A[static_cast<std::size_t>(i*n+j)] -= f * A[static_cast<std::size_t>(k*n+j)];
            b[static_cast<std::size_t>(i)] -= f * b[static_cast<std::size_t>(k)];
        }
    }
    std::vector<double> x(static_cast<std::size_t>(n));
    for (int i = n - 1; i >= 0; --i) {
        double s = b[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j)
            s -= A[static_cast<std::size_t>(i*n+j)] * x[static_cast<std::size_t>(j)];
        x[static_cast<std::size_t>(i)] = s / A[static_cast<std::size_t>(i*n+i)];
    }
    if(diagnostics){diagnostics->min_abs_pivot=std::isfinite(minp)?minp:0.0;diagnostics->max_abs_pivot=maxp;diagnostics->pivot_ratio=maxp>0.0?diagnostics->min_abs_pivot/maxp:1.0;}
    return x;
}

LowRankWoodburySolver::LowRankWoodburySolver(const SparseMatrixCSC& baseline,
                                             SparseUpdateBasis basis)
    : n_(baseline.rows()), m_(basis.cols()), factor_(baseline),
      basis_(std::move(basis)) {
    if (baseline.rows() != baseline.cols() || basis_.rows() != n_)
        throw std::invalid_argument("LowRankWoodburySolver dimension mismatch");
    const auto dense_B=basis_.dense_column_major();
    W_ = factor_.solve_multiple(dense_B, m_);
    G_.assign(static_cast<std::size_t>(m_ * m_), 0.0);
    for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < m_; ++j) {
            double g=0.0;
            const int pb=basis_.col_ptr()[static_cast<std::size_t>(i)];
            const int pe=basis_.col_ptr()[static_cast<std::size_t>(i+1)];
            for(int p=pb;p<pe;++p){
                const int row=basis_.row_ind()[static_cast<std::size_t>(p)];
                g += basis_.values()[static_cast<std::size_t>(p)] * W_[static_cast<std::size_t>(j*n_+row)];
            }
            G_[static_cast<std::size_t>(i*m_+j)] = g;
        }
    }
}

LowRankWoodburySolver::LowRankWoodburySolver(const SparseMatrixCSC& baseline,
                                             std::vector<double> B_column_major,
                                             int update_count)
    : LowRankWoodburySolver(baseline,
          SparseUpdateBasis::from_dense(baseline.rows(),update_count,B_column_major)) {}

std::vector<double> LowRankWoodburySolver::solve(const std::vector<double>& rhs,
                                                 const std::vector<double>& delta_k,
                                                 double active_tol) {
    if (static_cast<int>(rhs.size()) != n_ || static_cast<int>(delta_k.size()) != m_)
        throw std::invalid_argument("Woodbury solve dimension mismatch");

    auto y = factor_.solve(rhs);
    std::vector<int> active;
    active.reserve(static_cast<std::size_t>(m_));
    for (int i = 0; i < m_; ++i)
        if (std::abs(delta_k[static_cast<std::size_t>(i)]) > active_tol)
            active.push_back(i);
    if (active.empty()) return y;

    const int r = static_cast<int>(active.size());
    std::vector<double> z(static_cast<std::size_t>(r), 0.0);
    for (int a = 0; a < r; ++a) {
        const int i = active[static_cast<std::size_t>(a)];
        z[static_cast<std::size_t>(a)] = basis_.column_dot(i,y);
    }

    // S = I + D * G_active, rhs_small = D * B^T y.
    std::vector<double> S(static_cast<std::size_t>(r*r), 0.0);
    std::vector<double> q(static_cast<std::size_t>(r), 0.0);
    for (int a = 0; a < r; ++a) {
        const int i = active[static_cast<std::size_t>(a)];
        const double di = delta_k[static_cast<std::size_t>(i)];
        q[static_cast<std::size_t>(a)] = di * z[static_cast<std::size_t>(a)];
        for (int b = 0; b < r; ++b) {
            const int j = active[static_cast<std::size_t>(b)];
            S[static_cast<std::size_t>(a*r+b)] = (a == b ? 1.0 : 0.0) +
                di * G_[static_cast<std::size_t>(i*m_+j)];
        }
    }
    DenseSolveDiagnostics diag; const auto alpha = dense_solve(std::move(S), std::move(q), r, &diag); last_pivot_ratio_=diag.pivot_ratio;
    for (int k = 0; k < n_; ++k) {
        double corr = 0.0;
        for (int a = 0; a < r; ++a) {
            const int i = active[static_cast<std::size_t>(a)];
            corr += W_[static_cast<std::size_t>(i*n_+k)] * alpha[static_cast<std::size_t>(a)];
        }
        y[static_cast<std::size_t>(k)] -= corr;
    }
    return y;
}


LazyLowRankWoodburySolver::LazyLowRankWoodburySolver(const SparseMatrixCSC& baseline,
                                                     SparseUpdateBasis basis)
    : n_(baseline.rows()), m_(basis.cols()), factor_(baseline),
      basis_(std::move(basis)),
      index_to_cache_(static_cast<std::size_t>(m_), -1) {
    if (baseline.rows() != baseline.cols() || basis_.rows() != n_)
        throw std::invalid_argument("LazyLowRankWoodburySolver dimension mismatch");
}

LazyLowRankWoodburySolver::LazyLowRankWoodburySolver(const SparseMatrixCSC& baseline,
                                                     std::vector<double> B_column_major,
                                                     int update_count)
    : LazyLowRankWoodburySolver(baseline,
          SparseUpdateBasis::from_dense(baseline.rows(),update_count,B_column_major)) {}

void LazyLowRankWoodburySolver::cache_update(int update_index) {
    if (update_index < 0 || update_index >= m_) throw std::out_of_range("low-rank update index");
    if (index_to_cache_[static_cast<std::size_t>(update_index)] >= 0) return;

    auto b = basis_.dense_column(update_index);
    auto w = factor_.solve(b);

    const int old_k = static_cast<int>(cached_indices_.size());
    const int new_k = old_k + 1;
    std::vector<double> newG(static_cast<std::size_t>(new_k*new_k),0.0);
    for(int i=0;i<old_k;++i)
        for(int j=0;j<old_k;++j)
            newG[static_cast<std::size_t>(i*new_k+j)] = G_cache_[static_cast<std::size_t>(i*old_k+j)];

    for(int c=0;c<old_k;++c) {
        const int candidate = cached_indices_[static_cast<std::size_t>(c)];
        const double g_old_new=basis_.column_dot(candidate,w);
        double g_new_old=0.0;
        const int pb=basis_.col_ptr()[static_cast<std::size_t>(update_index)];
        const int pe=basis_.col_ptr()[static_cast<std::size_t>(update_index+1)];
        for(int p=pb;p<pe;++p){const int row=basis_.row_ind()[static_cast<std::size_t>(p)];g_new_old+=basis_.values()[static_cast<std::size_t>(p)]*W_cache_[static_cast<std::size_t>(c*n_+row)];}
        newG[static_cast<std::size_t>(c*new_k+old_k)] = g_old_new;
        newG[static_cast<std::size_t>(old_k*new_k+c)] = g_new_old;
    }
    const double gnn=basis_.column_dot(update_index,w);
    newG[static_cast<std::size_t>(old_k*new_k+old_k)] = gnn;

    W_cache_.insert(W_cache_.end(),w.begin(),w.end());
    G_cache_=std::move(newG);
    index_to_cache_[static_cast<std::size_t>(update_index)] = old_k;
    cached_indices_.push_back(update_index);
}

std::vector<double> LazyLowRankWoodburySolver::solve(const std::vector<double>& rhs,
                                                     const std::vector<double>& delta_k,
                                                     double active_tol) {
    if (static_cast<int>(rhs.size()) != n_ || static_cast<int>(delta_k.size()) != m_)
        throw std::invalid_argument("Lazy Woodbury solve dimension mismatch");

    auto y = factor_.solve(rhs);
    std::vector<int> active;
    active.reserve(static_cast<std::size_t>(m_));
    for(int i=0;i<m_;++i) {
        if(std::abs(delta_k[static_cast<std::size_t>(i)]) > active_tol) {
            cache_update(i);
            active.push_back(i);
        }
    }
    if(active.empty()){last_pivot_ratio_=1.0;return y;}

    const int r=static_cast<int>(active.size());
    const int kcache=static_cast<int>(cached_indices_.size());
    std::vector<double> S(static_cast<std::size_t>(r*r),0.0);
    std::vector<double> q(static_cast<std::size_t>(r),0.0);
    for(int a=0;a<r;++a) {
        const int ia=active[static_cast<std::size_t>(a)];
        const int ca=index_to_cache_[static_cast<std::size_t>(ia)];
        const double z=basis_.column_dot(ia,y);
        const double di=delta_k[static_cast<std::size_t>(ia)];
        q[static_cast<std::size_t>(a)] = di*z;
        for(int b=0;b<r;++b) {
            const int ib=active[static_cast<std::size_t>(b)];
            const int cb=index_to_cache_[static_cast<std::size_t>(ib)];
            S[static_cast<std::size_t>(a*r+b)] = (a==b?1.0:0.0) + di*G_cache_[static_cast<std::size_t>(ca*kcache+cb)];
        }
    }
    DenseSolveDiagnostics diag; const auto alpha=dense_solve(std::move(S),std::move(q),r,&diag); last_pivot_ratio_=diag.pivot_ratio;
    for(int i=0;i<n_;++i) {
        double corr=0.0;
        for(int a=0;a<r;++a) {
            const int ia=active[static_cast<std::size_t>(a)];
            const int ca=index_to_cache_[static_cast<std::size_t>(ia)];
            corr += W_cache_[static_cast<std::size_t>(ca*n_+i)] * alpha[static_cast<std::size_t>(a)];
        }
        y[static_cast<std::size_t>(i)] -= corr;
    }
    return y;
}

} // namespace quake
