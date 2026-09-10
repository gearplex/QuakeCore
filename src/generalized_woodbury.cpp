#include "quake/generalized_woodbury.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {

GeneralizedWoodburySolver::GeneralizedWoodburySolver(const SparseMatrixCSC& baseline,
                                                     SparseUpdateBasis basis)
    : n_(baseline.rows()), r_(basis.cols()), factor_(baseline), basis_(std::move(basis)) {
    if(baseline.rows()!=baseline.cols()||basis_.rows()!=n_)throw std::invalid_argument("GeneralizedWoodburySolver dimension mismatch");
    if(r_==0)return;
    W_=factor_.solve_multiple(basis_.dense_column_major(),r_);
    H_.assign(static_cast<std::size_t>(r_*r_),0.0);
    for(int i=0;i<r_;++i)for(int j=0;j<r_;++j){
        double h=0.0;const int b=basis_.col_ptr()[static_cast<std::size_t>(i)],e=basis_.col_ptr()[static_cast<std::size_t>(i+1)];
        for(int p=b;p<e;++p){const int row=basis_.row_ind()[static_cast<std::size_t>(p)];h+=basis_.values()[static_cast<std::size_t>(p)]*W_[static_cast<std::size_t>(j*n_+row)];}
        H_[static_cast<std::size_t>(i*r_+j)]=h;
    }
}

std::vector<double> GeneralizedWoodburySolver::solve(const std::vector<double>& rhs,
                                                     const std::vector<double>& C,
                                                     double active_tol){
    if(static_cast<int>(rhs.size())!=n_||static_cast<int>(C.size())!=r_*r_)throw std::invalid_argument("generalized Woodbury solve dimension mismatch");
    auto y=factor_.solve(rhs);if(r_==0){last_active_dimension_=0;last_pivot_ratio_=1.0;return y;}
    std::vector<int> active;active.reserve(static_cast<std::size_t>(r_));
    for(int i=0;i<r_;++i){double mx=0.0;for(int j=0;j<r_;++j)mx=std::max({mx,std::abs(C[static_cast<std::size_t>(i*r_+j)]),std::abs(C[static_cast<std::size_t>(j*r_+i)])});if(mx>active_tol)active.push_back(i);}
    last_active_dimension_=active.size();if(active.empty()){last_pivot_ratio_=1.0;return y;}
    const int aN=static_cast<int>(active.size());
    std::vector<double> z(static_cast<std::size_t>(aN),0.0);
    for(int a=0;a<aN;++a)z[static_cast<std::size_t>(a)]=basis_.column_dot(active[static_cast<std::size_t>(a)],y);
    std::vector<double> Ca(static_cast<std::size_t>(aN*aN),0.0);
    for(int a=0;a<aN;++a)for(int b=0;b<aN;++b)Ca[static_cast<std::size_t>(a*aN+b)]=C[static_cast<std::size_t>(active[static_cast<std::size_t>(a)]*r_+active[static_cast<std::size_t>(b)])];
    std::vector<double> q(static_cast<std::size_t>(aN),0.0);
    for(int a=0;a<aN;++a)for(int b=0;b<aN;++b)q[static_cast<std::size_t>(a)]+=Ca[static_cast<std::size_t>(a*aN+b)]*z[static_cast<std::size_t>(b)];
    // S = I + C * H_active. This form does not require C^{-1}, which matters
    // because individual geometric blocks may be rank deficient at some states.
    std::vector<double> S(static_cast<std::size_t>(aN*aN),0.0);
    for(int a=0;a<aN;++a)for(int b=0;b<aN;++b){double v=(a==b?1.0:0.0);for(int k=0;k<aN;++k)v+=Ca[static_cast<std::size_t>(a*aN+k)]*H_[static_cast<std::size_t>(active[static_cast<std::size_t>(k)]*r_+active[static_cast<std::size_t>(b)])];S[static_cast<std::size_t>(a*aN+b)]=v;}
    DenseSolveDiagnostics diag;const auto alpha=dense_solve(std::move(S),std::move(q),aN,&diag);last_pivot_ratio_=diag.pivot_ratio;
    for(int i=0;i<n_;++i){double corr=0.0;for(int a=0;a<aN;++a)corr+=W_[static_cast<std::size_t>(active[static_cast<std::size_t>(a)]*n_+i)]*alpha[static_cast<std::size_t>(a)];y[static_cast<std::size_t>(i)]-=corr;}
    return y;
}

} // namespace quake
