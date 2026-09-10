#include "quake/stability.hpp"
#include "quake/superlu_solver.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>

extern "C" {
void dsyev_(const char* jobz,const char* uplo,const int* n,double* a,const int* lda,
            double* w,double* work,const int* lwork,int* info);
}

namespace quake {
namespace {
double dot(const std::vector<double>& a,const std::vector<double>& b){double s=0;for(std::size_t i=0;i<a.size();++i)s+=a[i]*b[i];return s;}
double norm2(const std::vector<double>& a){return std::sqrt(std::max(0.0,dot(a,a)));}

struct SymEntry { double lower_sum{0.0},upper_sum{0.0}; int lower_count{0},upper_count{0}; };

using Pair = std::pair<int,int>;

struct SymmetricCanonical {
    std::map<Pair,double> entries; // key=(min,max)
    double relative_error{0.0};
    double max_abs_diagonal{0.0};
};

SymmetricCanonical canonicalize_symmetric(const SparseMatrixCSC& a){
    if(a.rows()!=a.cols()) throw std::invalid_argument("stability matrix must be square");
    const int n=a.rows();
    std::map<Pair,SymEntry> acc;
    double max_scale=0.0,max_diff=0.0;
    for(int c=0;c<n;++c){
        for(int p=a.col_ptr()[static_cast<std::size_t>(c)];p<a.col_ptr()[static_cast<std::size_t>(c+1)];++p){
            const int r=a.row_ind()[static_cast<std::size_t>(p)];const double v=a.values()[static_cast<std::size_t>(p)];
            const int lo=std::min(r,c),hi=std::max(r,c);auto& e=acc[{lo,hi}];
            if(r>=c){e.lower_sum+=v;++e.lower_count;}else{e.upper_sum+=v;++e.upper_count;}
        }
    }
    SymmetricCanonical out;
    for(const auto& [ij,e]:acc){
        const int i=ij.first,j=ij.second;double v=0.0;
        if(i==j){
            const int cnt=e.lower_count+e.upper_count;v=(e.lower_sum+e.upper_sum)/std::max(1,cnt);out.max_abs_diagonal=std::max(out.max_abs_diagonal,std::abs(v));
        }else if(e.lower_count>0&&e.upper_count>0){
            const double lo=e.lower_sum/e.lower_count,up=e.upper_sum/e.upper_count;v=0.5*(lo+up);max_diff=std::max(max_diff,std::abs(lo-up));max_scale=std::max({max_scale,std::abs(lo),std::abs(up)});
        }else if(e.lower_count>0){v=e.lower_sum/e.lower_count;}
        else{v=e.upper_sum/e.upper_count;}
        if(v!=0.0) out.entries[ij]=v;
    }
    out.relative_error=max_diff/std::max(1.0,max_scale);
    return out;
}

std::vector<int> reverse_cuthill_mckee(int n,const std::map<Pair,double>& entries){
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n));
    for(const auto& [ij,v]:entries){
        if(ij.first==ij.second||v==0.0) continue;
        adj[static_cast<std::size_t>(ij.first)].push_back(ij.second);
        adj[static_cast<std::size_t>(ij.second)].push_back(ij.first);
    }
    for(auto& v:adj){std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}
    std::vector<int> deg(static_cast<std::size_t>(n));for(int i=0;i<n;++i)deg[static_cast<std::size_t>(i)]=static_cast<int>(adj[static_cast<std::size_t>(i)].size());
    std::vector<char> seen(static_cast<std::size_t>(n),0);std::vector<int> order;order.reserve(static_cast<std::size_t>(n));
    while(static_cast<int>(order.size())<n){
        int root=-1;for(int i=0;i<n;++i)if(!seen[static_cast<std::size_t>(i)]&&(root<0||deg[static_cast<std::size_t>(i)]<deg[static_cast<std::size_t>(root)]))root=i;
        if(root<0) break;
        std::deque<int> q;
        q.push_back(root);
        seen[static_cast<std::size_t>(root)]=1;
        while(!q.empty()){
            const int u=q.front();q.pop_front();order.push_back(u);auto nbr=adj[static_cast<std::size_t>(u)];
            std::sort(nbr.begin(),nbr.end(),[&](int x,int y){if(deg[static_cast<std::size_t>(x)]!=deg[static_cast<std::size_t>(y)])return deg[static_cast<std::size_t>(x)]<deg[static_cast<std::size_t>(y)];return x<y;});
            for(int v:nbr)if(!seen[static_cast<std::size_t>(v)]){seen[static_cast<std::size_t>(v)]=1;q.push_back(v);}
        }
    }
    std::reverse(order.begin(),order.end());return order; // new index -> old index
}

void dense_eigen_cross_check(const SymmetricCanonical& canon,int n,double tol,InitialStabilityAssessment& out){
    std::vector<double> dense(static_cast<std::size_t>(n*n),0.0);
    for(const auto& [ij,v]:canon.entries){const int i=ij.first,j=ij.second;dense[static_cast<std::size_t>(j*n+i)]=v;dense[static_cast<std::size_t>(i*n+j)]=v;}
    std::vector<double> w(static_cast<std::size_t>(n));const char job='N',uplo='L';const int lda=n;int info=0,lwork=-1;double query=0.0;
    dsyev_(&job,&uplo,&n,dense.data(),&lda,w.data(),&query,&lwork,&info);if(info!=0)return;lwork=std::max(3*n-1,static_cast<int>(std::ceil(query)));std::vector<double> work(static_cast<std::size_t>(lwork));
    // Rebuild because workspace query may inspect A.
    std::fill(dense.begin(),dense.end(),0.0);for(const auto& [ij,v]:canon.entries){const int i=ij.first,j=ij.second;dense[static_cast<std::size_t>(j*n+i)]=v;dense[static_cast<std::size_t>(i*n+j)]=v;}
    dsyev_(&job,&uplo,&n,dense.data(),&lda,w.data(),work.data(),&lwork,&info);if(info!=0)return;
    out.dense_check_performed=true;out.dense_minimum_eigenvalue=w.empty()?std::numeric_limits<double>::quiet_NaN():w.front();
    const double scale=std::max(1.0,canon.max_abs_diagonal);for(double x:w){if(x<-tol*scale)++out.dense_negative_eigenvalues;else if(std::abs(x)<=tol*scale)++out.dense_near_zero_eigenvalues;}
}
}

InitialStabilityAssessment assess_positive_definiteness(const SparseMatrixCSC& matrix,double rel_tol,int dense_limit){
    if(rel_tol<=0.0)throw std::invalid_argument("positive-definiteness tolerance must be positive");
    InitialStabilityAssessment out;out.dof=matrix.rows();
    if(matrix.rows()!=matrix.cols()||matrix.rows()<=0){out.status=InitialStabilityStatus::NumericalFailure;return out;}
    const int n=matrix.rows();const auto canon=canonicalize_symmetric(matrix);out.symmetry_relative_error=canon.relative_error;
    if(canon.relative_error>1e-8){out.status=InitialStabilityStatus::NumericalFailure;return out;}
    if(dense_limit>0&&n<=dense_limit)dense_eigen_cross_check(canon,n,rel_tol,out);
    const auto new_to_old=reverse_cuthill_mckee(n,canon.entries);if(static_cast<int>(new_to_old.size())!=n){out.status=InitialStabilityStatus::NumericalFailure;return out;}
    std::vector<int> old_to_new(static_cast<std::size_t>(n));for(int ni=0;ni<n;++ni)old_to_new[static_cast<std::size_t>(new_to_old[static_cast<std::size_t>(ni)])]=ni;
    std::vector<std::map<int,double>> col(static_cast<std::size_t>(n));
    for(const auto& [ij,v]:canon.entries){int i=old_to_new[static_cast<std::size_t>(ij.first)],j=old_to_new[static_cast<std::size_t>(ij.second)];if(i<j)std::swap(i,j);col[static_cast<std::size_t>(j)][i]+=v;}
    const double scale=std::max(1.0,canon.max_abs_diagonal);const double pivot_tol=rel_tol*scale;out.minimum_pivot=std::numeric_limits<double>::infinity();
    for(int k=0;k<n;++k){
        auto& ck=col[static_cast<std::size_t>(k)];auto dit=ck.find(k);const double pivot=dit==ck.end()?0.0:dit->second;out.minimum_pivot=std::min(out.minimum_pivot,pivot);out.minimum_pivot_ratio=out.minimum_pivot/scale;
        if(!std::isfinite(pivot)||pivot<=pivot_tol){out.status=InitialStabilityStatus::NotPositiveDefinite;out.positive_definite=false;out.failing_pivot=k;return out;}
        std::vector<std::pair<int,double>> nbr;nbr.reserve(ck.size());for(auto it=ck.upper_bound(k);it!=ck.end();++it)if(std::abs(it->second)>1e-30)nbr.push_back(*it);
        out.factor_nonzeros+=1+nbr.size();
        for(std::size_t a=0;a<nbr.size();++a){const int j=nbr[a].first;const double ajk=nbr[a].second;auto& cj=col[static_cast<std::size_t>(j)];
            for(std::size_t b=a;b<nbr.size();++b){const int i=nbr[b].first;const double aik=nbr[b].second;const double dv=-(aik*ajk)/pivot;auto it=cj.find(i);if(it==cj.end()){if(std::abs(dv)>1e-30)cj.emplace(i,dv);}else{it->second+=dv;if(std::abs(it->second)<1e-30&&i!=j)cj.erase(it);}}
        }
        std::map<int,double>().swap(ck);
    }
    out.status=InitialStabilityStatus::PositiveDefinite;out.positive_definite=true;out.failing_pivot=-1;return out;
}

InitialStabilityAssessment assess_initial_stability_from_tangents(const NonlinearDynamicModel& model,const std::vector<double>& u,const std::vector<double>& tang,double rel_tol,int dense_limit){
    if(static_cast<int>(u.size())!=model.dof()||static_cast<int>(tang.size())!=model.nonlinear_count())throw std::invalid_argument("invalid initial stability input");
    return assess_positive_definiteness(model.effective_state_tangent_matrix(u,tang,0.0,0.0),rel_tol,dense_limit);
}

namespace {
TangentStabilityEstimate estimate_matrix_stability(const SparseMatrixCSC& K,int iters){
    TangentStabilityEstimate out;
    try{
        SuperLUFactor factor(K);out.factorization_ok=true;
        std::vector<double> x(static_cast<std::size_t>(K.rows()),1.0);
        for(std::size_t i=0;i<x.size();++i)x[i]=1.0+0.013*static_cast<double>((i%17)+1);
        double n=norm2(x);for(double& v:x)v/=n;
        for(int k=0;k<iters;++k){auto y=factor.solve(x);n=norm2(y);if(!std::isfinite(n)||n<1e-30){out.factorization_ok=false;return out;}for(double& v:y)v/=n;x=std::move(y);}
        auto kx=K.multiply(x);out.near_zero_eigenvalue=dot(x,kx);out.rayleigh_eigenvalue=out.near_zero_eigenvalue;out.absolute_eigenvalue=std::abs(out.near_zero_eigenvalue);
    }catch(const std::exception&){out.factorization_ok=false;}
    return out;
}
}
TangentStabilityEstimate estimate_tangent_stability_from_tangents(const NonlinearDynamicModel& model,const std::vector<double>& u,const std::vector<double>& t,int iters){
    if(static_cast<int>(u.size())!=model.dof()||static_cast<int>(t.size())!=model.nonlinear_count()||iters<1)
        throw std::invalid_argument("invalid tangent-stability input");
    return estimate_matrix_stability(model.effective_state_tangent_matrix(u,t,0.0,0.0),iters);
}
TangentStabilityEstimate estimate_tangent_stability(const NonlinearDynamicModel& model,
                                                     const std::vector<double>& u,
                                                     const std::vector<double>& committed,
                                                     int iters){
    if(static_cast<int>(committed.size())!=model.nonlinear_state_size()) throw std::invalid_argument("invalid tangent-stability state");
    std::vector<double> f,t,trial;model.internal_force_and_tangent(u,committed,f,t,trial);
    if(iters<1)throw std::invalid_argument("invalid tangent-stability iterations");
    return estimate_matrix_stability(model.effective_state_tangent_matrix_with_state(u,t,committed,0.0,0.0),iters);
}
} // namespace quake
