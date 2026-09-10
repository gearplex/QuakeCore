#include "quake/reduction.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

extern "C" {
void dggev_(const char* jobvl, const char* jobvr, const int* n,
            double* a, const int* lda, double* b, const int* ldb,
            double* alphar, double* alphai, double* beta,
            double* vl, const int* ldvl, double* vr, const int* ldvr,
            double* work, const int* lwork, int* info);
}

namespace quake {
namespace {

std::vector<double> dense_from_csc(const SparseMatrixCSC& a) {
    std::vector<double> out(static_cast<std::size_t>(a.rows()*a.cols()),0.0); // col-major
    for(int c=0;c<a.cols();++c)
        for(int p=a.col_ptr()[static_cast<std::size_t>(c)];p<a.col_ptr()[static_cast<std::size_t>(c+1)];++p)
            out[static_cast<std::size_t>(c*a.rows()+a.row_ind()[static_cast<std::size_t>(p)])]+=a.values()[static_cast<std::size_t>(p)];
    return out;
}

SparseMatrixCSC dense_to_csc(int rows,int cols,const std::vector<double>& a,double tol=1e-14){
    std::vector<Triplet> t;
    for(int c=0;c<cols;++c)for(int r=0;r<rows;++r){
        const double v=a[static_cast<std::size_t>(c*rows+r)];
        if(std::abs(v)>tol || (rows==cols && r==c)) t.push_back({r,c,v});
    }
    return SparseMatrixCSC::from_triplets(rows,cols,t,-1.0);
}

std::vector<double> dense_mass(const NonlinearDynamicModel& model){
    const int n=model.dof();std::vector<double> M(static_cast<std::size_t>(n*n),0.0),e(static_cast<std::size_t>(n),0.0);
    for(int c=0;c<n;++c){std::fill(e.begin(),e.end(),0.0);e[static_cast<std::size_t>(c)]=1.0;auto me=model.mass_multiply(e);for(int r=0;r<n;++r)M[static_cast<std::size_t>(c*n+r)]=me[static_cast<std::size_t>(r)];}
    return M;
}

std::vector<double> dense_damping(const NonlinearDynamicModel& model){
    const int n=model.dof();std::vector<double> C(static_cast<std::size_t>(n*n),0.0),e(static_cast<std::size_t>(n),0.0);
    for(int c=0;c<n;++c){std::fill(e.begin(),e.end(),0.0);e[static_cast<std::size_t>(c)]=1.0;auto ce=model.damping_multiply(e);for(int r=0;r<n;++r)C[static_cast<std::size_t>(c*n+r)]=ce[static_cast<std::size_t>(r)];}
    return C;
}

// T^T A T, all matrices column-major.
std::vector<double> project_dense(const std::vector<double>& A,int nf,
                                  const std::vector<double>& T,int nr){
    std::vector<double> AT(static_cast<std::size_t>(nf*nr),0.0),R(static_cast<std::size_t>(nr*nr),0.0);
    for(int j=0;j<nr;++j)for(int k=0;k<nf;++k){
        const double tk=T[static_cast<std::size_t>(j*nf+k)];if(tk==0.0)continue;
        for(int i=0;i<nf;++i)AT[static_cast<std::size_t>(j*nf+i)]+=A[static_cast<std::size_t>(k*nf+i)]*tk;
    }
    for(int j=0;j<nr;++j)for(int i=0;i<nr;++i){double s=0.0;for(int k=0;k<nf;++k)s+=T[static_cast<std::size_t>(i*nf+k)]*AT[static_cast<std::size_t>(j*nf+k)];R[static_cast<std::size_t>(j*nr+i)]=s;}
    return R;
}

std::vector<double> finite_fixed_interface_modes(const SparseMatrixCSC& Kii,
                                                 const std::vector<double>& Mii,
                                                 int requested){
    const int n=Kii.rows(); if(requested<=0||n<=0)return {};
    auto A=dense_from_csc(Kii);auto B=Mii;
    std::vector<double> ar(static_cast<std::size_t>(n)),ai(static_cast<std::size_t>(n)),be(static_cast<std::size_t>(n)),vr(static_cast<std::size_t>(n*n));
    double vl=0.0;const int one=1,ld=n;const char N='N',V='V';int info=0,lwork=-1;double q=0.0;
    dggev_(&N,&V,&n,A.data(),&ld,B.data(),&ld,ar.data(),ai.data(),be.data(),&vl,&one,vr.data(),&ld,&q,&lwork,&info);
    if(info!=0)throw std::runtime_error("Craig-Bampton eigen workspace query failed");
    lwork=std::max(8*n,static_cast<int>(std::ceil(q)));std::vector<double> work(static_cast<std::size_t>(lwork));A=dense_from_csc(Kii);B=Mii;
    dggev_(&N,&V,&n,A.data(),&ld,B.data(),&ld,ar.data(),ai.data(),be.data(),&vl,&one,vr.data(),&ld,work.data(),&lwork,&info);
    if(info!=0)throw std::runtime_error("Craig-Bampton fixed-interface eigensolve failed");
    struct Cand{double lam;int idx;};std::vector<Cand> cs;
    for(int i=0;i<n;++i){const double sc=std::max({1.0,std::abs(ar[static_cast<std::size_t>(i)]),std::abs(be[static_cast<std::size_t>(i)])});if(std::abs(be[static_cast<std::size_t>(i)])<=1e-12*sc||std::abs(ai[static_cast<std::size_t>(i)])>1e-8*sc)continue;double l=ar[static_cast<std::size_t>(i)]/be[static_cast<std::size_t>(i)];if(std::isfinite(l)&&l>0)cs.push_back({l,i});}
    std::sort(cs.begin(),cs.end(),[](auto&a,auto&b){return a.lam<b.lam;});if(static_cast<int>(cs.size())>requested)cs.resize(static_cast<std::size_t>(requested));
    std::vector<double> out(static_cast<std::size_t>(n*cs.size()),0.0);
    for(std::size_t c=0;c<cs.size();++c){double mx=0.0;for(int r=0;r<n;++r)mx=std::max(mx,std::abs(vr[static_cast<std::size_t>(cs[c].idx*n+r)]));if(mx==0.0)continue;for(int r=0;r<n;++r)out[c*static_cast<std::size_t>(n)+static_cast<std::size_t>(r)]=vr[static_cast<std::size_t>(cs[c].idx*n+r)]/mx;}
    return out;
}

} // namespace

ReducedDynamicModel::ReducedDynamicModel(const NonlinearDynamicModel& parent,
                                         std::vector<double> transform_column_major,
                                         int reduced_dof)
    : parent_(&parent),nf_(parent.dof()),nr_(reduced_dof),transform_(std::move(transform_column_major)){
    if(nr_<=0||static_cast<int>(transform_.size())!=nf_*nr_)throw std::invalid_argument("reduction transform dimension mismatch");
    if(parent.has_state_dependent_global_tangent())throw std::invalid_argument("scalar-bank reduction does not support state-dependent/coupled element tangents");
    k_initial_=project_matrix(parent.K_initial());
    auto M=project_dense(dense_mass(parent),nf_,transform_,nr_);mass_matrix_=dense_to_csc(nr_,nr_,M,1e-13);mass_diag_=mass_matrix_.diagonal();
    auto C=project_dense(dense_damping(parent),nf_,transform_,nr_);damping_matrix_=dense_to_csc(nr_,nr_,C,1e-13);

    std::vector<std::vector<std::pair<int,double>>> cols(static_cast<std::size_t>(parent.nonlinear_count()));
    const auto& B=parent.nonlinear_basis();
    for(int j=0;j<B.cols();++j){for(int r=0;r<nr_;++r){double s=0.0;for(int p=B.col_ptr()[static_cast<std::size_t>(j)];p<B.col_ptr()[static_cast<std::size_t>(j+1)];++p){const int i=B.row_ind()[static_cast<std::size_t>(p)];s+=transform_[static_cast<std::size_t>(r*nf_+i)]*B.values()[static_cast<std::size_t>(p)];}if(std::abs(s)>1e-13)cols[static_cast<std::size_t>(j)].push_back({r,s});}}
    basis_=SparseUpdateBasis::from_columns(nr_,cols,1e-13);
    // Recover the reduced linear skeleton by subtracting the initial component
    // tangents from K_initial. This lets residual evaluation remain entirely in
    // reduced coordinates; the full parent model is touched only for its
    // constitutive bank and output reducers.
    std::vector<Triplet> lt;
    for(int c=0;c<k_initial_.cols();++c) for(int p=k_initial_.col_ptr()[static_cast<std::size_t>(c)];p<k_initial_.col_ptr()[static_cast<std::size_t>(c+1)];++p)
        lt.push_back({k_initial_.row_ind()[static_cast<std::size_t>(p)],c,k_initial_.values()[static_cast<std::size_t>(p)]});
    const auto& init=parent_->initial_nonlinear_tangents();
    for(int j=0;j<nonlinear_count();++j){const double k=init[static_cast<std::size_t>(j)];const int b=basis_.col_ptr()[static_cast<std::size_t>(j)],e=basis_.col_ptr()[static_cast<std::size_t>(j+1)];for(int pi=b;pi<e;++pi)for(int pj=b;pj<e;++pj)lt.push_back({basis_.row_ind()[static_cast<std::size_t>(pi)],basis_.row_ind()[static_cast<std::size_t>(pj)],-k*basis_.values()[static_cast<std::size_t>(pi)]*basis_.values()[static_cast<std::size_t>(pj)]});}
    k_linear_=SparseMatrixCSC::from_triplets(nr_,nr_,lt,1e-13);
}

std::vector<double> ReducedDynamicModel::expand(const std::vector<double>& q) const {
    if(static_cast<int>(q.size())!=nr_) throw std::invalid_argument("reduced vector size");
    std::vector<double> u(static_cast<std::size_t>(nf_),0.0);
    for(int c=0;c<nr_;++c){const double x=q[static_cast<std::size_t>(c)];for(int r=0;r<nf_;++r)u[static_cast<std::size_t>(r)]+=transform_[static_cast<std::size_t>(c*nf_+r)]*x;}
    return u;
}
std::vector<double> ReducedDynamicModel::project(const std::vector<double>& f) const {
    if(static_cast<int>(f.size())!=nf_) throw std::invalid_argument("full vector size");
    std::vector<double> r(static_cast<std::size_t>(nr_),0.0);
    for(int c=0;c<nr_;++c)for(int i=0;i<nf_;++i)r[static_cast<std::size_t>(c)]+=transform_[static_cast<std::size_t>(c*nf_+i)]*f[static_cast<std::size_t>(i)];
    return r;
}
SparseMatrixCSC ReducedDynamicModel::project_matrix(const SparseMatrixCSC& full) const{return dense_to_csc(nr_,nr_,project_dense(dense_from_csc(full),nf_,transform_,nr_),1e-13);}
std::vector<double> ReducedDynamicModel::mass_multiply(const std::vector<double>& a) const{return mass_matrix_.multiply(a);}
std::vector<double> ReducedDynamicModel::damping_multiply(const std::vector<double>& v) const{return damping_matrix_.multiply(v);}
SparseMatrixCSC ReducedDynamicModel::effective_initial_matrix(double a0,double a1) const{return add(add(k_initial_,mass_matrix_,1.0,a0),damping_matrix_,1.0,a1);}
SparseMatrixCSC ReducedDynamicModel::effective_tangent_matrix(const std::vector<double>& tangents,double a0,double a1) const{
    if(static_cast<int>(tangents.size())!=nonlinear_count()) throw std::invalid_argument("reduced tangent size");
    std::vector<Triplet> t;
    auto base=effective_initial_matrix(a0,a1);for(int c=0;c<base.cols();++c)for(int p=base.col_ptr()[static_cast<std::size_t>(c)];p<base.col_ptr()[static_cast<std::size_t>(c+1)];++p)t.push_back({base.row_ind()[static_cast<std::size_t>(p)],c,base.values()[static_cast<std::size_t>(p)]});
    const auto& init=initial_nonlinear_tangents();for(int j=0;j<nonlinear_count();++j){const double dk=tangents[static_cast<std::size_t>(j)]-init[static_cast<std::size_t>(j)];if(std::abs(dk)<1e-18)continue;const int b=basis_.col_ptr()[static_cast<std::size_t>(j)],e=basis_.col_ptr()[static_cast<std::size_t>(j+1)];for(int pi=b;pi<e;++pi)for(int pj=b;pj<e;++pj)t.push_back({basis_.row_ind()[static_cast<std::size_t>(pi)],basis_.row_ind()[static_cast<std::size_t>(pj)],dk*basis_.values()[static_cast<std::size_t>(pi)]*basis_.values()[static_cast<std::size_t>(pj)]});}
    return SparseMatrixCSC::from_triplets(nr_,nr_,t,-1.0);
}
void ReducedDynamicModel::evaluate_nonlinear_deformations(const std::vector<double>& deformations,const std::vector<double>& committed,std::vector<double>& component_forces,std::vector<double>& tangents,std::vector<double>& trial) const{parent_->evaluate_nonlinear_deformations(deformations,committed,component_forces,tangents,trial);}
void ReducedDynamicModel::internal_force_and_tangent(const std::vector<double>& q,const std::vector<double>& committed,std::vector<double>& force,std::vector<double>& tangents,std::vector<double>& trial) const{if(static_cast<int>(q.size())!=nr_)throw std::invalid_argument("reduced state size");force=k_linear_.multiply(q);std::vector<double> deformations(static_cast<std::size_t>(nonlinear_count())),cf;for(int j=0;j<nonlinear_count();++j)deformations[static_cast<std::size_t>(j)]=basis_.column_dot(j,q);evaluate_nonlinear_deformations(deformations,committed,cf,tangents,trial);for(int j=0;j<nonlinear_count();++j)basis_.axpy_column(j,cf[static_cast<std::size_t>(j)],force);}
void ReducedDynamicModel::internal_force_and_tangent_diagnostics(const std::vector<double>& q,const std::vector<double>& committed,std::vector<double>& force,std::vector<double>& tangents,std::vector<double>& trial,NonlinearEvalDiagnostics* diagnostics) const{
    // The generic reduced wrapper has no material-specific fast path; retain
    // exact reduced constitutive evaluation and report tangent activity.
    internal_force_and_tangent(q,committed,force,tangents,trial);
    if(!diagnostics) return;
    diagnostics->component_evaluations += static_cast<std::size_t>(nonlinear_count());
    const auto& init=initial_nonlinear_tangents();for(int j=0;j<nonlinear_count();++j){const bool a=std::abs(tangents[static_cast<std::size_t>(j)]-init[static_cast<std::size_t>(j)])>1e-14;diagnostics->active_tangent_evaluations+=a;diagnostics->fast_path_evaluations+=!a;diagnostics->full_state_evaluations+=a;}
}
std::vector<double> ReducedDynamicModel::base_excitation(double ag) const{return project(parent_->base_excitation(ag));}
double ReducedDynamicModel::response_value(const std::vector<double>& q) const{return parent_->response_value(expand(q));}
double ReducedDynamicModel::max_drift_measure(const std::vector<double>& q) const{return parent_->max_drift_measure(expand(q));}

std::vector<int> nonlinear_support_dofs(const NonlinearDynamicModel& model){std::vector<int> d;const auto& B=model.nonlinear_basis();for(int p:B.row_ind())d.push_back(p);std::sort(d.begin(),d.end());d.erase(std::unique(d.begin(),d.end()),d.end());return d;}

ReducedDynamicModel craig_bampton_reduce(const NonlinearDynamicModel& model,std::vector<int> retained,int fixed_modes,ReductionBuildInfo* info){
    const int n=model.dof();if(n<=0)throw std::invalid_argument("empty model");if(fixed_modes<0)throw std::invalid_argument("negative mode count");std::sort(retained.begin(),retained.end());retained.erase(std::unique(retained.begin(),retained.end()),retained.end());for(int d:retained)if(d<0||d>=n)throw std::out_of_range("retained DOF");if(retained.empty())throw std::invalid_argument("Craig-Bampton requires retained DOFs");
    std::vector<char> keep(static_cast<std::size_t>(n),0);for(int d:retained)keep[static_cast<std::size_t>(d)]=1;std::vector<int> interior;for(int i=0;i<n;++i)if(!keep[static_cast<std::size_t>(i)])interior.push_back(i);
    const int nb=static_cast<int>(retained.size()),ni=static_cast<int>(interior.size());
    if(ni==0){std::vector<double>T(static_cast<std::size_t>(n*n),0.0);for(int i=0;i<n;++i)T[static_cast<std::size_t>(i*n+i)]=1.0;if(info)*info={n,nb,0,n};return ReducedDynamicModel(model,std::move(T),n);}
    std::vector<int> imap(static_cast<std::size_t>(n),-1),bmap(static_cast<std::size_t>(n),-1);for(int i=0;i<ni;++i)imap[static_cast<std::size_t>(interior[static_cast<std::size_t>(i)])]=i;for(int i=0;i<nb;++i)bmap[static_cast<std::size_t>(retained[static_cast<std::size_t>(i)])]=i;
    std::vector<Triplet> kii_t;std::vector<double> kib(static_cast<std::size_t>(ni*nb),0.0);const auto& K=model.K_initial();
    for(int c=0;c<n;++c)for(int p=K.col_ptr()[static_cast<std::size_t>(c)];p<K.col_ptr()[static_cast<std::size_t>(c+1)];++p){const int r=K.row_ind()[static_cast<std::size_t>(p)];const double v=K.values()[static_cast<std::size_t>(p)];const int ic=imap[static_cast<std::size_t>(c)],ir=imap[static_cast<std::size_t>(r)],bc=bmap[static_cast<std::size_t>(c)];if(ic>=0&&ir>=0)kii_t.push_back({ir,ic,v});else if(ir>=0&&bc>=0)kib[static_cast<std::size_t>(bc*ni+ir)]+=v;}
    auto Kii=SparseMatrixCSC::from_triplets(ni,ni,kii_t,1e-18);SuperLUFactor f(Kii);for(double& v:kib)v=-v;auto psi=f.solve_multiple(kib,nb);

    auto Mfull=dense_mass(model);std::vector<double>Mii(static_cast<std::size_t>(ni*ni),0.0);for(int c=0;c<ni;++c)for(int r=0;r<ni;++r)Mii[static_cast<std::size_t>(c*ni+r)]=Mfull[static_cast<std::size_t>(interior[static_cast<std::size_t>(c)]*n+interior[static_cast<std::size_t>(r)])];
    auto phi=finite_fixed_interface_modes(Kii,Mii,std::min(fixed_modes,ni));const int nm=ni?static_cast<int>(phi.size()/static_cast<std::size_t>(ni)):0,nr=nb+nm;std::vector<double>T(static_cast<std::size_t>(n*nr),0.0);
    for(int b=0;b<nb;++b){T[static_cast<std::size_t>(b*n+retained[static_cast<std::size_t>(b)])]=1.0;for(int i=0;i<ni;++i)T[static_cast<std::size_t>(b*n+interior[static_cast<std::size_t>(i)])]=psi[static_cast<std::size_t>(b*ni+i)];}
    for(int m=0;m<nm;++m)for(int i=0;i<ni;++i)T[static_cast<std::size_t>((nb+m)*n+interior[static_cast<std::size_t>(i)])]=phi[static_cast<std::size_t>(m*ni+i)];
    if(info) *info={n,nb,nm,nr};
    return ReducedDynamicModel(model,std::move(T),nr);
}

} // namespace quake
