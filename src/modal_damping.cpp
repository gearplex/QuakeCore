#include "quake/modal_damping.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {
double dot(const std::vector<double>& a,const std::vector<double>& b){
    if(a.size()!=b.size()) throw std::invalid_argument("modal damping dot size");
    double s=0.0;for(std::size_t i=0;i<a.size();++i)s+=a[i]*b[i];return s;
}
}

FixedModalDampingModel::FixedModalDampingModel(const NonlinearDynamicModel& base,double zeta,int requested_modes,bool exact_low_rank_newton_tangent)
    : FixedModalDampingModel(base,std::vector<double>{zeta},requested_modes,exact_low_rank_newton_tangent) {}

FixedModalDampingModel::FixedModalDampingModel(const NonlinearDynamicModel& base,std::vector<double> zeta,int requested_modes,bool exact_low_rank_newton_tangent)
    : base_(&base),exact_low_rank_newton_tangent_(exact_low_rank_newton_tangent) {
    if(requested_modes<=0||zeta.empty()) throw std::invalid_argument("invalid modal damping setup");
    for(double x:zeta) if(x<0.0||!std::isfinite(x)) throw std::invalid_argument("invalid modal damping ratio");
    raw_modes_=modal_analysis(base,requested_modes);
    if(raw_modes_.empty()) throw std::invalid_argument("modal damping requires finite positive modes");
    zeta_.resize(raw_modes_.size());
    for(std::size_t i=0;i<zeta_.size();++i) zeta_[i]=zeta.size()==1?zeta[0]:(i<zeta.size()?zeta[i]:zeta.back());
    modes_.reserve(raw_modes_.size());
    for(const auto& md:raw_modes_){
        auto phi=md.shape;auto mphi=base.mass_multiply(phi);const double mm=dot(phi,mphi);
        if(!(mm>0.0)||!std::isfinite(mm)) continue;
        const double inv=1.0/std::sqrt(mm);for(double& x:phi)x*=inv;for(double& x:mphi)x*=inv;
        modes_.push_back({md.omega,std::move(phi),std::move(mphi)});
    }
    if(modes_.empty()) throw std::invalid_argument("modal damping has no mass-normalizable modes");
    if(modes_.size()!=zeta_.size()) zeta_.resize(modes_.size());
    std::vector<double> dense(static_cast<std::size_t>(dof()*static_cast<int>(modes_.size())),0.0);
    for(std::size_t c=0;c<modes_.size();++c)
        for(int r=0;r<dof();++r)
            dense[static_cast<std::size_t>(c*dof()+r)]=modes_[c].mphi[static_cast<std::size_t>(r)];
    modal_tangent_basis_=SparseUpdateBasis::from_dense(dof(),static_cast<int>(modes_.size()),dense,1e-15);
}

std::vector<double> FixedModalDampingModel::additional_effective_low_rank_coefficients(double,double a1) const{
    const int r=static_cast<int>(modes_.size());
    std::vector<double> C(static_cast<std::size_t>(r*r),0.0);
    for(int i=0;i<r;++i) C[static_cast<std::size_t>(i*r+i)]=a1*2.0*zeta_[static_cast<std::size_t>(i)]*modes_[static_cast<std::size_t>(i)].omega;
    return C;
}

std::vector<double> FixedModalDampingModel::damping_multiply(const std::vector<double>& v) const{
    if(static_cast<int>(v.size())!=dof()) throw std::invalid_argument("modal damping velocity size");
    auto out=base_->damping_multiply(v);
    for(std::size_t i=0;i<modes_.size();++i){
        const auto& m=modes_[i];const double qdot=dot(m.mphi,v);const double c=2.0*zeta_[i]*m.omega*qdot;
        for(std::size_t j=0;j<out.size();++j)out[j]+=c*m.mphi[j];
    }
    return out;
}

} // namespace quake
