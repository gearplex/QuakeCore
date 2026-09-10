#include "quake/steel2d.hpp"
#include "quake/dynamic_model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace quake {
namespace {

double dot6(const std::array<double,6>& a,const std::array<double,6>& b){
    double r=0.0;for(int i=0;i<6;++i)r+=a[i]*b[i];return r;
}
void outer6(std::array<double,36>& K,const std::array<double,6>& a,
            const std::array<double,6>& b,double scale){
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)K[6*i+j]+=scale*a[i]*b[j];
}
std::array<double,4> inverse2(const std::array<double,4>& a){
    const double scale=std::max({1.0,std::abs(a[0]),std::abs(a[1]),std::abs(a[2]),std::abs(a[3])});
    const double det=a[0]*a[3]-a[1]*a[2];
    if(!std::isfinite(det)||std::abs(det)<=1e-14*scale*scale)
        throw ConstitutiveIntegrationError("steel member local hinge Jacobian is singular");
    return {a[3]/det,-a[1]/det,-a[2]/det,a[0]/det};
}
std::array<double,2> multiply2(const std::array<double,4>& a,
                               const std::array<double,2>& x){
    return {a[0]*x[0]+a[1]*x[1],a[2]*x[0]+a[3]*x[1]};
}
std::array<double,4> matmul2(const std::array<double,4>& a,
                             const std::array<double,4>& b){
    return {a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]};
}
double norm2(const std::array<double,2>& a){return std::max(std::abs(a[0]),std::abs(a[1]));}

std::array<double,36> geometric_stiffness(double L,double c,double s,double compression){
    std::array<double,36> out{};
    if(compression==0.0)return out;
    const double q=-compression/(30.0*L); // compression subtracts stiffness
    const double g[4][4]={{36.0,3.0*L,-36.0,3.0*L},
                          {3.0*L,4.0*L*L,-3.0*L,-L*L},
                          {-36.0,-3.0*L,36.0,-3.0*L},
                          {3.0*L,-L*L,-3.0*L,4.0*L*L}};
    // q_local = T q_global; retain [v_i,rz_i,v_j,rz_j].
    const double T[4][6]={{-s,c,0,0,0,0},{0,0,1,0,0,0},
                          {0,0,0,-s,c,0},{0,0,0,0,0,1}};
    for(int a=0;a<6;++a)for(int b=0;b<6;++b)
        for(int i=0;i<4;++i)for(int j=0;j<4;++j)
            out[6*a+b]+=q*T[i][a]*g[i][j]*T[j][b];
    return out;
}

} // namespace

SteelMember2D::SteelMember2D(double xi,double yi,double xj,double yj,
                             SteelMember2DProperties p):properties_(std::move(p)){
    const double dx=xj-xi,dy=yj-yi;length_=std::hypot(dx,dy);
    if(!std::isfinite(length_)||length_<=0.0||!std::isfinite(properties_.E)||properties_.E<=0.0||
       !std::isfinite(properties_.A)||properties_.A<=0.0||!std::isfinite(properties_.I)||properties_.I<=0.0||
       !std::isfinite(properties_.axial_compression)||properties_.axial_compression<0.0)
        throw std::invalid_argument("invalid steel member geometry or section property");
    if(properties_.local_max_iterations<1||properties_.local_max_iterations>1000||
       !std::isfinite(properties_.local_relative_tolerance)||properties_.local_relative_tolerance<=0.0||
       properties_.local_relative_tolerance>1e-3)
        throw std::invalid_argument("invalid steel member local solver controls");
    c_=dx/length_;s_=dy/length_;
    axial_B_={-c_,-s_,0.0,c_,s_,0.0};
    const std::array<double,6> chord_B={-s_/length_,c_/length_,0.0,
                                        s_/length_,-c_/length_,0.0};
    for(int j=0;j<6;++j){
        bending_B_[j]=chord_B[j];
        bending_B_[6+j]=chord_B[j];
    }
    bending_B_[2]+=1.0;bending_B_[11]+=1.0;
    hinge_j_offset_=hinge_i_offset_+properties_.hinge_i.state_size();
    state_size_=hinge_j_offset_+properties_.hinge_j.state_size();
    geometric_K_=geometric_stiffness(length_,c_,s_,properties_.axial_compression);
}

std::array<double,4> SteelMember2D::elastic_bending_stiffness() const{
    const double f=properties_.E*properties_.I/length_;
    return {4.0*f,2.0*f,2.0*f,4.0*f};
}

std::array<double,4> SteelMember2D::condensed_bending_stiffness(double ki,double kj) const{
    const auto ke=elastic_bending_stiffness();
    const std::array<double,4> a={ke[0]+ki,ke[1],ke[2],ke[3]+kj};
    const auto correction=matmul2(matmul2(ke,inverse2(a)),ke);
    return {ke[0]-correction[0],ke[1]-correction[1],
            ke[2]-correction[2],ke[3]-correction[3]};
}

std::vector<double> SteelMember2D::initial_state() const{
    std::vector<double> s(static_cast<std::size_t>(state_size_),0.0);
    properties_.hinge_i.initialize_state(s.data()+hinge_i_offset_);
    properties_.hinge_j.initialize_state(s.data()+hinge_j_offset_);
    return s;
}

std::array<double,36> SteelMember2D::initial_tangent() const{
    std::array<double,36> K=geometric_K_;
    outer6(K,axial_B_,axial_B_,properties_.E*properties_.A/length_);
    const auto kb=condensed_bending_stiffness(properties_.hinge_i.initial_stiffness(),
                                               properties_.hinge_j.initial_stiffness());
    std::array<double,6> bi{},bj{};
    for(int k=0;k<6;++k){bi[k]=bending_B_[k];bj[k]=bending_B_[6+k];}
    outer6(K,bi,bi,kb[0]);outer6(K,bi,bj,kb[1]);
    outer6(K,bj,bi,kb[2]);outer6(K,bj,bj,kb[3]);
    return K;
}

SteelMember2DResponse SteelMember2D::trial(const std::array<double,6>& u,
                                            const double* committed) const{
    for(double x:u)if(!std::isfinite(x))throw std::invalid_argument("nonfinite steel member displacement");
    SteelMember2DResponse out;out.state.resize(static_cast<std::size_t>(state_size_));
    out.axial_deformation=dot6(axial_B_,u);
    out.axial_force=properties_.E*properties_.A/length_*out.axial_deformation;
    const double vi=-s_*u[0]+c_*u[1],vj=-s_*u[3]+c_*u[4];
    out.chord_rotation=(vj-vi)/length_;
    for(int row=0;row<2;++row){
        double q=0.0;for(int j=0;j<6;++j)q+=bending_B_[6*row+j]*u[j];
        out.end_rotation[row]=q;
        out.hinge_rotation[row]=committed[row];
    }
    const auto ke=elastic_bending_stiffness();
    MaterialTrialResult hi{},hj{};
    std::vector<double> si(static_cast<std::size_t>(properties_.hinge_i.state_size()));
    std::vector<double> sj(static_cast<std::size_t>(properties_.hinge_j.state_size()));
    auto evaluate=[&](const std::array<double,2>& p,MaterialTrialResult& ti,MaterialTrialResult& tj,
                      std::vector<double>& tsi,std::vector<double>& tsj){
        ti=properties_.hinge_i.trial(p[0],committed+hinge_i_offset_,tsi.data());
        tj=properties_.hinge_j.trial(p[1],committed+hinge_j_offset_,tsj.data());
        const std::array<double,2> elastic={out.end_rotation[0]-p[0],out.end_rotation[1]-p[1]};
        const auto me=multiply2(ke,elastic);
        return std::array<double,2>{me[0]-ti.force,me[1]-tj.force};
    };
    bool converged=false;
    auto residual=evaluate(out.hinge_rotation,hi,hj,si,sj);
    for(int iter=0;iter<properties_.local_max_iterations;++iter){
        ++out.local_iterations;
        const double scale=std::max({1.0,std::abs(hi.force),std::abs(hj.force)});
        if(norm2(residual)<=properties_.local_relative_tolerance*scale){converged=true;break;}
        const std::array<double,4> jac={ke[0]+hi.tangent,ke[1],ke[2],ke[3]+hj.tangent};
        const auto dp=multiply2(inverse2(jac),residual);
        double alpha=1.0;bool accepted=false;
        for(int ls=0;ls<18;++ls){
            const std::array<double,2> candidate={out.hinge_rotation[0]+alpha*dp[0],
                                                  out.hinge_rotation[1]+alpha*dp[1]};
            MaterialTrialResult ci{},cj{};std::vector<double> csi(si.size()),csj(sj.size());
            const auto cr=evaluate(candidate,ci,cj,csi,csj);
            if(norm2(cr)<norm2(residual)){
                out.hinge_rotation=candidate;hi=ci;hj=cj;si=std::move(csi);sj=std::move(csj);
                residual=cr;accepted=true;break;
            }
            alpha*=0.5;
        }
        if(!accepted)break;
    }
    if(!converged){
        const double scale=std::max({1.0,std::abs(hi.force),std::abs(hj.force)});
        if(norm2(residual)>properties_.local_relative_tolerance*scale)
            throw ConstitutiveIntegrationError("steel member end-hinge equilibrium failed");
    }
    out.state[0]=out.hinge_rotation[0];out.state[1]=out.hinge_rotation[1];
    std::copy(si.begin(),si.end(),out.state.begin()+hinge_i_offset_);
    std::copy(sj.begin(),sj.end(),out.state.begin()+hinge_j_offset_);
    out.hinge_diagnostics={hi.diagnostics,hj.diagnostics};
    const std::array<double,2> elastic={out.end_rotation[0]-out.hinge_rotation[0],
                                        out.end_rotation[1]-out.hinge_rotation[1]};
    out.end_moment=multiply2(ke,elastic);
    for(int a=0;a<6;++a){
        out.force[a]=axial_B_[a]*out.axial_force;
        for(int b=0;b<6;++b)out.force[a]+=geometric_K_[6*a+b]*u[b];
        out.force[a]+=bending_B_[a]*out.end_moment[0]+bending_B_[6+a]*out.end_moment[1];
    }
    out.tangent=geometric_K_;
    outer6(out.tangent,axial_B_,axial_B_,properties_.E*properties_.A/length_);
    const auto kb=condensed_bending_stiffness(hi.tangent,hj.tangent);
    std::array<double,6> bi{},bj{};for(int k=0;k<6;++k){bi[k]=bending_B_[k];bj[k]=bending_B_[6+k];}
    outer6(out.tangent,bi,bi,kb[0]);outer6(out.tangent,bi,bj,kb[1]);
    outer6(out.tangent,bj,bi,kb[2]);outer6(out.tangent,bj,bj,kb[3]);
    return out;
}

ViscousDamper2D::ViscousDamper2D(ViscousDamper2DProperties p):properties_(p){
    if(!std::isfinite(p.coefficient)||p.coefficient<=0.0||!std::isfinite(p.alpha)||p.alpha<=0.0||p.alpha>2.0||
       !std::isfinite(p.regularization_velocity)||p.regularization_velocity<0.0||
       (p.alpha<1.0&&p.regularization_velocity<=0.0))
        throw std::invalid_argument("invalid viscous damper coefficient, exponent, or regularization velocity");
}

ViscousDamper2DTrial ViscousDamper2D::trial(double v) const{
    if(!std::isfinite(v))throw std::invalid_argument("nonfinite viscous damper velocity");
    const auto& p=properties_;ViscousDamper2DTrial r;r.deformation_rate=v;
    if(p.regularization_velocity>0.0){
        const double z=v*v+p.regularization_velocity*p.regularization_velocity;
        const double power=std::pow(z,0.5*(p.alpha-1.0));
        r.force=p.coefficient*v*power;
        r.tangent=p.coefficient*power*(1.0+(p.alpha-1.0)*v*v/z);
    }else if(v==0.0){
        r.force=0.0;r.tangent=p.alpha==1.0?p.coefficient:0.0;
    }else{
        const double a=std::abs(v);
        r.force=p.coefficient*std::copysign(std::pow(a,p.alpha),v);
        r.tangent=p.coefficient*p.alpha*std::pow(a,p.alpha-1.0);
    }
    if(!std::isfinite(r.force)||!std::isfinite(r.tangent)||r.tangent<0.0)
        throw ConstitutiveIntegrationError("nonfinite viscous damper response");
    return r;
}

} // namespace quake
