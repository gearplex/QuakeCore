#include "quake/wall2d.hpp"
#include "quake/dynamic_model.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {
void geometry(double h,double c,double rho,std::size_t count){
    if(!std::isfinite(h)||h<=0||!std::isfinite(c)||c<0||c>1||!std::isfinite(rho)||rho<0||count<2)
        throw std::invalid_argument("wall requires positive height, >=2 fibers, c in [0,1], and nonnegative density");
}
void dimensions(double b,double t){if(!std::isfinite(b)||b<=0||!std::isfinite(t)||t<=0)throw std::invalid_argument("wall fiber width/thickness must be positive");}
double dot(const std::array<double,6>& a,const std::array<double,6>& b){double r=0;for(int i=0;i<6;++i)r+=a[i]*b[i];return r;}
void outer(std::array<double,36>& K,const std::array<double,6>& a,const std::array<double,6>& b,double s){for(int i=0;i<6;++i)for(int j=0;j<6;++j)K[6*i+j]+=s*a[i]*b[j];}
std::array<double,4> condense(const std::array<double,9>& D,double ref){
    if(!std::isfinite(D[0])||std::abs(D[0])<=1e-14*ref)throw ConstitutiveIntegrationError("SFI-MVLEM singular transverse panel tangent");
    return {D[4]-D[3]*D[1]/D[0],D[5]-D[3]*D[2]/D[0],D[7]-D[6]*D[1]/D[0],D[8]-D[6]*D[2]/D[0]};
}
void panel_stiffness(std::array<double,36>& K,const std::array<double,6>& y,const std::array<double,6>& g,const std::array<double,4>& d,double v){
    outer(K,y,y,v*d[0]);outer(K,y,g,v*d[1]);outer(K,g,y,v*d[2]);outer(K,g,g,v*d[3]);
}
}
Wall2D::Wall2D(double h,MVLEMProperties p):h_(h),c_(p.c),properties_(std::move(p)) {
    const auto& a=std::get<MVLEMProperties>(properties_);geometry(h,c_,a.density,a.fibers.size());
    double width=0;for(const auto& f:a.fibers){dimensions(f.width,f.thickness);if(!std::isfinite(f.reinforcement_ratio)||f.reinforcement_ratio<0||f.reinforcement_ratio>=1)throw std::invalid_argument("MVLEM reinforcement ratio must be in [0,1)");x_.push_back(width+f.width/2);width+=f.width;total_mass_+=h*f.width*f.thickness*a.density;offsets_.push_back(state_size_);state_size_+=12;}
    for(auto& x:x_)x-=width/2;
    state_size_+=6;
}
Wall2D::Wall2D(double h,SFIMVLEMProperties p):h_(h),c_(p.c),properties_(std::move(p)) {
    const auto& a=std::get<SFIMVLEMProperties>(properties_);geometry(h,c_,a.density,a.panels.size());
    if(a.local_max_iterations<1||a.local_max_iterations>1000||!std::isfinite(a.local_relative_tolerance)||a.local_relative_tolerance<=0||a.local_relative_tolerance>1e-3)throw std::invalid_argument("invalid SFI panel solver controls");
    double width=0;for(const auto& f:a.panels){dimensions(f.width,f.thickness);x_.push_back(width+f.width/2);width+=f.width;total_mass_+=h*f.width*f.thickness*a.density;offsets_.push_back(state_size_);state_size_+=1+f.material.state_size();auto D=f.material.initial_tangent();if(D[0]<=0)throw std::invalid_argument("SFI panel requires positive initial transverse stiffness");}
    for(auto& x:x_)x-=width/2;
}
std::array<double,6> Wall2D::shear_B() const {return {-1/h_,0,c_,1/h_,0,1-c_};}
std::array<double,6> Wall2D::axial_B(int i) const {double x=x_.at(i);return {0,-1/h_,-x/h_,0,1/h_,x/h_};}
std::vector<double> Wall2D::initial_state() const {
    std::vector<double> s(state_size_,0);
    if(!is_sfi()){
        const auto& p=std::get<MVLEMProperties>(properties_);for(std::size_t i=0;i<p.fibers.size();++i){p.fibers[i].concrete.initialize(s.data()+offsets_[i]);p.fibers[i].steel.initialize(s.data()+offsets_[i]+6);}p.shear.initialize(s.data()+state_size_-6);
    }else{const auto& p=std::get<SFIMVLEMProperties>(properties_);for(std::size_t i=0;i<p.panels.size();++i){auto t=p.panels[i].material.initial_state();std::copy(t.begin(),t.end(),s.begin()+offsets_[i]+1);}}
    return s;
}
std::array<double,36> Wall2D::initial_tangent() const {
    std::array<double,36> K{};auto g=shear_B();
    if(!is_sfi()){
        const auto& p=std::get<MVLEMProperties>(properties_);
        for(std::size_t i=0;i<p.fibers.size();++i){const auto& f=p.fibers[i];auto y=axial_B(i);double E=(1-f.reinforcement_ratio)*f.concrete.initial_tangent()+f.reinforcement_ratio*f.steel.initial_tangent();outer(K,y,y,h_*f.width*f.thickness*E);}
        outer(K,g,g,h_*h_*p.shear.initial_tangent());
    }else{
        const auto& p=std::get<SFIMVLEMProperties>(properties_);for(std::size_t i=0;i<p.panels.size();++i){const auto& f=p.panels[i];auto D=f.material.initial_tangent();panel_stiffness(K,axial_B(i),g,condense(D,D[0]),h_*f.width*f.thickness);}
    }
    return K;
}
Wall2DResponse Wall2D::trial(const std::array<double,6>& u,const double* s) const {
    for(double v:u)if(!std::isfinite(v))throw std::invalid_argument("nonfinite wall displacement");
    Wall2DResponse r;r.state.resize(state_size_);auto g=shear_B();double gamma=dot(g,u);r.shear_deformation=h_*gamma;r.curvature=(u[5]-u[2])/h_;
    if(!is_sfi()){
        const auto& p=std::get<MVLEMProperties>(properties_);
        for(std::size_t i=0;i<p.fibers.size();++i){
            const auto& f=p.fibers[i];auto y=axial_B(i);double e=dot(y,u);int o=offsets_[i];
            auto a=f.concrete.trial(e,s+o,r.state.data()+o),b=f.steel.trial(e,s+o+6,r.state.data()+o+6);
            double A=f.width*f.thickness,V=h_*A,rho=f.reinforcement_ratio;
            for(int j=0;j<6;++j)r.force[j]+=V*((1-rho)*a.stress+rho*b.stress)*y[j];
            outer(r.tangent,y,y,V*((1-rho)*a.tangent+rho*b.tangent));
            r.fiber_strain.push_back(e);r.concrete_stress.push_back(a.stress);r.steel_stress.push_back(b.stress);
        }
        auto sh=p.shear.trial(r.shear_deformation,s+state_size_-6,r.state.data()+state_size_-6);
        for(int j=0;j<6;++j)r.force[j]+=h_*sh.stress*g[j];
        outer(r.tangent,g,g,h_*h_*sh.tangent);
    }else{
        const auto& p=std::get<SFIMVLEMProperties>(properties_);
        for(std::size_t i=0;i<p.panels.size();++i){
            const auto& f=p.panels[i];auto y=axial_B(i);double ey=dot(y,u);int o=offsets_[i];double ex=s[o];
            if(!std::isfinite(ex))throw std::invalid_argument("nonfinite transverse panel history");
            const auto D0=f.material.initial_tangent();double ref=D0[0];
            double tol=p.local_relative_tolerance*ref*std::max({1e-6,std::abs(ey),std::abs(gamma)});
            auto v=f.material.trial({ex,ey,gamma},s+o+1);bool done=false;
            for(int it=0;it<p.local_max_iterations;++it){
                ++r.local_iterations;
                if(!std::isfinite(v.stress[0]))throw ConstitutiveIntegrationError("SFI-MVLEM nonfinite transverse stress");
                if(std::abs(v.stress[0])<=tol){done=true;break;}
                (void)condense(v.tangent,ref);
                double dx=-v.stress[0]/v.tangent[0],alpha=1;bool accepted=false;
                for(int ls=0;ls<18;++ls){auto w=f.material.trial({ex+alpha*dx,ey,gamma},s+o+1);
                    if(std::isfinite(w.stress[0])&&(std::abs(w.stress[0])<=tol||std::abs(w.stress[0])<(1-1e-4*alpha)*std::abs(v.stress[0]))){ex+=alpha*dx;v=std::move(w);accepted=true;break;}alpha*=.5;}
                if(!accepted)break;
            }
            if(!done && std::abs(v.stress[0])>tol)throw ConstitutiveIntegrationError("SFI-MVLEM transverse equilibrium failed in panel "+std::to_string(i));
            const auto d=condense(v.tangent,ref);double V=h_*f.width*f.thickness;
            for(int j=0;j<6;++j)r.force[j]+=V*(v.stress[1]*y[j]+v.stress[2]*g[j]);
            panel_stiffness(r.tangent,y,g,d,V);
            r.state[o]=ex;std::copy(v.state.begin(),v.state.end(),r.state.begin()+o+1);
            r.panel_strain.push_back({ex,ey,gamma});r.panel_stress.push_back(v.stress);
            r.maximum_transverse_stress_residual=std::max(r.maximum_transverse_stress_residual,std::abs(v.stress[0]));
        }
    }
    return r;
}
} // namespace quake
