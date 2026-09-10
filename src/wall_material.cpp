#include "quake/wall_material.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {
void positive(double x) { if(!std::isfinite(x)||x<=0)throw std::invalid_argument("wall material modulus/strength must be finite and positive"); }
}
WallUniaxial WallUniaxial::elastic(double E) { positive(E);WallUniaxial m;m.E_=E;return m; }
WallUniaxial WallUniaxial::steel(double E,double fy,double b) {
    positive(E);positive(fy);if(!std::isfinite(b)||b<0||b>=1)throw std::invalid_argument("wall steel hardening ratio must be in [0,1)");
    WallUniaxial m;m.kind_=Kind::SteelBilinear;m.E_=E;m.fy_=fy;m.b_=b;return m;
}
WallUniaxial WallUniaxial::concrete01(double fc,double ec,double fu,double eu) {
    if(!std::isfinite(fc)||!std::isfinite(ec)||!std::isfinite(fu)||!std::isfinite(eu)||fc>=0||ec>=0||fu>0||fu<fc||eu>=ec)
        throw std::invalid_argument("concrete01 requires fc<0, epsc<0, fc<=fcu<=0, epsu<epsc");
    WallUniaxial m;m.kind_=Kind::Concrete01;m.fc_=fc;m.ec_=ec;m.fu_=fu;m.eu_=eu;m.E_=2*fc/ec;positive(m.E_);return m;
}
void WallUniaxial::initialize(double* s) const { std::fill(s,s+6,0.0);if(kind_==Kind::Concrete01)s[2]=E_; }
WallUniaxial::Result WallUniaxial::trial(double e,const double* s,double* t) const {
    if(!std::isfinite(e))throw std::invalid_argument("nonfinite wall material strain");
    for(int i=0;i<6;++i)if(!std::isfinite(s[i]))throw std::invalid_argument("nonfinite wall material state");
    std::copy(s,s+6,t);
    if(kind_==Kind::Elastic)return {E_*e,E_};
    if(kind_==Kind::SteelBilinear){auto r=BilinearSpring(E_,fy_,b_).trial(e,{s[0],s[1]});t[0]=r.state.plastic;t[1]=r.state.backstress;return {r.force,r.tangent};}
    // Concrete01 history: minimum strain, zero-stress intercept, unloading E,
    // last strain, last stress, last tangent. No tension; residual compression
    // plateau is explicitly specified by fcu (zero is permitted).
    double stress=0,k=0;
    if(e==s[3])return {s[4],s[5]==0 && e==0 && s[0]==0 ? E_ : s[5]};
    if(e<0){
        if(e<=s[0]){
            if(e>ec_){double r=e/ec_;stress=fc_*r*(2-r);k=E_*(1-r);}
            else if(e>eu_){k=(fu_-fc_)/(eu_-ec_);stress=fc_+k*(e-ec_);}
            else {stress=fu_;k=0;}
            t[0]=e;
            double eta=std::max(e,eu_)/ec_;
            double ratio=eta<2 ? .145*eta*eta+.13*eta : .707*(eta-2)+.834;
            double zero=ratio*ec_;
            double span=e-zero;
            double elastic_span=stress/E_;
            if(span<elastic_span){t[1]=zero;t[2]=stress/span;}
            else {t[1]=e-elastic_span;t[2]=E_;}
        }else if(e<s[1]){stress=s[2]*(e-s[1]);k=s[2];}
    }
    // Do not cross the unloading/reloading straight line from the last point
    // when a descending envelope is reached during a reversal.
    double line=s[4]+s[2]*(e-s[3]);
    if(e<s[3] && line>stress){stress=line;k=s[2];}
    t[3]=e;t[4]=stress;t[5]=k;
    return {stress,k};
}
WallPanel::WallPanel(double E,double nu,std::vector<WallPanelLayer> layers):layers_(std::move(layers)) {
    if(!std::isfinite(E)||E<0||!std::isfinite(nu)||nu<=-1||nu>=.5)throw std::invalid_argument("invalid wall panel elastic background");
    if(E==0 && layers_.empty())throw std::invalid_argument("empty wall panel");
    double d=E/(1-nu*nu);background_={d,nu*d,0,nu*d,d,0,0,0,E/(2*(1+nu))};
    for(const auto& l:layers_){
        if(!std::isfinite(l.angle_rad)||!std::isfinite(l.weight)||l.weight<=0)throw std::invalid_argument("invalid wall panel layer");
        double c=std::cos(l.angle_rad),s=std::sin(l.angle_rad);directions_.push_back({c*c,s*s,c*s});
    }
}
std::vector<double> WallPanel::initial_state() const {
    std::vector<double> s(state_size());for(std::size_t i=0;i<layers_.size();++i)layers_[i].material.initialize(s.data()+6*i);return s;
}
std::array<double,9> WallPanel::initial_tangent() const {
    auto D=background_;for(std::size_t i=0;i<layers_.size();++i)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        D[3*a+b]+=layers_[i].weight*layers_[i].material.initial_tangent()*directions_[i][a]*directions_[i][b];
    return D;
}
WallPanelResult WallPanel::trial(const std::array<double,3>& e,const double* s) const {
    for(double v:e)if(!std::isfinite(v))throw std::invalid_argument("nonfinite panel strain");
    WallPanelResult r;r.tangent=background_;r.state.resize(state_size());
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)r.stress[a]+=background_[3*a+b]*e[b];
    for(std::size_t i=0;i<layers_.size();++i){
        const auto& n=directions_[i];double q=n[0]*e[0]+n[1]*e[1]+n[2]*e[2];
        auto v=layers_[i].material.trial(q,s+6*i,r.state.data()+6*i);
        for(int a=0;a<3;++a){r.stress[a]+=layers_[i].weight*v.stress*n[a];for(int b=0;b<3;++b)r.tangent[3*a+b]+=layers_[i].weight*v.tangent*n[a]*n[b];}
    }
    return r;
}
} // namespace quake
