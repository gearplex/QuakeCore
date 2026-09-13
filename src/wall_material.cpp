#include "quake/wall_material.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace quake {
namespace {
void positive(double x) { if(!std::isfinite(x)||x<=0)throw std::invalid_argument("wall material modulus/strength must be finite and positive"); }
constexpr int concrete_cm_state_size=28;
constexpr int pinching4_state_size=17;
void encode_concrete_cm_state(const ConcreteCMState& s,double* v){
    v[0]=static_cast<int>(s.rule);v[1]=s.strain;v[2]=s.stress;v[3]=s.tangent;v[4]=s.increment;
    v[5]=s.unloading_strain;v[6]=s.unloading_stress;v[7]=s.zero_stress_strain;v[8]=s.zero_stress_tangent;
    v[9]=s.tension_zero_strain;v[10]=s.tension_peak_strain;v[11]=s.tension_peak_stress;v[12]=s.tension_new_stress;v[13]=s.tension_new_tangent;
    v[14]=s.tension_rejoin_strain;v[15]=s.tension_rejoin_stress;v[16]=s.tension_rejoin_tangent;
    v[17]=s.has_positive_to_negative_reversal?1.0:0.0;v[18]=s.positive_reversal_strain;v[19]=s.positive_reversal_stress;
    v[20]=s.positive_zero_stress_strain;v[21]=s.positive_zero_stress_tangent;v[22]=s.compression_new_stress;v[23]=s.compression_new_tangent;
    v[24]=s.compression_rejoin_strain;v[25]=s.compression_rejoin_stress;v[26]=s.compression_rejoin_tangent;
    v[27]=s.has_second_negative_to_positive_reversal?1.0:0.0;
}
ConcreteCMState decode_concrete_cm_state(const double* v){
    ConcreteCMState s;s.rule=static_cast<ConcreteCMRule>(static_cast<int>(v[0]));s.strain=v[1];s.stress=v[2];s.tangent=v[3];s.increment=v[4];
    s.unloading_strain=v[5];s.unloading_stress=v[6];s.zero_stress_strain=v[7];s.zero_stress_tangent=v[8];
    s.tension_zero_strain=v[9];s.tension_peak_strain=v[10];s.tension_peak_stress=v[11];s.tension_new_stress=v[12];s.tension_new_tangent=v[13];
    s.tension_rejoin_strain=v[14];s.tension_rejoin_stress=v[15];s.tension_rejoin_tangent=v[16];
    s.has_positive_to_negative_reversal=v[17]!=0.0;s.positive_reversal_strain=v[18];s.positive_reversal_stress=v[19];
    s.positive_zero_stress_strain=v[20];s.positive_zero_stress_tangent=v[21];s.compression_new_stress=v[22];s.compression_new_tangent=v[23];
    s.compression_rejoin_strain=v[24];s.compression_rejoin_stress=v[25];s.compression_rejoin_tangent=v[26];
    s.has_second_negative_to_positive_reversal=v[27]!=0.0;return s;
}
void encode_pinching4_state(const Pinching4State& s,double* v){
    v[0]=static_cast<int>(s.kind);v[1]=s.deformation;v[2]=s.force;v[3]=s.tangent;v[4]=s.deformation_rate;
    v[5]=s.low_state_deformation;v[6]=s.low_state_force;v[7]=s.high_state_deformation;v[8]=s.high_state_force;
    v[9]=s.min_demand;v[10]=s.max_demand;v[11]=s.energy;v[12]=s.gamma_d;v[13]=s.cycles;v[14]=s.damaged_min;v[15]=s.damaged_max;
    v[16]=static_cast<double>(s.reversal_count);
}
Pinching4State decode_pinching4_state(const double* v){
    Pinching4State s;s.kind=static_cast<Pinching4StateKind>(static_cast<int>(v[0]));s.deformation=v[1];s.force=v[2];s.tangent=v[3];s.deformation_rate=v[4];
    s.low_state_deformation=v[5];s.low_state_force=v[6];s.high_state_deformation=v[7];s.high_state_force=v[8];
    s.min_demand=v[9];s.max_demand=v[10];s.energy=v[11];s.gamma_d=v[12];s.cycles=v[13];s.damaged_min=v[14];s.damaged_max=v[15];
    s.reversal_count=static_cast<unsigned>(v[16]);return s;
}
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
WallUniaxial WallUniaxial::concrete_cm(ConcreteCMParameters parameters) {
    ConcreteCM material(parameters);
    WallUniaxial m;m.kind_=Kind::ConcreteCM;m.concrete_cm_parameters_=std::make_shared<const ConcreteCMParameters>(std::move(parameters));
    m.E_=material.initial_state().tangent;return m;
}
WallUniaxial WallUniaxial::pinching4(Pinching4CyclicParameters parameters) {
    Pinching4 material(parameters);
    WallUniaxial m;m.kind_=Kind::Pinching4;m.pinching4_parameters_=std::make_shared<const Pinching4CyclicParameters>(std::move(parameters));
    m.E_=material.initial_state().tangent;return m;
}
WallUniaxial WallUniaxial::minmax(WallUniaxial material,double min_strain,double max_strain) {
    if(!std::isfinite(min_strain)||!std::isfinite(max_strain)||min_strain>=max_strain)
        throw std::invalid_argument("MinMax requires finite min strain < max strain");
    WallUniaxial m;m.kind_=Kind::MinMax;m.min_strain_=min_strain;m.max_strain_=max_strain;
    m.children_=std::make_shared<const std::vector<WallUniaxial>>(std::vector<WallUniaxial>{std::move(material)});return m;
}
WallUniaxial WallUniaxial::parallel(std::vector<WallUniaxial> materials) {
    if(materials.empty())throw std::invalid_argument("Parallel requires at least one material");
    WallUniaxial m;m.kind_=Kind::Parallel;m.children_=std::make_shared<const std::vector<WallUniaxial>>(std::move(materials));return m;
}
double WallUniaxial::initial_tangent() const {
    if(kind_==Kind::MinMax)return children_->front().initial_tangent();
    if(kind_==Kind::Parallel){double k=0;for(const auto& child:*children_)k+=child.initial_tangent();return k;}
    return E_;
}
int WallUniaxial::state_size() const {
    if(kind_==Kind::ConcreteCM)return concrete_cm_state_size;
    if(kind_==Kind::Pinching4)return pinching4_state_size;
    if(kind_==Kind::MinMax){const int n=children_->front().state_size();if(n==std::numeric_limits<int>::max())throw std::overflow_error("wall material state size overflow");return n+1;}
    if(kind_==Kind::Parallel){int n=0;for(const auto& child:*children_){const int c=child.state_size();if(c>std::numeric_limits<int>::max()-n)throw std::overflow_error("wall material state size overflow");n+=c;}return n;}
    return 6;
}
void WallUniaxial::initialize(double* s) const {
    if(kind_==Kind::ConcreteCM){encode_concrete_cm_state(ConcreteCM(*concrete_cm_parameters_).initial_state(),s);return;}
    if(kind_==Kind::Pinching4){encode_pinching4_state(Pinching4(*pinching4_parameters_).initial_state(),s);return;}
    if(kind_==Kind::MinMax){const auto& child=children_->front();child.initialize(s);s[child.state_size()]=0;return;}
    if(kind_==Kind::Parallel){int o=0;for(const auto& child:*children_){child.initialize(s+o);o+=child.state_size();}return;}
    std::fill(s,s+6,0.0);if(kind_==Kind::Concrete01)s[2]=E_;
}
WallUniaxial::Result WallUniaxial::trial(double e,const double* s,double* t) const {
    if(!std::isfinite(e))throw std::invalid_argument("nonfinite wall material strain");
    const int nstate=state_size();for(int i=0;i<nstate;++i)if(!std::isfinite(s[i]))throw std::invalid_argument("nonfinite wall material state");
    if(kind_==Kind::ConcreteCM){const auto r=ConcreteCM(*concrete_cm_parameters_).trial(e,decode_concrete_cm_state(s));encode_concrete_cm_state(r.state,t);return {r.response.stress,r.response.tangent};}
    if(kind_==Kind::Pinching4){const auto r=Pinching4(*pinching4_parameters_).trial(e,decode_pinching4_state(s));encode_pinching4_state(r.state,t);return {r.response.force,r.response.tangent};}
    if(kind_==Kind::MinMax){
        const auto& child=children_->front();const int n=child.state_size();std::copy(s,s+n+1,t);
        // Match OpenSees MinMaxMaterial exactly: the limits themselves fail,
        // a failed wrapper returns zero stress, and the tangent is a tiny
        // residual stiffness based on the child's initial tangent.
        if(s[n]!=0.0||e<=min_strain_||e>=max_strain_){t[n]=1.0;return {0,1.0e-8*child.initial_tangent()};}
        auto r=child.trial(e,s,t);t[n]=0.0;return r;
    }
    if(kind_==Kind::Parallel){
        double stress=0,tangent=0;int o=0;
        for(const auto& child:*children_){auto r=child.trial(e,s+o,t+o);stress+=r.stress;tangent+=r.tangent;o+=child.state_size();}
        return {stress,tangent};
    }
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
        double c=std::cos(l.angle_rad),s=std::sin(l.angle_rad);directions_.push_back({c*c,s*s,c*s});offsets_.push_back(state_size_);state_size_+=l.material.state_size();
    }
}
std::vector<double> WallPanel::initial_state() const {
    std::vector<double> s(state_size());for(std::size_t i=0;i<layers_.size();++i)layers_[i].material.initialize(s.data()+offsets_[i]);return s;
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
        auto v=layers_[i].material.trial(q,s+offsets_[i],r.state.data()+offsets_[i]);
        for(int a=0;a<3;++a){r.stress[a]+=layers_[i].weight*v.stress*n[a];for(int b=0;b<3;++b)r.tangent[3*a+b]+=layers_[i].weight*v.tangent*n[a]*n[b];}
    }
    return r;
}
} // namespace quake
