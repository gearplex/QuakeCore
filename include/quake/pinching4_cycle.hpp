#pragma once
#include "quake/pinching4.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace quake {

enum class Pinching4StateKind:int{Initial=0,PositiveEnvelope=1,NegativeEnvelope=2,PositiveToNegative=3,NegativeToPositive=4};
enum class Pinching4DamageMode:int{Energy=0,Cycle=1};
struct Pinching4CyclicParameters{
    Pinching4EnvelopeParameters envelope{};
    double r_disp_positive{},r_force_positive{},u_force_positive{};
    double r_disp_negative{},r_force_negative{},u_force_negative{};
    std::array<double,4> gamma_k{};double gamma_k_limit{};
    std::array<double,4> gamma_d{};double gamma_d_limit{};
    std::array<double,4> gamma_f{};double gamma_f_limit{};
    double gamma_e{};Pinching4DamageMode damage_mode{Pinching4DamageMode::Energy};
    unsigned admitted_reversal_count{2};
};
struct Pinching4State{
    Pinching4StateKind kind{Pinching4StateKind::Initial};
    double deformation{},force{},tangent{},deformation_rate{};
    double low_state_deformation{},low_state_force{},high_state_deformation{},high_state_force{};
    double min_demand{},max_demand{},energy{},gamma_d{},cycles{},damaged_min{},damaged_max{};
    unsigned reversal_count{};
};
struct Pinching4Trial{Pinching4Response response{};Pinching4State state{};};

// Gate 4 stateful slice: reversal admission is explicit and protocol-scoped.
// The default preserves the frozen steel_raw step-120 boundary; focused oracle
// tests may opt into later admitted reversals as each milestone is verified.
class Pinching4{
public:
    explicit Pinching4(Pinching4CyclicParameters p):p_(p),env_(p.envelope){
        if(!(p_.gamma_e>0.0))throw std::invalid_argument("Pinching4 gammaE must be positive");
        if(p_.admitted_reversal_count<2)throw std::invalid_argument("Pinching4 Gate 4 slice requires at least two admitted reversals");
        for(double x:p_.gamma_k)if(x!=0.0)throw std::invalid_argument("Pinching4 Gate 4 slice requires zero gammaK");
        for(double x:p_.gamma_f)if(x!=0.0)throw std::invalid_argument("Pinching4 Gate 4 slice requires zero gammaF");
    }
    Pinching4State initial_state()const{
        Pinching4State s;const auto p=positive_points();const auto n=negative_points();
        s.tangent=env_.positive_elastic_tangent();s.low_state_deformation=n[0].deformation;s.low_state_force=n[0].force;
        s.high_state_deformation=p[0].deformation;s.high_state_force=p[0].force;s.min_demand=n[1].deformation;s.max_demand=p[1].deformation;
        s.damaged_min=s.min_demand;s.damaged_max=s.max_demand;return s;
    }
    Pinching4Trial trial(double u,const Pinching4State&c)const{
        Pinching4State t=c;double du=u-c.deformation;if(du<1e-12&&du>-1e-12)du=0.0;
        if(u<t.low_state_deformation||u>t.high_state_deformation||du*c.deformation_rate<=0.0)update_state(u,du,c,t);
        Pinching4Response r{};
        switch(t.kind){
        case Pinching4StateKind::Initial:r={env_.positive_elastic_tangent()*u,env_.positive_elastic_tangent()};break;
        case Pinching4StateKind::PositiveEnvelope:r=env_.positive(u);break;
        case Pinching4StateKind::NegativeEnvelope:r=env_.negative(u);break;
        case Pinching4StateKind::PositiveToNegative:r=state3(t,u);break;
        case Pinching4StateKind::NegativeToPositive:r=state4(t,u);break;}
        t.energy=c.energy+0.5*(r.force+c.force)*du;
        const double ke=u>0.0?env_.positive_elastic_tangent():env_.negative_elastic_tangent();
        update_damage(t,u,du,0.5*r.force/ke*r.force);
        t.deformation=u;t.force=r.force;t.tangent=r.tangent;if(du>1e-12||du<-1e-12)t.deformation_rate=du;
        t.damaged_max=t.max_demand*(1.0+t.gamma_d);t.damaged_min=t.min_demand*(1.0+t.gamma_d);return{r,t};
    }
private:
    std::array<Pinching4Point,6> positive_points()const{return envelope_points(true);}
    std::array<Pinching4Point,6> negative_points()const{return envelope_points(false);}
    std::array<Pinching4Point,6> envelope_points(bool positive)const{
        std::array<Pinching4Point,6> q{};const auto& src=positive?p_.envelope.positive:p_.envelope.negative;
        const double kp=p_.envelope.positive[0].force/p_.envelope.positive[0].deformation;
        const double kn=p_.envelope.negative[0].force/p_.envelope.negative[0].deformation;
        const double k=std::max(kp,kn),u=std::max(p_.envelope.positive[0].deformation,-p_.envelope.negative[0].deformation)*1.0e-4;
        q[0]=positive?Pinching4Point{u,u*k}:Pinching4Point{-u,-u*k};for(std::size_t i=0;i<4;++i)q[i+1]=src[i];
        const double last=slope(src[2],src[3]);q[5].deformation=1.0e6*src[3].deformation;
        q[5].force=last>0.0?src[3].force+last*(q[5].deformation-src[3].deformation):1.1*src[3].force;return q;
    }
    void update_state(double u,double du,const Pinching4State&c,Pinching4State&t)const{
        if(t.kind==Pinching4StateKind::Initial){
            if(u>t.high_state_deformation){to_positive_envelope(t);}
            else if(u<t.low_state_deformation){to_negative_envelope(t);}return;}
        if(t.kind==Pinching4StateKind::PositiveEnvelope&&du<0.0){
            if(t.reversal_count>=p_.admitted_reversal_count)throw std::logic_error("Pinching4 additional positive-to-negative reversal is not yet admitted in Gate 4");
            ++t.reversal_count;
            if(c.deformation>t.max_demand)t.max_demand=c.deformation;
            if(t.max_demand<t.damaged_max)t.max_demand=t.damaged_max;
            if(u<t.damaged_min){to_negative_envelope(t);}else{t.kind=Pinching4StateKind::PositiveToNegative;t.low_state_deformation=t.damaged_min;t.low_state_force=env_.negative(t.damaged_min).force;t.high_state_deformation=c.deformation;t.high_state_force=c.force;}return;}
        if(t.kind==Pinching4StateKind::PositiveToNegative){if(u<t.low_state_deformation)to_negative_envelope(t);else if(du>0.0)throw std::logic_error("Pinching4 reversal from state 3 is not yet admitted in Gate 4");return;}
        if(t.kind==Pinching4StateKind::NegativeEnvelope&&du>0.0){
            if(t.reversal_count>=p_.admitted_reversal_count)throw std::logic_error("Pinching4 additional negative-to-positive reversal is not yet admitted in Gate 4");
            ++t.reversal_count;
            if(c.deformation<t.min_demand)t.min_demand=c.deformation;
            if(t.min_demand>t.damaged_min)t.min_demand=t.damaged_min;
            if(u>t.damaged_max){to_positive_envelope(t);}else{t.kind=Pinching4StateKind::NegativeToPositive;t.low_state_deformation=c.deformation;t.low_state_force=c.force;t.high_state_deformation=t.damaged_max;t.high_state_force=env_.positive(t.damaged_max).force;}return;}
        if(t.kind==Pinching4StateKind::NegativeToPositive){
            if(u>t.high_state_deformation)to_positive_envelope(t);
            else if(du<0.0)throw std::logic_error("Pinching4 reversal from state 4 is not yet admitted in Gate 4");
        }
    }
    void to_positive_envelope(Pinching4State&t)const{const auto p=positive_points();t.kind=Pinching4StateKind::PositiveEnvelope;t.low_state_deformation=p[0].deformation;t.low_state_force=p[0].force;t.high_state_deformation=p[5].deformation;t.high_state_force=p[5].force;}
    void to_negative_envelope(Pinching4State&t)const{const auto n=negative_points();t.kind=Pinching4StateKind::NegativeEnvelope;t.low_state_deformation=n[5].deformation;t.low_state_force=n[5].force;t.high_state_deformation=n[0].deformation;t.high_state_force=n[0].force;}
    Pinching4Response state3(const Pinching4State&s,double u)const{
        std::array<Pinching4Point,4>q{};q[0]={s.low_state_deformation,s.low_state_force};q[3]={s.high_state_deformation,s.high_state_force};
        const double kunload=q[3].deformation<0.0?env_.negative_elastic_tangent():env_.positive_elastic_tangent();const double kmax=std::max(kunload,env_.negative_elastic_tangent());
        if(q[0].deformation*q[3].deformation>=0.0)linear(q);else{
            q[1].deformation=q[0].deformation*p_.r_disp_negative;if(p_.r_force_negative-p_.u_force_negative<=1e-8)throw std::logic_error("Pinching4 alternate state-3 force branch is not yet admitted");
            q[1].force=q[0].force*p_.r_force_negative;if(slope(q[0],q[1])>env_.negative_elastic_tangent())q[1].deformation=q[0].deformation+(q[1].force-q[0].force)/env_.negative_elastic_tangent();
            if(q[1].deformation>q[3].deformation)linear(q);else{const auto n=negative_points();q[2].force=p_.u_force_negative*(s.min_demand<n[3].deformation?n[4].force:n[3].force);q[2].deformation=q[3].deformation-(q[3].force-q[2].force)/kunload;
                if(q[2].deformation>q[3].deformation)midpoint(q);else if(slope(q[1],q[2])>kmax)linear(q);else if(q[2].deformation<q[1].deformation||slope(q[1],q[2])<0.0)throw std::logic_error("Pinching4 alternate state-3 ordering branch is not yet admitted");}}
        for(std::size_t i=0;i<3;++i){if(u>=q[i].deformation&&(i==2||u<q[i+1].deformation))return interp(u,q[i],q[i+1]);}return interp(u,q[0],q[1]);
    }
    Pinching4Response state4(const Pinching4State&s,double u)const{
        std::array<Pinching4Point,4>q{};q[0]={s.low_state_deformation,s.low_state_force};q[3]={s.high_state_deformation,s.high_state_force};
        const double kunload=q[0].deformation<0.0?env_.negative_elastic_tangent():env_.positive_elastic_tangent();const double kmax=std::max(kunload,env_.positive_elastic_tangent());
        if(q[0].deformation*q[3].deformation>=0.0)linear(q);else{
            q[2].deformation=q[3].deformation*p_.r_disp_positive;if(p_.u_force_positive==0.0||p_.r_force_positive-p_.u_force_positive>1e-8)q[2].force=q[3].force*p_.r_force_positive;else throw std::logic_error("Pinching4 alternate state-4 force branch is not yet admitted");
            if(slope(q[2],q[3])>env_.positive_elastic_tangent())q[2].deformation=q[3].deformation-(q[3].force-q[2].force)/env_.positive_elastic_tangent();
            if(q[2].deformation<q[0].deformation)linear(q);else{const auto p=positive_points();q[1].force=p_.u_force_positive*(s.max_demand>p[3].deformation?p[4].force:p[3].force);q[1].deformation=q[0].deformation+(-q[0].force+q[1].force)/kunload;
                if(q[1].deformation<q[0].deformation)midpoint_left(q);else if(slope(q[1],q[2])>kmax)linear(q);else if(q[2].deformation<q[1].deformation||slope(q[1],q[2])<0.0)throw std::logic_error("Pinching4 alternate state-4 ordering branch is not yet admitted");}}
        for(std::size_t i=0;i<3;++i){if(u>=q[i].deformation&&(i==2||u<q[i+1].deformation))return interp(u,q[i],q[i+1]);}return interp(u,q[0],q[1]);
    }
    void update_damage(Pinching4State&s,double u,double du,double elastic)const{
        const auto p=positive_points();const auto n=negative_points();const double umax=std::max(s.max_demand,-s.min_demand),uult=std::max(p[4].deformation,-n[4].deformation);s.cycles+=std::abs(du)/(4.0*umax);
        if(u<uult&&u>-uult&&s.energy<energy_capacity()){s.gamma_d=p_.gamma_d[0]*std::pow(umax/uult,p_.gamma_d[2]);if(s.energy>elastic&&p_.damage_mode==Pinching4DamageMode::Energy){double e=(s.energy-elastic)/energy_capacity();s.gamma_d+=p_.gamma_d[1]*std::pow(e,p_.gamma_d[3]);}else if(p_.damage_mode==Pinching4DamageMode::Cycle)s.gamma_d+=p_.gamma_d[1]*std::pow(s.cycles,p_.gamma_d[3]);s.gamma_d=std::min(s.gamma_d,p_.gamma_d_limit);}else if(u<uult&&u>-uult)s.gamma_d=p_.gamma_d_limit;
    }
    double energy_capacity()const{const auto p=positive_points();const auto n=negative_points();double ep=.5*p[0].deformation*p[0].force;double en=.5*n[0].deformation*n[0].force;for(std::size_t i=0;i<4;++i){ep+=.5*(p[i].force+p[i+1].force)*(p[i+1].deformation-p[i].deformation);en+=.5*(n[i].force+n[i+1].force)*(n[i+1].deformation-n[i].deformation);}return p_.gamma_e*std::max(ep,en);}
    static double slope(const Pinching4Point&a,const Pinching4Point&b){return(b.force-a.force)/(b.deformation-a.deformation);}static Pinching4Response interp(double u,const Pinching4Point&a,const Pinching4Point&b){double k=slope(a,b);return{a.force+(u-a.deformation)*k,k};}
    static void linear(std::array<Pinching4Point,4>&q){double du=q[3].deformation-q[0].deformation,df=q[3].force-q[0].force;q[1]={q[0].deformation+.33*du,q[0].force+.33*df};q[2]={q[0].deformation+.67*du,q[0].force+.67*df};}
    static void midpoint(std::array<Pinching4Point,4>&q){double du=q[3].deformation-q[1].deformation,df=q[3].force-q[1].force;q[2]={q[1].deformation+.5*du,q[1].force+.5*df};}
    static void midpoint_left(std::array<Pinching4Point,4>&q){double du=q[2].deformation-q[0].deformation,df=q[2].force-q[0].force;q[1]={q[0].deformation+.5*du,q[0].force+.5*df};}
    Pinching4CyclicParameters p_;Pinching4Envelope env_;
};
} // namespace quake
