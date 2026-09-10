#include "quake/asce41_hinge.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {
constexpr double eps=1e-14;
inline double signum(double x){ return x>=0.0?1.0:-1.0; }
inline double clamp01(double x){ return std::clamp(x,0.0,1.0); }
}

ASCE41HingeMaterial::ASCE41HingeMaterial(ASCE41HingeParams p):p_(p){
    if(p_.Ke<=0.0||p_.posFy<=0.0||p_.negFy<=0.0||
       p_.hardening_ratio<0.0||p_.pos_a<=0.0||p_.neg_a<=0.0||
       p_.pos_b<=p_.pos_a||p_.neg_b<=p_.neg_a||
       (p_.pos_f>0.0&&p_.pos_f<p_.pos_b)||(p_.neg_f>0.0&&p_.neg_f<p_.neg_b)||
       p_.posMc<0.0||p_.negMc<0.0||
       p_.pos_c<0.0||p_.pos_c>1.0||p_.neg_c<0.0||p_.neg_c>1.0||
       p_.pos_drop_span<0.0||p_.neg_drop_span<0.0||p_.pos_e_drop_span<0.0||p_.neg_e_drop_span<0.0||
       p_.lambda_strength<=0.0||p_.lambda_unloading<=0.0||p_.cyclic_exponent<=0.0)
        throw std::invalid_argument("invalid ASCE41 hinge parameters");
    if(p_.backbone_shape==ASCE41BackboneShape::ResearchExtendedCDE){
        const double ps=p_.pos_drop_span>0.0?p_.pos_drop_span:0.05*(p_.pos_b-p_.pos_a);
        const double ns=p_.neg_drop_span>0.0?p_.neg_drop_span:0.05*(p_.neg_b-p_.neg_a);
        const double pe=p_.pos_e_drop_span>0.0?p_.pos_e_drop_span:0.05*(p_.pos_b-p_.pos_a-ps);
        const double ne=p_.neg_e_drop_span>0.0?p_.neg_e_drop_span:0.05*(p_.neg_b-p_.neg_a-ns);
        if(ps>=p_.pos_b-p_.pos_a||ns>=p_.neg_b-p_.neg_a||pe<=0.0||ne<=0.0||
           ps+pe>=p_.pos_b-p_.pos_a||ns+ne>=p_.neg_b-p_.neg_a)
            throw std::invalid_argument("ASCE41 degradation spans must fit between a and b");
    }
}

void ASCE41HingeMaterial::initialize_state(double* s) const{
    std::fill(s,s+kStateSize,0.0);
    s[4]=p_.Ke;      // unloading stiffness
    s[10]=1.0;       // positive strength scale
    s[11]=1.0;       // negative strength scale
    s[18]=p_.Ke;      // current tangent
}

double ASCE41HingeMaterial::degradation_beta(double excursion,double cumulative,double lambda) const{
    if(excursion<=0.0) return 0.0;
    const double ref=lambda*std::max(p_.posFy,p_.negFy);
    const double remaining=ref-cumulative;
    if(remaining<=eps) return 1.0;
    return clamp01(std::pow(std::max(0.0,excursion/remaining),p_.cyclic_exponent));
}

ASCE41HingeMaterial::EnvelopePoint ASCE41HingeMaterial::envelope(double u,double strength_scale) const{
    EnvelopePoint r;
    const bool pos=u>=0.0;
    const double sg=pos?1.0:-1.0;
    const double x=std::abs(u);
    const double fy0=pos?p_.posFy:p_.negFy;
    const double a=pos?p_.pos_a:p_.neg_a;
    const double b=pos?p_.pos_b:p_.neg_b;
    const double fpar=pos?p_.pos_f:p_.neg_f;
    const double f=(fpar>0.0)?fpar:b;
    const double c=pos?p_.pos_c:p_.neg_c;
    const double drop_requested=pos?p_.pos_drop_span:p_.neg_drop_span;
    const double drop=drop_requested>0.0?drop_requested:0.05*(b-a);
    const double e_requested=pos?p_.pos_e_drop_span:p_.neg_e_drop_span;
    const double e_drop=e_requested>0.0?e_requested:0.05*(b-a-drop);
    const double fy=std::max(1e-12*fy0,fy0*strength_scale);
    const double uy=fy/p_.Ke;
    const double uc=uy+a;
    const double ud=uc+drop;
    const double ue=uy+b;
    const double ue0=ue-e_drop;
    const double uf=uy+f;
    const double mc0=pos?p_.posMc:p_.negMc;
    const double kh_default=(p_.hardening_stiffness>0.0)?p_.hardening_stiffness:p_.hardening_ratio*p_.Ke;
    const double fc=(mc0>0.0)?mc0*strength_scale:(fy+kh_default*a);
    const double kh=(mc0>0.0)?(fc-fy)/a:kh_default;
    const double fr=fy0*c*strength_scale;

    if(x>=uf){
        r.force=0.0; r.tangent=0.0; r.events|=ASCE41_EVENT_LATERAL_LOSS|ASCE41_EVENT_FAILURE; r.failed=true; return r;
    }
    if(x>=ue){
        r.force=0.0; r.tangent=0.0; r.events|=ASCE41_EVENT_LATERAL_LOSS; return r;
    }
    if(x<=uy){r.force=sg*p_.Ke*x;r.tangent=p_.Ke;return r;}
    if(x<=uc){
        r.force=sg*(fy+kh*(x-uy));
        r.tangent=kh;
        r.events|=ASCE41_EVENT_YIELD;return r;
    }
    if(p_.backbone_shape==ASCE41BackboneShape::StraightCE){
        const double kde=(0.0-fc)/(ue-uc);
        const double mag=std::max(0.0,fc+kde*(x-uc));
        r.force=sg*mag;
        r.tangent=(mag>0.0)?kde:0.0;
        r.events|=ASCE41_EVENT_YIELD|ASCE41_EVENT_CAP|ASCE41_EVENT_STRENGTH_DROP;
        return r;
    }
    if(x<=ud){
        const double kd=(fr-fc)/drop;
        r.force=sg*(fc+kd*(x-uc));r.tangent=kd;
        r.events|=ASCE41_EVENT_YIELD|ASCE41_EVENT_CAP|ASCE41_EVENT_STRENGTH_DROP;return r;
    }
    if(x<=ue0){
        r.force=sg*fr;r.tangent=0.0;
        r.events|=ASCE41_EVENT_YIELD|ASCE41_EVENT_CAP|ASCE41_EVENT_RESIDUAL;
        return r;
    }
    const double kE=-fr/e_drop;
    const double mag=std::max(0.0,fr+kE*(x-ue0));
    r.force=sg*mag;r.tangent=(mag>0.0)?kE:0.0;
    r.events|=ASCE41_EVENT_YIELD|ASCE41_EVENT_CAP|ASCE41_EVENT_RESIDUAL;
    return r;
}

ASCE41TrialResult ASCE41HingeMaterial::trial(double u,const double* c,double* s) const{
    enum : int { Q=0,F=1,CUM=2,EXC=3,KUN=4,DIR=5,POSQ=6,POSF=7,NEGQ=8,NEGF=9,
                 POSS=10,NEGS=11,FAIL=12,LATLOSS=13,BRANCH=14,UZERO=15,TARGETQ=16,TARGETF=17,KT=18 };
    std::copy(c,c+kStateSize,s);
    ASCE41TrialResult out;

    const auto side_f=[&](double q){
        const bool pos=q>=0.0;
        const double fy0=pos?p_.posFy:p_.negFy;
        const double fpar=pos?p_.pos_f:p_.neg_f;
        const double b=pos?p_.pos_b:p_.neg_b;
        return fy0/p_.Ke + ((fpar>0.0)?fpar:b);
    };

    // F is an irreversible effective/gravity component failure. E is an
    // irreversible loss of lateral resistance, but the component can remain
    // present for gravity/effective-capacity bookkeeping until F is reached.
    if(c[FAIL]>0.5){
        s[Q]=u;s[F]=0.0;s[KT]=0.0;
        out.events=ASCE41_EVENT_LATERAL_LOSS|ASCE41_EVENT_FAILURE;
        return out;
    }
    if(c[LATLOSS]>0.5){
        s[Q]=u;s[F]=0.0;s[KT]=0.0;
        out.events=ASCE41_EVENT_LATERAL_LOSS;
        if(std::abs(u)>=side_f(u)-eps){s[FAIL]=1.0;out.events|=ASCE41_EVENT_FAILURE;}
        return out;
    }

    const double up=c[Q],fp=c[F],du=u-up;
    if(std::abs(du)<=eps){
        out.force=fp;out.tangent=(c[KT]!=0.0?c[KT]:p_.Ke);
        const bool virgin=c[CUM]<=eps&&std::abs(fp-p_.Ke*up)<=1e-10*std::max(1.0,p_.posFy);
        if(virgin)out.events|=ASCE41_EVENT_FAST_ELASTIC;
        return out;
    }

    const double dir=signum(du),last=c[DIR];
    const bool reversal=std::abs(last)>0.5&&dir*last<0.0;
    double unload=std::clamp(c[KUN]>0.0?c[KUN]:p_.Ke,1e-9*p_.Ke,p_.Ke);
    double pos_scale=c[POSS]>0.0?c[POSS]:1.0;
    double neg_scale=c[NEGS]>0.0?c[NEGS]:1.0;
    double excursion=c[EXC];
    unsigned events=ASCE41_EVENT_NONE;
    bool branch=c[BRANCH]>0.5;
    double u0=c[UZERO],target_q=c[TARGETQ],target_f=c[TARGETF];

    // Virgin elastic evaluation remains the cheapest active-front path.
    const bool pos_trial=u>=0.0;
    const double fy_trial=(pos_trial?p_.posFy:p_.negFy)*(pos_trial?pos_scale:neg_scale);
    const double uy_trial=fy_trial/p_.Ke;
    const bool virgin=c[CUM]<=eps && !branch && std::abs(fp-p_.Ke*up)<=1e-10*std::max(1.0,fy_trial)
                      && std::abs(u)<=uy_trial+1e-14;
    if(virgin){
        out.force=p_.Ke*u;out.tangent=p_.Ke;out.events=ASCE41_EVENT_FAST_ELASTIC;
        s[Q]=u;s[F]=out.force;s[KUN]=p_.Ke;s[DIR]=dir;s[KT]=out.tangent;
        if(u>s[POSQ]){s[POSQ]=u;s[POSF]=out.force;}
        if(u<s[NEGQ]){s[NEGQ]=u;s[NEGF]=out.force;}
        return out;
    }

    if(reversal){
        events|=ASCE41_EVENT_REVERSAL;
        const double bK=degradation_beta(excursion,c[CUM],p_.lambda_unloading);
        const double bS=degradation_beta(excursion,c[CUM],p_.lambda_strength);
        unload=std::max(1e-9*p_.Ke,unload*(1.0-bK));
        if(dir>0.0)pos_scale=std::max(1e-6,pos_scale*(1.0-bS));
        else neg_scale=std::max(1e-6,neg_scale*(1.0-bS));
        if(bK>0.0||bS>0.0)events|=ASCE41_EVENT_DETERIORATION;
        excursion=0.0;

        // Explicit continuous reversal path. The unloading line passes exactly
        // through the committed point and zero force at u0. Reloading then
        // targets the current degraded envelope at the historical opposite
        // peak (or a fallback target beyond u0). No branch is inferred from
        // the Newton trial endpoint.
        u0=up-fp/unload;
        const bool to_pos=dir>0.0;
        const double scale=to_pos?pos_scale:neg_scale;
        const double fy0=to_pos?p_.posFy:p_.negFy;
        const double uy=fy0*scale/p_.Ke;
        const double a=to_pos?p_.pos_a:p_.neg_a;
        const double b=to_pos?p_.pos_b:p_.neg_b;
        const double ed_req=to_pos?p_.pos_e_drop_span:p_.neg_e_drop_span;
        const double drop_req=to_pos?p_.pos_drop_span:p_.neg_drop_span;
        const double drop=drop_req>0.0?drop_req:0.05*(b-a);
        const double ed=ed_req>0.0?ed_req:0.05*(b-a-drop);
        const double hist=to_pos?c[POSQ]:c[NEGQ];
        const double cap_q=(to_pos?1.0:-1.0)*(uy+a);
        const double pre_e=(to_pos?1.0:-1.0)*(uy+b-ed);
        auto beyond=[&](double q,double ref){return to_pos?q>ref+eps:q<ref-eps;};
        if(beyond(hist,u0)) target_q=hist;
        else if(beyond(cap_q,u0)) target_q=cap_q;
        else if(beyond(pre_e,u0)) target_q=pre_e;
        else target_q=(to_pos?1.0:-1.0)*(std::abs(u0)+std::max(uy,1e-8));
        auto tenv=envelope(target_q,scale);
        if((tenv.events&ASCE41_EVENT_LATERAL_LOSS)!=0 || std::abs(tenv.force)<=eps){
            // If no positive-strength envelope target remains beyond the zero
            // intercept, retain a continuous zero-force approach to E rather
            // than introduce a force jump.
            target_f=0.0;
        }else target_f=tenv.force;
        branch=true;
    }

    bool on_envelope=!branch;
    if(branch){
        const bool to_pos=(target_q>=u0);
        const bool before_zero=to_pos?(u<=u0):(u>=u0);
        const bool before_target=to_pos?(u<target_q):(u>target_q);
        if(before_zero){
            out.force=unload*(u-u0);out.tangent=unload;events|=ASCE41_EVENT_RELOAD;
        }else if(before_target){
            const double den=target_q-u0;
            const double kr=(std::abs(den)>eps)?target_f/den:0.0;
            out.force=kr*(u-u0);out.tangent=kr;events|=ASCE41_EVENT_RELOAD;
        }else{
            on_envelope=true;branch=false;
        }
    }

    if(on_envelope){
        const bool pos=u>=0.0;
        auto env=envelope(u,pos?pos_scale:neg_scale);
        out.force=env.force;out.tangent=env.tangent;events|=env.events;
        if((env.events&ASCE41_EVENT_LATERAL_LOSS)!=0)s[LATLOSS]=1.0;
        if(env.failed){s[FAIL]=1.0;s[LATLOSS]=1.0;}
    }

    const double dE=0.5*(fp+out.force)*du;
    s[CUM]=c[CUM]+std::abs(dE);s[EXC]=excursion+std::abs(dE);
    s[Q]=u;s[F]=out.force;s[KUN]=unload;s[DIR]=dir;s[POSS]=pos_scale;s[NEGS]=neg_scale;
    s[BRANCH]=branch?1.0:0.0;s[UZERO]=u0;s[TARGETQ]=target_q;s[TARGETF]=target_f;s[KT]=out.tangent;
    if(on_envelope && u>s[POSQ]&&out.force>=0.0){s[POSQ]=u;s[POSF]=out.force;}
    if(on_envelope && u<s[NEGQ]&&out.force<=0.0){s[NEGQ]=u;s[NEGF]=out.force;}
    out.events|=events;
    return out;
}

ASCE41PerformanceLevel ASCE41HingeMaterial::performance_level(double u) const{
    const bool pos=u>=0.0;const double fy=pos?p_.posFy:p_.negFy;
    const double plastic=std::max(0.0,std::abs(u)-fy/p_.Ke);
    const double io=pos?p_.pos_io:p_.neg_io,ls=pos?p_.pos_ls:p_.neg_ls,cp=pos?p_.pos_cp:p_.neg_cp;
    const double b=pos?p_.pos_b:p_.neg_b;
    const double fpar=pos?p_.pos_f:p_.neg_f;
    const double f=(fpar>0.0)?fpar:b;
    if(plastic>=f)return ASCE41PerformanceLevel::Failed;
    if(cp>0.0&&plastic>cp)return ASCE41PerformanceLevel::BeyondCP;
    if(cp>0.0&&plastic>=cp)return ASCE41PerformanceLevel::CP;
    if(ls>0.0&&plastic>=ls)return ASCE41PerformanceLevel::LS;
    if(io>0.0&&plastic>=io)return ASCE41PerformanceLevel::IO;
    if(plastic>0.0)return ASCE41PerformanceLevel::InelasticBelowIO;
    return ASCE41PerformanceLevel::Elastic;
}

} // namespace quake
