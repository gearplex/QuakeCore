#include "quake/pm_interaction_hinge.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace quake {
namespace {
constexpr double tiny=1e-14;

std::vector<double> monotone_slopes(const std::vector<double>& x,const std::vector<double>& y){
    if(x.size()<2||x.size()!=y.size())throw std::invalid_argument("P-M capacity table size");
    const std::size_t n=x.size();std::vector<double> d(n-1),h(n-1),m(n,0.0);
    for(std::size_t i=0;i+1<n;++i){
        if(!(x[i+1]>x[i]))throw std::invalid_argument("P-M axial force points must increase");
        h[i]=x[i+1]-x[i];d[i]=(y[i+1]-y[i])/h[i];
    }
    m[0]=d[0];m[n-1]=d[n-2];
    for(std::size_t i=1;i+1<n;++i){
        if(d[i-1]*d[i]<=0.0)m[i]=0.0;
        else{const double w1=2*h[i]+h[i-1],w2=h[i]+2*h[i-1];m[i]=(w1+w2)/(w1/d[i-1]+w2/d[i]);}
    }
    if(m[0]*d[0]<=0.0)m[0]=0.0;
    if(m[n-1]*d[n-2]<=0.0)m[n-1]=0.0;
    return m;
}

// A convex P-M surface has a unique elastic-metric projection. When the
// trial moment is below MyB, only a narrow interval near an axial intercept
// has a nonnegative plastic multiplier. A uniform P scan can miss that entire
// interval. Bound it analytically by My(P)=|Mtrial|, then bisect the return
// equation. No arbitrary scan spacing or response calibration is involved.
template<class Eval>
bool axial_coordinate_projection(double ptr,double mabs,double pb,double pt,double pc,
                        double my,double at,double ac,double beta,double tol,Eval eval,double& root){
    if(mabs<=0.0)return false;
    double lo=pt,hi=pc;
    if(mabs<my){const double v=std::max(0.0,1.0-std::pow(mabs/my,beta));
        if(ptr<pb)hi=pb+(pt-pb)*std::pow(v,1.0/at);else lo=pb+(pc-pb)*std::pow(v,1.0/ac);}
    for(int i=0;i<120;++i){const double mid=lo+0.5*(hi-lo);if(mid==lo||mid==hi)break;
        root=mid;const double g=eval(mid);if(!std::isfinite(g))return false;
        if(std::abs(g)<=tol)return true;
        if(g>0.0)lo=mid;else hi=mid;
    }
    root=lo+0.5*(hi-lo);return root>pt&&root<pc&&std::abs(eval(root))<=10.0*tol;
}
struct TipProjection { bool ok{}; double P{},M{},dpdm{},d2pdm2{}; };
TipProjection bracket_projection(double ptr,double mabs,double pb,double pt,double pc,
                        double my,double at,double ac,double beta,double ka,double kr,double tol){
    TipProjection out;
    if(mabs<=0.0)return out; // exact axial tips are handled separately
    // P can round to its intercept while the corresponding moment is still
    // nonzero. Solve in q=(M/My)^(beta-1) instead of P: the normal direction
    // remains resolvable even in that limit (beta is greater than one).
    const double den=(ptr<pb?pt:pc)-pb,alpha=ptr<pb?at:ac;
    double lo=0.0,hi=mabs<my?std::pow(mabs/my,beta-1.0):1.0;
    for(int i=0;i<160;++i){
        const double q=lo+0.5*(hi-lo);
        const double t=std::pow(q,1.0/(beta-1.0));
        const double v=1.0-std::pow(t,beta);
        if(v<=0.0){hi=q;continue;}
        out.P=pb+den*std::pow(v,1.0/alpha);out.M=my*t;
        out.dpdm=-den*beta/(alpha*my)*q*std::pow(v,1.0/alpha-1.0);
        const double g=(out.P-ptr)*(kr/ka)*out.dpdm+out.M-mabs;
        if(std::isfinite(g)&&std::abs(g)<=tol*std::max(my,mabs)){
            out.ok=true;
            out.d2pdm2=-den*beta/(alpha*my*my)*(
                (beta-1.0)*std::pow(t,beta-2.0)*std::pow(v,1.0/alpha-1.0)
                -beta*(1.0/alpha-1.0)*q*q*std::pow(v,1.0/alpha-2.0));
            return out;
        }
        if(q==lo||q==hi)break;
        if(g>0.0)hi=q;else lo=q;
    }
    return out;
}

inline double sgn(double x){return x>=0.0?1.0:-1.0;}
}

PMInteractionHinge2D::PMInteractionHinge2D(PMInteractionHingeParams p):p_(std::move(p)){
    for(double v:{p_.axial_stiffness,p_.hinge.Ke,p_.axial_preload,p_.return_tolerance,p_.mroz_outer_scale})
        if(!std::isfinite(v))throw std::invalid_argument("non-finite P-M hinge parameter");
    if(!(p_.axial_stiffness>0.0)||!(p_.hinge.Ke>0.0)||p_.max_return_iterations<3||
       !(p_.return_tolerance>0.0))throw std::invalid_argument("invalid P-M interaction hinge parameters");
    if(p_.surface_shape==PMInteractionSurfaceShape::TabulatedMomentCapacity){
        if(p_.axial_force_points.size()<2||p_.axial_force_points.size()!=p_.moment_capacity_points.size())
            throw std::invalid_argument("invalid tabulated P-M interaction surface");
        for(double q:p_.axial_force_points)if(!std::isfinite(q))throw std::invalid_argument("non-finite P-M axial coordinate");
        for(double m:p_.moment_capacity_points)if(!std::isfinite(m)||!(m>0.0))throw std::invalid_argument("P-M moment capacities must be finite and positive");
        capacity_slopes_=monotone_slopes(p_.axial_force_points,p_.moment_capacity_points);
    }else{
        for(double v:{p_.py_tension,p_.p_balance,p_.py_compression,p_.my_balance,p_.alpha_tension,p_.alpha_compression,p_.beta_pm})
            if(!std::isfinite(v))throw std::invalid_argument("non-finite P-M surface parameter");
        if(!(p_.py_tension<p_.p_balance&&p_.p_balance<p_.py_compression)||!(p_.my_balance>0.0)||
           !(p_.alpha_tension>=1.0)||!(p_.alpha_compression>=1.0)||!(p_.beta_pm>1.0))
            throw std::invalid_argument("PERFORM concrete projection requires finite convex surface: alpha >= 1, beta > 1, ordered axial intercepts");
    }
    // Reuse ASCE41 validation for the directional deformation parameters.
    (void)ASCE41HingeMaterial(p_.hinge);
}

void PMInteractionHinge2D::initialize_state(double* s) const{
    std::fill(s,s+kStateSize,0.0);
    s[8]=base_capacity(p_.axial_preload).value;
}

PMInteractionHinge2D::CapacityEval PMInteractionHinge2D::base_capacity(double q) const{
    if(p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete){
        const bool compression_branch=q>=p_.p_balance;
        const double py=compression_branch?p_.py_compression:p_.py_tension;
        const double alpha=compression_branch?p_.alpha_compression:p_.alpha_tension;
        const double den=py-p_.p_balance;
        const double z=(q-p_.p_balance)/den; // nonnegative within either active branch
        const double r=std::abs(z);
        if(r>=1.0)return {1e-12,0.0,0.0};
        const double beta=p_.beta_pm;
        const double rr=std::max(0.0,1.0-std::pow(r,alpha));
        if(rr<=1e-16)return {1e-12,0.0,0.0};
        const double n=1.0/beta;
        const double v=p_.my_balance*std::pow(rr,n);
        // q lies between PB and the selected intercept, so r=(q-PB)/den
        // is positive despite den changing sign on the tension branch.
        const double drdq=1.0/den;
        const double A=-alpha*std::pow(std::max(r,1e-14),alpha-1.0)*drdq; // d(rr)/dq
        double d1=p_.my_balance*n*std::pow(rr,n-1.0)*A;
        double dA=-alpha*(alpha-1.0)*std::pow(std::max(r,1e-14),alpha-2.0)*drdq*drdq;
        double d2=p_.my_balance*n*((n-1.0)*std::pow(rr,n-2.0)*A*A+std::pow(rr,n-1.0)*dA);
        return {std::max(1e-12,v),d1,d2};
    }
    const auto& x=p_.axial_force_points;const auto& y=p_.moment_capacity_points;const auto& m=capacity_slopes_;
    if(q<=x.front())return {y.front(),0.0,0.0};
    if(q>=x.back())return {y.back(),0.0,0.0};
    auto it=std::upper_bound(x.begin(),x.end(),q);const std::size_t i=static_cast<std::size_t>(it-x.begin()-1);
    const double h=x[i+1]-x[i],t=(q-x[i])/h,t2=t*t,t3=t2*t;
    const double h00=2*t3-3*t2+1,h10=t3-2*t2+t,h01=-2*t3+3*t2,h11=t3-t2;
    const double v=h00*y[i]+h10*h*m[i]+h01*y[i+1]+h11*h*m[i+1];
    const double d1=((6*t2-6*t)/h)*y[i]+(3*t2-4*t+1)*m[i]+((-6*t2+6*t)/h)*y[i+1]+(3*t2-2*t)*m[i+1];
    const double d2=((12*t-6)/(h*h))*y[i]+((6*t-4)/h)*m[i]+((-12*t+6)/(h*h))*y[i+1]+((6*t-2)/h)*m[i+1];
    return {std::max(1e-12,v),d1,d2};
}


PMInteractionHinge2D::CapacityEval PMInteractionHinge2D::base_capacity_scaled(double q,double scale) const{
    if(!(scale>0.0)) throw std::invalid_argument("P-M surface scale must be positive");
    if(p_.surface_shape!=PMInteractionSurfaceShape::PerformConcrete){
        const auto b=base_capacity(q);
        return {scale*b.value,scale*b.d1,scale*b.d2};
    }
    const bool compression_branch=q>=p_.p_balance;
    const double py=compression_branch?p_.py_compression:p_.py_tension;
    const double alpha=compression_branch?p_.alpha_compression:p_.alpha_tension;
    const double den=scale*(py-p_.p_balance);
    const double z=(q-p_.p_balance)/den;
    const double r=std::abs(z);
    if(r>=1.0)return {1e-12,0.0,0.0};
    const double beta=p_.beta_pm;
    const double rr=std::max(0.0,1.0-std::pow(r,alpha));
    if(rr<=1e-16)return {1e-12,0.0,0.0};
    const double n=1.0/beta;
    const double v=scale*p_.my_balance*std::pow(rr,n);
    const double drdq=1.0/den;
    const double A=-alpha*std::pow(std::max(r,1e-14),alpha-1.0)*drdq;
    const double d1=scale*p_.my_balance*n*std::pow(rr,n-1.0)*A;
    const double dA=-alpha*(alpha-1.0)*std::pow(std::max(r,1e-14),alpha-2.0)*drdq*drdq;
    const double d2=scale*p_.my_balance*n*((n-1.0)*std::pow(rr,n-2.0)*A*A+std::pow(rr,n-1.0)*dA);
    return {std::max(1e-12,v),d1,d2};
}

PMInteractionHinge2D::BackboneEval PMInteractionHinge2D::directional_capacity(double P,double kappa,bool pos) const{
    const auto bcap=base_capacity(P);
    const double fyref=pos?p_.hinge.posFy:p_.hinge.negFy;
    const double mcref=pos?p_.hinge.posMc:p_.hinge.negMc;
    const double a=pos?p_.hinge.pos_a:p_.hinge.neg_a;
    const double b=pos?p_.hinge.pos_b:p_.hinge.neg_b;
    const double c=pos?p_.hinge.pos_c:p_.hinge.neg_c;
    const double drop_req=pos?p_.hinge.pos_drop_span:p_.hinge.neg_drop_span;
    const double drop=drop_req>0.0?drop_req:0.05*(b-a);
    const double e_req=pos?p_.hinge.pos_e_drop_span:p_.hinge.neg_e_drop_span;
    const double e_drop=e_req>0.0?e_req:0.05*(b-a-drop);
    const double kh=(p_.hinge.hardening_stiffness>0.0)?p_.hinge.hardening_stiffness:p_.hinge.hardening_ratio*p_.hinge.Ke;

    BackboneEval out;
    auto scale_base=[&](double r,double drdk){
        out.value=bcap.value*r;out.dP=bcap.d1*r;out.dPP=bcap.d2*r;out.dKappa=bcap.value*drdk;out.dPdKappa=bcap.d1*drdk;
    };

    const double x=std::max(0.0,kappa);
    if(p_.surface_evolution==PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic){
        out.value=bcap.value;out.dP=bcap.d1;out.dPP=bcap.d2;
        out.dKappa=0.0;out.dPdKappa=0.0;
        if(x>tiny)out.events|=PM_EVENT_YIELD;
        if(x>=a)out.events|=PM_EVENT_CAP;
        if(p_.enforce_deformation_capacity && x>=b){
            out.events|=PM_EVENT_LATERAL_LOSS;out.lateral_loss=true;out.value=0.0;
            out.dP=out.dPP=out.dKappa=out.dPdKappa=0.0;
        }
        return out;
    }
    if(x>=b){out.events|=PM_EVENT_LATERAL_LOSS;out.lateral_loss=true;return out;}
    if(x<=a){
        out.events|=(x>tiny?PM_EVENT_YIELD:PM_EVENT_NONE);
        if(mcref>0.0&&fyref>0.0){
            const double rc=mcref/fyref;const double r=1.0+(rc-1.0)*(x/a);scale_base(r,(rc-1.0)/a);
        }else{
            out.value=bcap.value+kh*x;out.dP=bcap.d1;out.dPP=bcap.d2;out.dKappa=kh;out.dPdKappa=0.0;
        }
        return out;
    }

    const double mc_scale=(mcref>0.0&&fyref>0.0)?(mcref/fyref):1.0;
    const double mc=bcap.value*mc_scale + ((mcref>0.0)?0.0:kh*a);
    const double mcP=bcap.d1*mc_scale,mcPP=bcap.d2*mc_scale;
    out.events|=PM_EVENT_YIELD|PM_EVENT_CAP;
    if(p_.hinge.backbone_shape==ASCE41BackboneShape::StraightCE){
        const double r=(b-x)/(b-a);out.value=std::max(0.0,mc*r);out.dP=mcP*r;out.dPP=mcPP*r;
        out.dKappa=-mc/(b-a);out.dPdKappa=-mcP/(b-a);out.events|=PM_EVENT_STRENGTH_DROP;return out;
    }
    const double kd_end=a+drop;const double e0=b-e_drop;const double fr=c*bcap.value;
    if(x<=kd_end){
        const double t=(x-a)/drop;out.value=mc+(fr-mc)*t;out.dP=mcP+(c*bcap.d1-mcP)*t;out.dPP=mcPP+(c*bcap.d2-mcPP)*t;
        out.dKappa=(fr-mc)/drop;out.dPdKappa=(c*bcap.d1-mcP)/drop;out.events|=PM_EVENT_STRENGTH_DROP;return out;
    }
    if(x<=e0){scale_base(c,0.0);out.events|=PM_EVENT_STRENGTH_DROP;return out;}
    const double r=c*(b-x)/e_drop;scale_base(std::max(0.0,r),-c/e_drop);out.events|=PM_EVENT_STRENGTH_DROP;return out;
}


PMInteractionTrialResult PMInteractionHinge2D::trial_mroz(double da,double th,const double* c,double* s,bool compute_tangent) const{
    enum : int { EPA=0,EPR=1,KPOS=2,KNEG=3,LN=4,LM=5,LAT=6,FAIL=7,CAP=8,CP=9,CM=10,OUTER=11 };
    std::copy(c,c+kStateSize,s);
    PMInteractionTrialResult out;
    const double Ka=p_.axial_stiffness,Kr=p_.hinge.Ke;
    const double Ntr=Ka*(da-c[EPA]);
    const double Mtr=Kr*(th-c[EPR]);
    const double Ptr=p_.axial_preload-Ntr;
    if(c[FAIL]>0.5){out.axial_force=Ntr;out.compression=Ptr;out.moment=0.0;out.tangent[0]=Ka;out.events=PM_EVENT_LATERAL_LOSS|PM_EVENT_FAILURE;s[LN]=Ntr;s[LM]=0.0;return out;}
    if(c[LAT]>0.5){out.axial_force=Ntr;out.compression=Ptr;out.moment=0.0;out.tangent[0]=Ka;out.events=PM_EVENT_LATERAL_LOSS;s[LN]=Ntr;s[LM]=0.0;return out;}
    if(p_.surface_shape!=PMInteractionSurfaceShape::PerformConcrete){out.converged=false;out.events=PM_EVENT_RETURN_FAIL;return out;}

    auto outer_scale=[&](bool pos){
        if(p_.mroz_outer_scale>1.0) return p_.mroz_outer_scale;
        const double fy=pos?p_.hinge.posFy:p_.hinge.negFy;
        const double mc=pos?p_.hinge.posMc:p_.hinge.negMc;
        const double a=pos?p_.hinge.pos_a:p_.hinge.neg_a;
        const double kh=(p_.hinge.hardening_stiffness>0.0)?p_.hinge.hardening_stiffness:p_.hinge.hardening_ratio*p_.hinge.Ke;
        if(fy<=0.0) return 1.0;
        const double u=(mc>0.0)?mc:(fy+std::max(0.0,kh)*std::max(0.0,a));
        return std::max(1.0,u/fy);
    };
    struct FR { bool ok{}; double P{},lam{},N{},M{},cap{},d1{},d2{},sm{}; bool pos{}; bool axial_tip{}; };
    auto fixed_return=[&](double centerP,double centerM,double scale)->FR{
        FR r;
        const double mrel=Mtr-centerM;
        r.sm=(std::abs(mrel)>1e-14)?sgn(mrel):sgn(th-c[EPR]);r.pos=r.sm>0.0;
        const double mabs=r.sm*mrel;
        auto ce0=base_capacity_scaled(Ptr-centerP,scale);
        const double ftol=p_.return_tolerance*std::max({1.0,std::abs(mrel),ce0.value});
        const double lower=p_.p_balance+scale*(p_.py_tension-p_.p_balance)+centerP;
        const double upper=p_.p_balance+scale*(p_.py_compression-p_.p_balance)+centerP;
        if(mabs-ce0.value<=ftol && Ptr>lower && Ptr<upper){r.ok=true;r.P=Ptr;r.lam=0.0;r.N=Ntr;r.M=Mtr;r.cap=ce0.value;r.d1=ce0.d1;r.d2=ce0.d2;return r;}
        auto eval=[&](double pp,double& lam,double& g,double& dg,CapacityEval& ce){
            ce=base_capacity_scaled(pp-centerP,scale);
            lam=(mabs-ce.value)/Kr;
            g=p_.axial_preload-pp-Ntr+Ka*lam*ce.d1;
            dg=-1.0+Ka*((-ce.d1/Kr)*ce.d1+lam*ce.d2);
        };
        double P=Ptr,lam=0,g=0,dg=0;CapacityEval ce;eval(P,lam,g,dg,ce);
        const double pyT=p_.p_balance+scale*(p_.py_tension-p_.p_balance)+centerP;
        const double pyC=p_.p_balance+scale*(p_.py_compression-p_.p_balance)+centerP;
        const double fscale=std::max({pyC-pyT,std::abs(Ptr),std::abs(p_.axial_preload)});
        if(mabs==0.0 && (Ptr<=pyT||Ptr>=pyC)){
            r.ok=true;r.axial_tip=true;r.P=Ptr<=pyT?pyT:pyC;r.N=p_.axial_preload-r.P;r.M=Mtr;r.cap=0.0;return r;
        }
        bool conv=false;
        for(int it=0;it<p_.max_return_iterations;++it){
            if(P>pyT&&P<pyC&&lam>=-p_.return_tolerance&&std::abs(g)<=p_.return_tolerance*fscale){lam=std::max(0.0,lam);conv=true;break;}
            if(!std::isfinite(dg)||std::abs(dg)<1e-16)break;
            const double dP=-g/dg;bool accepted=false;double fac=1.0;
            for(int ls=0;ls<18;++ls){
                const double pc=std::clamp(P+fac*dP,pyT+1e-10,pyC-1e-10);double lc,gc,dgc;CapacityEval cc;eval(pc,lc,gc,dgc,cc);
                if(lc>=-1e-12&&std::isfinite(gc)&&std::abs(gc)<std::abs(g)){P=pc;lam=std::max(0.0,lc);g=gc;dg=dgc;ce=cc;accepted=true;break;}fac*=0.5;
            }
            if(!accepted)break;
        }
        if(!conv){
            auto residual=[&](double pp){eval(pp,lam,g,dg,ce);return g;};
            conv=axial_coordinate_projection(Ptr,mabs,p_.p_balance+centerP,pyT,pyC,
                scale*p_.my_balance,p_.alpha_tension,p_.alpha_compression,p_.beta_pm,
                p_.return_tolerance*fscale,residual,P);
            if(conv){eval(P,lam,g,dg,ce);lam=std::max(0.0,lam);}
        }
        if(!conv){
            const auto tip=bracket_projection(Ptr,mabs,p_.p_balance+centerP,pyT,pyC,
                scale*p_.my_balance,p_.alpha_tension,p_.alpha_compression,p_.beta_pm,
                Ka,Kr,p_.return_tolerance);
            if(tip.ok){r.ok=true;r.P=tip.P;r.N=p_.axial_preload-tip.P;r.M=centerM+r.sm*tip.M;
                r.cap=tip.M;r.lam=std::max(0.0,(mabs-tip.M)/Kr);return r;}
        }
        if(!conv)return r;
        r.ok=true;r.P=P;r.lam=lam;r.N=p_.axial_preload-P;r.cap=ce.value;r.d1=ce.d1;r.d2=ce.d2;r.M=centerM+r.sm*ce.value;return r;
    };

    const bool was_outer=c[OUTER]>0.5;
    double cp=was_outer?0.0:c[CP], cm=was_outer?0.0:c[CM];
    const double relM=Mtr-cm;const bool pred_pos=((std::abs(relM)>1e-14)?relM:(th-c[EPR]))>=0.0;
    const double us=outer_scale(pred_pos);
    if(was_outer||us<=1.0+1e-10){
        auto fr=fixed_return(0.0,0.0,us);if(!fr.ok){out.converged=false;out.events=PM_EVENT_RETURN_FAIL;return out;}
        const double epA=da-fr.N/Ka,epR=th-fr.M/Kr;
        s[EPA]=epA;s[EPR]=epR;if(fr.pos)s[KPOS]=c[KPOS]+fr.lam;else s[KNEG]=c[KNEG]+fr.lam;s[LN]=fr.N;s[LM]=fr.M;s[CAP]=fr.cap;s[OUTER]=1.0;s[CP]=s[CM]=0.0;
        out.axial_force=fr.N;out.compression=fr.P;out.moment=fr.M;out.active_capacity=fr.cap;out.events=((fr.lam>0||fr.axial_tip)?PM_EVENT_YIELD:PM_EVENT_ELASTIC)|PM_EVENT_CAP;out.plastic_axial_deformation=epA;out.plastic_rotation=epR;out.positive_plastic_rotation=s[KPOS];out.negative_plastic_rotation=s[KNEG];
    }else{
        FR fr;double cpnew=cp,cmnew=cm;bool hit=false;
        const double pmetric=p_.legacy_mroz_mixed_unit_metric?1.0:(p_.py_compression-p_.py_tension);
        const double mmetric=p_.legacy_mroz_mixed_unit_metric?1.0:p_.my_balance;
        auto distance=[&](double dp,double dm){return std::hypot(dp/pmetric,dm/mmetric);};
        for(int hi=0;hi<8;++hi){
            fr=fixed_return(cpnew,cmnew,1.0);if(!fr.ok){out.converged=false;out.events=PM_EVENT_RETURN_FAIL;return out;}if(fr.lam<=p_.return_tolerance)break;
            const double pl=fr.P-cpnew;const auto cy=base_capacity_scaled(pl,1.0);
            const double targetP=(us-1.0)*(pl-p_.p_balance);
            const double targetM=fr.sm*(us-1.0)*cy.value;
            const double vp=targetP-c[CP],vm=targetM-c[CM];const double dist=distance(vp,vm);const double dref=std::max(1e-12,distance(targetP,targetM));
            const double a=fr.pos?p_.hinge.pos_a:p_.hinge.neg_a;const double move=(a>1e-12)?fr.lam/a*dref:dref;
            const double fac=(dist>1e-14)?std::min(1.0,move/dist):1.0;
            const double np=c[CP]+fac*vp,nm=c[CM]+fac*vm;
            hit=fac>=1.0-1e-8;
            if(distance(np-cpnew,nm-cmnew)<=1e-10*std::max(1.0,dref)){cpnew=np;cmnew=nm;break;}
            cpnew=np;cmnew=nm;
        }
        fr=fixed_return(cpnew,cmnew,1.0);if(!fr.ok){out.converged=false;out.events=PM_EVENT_RETURN_FAIL;return out;}
        if(hit){auto fu=fixed_return(0.0,0.0,us);if(fu.ok){fr=fu;cpnew=cmnew=0.0;s[OUTER]=1.0;}}
        const double epA=da-fr.N/Ka,epR=th-fr.M/Kr;
        s[EPA]=epA;s[EPR]=epR;if(fr.pos)s[KPOS]=c[KPOS]+fr.lam;else s[KNEG]=c[KNEG]+fr.lam;s[LN]=fr.N;s[LM]=fr.M;s[CAP]=fr.cap;s[CP]=cpnew;s[CM]=cmnew;
        out.axial_force=fr.N;out.compression=fr.P;out.moment=fr.M;out.active_capacity=fr.cap;out.events=((fr.lam>0||fr.axial_tip)?PM_EVENT_YIELD:PM_EVENT_ELASTIC)|(s[OUTER]>0.5?PM_EVENT_CAP:PM_EVENT_NONE);out.plastic_axial_deformation=epA;out.plastic_rotation=epR;out.positive_plastic_rotation=s[KPOS];out.negative_plastic_rotation=s[KNEG];
    }
    if(p_.enforce_deformation_capacity){const bool pos=out.moment>=0.0;const double b=pos?p_.hinge.pos_b:p_.hinge.neg_b;if((pos?s[KPOS]:s[KNEG])>=b-1e-12){s[LAT]=1.0;out.events|=PM_EVENT_LATERAL_LOSS;}}
    if(!compute_tangent){out.tangent[0]=Ka;out.tangent[3]=Kr;return out;}
    if((out.events&PM_EVENT_ELASTIC)!=0){out.tangent[0]=Ka;out.tangent[3]=Kr;return out;}
    // Perturbations must scale with physical force/deformation scales, rather
    // than a hidden assumption that the unit of length is one inch.
    const double ha=p_.legacy_mroz_mixed_unit_metric?std::max(1e-8,2e-6*std::max(1.0,std::abs(da))):
        2e-6*std::max(std::abs(da),(p_.py_compression-p_.py_tension)/Ka);
    const double ht=p_.legacy_mroz_mixed_unit_metric?std::max(1e-8,2e-6*std::max(1.0,std::abs(th))):
        2e-6*std::max(std::abs(th),p_.my_balance/Kr);
    auto fd=[&](double dap,double thp){std::vector<double> st(static_cast<std::size_t>(kStateSize));return trial_mroz(dap,thp,c,st.data(),false);};
    auto ap=fd(da+ha,th),am=fd(da-ha,th),tp=fd(da,th+ht),tm=fd(da,th-ht);
    if(!ap.converged||!am.converged||!tp.converged||!tm.converged){out.converged=false;out.events|=PM_EVENT_RETURN_FAIL;return out;}
    out.tangent[0]=(ap.axial_force-am.axial_force)/(2*ha);out.tangent[2]=(ap.moment-am.moment)/(2*ha);out.tangent[1]=(tp.axial_force-tm.axial_force)/(2*ht);out.tangent[3]=(tp.moment-tm.moment)/(2*ht);
    return out;
}

PMInteractionTrialResult PMInteractionHinge2D::trial(double da,double th,const double* c,double* s) const{
    if(p_.surface_evolution==PMInteractionSurfaceEvolution::MrozTwoSurface)
        return trial_mroz(da,th,c,s,true);
    enum : int { EPA=0,EPR=1,KPOS=2,KNEG=3,LN=4,LM=5,LAT=6,FAIL=7,CAP=8 };
    std::copy(c,c+kStateSize,s);
    PMInteractionTrialResult out;
    const double Ka=p_.axial_stiffness,Kr=p_.hinge.Ke;
    const double Ntr=Ka*(da-c[EPA]);
    const double Mtr=Kr*(th-c[EPR]);
    const double Ptr=p_.axial_preload-Ntr;

    if(c[FAIL]>0.5){
        out.axial_force=Ntr;out.compression=Ptr;out.moment=0.0;out.tangent[0]=Ka;out.events=PM_EVENT_LATERAL_LOSS|PM_EVENT_FAILURE;
        s[LN]=out.axial_force;s[LM]=0.0;return out;
    }
    if(c[LAT]>0.5){
        out.axial_force=Ntr;out.compression=Ptr;out.moment=0.0;out.tangent[0]=Ka;out.events=PM_EVENT_LATERAL_LOSS;
        const double p_rot=th-c[EPR];
        const double fpos=p_.hinge.pos_f>0.0?p_.hinge.pos_f:p_.hinge.pos_b;
        const double fneg=p_.hinge.neg_f>0.0?p_.hinge.neg_f:p_.hinge.neg_b;
        if((p_rot>=0.0&&std::abs(p_rot)>=fpos)||(p_rot<0.0&&std::abs(p_rot)>=fneg)){s[FAIL]=1.0;out.events|=PM_EVENT_FAILURE;}
        s[LN]=out.axial_force;s[LM]=0.0;return out;
    }

    if(p_.surface_evolution==PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic &&
       p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete && Mtr==0.0 &&
       (Ptr<p_.py_tension||Ptr>p_.py_compression)){
        const bool compression=Ptr>p_.py_compression;
        const double P=compression?p_.py_compression:p_.py_tension;
        const double alpha=compression?p_.alpha_compression:p_.alpha_tension;
        const double normal=alpha/std::abs(P-p_.p_balance);
        const double multiplier=std::abs(Ptr-P)/(Ka*normal);
        out.axial_force=p_.axial_preload-P;out.compression=P;out.moment=0.0;
        out.plastic_axial_deformation=da-out.axial_force/Ka;out.plastic_rotation=c[EPR];
        out.positive_plastic_rotation=c[KPOS];out.negative_plastic_rotation=c[KNEG];out.events=PM_EVENT_YIELD;
        out.tangent[3]=p_.beta_pm<2.0?0.0:(p_.beta_pm==2.0?Kr/(1.0+Kr*multiplier*2.0/(p_.my_balance*p_.my_balance)):Kr);
        s[EPA]=out.plastic_axial_deformation;s[LN]=out.axial_force;s[LM]=0.0;s[CAP]=0.0;return out;
    }

    const double sign_hint=(std::abs(Mtr)>1e-14)?sgn(Mtr):sgn(th-c[EPR]);
    const bool pos=sign_hint>0.0;const double kn=pos?c[KPOS]:c[KNEG];
    auto bt=directional_capacity(Ptr,kn,pos);
    const double ftr=std::abs(Mtr)-bt.value;
    const double ftol=p_.return_tolerance*std::max({1.0,std::abs(Mtr),bt.value});
    const bool inside_axial_domain=p_.surface_shape!=PMInteractionSurfaceShape::PerformConcrete ||
        (Ptr>p_.py_tension&&Ptr<p_.py_compression);
    if(ftr<=ftol && !bt.lateral_loss && inside_axial_domain){
        out.axial_force=Ntr;out.compression=Ptr;out.moment=Mtr;out.tangent[0]=Ka;out.tangent[3]=Kr;out.active_capacity=bt.value;out.events=PM_EVENT_ELASTIC|bt.events;
        out.plastic_axial_deformation=c[EPA];out.plastic_rotation=c[EPR];out.positive_plastic_rotation=c[KPOS];out.negative_plastic_rotation=c[KNEG];
        s[LN]=Ntr;s[LM]=Mtr;s[CAP]=bt.value;return out;
    }

    const double sm=sign_hint;
    double lam=std::max(0.0,ftr/(Kr+std::max(0.0,std::abs(bt.dKappa))));
    double P=Ptr;bool conv=false;BackboneEval be=bt;
    double J00=0,J01=0,J10=0,J11=0;

    if(p_.surface_evolution==PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic){
        // For a fixed P-M surface, eliminate lambda using rotational
        // consistency and solve the closest-point condition as a scalar
        // equation in compression P. This is mathematically equivalent to
        // the two-equation return map but is much better conditioned when the
        // rigid-plastic axial hinge uses a very large penalty stiffness.
        const double mabs=sm*Mtr;
        auto eval_scalar=[&](double pp,BackboneEval& bb,double& ll,double& gg,double& dg){
            bb=directional_capacity(pp,kn,pos);
            ll=(mabs-bb.value)/Kr;
            gg=p_.axial_preload-pp-Ntr+Ka*ll*bb.dP;
            dg=-1.0+Ka*((-bb.dP/Kr)*bb.dP+ll*bb.dPP);
        };
        double g=0.0,dg=0.0;eval_scalar(P,be,lam,g,dg);
        const double axial_span=p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete?p_.py_compression-p_.py_tension:p_.axial_force_points.back()-p_.axial_force_points.front();
        const double force_scale=std::max({axial_span,std::abs(Ptr),std::abs(p_.axial_preload),1e-15});
        for(int it=0;it<p_.max_return_iterations;++it){
            const bool admissible=p_.surface_shape!=PMInteractionSurfaceShape::PerformConcrete ||
                (P>p_.py_tension&&P<p_.py_compression);
            if(admissible && lam>=-p_.return_tolerance && std::abs(g)<=p_.return_tolerance*force_scale){
                lam=std::max(0.0,lam);conv=true;break;
            }
            if(!std::isfinite(dg)||std::abs(dg)<1e-16)break;
            const double stepP=-g/dg;
            bool accepted=false;double fac=1.0;
            for(int ls=0;ls<18;++ls){
                const double pc=P+fac*stepP;BackboneEval bc;double lc=0,gc=0,dgc=0;
                eval_scalar(pc,bc,lc,gc,dgc);
                if(lc>=-1e-12 && std::isfinite(gc) && std::abs(gc)<std::abs(g)){
                    P=pc;be=bc;lam=std::max(0.0,lc);g=gc;dg=dgc;accepted=true;break;
                }
                fac*=0.5;
            }
            if(!accepted)break;
        }
        if(!conv && p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete){
            auto residual=[&](double pp){eval_scalar(pp,be,lam,g,dg);return g;};
            conv=axial_coordinate_projection(Ptr,mabs,p_.p_balance,p_.py_tension,p_.py_compression,
                p_.my_balance,p_.alpha_tension,p_.alpha_compression,p_.beta_pm,
                p_.return_tolerance*force_scale,residual,P);
            if(conv){eval_scalar(P,be,lam,g,dg);lam=std::max(0.0,lam);}
        }
        if(!conv && p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete){
            const auto tip=bracket_projection(Ptr,mabs,p_.p_balance,p_.py_tension,p_.py_compression,
                p_.my_balance,p_.alpha_tension,p_.alpha_compression,p_.beta_pm,
                Ka,Kr,p_.return_tolerance);
            if(tip.ok){
                lam=std::max(0.0,(mabs-tip.M)/Kr);
                out.axial_force=p_.axial_preload-tip.P;out.compression=tip.P;out.moment=sm*tip.M;
                out.plastic_axial_deformation=da-out.axial_force/Ka;out.plastic_rotation=th-out.moment/Kr;
                s[EPA]=out.plastic_axial_deformation;s[EPR]=out.plastic_rotation;
                if(pos)s[KPOS]=kn+lam;else s[KNEG]=kn+lam;
                s[LN]=out.axial_force;s[LM]=out.moment;s[CAP]=tip.M;
                out.positive_plastic_rotation=s[KPOS];out.negative_plastic_rotation=s[KNEG];
                out.active_capacity=tip.M;out.events=PM_EVENT_YIELD;
                const double H=(tip.dpdm*tip.dpdm+(tip.P-Ptr)*tip.d2pdm2)/Ka+1.0/Kr;
                if(std::isfinite(H)&&H>0.0){
                    out.tangent[0]=tip.dpdm*tip.dpdm/H;
                    out.tangent[1]=out.tangent[2]=-sm*tip.dpdm/H;out.tangent[3]=1.0/H;
                }
                const double a=pos?p_.hinge.pos_a:p_.hinge.neg_a,b=pos?p_.hinge.pos_b:p_.hinge.neg_b;
                if(kn+lam>=a)out.events|=PM_EVENT_CAP;
                if(p_.enforce_deformation_capacity&&kn+lam>=b-1e-12){s[LAT]=1.0;out.events|=PM_EVENT_LATERAL_LOSS;}
                return out;
            }
        }
        if(!conv && p_.surface_shape==PMInteractionSurfaceShape::TabulatedMomentCapacity){
            // Safeguarded fallback: locate a sign-changing stationary root
            // across the supplied interaction table, preferring the root
            // closest to the elastic trial point.
            const double lo=(p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete)?p_.py_tension:p_.axial_force_points.front();
            const double hi=(p_.surface_shape==PMInteractionSurfaceShape::PerformConcrete)?p_.py_compression:p_.axial_force_points.back();
            constexpr int nscan=256;
            bool have=false;double best_dist=std::numeric_limits<double>::infinity();
            double pa=lo;BackboneEval ba;double la=0,ga=0,dga=0;eval_scalar(pa,ba,la,ga,dga);
            for(int is=1;is<=nscan;++is){
                const double pb=lo+(hi-lo)*static_cast<double>(is)/nscan;
                BackboneEval bb;double lb=0,gb=0,dgb=0;eval_scalar(pb,bb,lb,gb,dgb);
                if(la>=-1e-12 && lb>=-1e-12 && ga*gb<=0.0){
                    double a0=pa,b0=pb,fa=ga;BackboneEval bm=bb;double lm=lb,gm=gb,dgm=dgb;
                    for(int bi=0;bi<70;++bi){
                        const double mid=0.5*(a0+b0);eval_scalar(mid,bm,lm,gm,dgm);
                        if(std::abs(gm)<=p_.return_tolerance*force_scale){a0=b0=mid;break;}
                        if(fa*gm<=0.0)b0=mid;else{a0=mid;fa=gm;}
                    }
                    const double root=0.5*(a0+b0);BackboneEval br;double lr=0,gr=0,dgr=0;eval_scalar(root,br,lr,gr,dgr);
                    const double dist=std::abs(root-Ptr);
                    if(lr>=-1e-10 && dist<best_dist){best_dist=dist;P=root;be=br;lam=std::max(0.0,lr);g=gr;dg=dgr;have=true;}
                }
                pa=pb;ba=bb;la=lb;ga=gb;dga=dgb;
            }
            if(have && std::abs(g)<=10.0*p_.return_tolerance*force_scale)conv=true;
        }
    }else{
        for(int it=0;it<p_.max_return_iterations;++it){
            const double kap=kn+lam;be=directional_capacity(P,kap,pos);
            if(be.lateral_loss){break;}
            const double r1=p_.axial_preload-P-Ntr+Ka*lam*be.dP;
            const double r2=sm*Mtr-Kr*lam-be.value;
            const double scale=std::max({1.0,std::abs(p_.axial_preload),std::abs(Ntr),std::abs(Mtr),be.value});
            if(std::max(std::abs(r1),std::abs(r2))<=p_.return_tolerance*scale){conv=true;break;}
            J00=-1.0+Ka*lam*be.dPP;
            J01=Ka*(be.dP+lam*be.dPdKappa);
            J10=-be.dP;
            J11=-Kr-be.dKappa;
            const double det=J00*J11-J01*J10;
            if(std::abs(det)<1e-18*std::max(1.0,std::abs(J00*J11)))break;
            const double dP=(-r1*J11+J01*r2)/det;
            const double dl=(-J00*r2+J10*r1)/det;
            double step=1.0;
            if(lam+dl<0.0)step=std::min(step,0.8*lam/std::max(1e-30,-dl));
            P+=step*dP;lam=std::max(0.0,lam+step*dl);
        }
    }

    if(!conv){
        const double kap=kn+lam;be=directional_capacity(P,kap,pos);
        if(be.lateral_loss){
            s[LAT]=1.0;out.events=PM_EVENT_LATERAL_LOSS|PM_EVENT_YIELD|PM_EVENT_CAP|PM_EVENT_STRENGTH_DROP;
            out.axial_force=Ntr;out.compression=Ptr;out.moment=0.0;out.tangent[0]=Ka;
            s[EPR]=th;s[LN]=Ntr;s[LM]=0.0;s[CAP]=0.0;return out;
        }
        out.converged=false;out.events=PM_EVENT_RETURN_FAIL;return out;
    }

    const double N=p_.axial_preload-P;const double M=sm*be.value;
    const double epA=c[EPA]+lam*be.dP;const double epR=c[EPR]+lam*sm;
    s[EPA]=epA;s[EPR]=epR;if(pos)s[KPOS]=kn+lam;else s[KNEG]=kn+lam;s[LN]=N;s[LM]=M;s[CAP]=be.value;
    out.axial_force=N;out.compression=P;out.moment=M;out.active_capacity=be.value;out.events=PM_EVENT_YIELD|be.events;
    out.plastic_axial_deformation=epA;out.plastic_rotation=epR;out.positive_plastic_rotation=s[KPOS];out.negative_plastic_rotation=s[KNEG];

    // Exact algorithmic 2x2 tangent from the converged local return equations.
    J00=-1.0+Ka*lam*be.dPP;J01=Ka*(be.dP+lam*be.dPdKappa);J10=-be.dP;J11=-Kr-be.dKappa;
    const double det=J00*J11-J01*J10;
    if(std::abs(det)<1e-20){out.converged=false;out.events|=PM_EVENT_RETURN_FAIL;return out;}
    auto column=[&](double b0,double b1,int j){
        // J dx + b = 0
        const double dP=(-b0*J11+J01*b1)/det;
        const double dl=(-J00*b1+J10*b0)/det;
        out.tangent[j]=-dP;
        out.tangent[2+j]=sm*(be.dP*dP+be.dKappa*dl);
    };
    column(-Ka,0.0,0);
    column(0.0,sm*Kr,1);

    const double bside=pos?p_.hinge.pos_b:p_.hinge.neg_b;
    if(p_.enforce_deformation_capacity && (pos?s[KPOS]:s[KNEG])>=bside-1e-12){s[LAT]=1.0;out.events|=PM_EVENT_LATERAL_LOSS;}
    return out;
}

} // namespace quake
