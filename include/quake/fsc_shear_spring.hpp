#pragma once

#include "quake/nonlinear_material.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {

// Research flexure-shear-critical (FSC) column shear spring.
//
// The spring is intended to sit in series with a flexural column model.  It is
// elastic before a remote rotation/force limit is reached.  After initiation,
// it unloads/reloads elastically but is capped by a symmetric degrading shear
// envelope.  The post-failure slope and residual strength are explicit inputs;
// no benchmark EDP is embedded in this law.
//
// Axial force and end rotations are supplied by the parent frame at each trial.
// Axial capacity loss is deliberately NOT represented here; that remains a
// separate component/state so an FSC shear event cannot silently redistribute
// gravity load.
struct FSCShearSpringParams {
    double Ke{};                             // kip/in
    double post_failure_stiffness{};         // kip/in, <= 0
    double residual_strength_ratio{0.20};
    double nominal_shear_limit_kip{-1.0};    // <=0 disables force trigger

    // Ghannoum-Moehle total-plastic rotation initiation inputs.
    double b_in{6.0};
    double d_in{4.8};
    double h_in{6.0};
    double clear_length_in{39.0};
    double tie_spacing_in{4.0};
    double longitudinal_steel_area_in2{0.88};
    double confined_concrete_area_in2{18.0};
    double fc_ksi{3.57};
    double fy_ksi{64.0};

    // Cyclic envelope loss coefficient. <0 requests the calibrated expression
    // used by the source-inspired diagnostic; >=0 is a direct explicit input.
    double cyclic_strength_coefficient{-1.0};
    double sign_crossing_fraction{0.05};
};

struct FSCShearSpringTrial {
    MaterialTrialResult material;
    bool initiated{};
    bool residual_reached{};
    double rotation_limit{};
    double retained_strength_ratio{1.0};
};

class FSCShearSpringLaw {
public:
    // q, f, initiated, |q| at initiation, |V| at initiation,
    // cumulative cyclic decrement, retained ratio, lobe sign, lobe peak,
    // residual-reached flag.
    static constexpr int kStateSize = 10;

    explicit FSCShearSpringLaw(FSCShearSpringParams p) : p_(p) {
        if (!(p_.Ke > 0.0)) throw std::invalid_argument("FSC shear Ke must be positive");
        if (p_.post_failure_stiffness > 0.0) throw std::invalid_argument("FSC post-failure stiffness must be nonpositive");
        if (!(p_.b_in>0.0 && p_.d_in>0.0 && p_.h_in>0.0 && p_.clear_length_in>0.0 && p_.fc_ksi>0.0))
            throw std::invalid_argument("invalid FSC shear geometry/material parameters");
        p_.residual_strength_ratio=std::clamp(p_.residual_strength_ratio,0.0,1.0);
        if(p_.cyclic_strength_coefficient<0.0)p_.cyclic_strength_coefficient=calibrated_cyclic_strength_coefficient(p_);
        p_.cyclic_strength_coefficient=std::max(0.0,p_.cyclic_strength_coefficient);
    }

    const FSCShearSpringParams& params() const noexcept { return p_; }
    double initial_stiffness() const noexcept { return p_.Ke; }
    void initialize_state(double* s) const {
        for(int i=0;i<kStateSize;++i)s[i]=0.0;
        s[6]=1.0;
    }

    static double plastic_rotation_limit(const FSCShearSpringParams& p,
                                         double axial_compression_kip,
                                         double shear_kip) {
        const double P=std::max(0.0,axial_compression_kip);
        const double V=std::abs(shear_kip);
        const double Ag=p.b_in*p.h_in;
        const double fc_psi=p.fc_ksi*1000.0;
        const double v_psi=V/(p.b_in*p.d_in)*1000.0;
        return std::max(0.0,0.032
            -0.014*(p.tie_spacing_in/p.d_in)
            -0.017*(P/(Ag*p.fc_ksi))
            -0.0016*(v_psi/std::sqrt(fc_psi)));
    }

    static double calibrated_cyclic_strength_coefficient(const FSCShearSpringParams& p) {
        const double Ag=p.b_in*p.h_in;
        const double a=0.5*p.clear_length_in;
        return std::max(0.0,0.037133
            +0.251204*(p.fy_ksi*p.longitudinal_steel_area_in2/(p.fc_ksi*Ag))
            -0.354989*(p.confined_concrete_area_in2/Ag)
            +0.056569*(a/p.d_in));
    }

    FSCShearSpringTrial trial(double q,double local_rotation,double axial_compression_kip,
                              const double* c,double* s) const {
        enum : int { Q=0,F=1,INIT=2,QFAIL=3,VFAIL=4,CUMDEC=5,RETAIN=6,LOBESIGN=7,LOBEPEAK=8,RESID=9 };
        std::copy(c,c+kStateSize,s);
        const double elastic_force=p_.Ke*q;
        const double rot_lim=plastic_rotation_limit(p_,axial_compression_kip,elastic_force);
        const bool rot_hit=rot_lim>0.0 && std::abs(local_rotation)>=rot_lim;
        const bool force_hit=p_.nominal_shear_limit_kip>0.0 && std::abs(elastic_force)>=p_.nominal_shear_limit_kip;
        const bool was_initiated=c[INIT]>0.5;

        MaterialEvalDiagnostics d;
        if(!was_initiated){
            s[Q]=q;s[F]=elastic_force;s[RETAIN]=1.0;
            if(rot_hit||force_hit){
                s[INIT]=1.0;
                s[QFAIL]=std::abs(q);
                s[VFAIL]=std::max({std::abs(elastic_force),
                                    p_.nominal_shear_limit_kip>0.0?p_.nominal_shear_limit_kip:0.0,
                                    1e-12});
                const int sg=sign_with_deadband(elastic_force,s[VFAIL]);
                s[LOBESIGN]=static_cast<double>(sg);
                s[LOBEPEAK]=std::abs(elastic_force);
                d.transition=true;
                d.deterioration=true;
                d.lateral_resistance_lost=true;
            }else d.fast_path=true;
            d.tangent_active=false;
            return {{elastic_force,p_.Ke,d},s[INIT]>0.5,false,rot_lim,s[RETAIN]};
        }

        const double vfail=std::max(c[VFAIL],1e-12);
        const double vres=p_.residual_strength_ratio*vfail;
        double cumdec=std::max(0.0,c[CUMDEC]);
        int lobe_sign=static_cast<int>(std::lround(c[LOBESIGN]));
        double lobe_peak=std::max(0.0,c[LOBEPEAK]);

        // Elastic predictor from the last accepted point.  Cyclic degradation
        // is applied only when a completed excursion enters the opposite force
        // sign, so repeated Newton trials from one committed state cannot count
        // the same reversal more than once.
        double predictor=c[F]+p_.Ke*(q-c[Q]);
        int sg=sign_with_deadband(predictor,vfail);
        if(sg!=0 && lobe_sign!=0 && sg!=lobe_sign){
            cumdec += lobe_peak*p_.cyclic_strength_coefficient;
            lobe_sign=sg;
            lobe_peak=std::abs(predictor);
            d.transition=true;
            d.deterioration=true;
        }else if(sg!=0){
            if(lobe_sign==0)lobe_sign=sg;
            lobe_peak=std::max(lobe_peak,std::abs(predictor));
        }

        const double strength_cap=std::max(vres,vfail-cumdec);
        double deformation_cap=vfail;
        bool softening_controls=false;
        if(std::abs(q)>c[QFAIL]+1e-14){
            deformation_cap=std::max(vres,vfail+p_.post_failure_stiffness*(std::abs(q)-c[QFAIL]));
            softening_controls=deformation_cap<=strength_cap+1e-12;
        }
        const double cap=std::min(strength_cap,deformation_cap);
        double force=predictor,tangent=p_.Ke;
        if(std::abs(predictor)>cap){
            const int sf=predictor>=0.0?1:-1;
            force=sf*cap;
            if(softening_controls && deformation_cap>vres+1e-12 && sf*(q>=0.0?1.0:-1.0)>0.0)
                tangent=p_.post_failure_stiffness;
            else tangent=0.0;
        }

        sg=sign_with_deadband(force,vfail);
        if(sg!=0){
            if(lobe_sign==0)lobe_sign=sg;
            if(sg==lobe_sign)lobe_peak=std::max(lobe_peak,std::abs(force));
        }
        const double retained=cap/vfail;
        const bool resid=cap<=vres+1e-10*vfail;
        s[Q]=q;s[F]=force;s[INIT]=1.0;s[CUMDEC]=cumdec;s[RETAIN]=retained;
        s[LOBESIGN]=static_cast<double>(lobe_sign);s[LOBEPEAK]=lobe_peak;s[RESID]=resid?1.0:0.0;

        d.fast_path=false;
        d.tangent_active=std::abs(tangent-p_.Ke)>1e-14;
        d.deterioration=true;
        d.lateral_resistance_lost=true;
        return {{force,tangent,d},true,resid,rot_lim,retained};
    }

private:
    int sign_with_deadband(double force,double reference) const {
        const double eps=std::max(1e-9,p_.sign_crossing_fraction*std::max(reference,1e-12));
        if(force>eps)return 1;
        if(force<-eps)return -1;
        return 0;
    }
    FSCShearSpringParams p_{};
};

} // namespace quake
