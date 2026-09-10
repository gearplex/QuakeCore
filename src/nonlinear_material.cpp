#include "quake/nonlinear_material.hpp"

#include <algorithm>
#include <cmath>

namespace quake {

NonlinearMaterial::NonlinearMaterial(BilinearSpring law) : law_(std::move(law)) {}
NonlinearMaterial::NonlinearMaterial(IMKPeakOrientedMaterial law) : law_(std::move(law)) {}
NonlinearMaterial::NonlinearMaterial(ASCE41HingeMaterial law) : law_(std::move(law)) {}
NonlinearMaterial::NonlinearMaterial(AsymmetricElasticPerfectlyPlasticSpring law) : law_(std::move(law)) {}

NonlinearMaterialKind NonlinearMaterial::kind() const {
    if (std::holds_alternative<BilinearSpring>(law_)) return NonlinearMaterialKind::Bilinear;
    if (std::holds_alternative<IMKPeakOrientedMaterial>(law_)) return NonlinearMaterialKind::IMKPeakOriented;
    if (std::holds_alternative<ASCE41HingeMaterial>(law_)) return NonlinearMaterialKind::ASCE41Hinge;
    return NonlinearMaterialKind::AsymmetricElasticPerfectlyPlastic;
}
int NonlinearMaterial::state_size() const {
    if (std::holds_alternative<BilinearSpring>(law_)) return 2;
    if (std::holds_alternative<IMKPeakOrientedMaterial>(law_)) return IMKPeakOrientedMaterial::kStateSize;
    if (std::holds_alternative<ASCE41HingeMaterial>(law_)) return ASCE41HingeMaterial::kStateSize;
    return 1;
}
double NonlinearMaterial::initial_stiffness() const {
    if (auto* b = std::get_if<BilinearSpring>(&law_)) return b->initial_stiffness();
    if (auto* imk = std::get_if<IMKPeakOrientedMaterial>(&law_)) return imk->initial_stiffness();
    if (auto* h=std::get_if<ASCE41HingeMaterial>(&law_)) return h->initial_stiffness();
    return std::get<AsymmetricElasticPerfectlyPlasticSpring>(law_).initial_stiffness();
}
const ASCE41HingeParams* NonlinearMaterial::asce41_params() const {
    const auto* h=std::get_if<ASCE41HingeMaterial>(&law_);
    return h?&h->params():nullptr;
}
void NonlinearMaterial::initialize_state(double* s) const {
    if (auto* b = std::get_if<BilinearSpring>(&law_)) {
        (void)b; s[0]=0.0; s[1]=0.0;
    } else if (auto* imk=std::get_if<IMKPeakOrientedMaterial>(&law_)) imk->initialize_state(s);
    else if (auto* h=std::get_if<ASCE41HingeMaterial>(&law_)) h->initialize_state(s);
    else s[0]=0.0;
}
MaterialTrialResult NonlinearMaterial::trial(double q, const double* c, double* s) const {
    if (auto* b = std::get_if<BilinearSpring>(&law_)) {
        const BilinearState bs{c[0],c[1]};
        const auto tr=b->trial(q,bs);
        s[0]=tr.state.plastic; s[1]=tr.state.backstress;
        MaterialEvalDiagnostics d;
        d.tangent_active=std::abs(tr.tangent-b->initial_stiffness())>1e-14;
        d.fast_path=!d.tangent_active;
        d.transition=d.tangent_active != (std::abs(c[0])>1e-14 || std::abs(c[1])>1e-14);
        return {tr.force,tr.tangent,d};
    }
    if (auto* imkptr=std::get_if<IMKPeakOrientedMaterial>(&law_)) {
        const auto& imk=*imkptr;
        const auto tr=imk.trial(q,c,s);
        MaterialEvalDiagnostics d;
        d.fast_path=(tr.events & IMK_EVENT_FAST_ELASTIC)!=0;
        d.tangent_active=std::abs(tr.tangent-imk.initial_stiffness())>1e-14;
        d.reversal=(tr.events & IMK_EVENT_REVERSAL)!=0;
        const auto& p=imk.params();
        const bool pos=q>=0.0; const double fy=pos?p.posFy:p.negFy; const double up=pos?p.posUp:p.negUp;
        const double uy=fy/p.Ke, ucap=uy+up;
        const double x0=std::abs(c[0]), x1=std::abs(q);
        const bool yield_cross=(x0<=uy+1e-14 && x1>uy+1e-14);
        const bool cap_cross=(x0<=ucap+1e-14 && x1>ucap+1e-14);
        d.transition=yield_cross||cap_cross||d.reversal||(tr.events&IMK_EVENT_DETERIORATION)||(tr.events&IMK_EVENT_FAILURE);
        d.deterioration=(tr.events & IMK_EVENT_DETERIORATION)!=0;
        d.failed=(tr.events & IMK_EVENT_FAILURE)!=0;
        return {tr.force,tr.tangent,d};
    }
    if (const auto* a=std::get_if<AsymmetricElasticPerfectlyPlasticSpring>(&law_)) {
        const auto tr=a->trial(q,c[0]);s[0]=tr.plastic_deformation;
        MaterialEvalDiagnostics d;d.tangent_active=tr.tangent==0.0;d.fast_path=!d.tangent_active;
        d.transition=d.tangent_active!=(std::abs(c[0])>1e-14);return {tr.force,tr.tangent,d};
    }
    const auto& h=std::get<ASCE41HingeMaterial>(law_);
    const auto tr=h.trial(q,c,s);
    MaterialEvalDiagnostics d;
    d.fast_path=(tr.events&ASCE41_EVENT_FAST_ELASTIC)!=0;
    d.tangent_active=std::abs(tr.tangent-h.initial_stiffness())>1e-14;
    d.reversal=(tr.events&ASCE41_EVENT_REVERSAL)!=0;
    d.transition=(tr.events&(ASCE41_EVENT_YIELD|ASCE41_EVENT_CAP|ASCE41_EVENT_STRENGTH_DROP|ASCE41_EVENT_REVERSAL|ASCE41_EVENT_DETERIORATION|ASCE41_EVENT_FAILURE))!=0;
    d.deterioration=(tr.events&ASCE41_EVENT_DETERIORATION)!=0;
    d.failed=(tr.events&ASCE41_EVENT_FAILURE)!=0;
    d.lateral_resistance_lost=(tr.events&ASCE41_EVENT_LATERAL_LOSS)!=0;
    const auto level=h.performance_level(q);
    d.at_or_beyond_io=static_cast<int>(level)>=static_cast<int>(ASCE41PerformanceLevel::IO);
    d.at_or_beyond_ls=static_cast<int>(level)>=static_cast<int>(ASCE41PerformanceLevel::LS);
    d.at_or_beyond_cp=static_cast<int>(level)>=static_cast<int>(ASCE41PerformanceLevel::CP);
    d.beyond_cp=static_cast<int>(level)>=static_cast<int>(ASCE41PerformanceLevel::BeyondCP);
    return {tr.force,tr.tangent,d};
}

} // namespace quake
