#include "quake/rc_column_asce41.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace quake {
namespace {
void validate_resolved(const RCColumnResolvedParameters& r){
    if(r.numerical_hinge_Ke<=0.0||r.posMy<=0.0||r.negMy<=0.0||
       r.posMc<=0.0||r.negMc<=0.0||r.pos_a<=0.0||r.neg_a<=0.0||
       r.pos_b<=r.pos_a||r.neg_b<=r.neg_a||
       r.pos_c<0.0||r.pos_c>1.0||r.neg_c<0.0||r.neg_c>1.0)
        throw std::invalid_argument("invalid resolved RC column ASCE41 parameters");
}

void validate_section(const RCColumnSectionInput& s){
    if(s.width_in<=0.0||s.depth_in<=0.0||s.clear_length_in<=0.0||
       s.expected_fc_ksi<=0.0||s.expected_fy_long_ksi<=0.0||
       s.longitudinal_bar_count<=0||s.longitudinal_bar_area_in2<=0.0)
        throw std::invalid_argument("invalid RC column section input");
    if(s.gravity_axial_compression_kip<0.0)
        throw std::invalid_argument("RC column compression must be positive");
}

void validate_resolution(const RCColumnRuleResolution& r){
    validate_resolved(r.resolved);
    if(r.provenance.empty())
        throw std::invalid_argument("RC column rule resolution requires explicit provenance");
}

double rel_change(double a,double b){
    const double scale=std::max({1.0,std::abs(a),std::abs(b)});
    return std::abs(a-b)/scale;
}

double demand_change(const RCColumnDemandState& a,const RCColumnDemandState& b){
    return std::max({rel_change(a.max_compression_kip,b.max_compression_kip),
                     rel_change(a.max_tension_kip,b.max_tension_kip),
                     rel_change(a.max_abs_shear_kip,b.max_abs_shear_kip)});
}

double parameter_change(const RCColumnResolvedParameters& a,
                        const RCColumnResolvedParameters& b){
    const double av[]={a.numerical_hinge_Ke,a.posMy,a.negMy,a.posMc,a.negMc,
        a.pos_a,a.neg_a,a.pos_b,a.neg_b,a.pos_c,a.neg_c,
        a.pos_io,a.pos_ls,a.pos_cp,a.neg_io,a.neg_ls,a.neg_cp,
        a.pos_f,a.neg_f,a.pos_drop_span,a.neg_drop_span,
        a.pos_e_drop_span,a.neg_e_drop_span};
    const double bv[]={b.numerical_hinge_Ke,b.posMy,b.negMy,b.posMc,b.negMc,
        b.pos_a,b.neg_a,b.pos_b,b.neg_b,b.pos_c,b.neg_c,
        b.pos_io,b.pos_ls,b.pos_cp,b.neg_io,b.neg_ls,b.neg_cp,
        b.pos_f,b.neg_f,b.pos_drop_span,b.neg_drop_span,
        b.pos_e_drop_span,b.neg_e_drop_span};
    double out=0.0;
    for(std::size_t i=0;i<sizeof(av)/sizeof(av[0]);++i) out=std::max(out,rel_change(av[i],bv[i]));
    return out;
}

double mix(double old_value,double new_value,double r){
    return old_value+r*(new_value-old_value);
}

RCColumnDemandState relaxed_demand(const RCColumnDemandState& old_demand,
                                   const RCColumnDemandState& observed,
                                   double relaxation,
                                   double gravity_compression){
    RCColumnDemandState d;
    d.max_compression_kip=std::max(gravity_compression,
        mix(old_demand.max_compression_kip,observed.max_compression_kip,relaxation));
    d.max_tension_kip=std::max(0.0,
        mix(old_demand.max_tension_kip,observed.max_tension_kip,relaxation));
    d.max_abs_shear_kip=std::max(0.0,
        mix(old_demand.max_abs_shear_kip,observed.max_abs_shear_kip,relaxation));
    return d;
}
}

CallbackRCColumnRulesResolver::CallbackRCColumnRulesResolver(Function fn):fn_(std::move(fn)){
    if(!fn_) throw std::invalid_argument("RC column rules resolver callback is empty");
}

RCColumnRuleResolution CallbackRCColumnRulesResolver::resolve(
    const RCColumnSectionInput& section,const RCColumnDemandState& demand) const{
    auto r=fn_(section,demand);
    validate_resolution(r);
    return r;
}

ASCE41HingeParams RCColumnASCE41Provider::common_params(const RCColumnResolvedParameters& r){
    validate_resolved(r);
    ASCE41HingeParams p;
    p.Ke=r.numerical_hinge_Ke;
    p.posFy=r.posMy; p.negFy=r.negMy;
    p.posMc=r.posMc; p.negMc=r.negMc;
    p.pos_a=r.pos_a; p.neg_a=r.neg_a;
    p.pos_b=r.pos_b; p.neg_b=r.neg_b;
    p.pos_f=r.pos_f; p.neg_f=r.neg_f;
    p.pos_c=r.pos_c; p.neg_c=r.neg_c;
    p.pos_drop_span=r.pos_drop_span; p.neg_drop_span=r.neg_drop_span;
    p.pos_e_drop_span=r.pos_e_drop_span; p.neg_e_drop_span=r.neg_e_drop_span;
    p.pos_io=r.pos_io; p.pos_ls=r.pos_ls; p.pos_cp=r.pos_cp;
    p.neg_io=r.neg_io; p.neg_ls=r.neg_ls; p.neg_cp=r.neg_cp;
    p.lambda_strength=r.lambda_strength;
    p.lambda_unloading=r.lambda_unloading;
    p.cyclic_exponent=r.cyclic_exponent;
    return p;
}

RCColumnModelSpec RCColumnASCE41Provider::nist_asce41_17_benchmark(
    const RCColumnResolvedParameters& r, std::string provenance){
    RCColumnModelSpec spec;
    spec.edition=RCColumnCodeEdition::NIST_ASCE41_17_Benchmark;
    spec.hinge=common_params(r);
    spec.hinge.backbone_shape=ASCE41BackboneShape::StraightCE;
    // NIST benchmark interpretation: E is the collapse deformation for the
    // component backbone, so no QuakeCore-only post-E gravity interval is added.
    spec.hinge.pos_f=spec.hinge.pos_b;
    spec.hinge.neg_f=spec.hinge.neg_b;
    spec.provenance=provenance.empty()?"NIST ASCE 41-17 benchmark: caller-resolved values":std::move(provenance);
    spec.code_coefficients_embedded=false;
    return spec;
}

RCColumnModelSpec RCColumnASCE41Provider::nist_asce41_17_benchmark(
    const RCColumnRuleResolution& r){
    validate_resolution(r);
    auto spec=nist_asce41_17_benchmark(r.resolved,r.provenance);
    spec.demand_used=r.demand_used;
    spec.audit=r.audit;
    return spec;
}

RCColumnModelSpec RCColumnASCE41Provider::asce41_23_aci369_1_22(
    const RCColumnResolvedParameters& r, ASCE41BackboneShape topology,
    std::string provenance){
    if(provenance.empty())
        throw std::invalid_argument("ASCE 41-23 / ACI 369.1-22 model requires explicit parameter provenance");
    RCColumnModelSpec spec;
    spec.edition=RCColumnCodeEdition::ASCE41_23_ACI369_1_22;
    spec.hinge=common_params(r);
    spec.hinge.backbone_shape=topology;
    spec.provenance=std::move(provenance);
    spec.code_coefficients_embedded=false;
    return spec;
}

RCColumnModelSpec RCColumnASCE41Provider::asce41_23_aci369_1_22(
    const RCColumnRuleResolution& r,ASCE41BackboneShape topology){
    validate_resolution(r);
    auto spec=asce41_23_aci369_1_22(r.resolved,topology,r.provenance);
    spec.demand_used=r.demand_used;
    spec.audit=r.audit;
    return spec;
}

RCColumnAxialIterationResult RCColumnAxialIteration::run_two_pass_asce41_23_aci369_1_22(
    const RCColumnSectionInput& section,
    const RCColumnRulesResolver& resolver,
    const RCColumnResponseRunner& response_runner,
    ASCE41BackboneShape topology){
    RCColumnAxialIterationOptions options;
    options.max_iterations=2;
    options.min_iterations=2;
    options.relaxation=1.0;
    return run_asce41_23_aci369_1_22(section,resolver,response_runner,topology,options);
}

RCColumnAxialIterationResult RCColumnAxialIteration::run_asce41_23_aci369_1_22(
    const RCColumnSectionInput& section,
    const RCColumnRulesResolver& resolver,
    const RCColumnResponseRunner& response_runner,
    ASCE41BackboneShape topology,
    const RCColumnAxialIterationOptions& options){
    validate_section(section);
    if(!response_runner) throw std::invalid_argument("RC column response runner is empty");
    if(options.max_iterations==0||options.min_iterations==0||
       options.min_iterations>options.max_iterations||
       options.demand_relative_tolerance<0.0||options.parameter_relative_tolerance<0.0||
       options.relaxation<=0.0||options.relaxation>1.0)
        throw std::invalid_argument("invalid RC column axial iteration options");

    RCColumnAxialIterationResult result;
    RCColumnDemandState demand;
    demand.max_compression_kip=section.gravity_axial_compression_kip;
    demand.max_tension_kip=0.0;
    demand.max_abs_shear_kip=std::max(0.0,section.initial_shear_demand_kip);

    for(std::size_t i=0;i<options.max_iterations;++i){
        auto resolution=resolver.resolve(section,demand);
        // The resolver must report exactly what demand state it used so the
        // generated parameter record remains auditable.
        resolution.demand_used=demand;
        auto model=RCColumnASCE41Provider::asce41_23_aci369_1_22(resolution,topology);
        auto observed=response_runner(model,i);
        if(observed.max_compression_kip<0.0||observed.max_tension_kip<0.0||observed.max_abs_shear_kip<0.0)
            throw std::runtime_error("RC column response runner returned a negative demand magnitude");

        RCColumnDemandState target=observed;
        target.max_compression_kip=std::max(section.gravity_axial_compression_kip,target.max_compression_kip);
        auto candidate_resolution=resolver.resolve(section,target);
        candidate_resolution.demand_used=target;

        RCColumnAxialIterationStep step;
        step.iteration=i+1;
        step.demand_input=demand;
        step.demand_observed=observed;
        step.model=model;
        step.demand_relative_change=demand_change(demand,target);
        step.parameter_relative_change=parameter_change(resolution.resolved,candidate_resolution.resolved);
        result.steps.push_back(step);
        result.final_model=model;
        result.final_observed_demand=observed;

        if(i+1>=options.min_iterations &&
           step.demand_relative_change<=options.demand_relative_tolerance &&
           step.parameter_relative_change<=options.parameter_relative_tolerance){
            result.converged=true;
            break;
        }
        demand=relaxed_demand(demand,target,options.relaxation,section.gravity_axial_compression_kip);
    }
    return result;
}

} // namespace quake
