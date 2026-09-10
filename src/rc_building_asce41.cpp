#include "quake/rc_building_asce41.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace quake {
namespace {

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
    for(std::size_t i=0;i<sizeof(av)/sizeof(av[0]);++i)
        out=std::max(out,rel_change(av[i],bv[i]));
    return out;
}

double mix(double old_value,double new_value,double r){
    return old_value+r*(new_value-old_value);
}

RCColumnDemandState initial_demand(const RCColumnSectionInput& section){
    RCColumnDemandState d;
    d.max_compression_kip=section.gravity_axial_compression_kip;
    d.max_tension_kip=0.0;
    d.max_abs_shear_kip=std::max(0.0,section.initial_shear_demand_kip);
    return d;
}

RCColumnDemandState target_demand(const RCColumnSectionInput& section,
                                  const RCColumnDemandState& observed){
    if(observed.max_compression_kip<0.0||observed.max_tension_kip<0.0||
       observed.max_abs_shear_kip<0.0)
        throw std::runtime_error("building response runner returned a negative RC column demand magnitude");
    RCColumnDemandState target=observed;
    target.max_compression_kip=std::max(section.gravity_axial_compression_kip,
                                        target.max_compression_kip);
    return target;
}

RCColumnDemandState relaxed_demand(const RCColumnSectionInput& section,
                                   const RCColumnDemandState& old_demand,
                                   const RCColumnDemandState& target,
                                   double relaxation){
    RCColumnDemandState d;
    d.max_compression_kip=std::max(section.gravity_axial_compression_kip,
        mix(old_demand.max_compression_kip,target.max_compression_kip,relaxation));
    d.max_tension_kip=std::max(0.0,
        mix(old_demand.max_tension_kip,target.max_tension_kip,relaxation));
    d.max_abs_shear_kip=std::max(0.0,
        mix(old_demand.max_abs_shear_kip,target.max_abs_shear_kip,relaxation));
    return d;
}

void validate_options(const RCColumnAxialIterationOptions& o){
    if(o.max_iterations==0||o.min_iterations==0||o.min_iterations>o.max_iterations||
       o.demand_relative_tolerance<0.0||o.parameter_relative_tolerance<0.0||
       o.relaxation<=0.0||o.relaxation>1.0)
        throw std::invalid_argument("invalid RC building ASCE41 iteration options");
}

std::vector<RCBuildingColumnDefinition> canonical_columns(
    const std::vector<RCBuildingColumnDefinition>& columns){
    if(columns.empty()) throw std::invalid_argument("RC building ASCE41 coordinator requires at least one column");
    auto out=columns;
    std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){
        return a.section.component_id<b.section.component_id;
    });
    std::string previous;
    for(const auto& c:out){
        const auto& s=c.section;
        if(s.component_id.empty())
            throw std::invalid_argument("RC building column requires a nonempty component_id");
        if(!previous.empty()&&previous==s.component_id)
            throw std::invalid_argument("duplicate RC building column component_id: "+s.component_id);
        previous=s.component_id;
        if(s.width_in<=0.0||s.depth_in<=0.0||s.clear_length_in<=0.0||
           s.expected_fc_ksi<=0.0||s.expected_fy_long_ksi<=0.0||
           s.longitudinal_bar_count<=0||s.longitudinal_bar_area_in2<=0.0||
           s.gravity_axial_compression_kip<0.0)
            throw std::invalid_argument("invalid RC building column section input: "+s.component_id);
    }
    return out;
}

std::map<std::string,RCColumnDemandState> observation_map(
    const RCBuildingAnalysisObservation& observation,
    const std::vector<RCBuildingColumnDefinition>& columns){
    std::map<std::string,RCColumnDemandState> out;
    for(const auto& item:observation.column_demands){
        if(item.component_id.empty())
            throw std::runtime_error("building response runner returned an empty component_id");
        if(!out.emplace(item.component_id,item.demand).second)
            throw std::runtime_error("building response runner returned duplicate demand for "+item.component_id);
    }
    if(out.size()!=columns.size())
        throw std::runtime_error("building response runner did not return exactly one demand per coordinated RC column");
    for(const auto& c:columns){
        if(out.find(c.section.component_id)==out.end())
            throw std::runtime_error("building response runner omitted demand for "+c.section.component_id);
    }
    std::set<std::string> known_ids;
    for(const auto& c:columns) known_ids.insert(c.section.component_id);
    for(const auto& kv:out){
        if(known_ids.find(kv.first)==known_ids.end())
            throw std::runtime_error("building response runner returned unknown RC column "+kv.first);
    }
    return out;
}

} // namespace

RCBuildingASCE41IterationResult RCBuildingASCE41Iteration::run_two_pass_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,
    const RCColumnRulesResolver& resolver,
    const RCBuildingResponseRunner& response_runner){
    RCColumnAxialIterationOptions options;
    options.max_iterations=2;
    options.min_iterations=2;
    options.relaxation=1.0;
    return run_asce41_23_aci369_1_22(columns,resolver,response_runner,options);
}

RCBuildingASCE41IterationResult RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,
    const RCColumnRulesResolver& resolver,
    const RCBuildingResponseRunner& response_runner,
    const RCColumnAxialIterationOptions& options){
    validate_options(options);
    if(!response_runner) throw std::invalid_argument("RC building response runner is empty");
    const auto defs=canonical_columns(columns);

    std::map<std::string,RCColumnDemandState> demands;
    for(const auto& c:defs) demands[c.section.component_id]=initial_demand(c.section);

    RCBuildingASCE41IterationResult result;
    for(std::size_t i=0;i<options.max_iterations;++i){
        std::vector<RCBuildingColumnModel> models;
        std::map<std::string,RCColumnRuleResolution> resolutions;
        models.reserve(defs.size());

        // Synchronous generation: every model uses the complete demand field
        // committed at the start of this iteration.
        for(const auto& c:defs){
            const auto& id=c.section.component_id;
            auto resolution=resolver.resolve(c.section,demands.at(id));
            resolution.demand_used=demands.at(id);
            auto model=RCColumnASCE41Provider::asce41_23_aci369_1_22(
                resolution,c.topology);
            resolutions.emplace(id,std::move(resolution));
            models.push_back({id,std::move(model)});
        }

        auto observation=response_runner(models,i);
        if(!observation.analysis_succeeded){
            result.converged=false;
            result.termination=observation.physical_collapse?
                RCBuildingASCE41Termination::PhysicalCollapse:
                RCBuildingASCE41Termination::AnalysisFailure;
            result.termination_detail=observation.status.empty()?
                (observation.physical_collapse?"building analysis reached physical collapse":"building analysis reported failure"):
                observation.status;
            result.final_models=models;
            return result;
        }
        const auto observed_by_id=observation_map(observation,defs);

        RCBuildingASCE41IterationStep step;
        step.iteration=i+1;
        step.models=models;
        step.analysis_status=observation.status;
        step.component_metrics.reserve(defs.size());

        std::map<std::string,RCColumnDemandState> targets;
        double global_demand_change=0.0;
        double global_parameter_change=0.0;

        // Evaluate every candidate update before committing any of them.
        for(const auto& c:defs){
            const auto& id=c.section.component_id;
            const auto target=target_demand(c.section,observed_by_id.at(id));
            auto candidate=resolver.resolve(c.section,target);
            candidate.demand_used=target;
            const double dc=demand_change(demands.at(id),target);
            const double pc=parameter_change(resolutions.at(id).resolved,candidate.resolved);
            step.component_metrics.push_back({id,demands.at(id),observed_by_id.at(id),dc,pc});
            global_demand_change=std::max(global_demand_change,dc);
            global_parameter_change=std::max(global_parameter_change,pc);
            targets.emplace(id,target);
        }
        step.global_demand_relative_change=global_demand_change;
        step.global_parameter_relative_change=global_parameter_change;
        result.steps.push_back(std::move(step));
        result.final_models=models;
        result.final_observed_demands.clear();
        for(const auto& c:defs)
            result.final_observed_demands.push_back({c.section.component_id,
                                                     observed_by_id.at(c.section.component_id)});

        if(i+1>=options.min_iterations&&
           global_demand_change<=options.demand_relative_tolerance&&
           global_parameter_change<=options.parameter_relative_tolerance){
            result.converged=true;
            result.termination=RCBuildingASCE41Termination::Converged;
            result.termination_detail="whole-building RC column parameter field converged";
            return result;
        }

        // Jacobi-style simultaneous commit after all candidates have been formed.
        std::map<std::string,RCColumnDemandState> next_demands;
        for(const auto& c:defs){
            const auto& id=c.section.component_id;
            next_demands.emplace(id,relaxed_demand(c.section,demands.at(id),
                                                   targets.at(id),options.relaxation));
        }
        demands=std::move(next_demands);
    }

    result.converged=false;
    result.termination=RCBuildingASCE41Termination::IterationLimit;
    result.termination_detail="whole-building RC column parameter field reached iteration limit";
    return result;
}

} // namespace quake
