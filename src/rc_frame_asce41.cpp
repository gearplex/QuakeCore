#include "quake/rc_frame_asce41.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace quake {
namespace {

template<class Model>
void canonicalize_and_validate_bindings(const Model& model,
                                        std::vector<RCColumnElementBinding>& bindings){
    if(bindings.empty()) throw std::invalid_argument("native RC frame analysis requires column bindings");
    std::sort(bindings.begin(),bindings.end(),[](const auto& a,const auto& b){return a.component_id<b.component_id;});
    std::string prev;
    std::set<int> bound_elements;
    for(const auto& b:bindings){
        if(b.component_id.empty()) throw std::invalid_argument("native RC column binding requires component_id");
        if(!prev.empty()&&prev==b.component_id) throw std::invalid_argument("duplicate native RC column binding: "+b.component_id);
        prev=b.component_id;
        if(b.elastic_element_ids.empty()) throw std::invalid_argument("native RC column binding requires at least one elastic element: "+b.component_id);
        for(int id:b.elastic_element_ids){
            if(model.elastic_element_index(id)<0) throw std::invalid_argument("native RC column binding references unknown elastic element "+std::to_string(id));
            if(!bound_elements.insert(id).second) throw std::invalid_argument("elastic element is bound to more than one RC column: "+std::to_string(id));
        }
        for(int id:b.hinge_component_ids)
            if(model.nonlinear_component_index(id)<0) throw std::invalid_argument("native RC column binding references unknown hinge "+std::to_string(id));
    }
}

void validate_binding_coverage(const std::vector<RCBuildingColumnModel>& models,
                               const std::vector<RCColumnElementBinding>& bindings){
    std::set<std::string> expected,actual;
    for(const auto& m:models){
        if(!expected.insert(m.component_id).second) throw std::runtime_error("duplicate coordinated RC model id: "+m.component_id);
    }
    for(const auto& b:bindings) actual.insert(b.component_id);
    if(expected!=actual) throw std::runtime_error("native frame bindings do not exactly cover coordinated RC column models");
}

template<class Model>
void apply_asce_material_field(Model& model,
                               const std::vector<RCColumnElementBinding>& bindings,
                               const std::vector<RCBuildingColumnModel>& models){
    validate_binding_coverage(models,bindings);
    std::map<std::string,const RCBuildingColumnModel*> by_id;
    for(const auto& m:models) by_id.emplace(m.component_id,&m);
    std::vector<std::pair<int,NonlinearMaterial>> replacements;
    std::set<int> hinge_ids;
    for(const auto& b:bindings){
        const auto it=by_id.find(b.component_id);
        if(it==by_id.end()) throw std::runtime_error("missing coordinated RC model for native material field: "+b.component_id);
        if(b.hinge_component_ids.empty()) throw std::runtime_error("compiled ASCE material-field binding has no hinge IDs: "+b.component_id);
        for(int hid:b.hinge_component_ids){
            if(!hinge_ids.insert(hid).second) throw std::runtime_error("nonlinear hinge is assigned to more than one RC column: "+std::to_string(hid));
            replacements.emplace_back(hid,NonlinearMaterial{ASCE41HingeMaterial(it->second->model.hinge)});
        }
    }
    model.replace_nonlinear_material_field(replacements);
}

void update_demand(RCColumnDemandState& d,double axial_compression,double shear){
    if(axial_compression>=0.0) d.max_compression_kip=std::max(d.max_compression_kip,axial_compression);
    else d.max_tension_kip=std::max(d.max_tension_kip,-axial_compression);
    d.max_abs_shear_kip=std::max(d.max_abs_shear_kip,std::abs(shear));
}

std::string termination_status(const AnalysisResult& r){
    if(!r.termination_reason.empty()) return r.termination_reason;
    switch(r.termination){
        case AnalysisTermination::Completed: return "completed";
        case AnalysisTermination::PhysicalCollapse: return "physical collapse";
        case AnalysisTermination::NumericalFailure: return "numerical failure";
        case AnalysisTermination::InitialInstability: return "initial instability";
    }
    return "unknown analysis termination";
}

RCBuildingAnalysisObservation analyze_native_frame2d(
    RCNativeFrame2DBuild build,const RCNativeFrameAnalysisSettings& settings,
    const std::vector<RCBuildingColumnModel>& models){
    if(settings.dt<=0.0||settings.ground_accel.empty()) throw std::invalid_argument("invalid native 2D RC analysis settings");
    canonicalize_and_validate_bindings(build.model,build.bindings);
    validate_binding_coverage(models,build.bindings);
    RCFrame2DColumnDemandRecorder recorder(build.model,build.bindings);
    auto options=settings.newmark_options;
    const auto chained=options.accepted_substep_state_observer;
    options.return_numerical_failure=true;
    options.accepted_substep_state_observer=[&](std::size_t step,std::size_t sub,int depth,double time,double ag,
        const std::vector<double>& u,const std::vector<double>& v,const std::vector<double>& a,const std::vector<double>& state){
        recorder.observe(u,state);
        if(chained) chained(step,sub,depth,time,ag,u,v,a,state);
    };
    const auto result=run_newmark_robust(build.model,settings.ground_accel,settings.dt,settings.strategy,options);
    RCBuildingAnalysisObservation out;
    out.status=termination_status(result);
    out.column_demands=recorder.observations();
    if(result.termination==AnalysisTermination::Completed) out.analysis_succeeded=true;
    else {out.analysis_succeeded=false;out.physical_collapse=result.termination==AnalysisTermination::PhysicalCollapse;}
    return out;
}

RCBuildingAnalysisObservation analyze_native_frame3d(
    RCNativeFrame3DBuild build,const RCNativeFrameAnalysisSettings& settings,
    const std::vector<RCBuildingColumnModel>& models){
    if(settings.dt<=0.0||settings.ground_accel.empty()) throw std::invalid_argument("invalid native 3D RC analysis settings");
    canonicalize_and_validate_bindings(build.model,build.bindings);
    validate_binding_coverage(models,build.bindings);
    RCFrame3DColumnDemandRecorder recorder(build.model,build.bindings);
    auto options=settings.newmark_options;
    const auto chained=options.accepted_substep_state_observer;
    options.return_numerical_failure=true;
    options.accepted_substep_state_observer=[&](std::size_t step,std::size_t sub,int depth,double time,double ag,
        const std::vector<double>& u,const std::vector<double>& v,const std::vector<double>& a,const std::vector<double>& state){
        recorder.observe(u,state);
        if(chained) chained(step,sub,depth,time,ag,u,v,a,state);
    };
    const auto result=run_newmark_robust(build.model,settings.ground_accel,settings.dt,settings.strategy,options);
    RCBuildingAnalysisObservation out;
    out.status=termination_status(result);
    out.column_demands=recorder.observations();
    if(result.termination==AnalysisTermination::Completed) out.analysis_succeeded=true;
    else {out.analysis_succeeded=false;out.physical_collapse=result.termination==AnalysisTermination::PhysicalCollapse;}
    return out;
}

RCBuildingAnalysisObservation analyze_compiled_frame2d(
    CompiledFrame2D& model,const std::vector<RCColumnElementBinding>& bindings,
    const RCNativeFrameAnalysisSettings& settings,const std::vector<RCBuildingColumnModel>& models){
    if(settings.dt<=0.0||settings.ground_accel.empty()) throw std::invalid_argument("invalid compiled 2D RC analysis settings");
    validate_binding_coverage(models,bindings);
    RCFrame2DColumnDemandRecorder recorder(model,bindings);
    auto options=settings.newmark_options;const auto chained=options.accepted_substep_state_observer;options.return_numerical_failure=true;
    options.accepted_substep_state_observer=[&](std::size_t step,std::size_t sub,int depth,double time,double ag,
        const std::vector<double>& u,const std::vector<double>& v,const std::vector<double>& a,const std::vector<double>& state){
        recorder.observe(u,state);if(chained) chained(step,sub,depth,time,ag,u,v,a,state);
    };
    const auto result=run_newmark_robust(model,settings.ground_accel,settings.dt,settings.strategy,options);
    RCBuildingAnalysisObservation out;out.status=termination_status(result);out.column_demands=recorder.observations();
    if(result.termination==AnalysisTermination::Completed) out.analysis_succeeded=true;
    else{out.analysis_succeeded=false;out.physical_collapse=result.termination==AnalysisTermination::PhysicalCollapse;}
    return out;
}

RCBuildingAnalysisObservation analyze_compiled_frame3d(
    CompiledFrame3D& model,const std::vector<RCColumnElementBinding>& bindings,
    const RCNativeFrameAnalysisSettings& settings,const std::vector<RCBuildingColumnModel>& models){
    if(settings.dt<=0.0||settings.ground_accel.empty()) throw std::invalid_argument("invalid compiled 3D RC analysis settings");
    validate_binding_coverage(models,bindings);
    RCFrame3DColumnDemandRecorder recorder(model,bindings);
    auto options=settings.newmark_options;const auto chained=options.accepted_substep_state_observer;options.return_numerical_failure=true;
    options.accepted_substep_state_observer=[&](std::size_t step,std::size_t sub,int depth,double time,double ag,
        const std::vector<double>& u,const std::vector<double>& v,const std::vector<double>& a,const std::vector<double>& state){
        recorder.observe(u,state);if(chained) chained(step,sub,depth,time,ag,u,v,a,state);
    };
    const auto result=run_newmark_robust(model,settings.ground_accel,settings.dt,settings.strategy,options);
    RCBuildingAnalysisObservation out;out.status=termination_status(result);out.column_demands=recorder.observations();
    if(result.termination==AnalysisTermination::Completed) out.analysis_succeeded=true;
    else{out.analysis_succeeded=false;out.physical_collapse=result.termination==AnalysisTermination::PhysicalCollapse;}
    return out;
}

} // namespace

RCFrame2DColumnDemandRecorder::RCFrame2DColumnDemandRecorder(
    const CompiledFrame2D& model,std::vector<RCColumnElementBinding> bindings)
    :model_(&model),bindings_(std::move(bindings)){
    canonicalize_and_validate_bindings(model,bindings_);
    demands_.resize(bindings_.size());
    // Include the gravity/preload state even if the dynamic analysis terminates
    // before the first accepted step.
    observe(std::vector<double>(static_cast<std::size_t>(model.dof()),0.0),model.initial_nonlinear_state());
}

void RCFrame2DColumnDemandRecorder::observe(const std::vector<double>& u,const std::vector<double>& committed_state){
    if(static_cast<int>(committed_state.size())!=model_->nonlinear_state_size()) throw std::invalid_argument("2D native recorder committed-state size");
    for(std::size_t i=0;i<bindings_.size();++i){
        for(int eid:bindings_[i].elastic_element_ids){
            const auto r=model_->elastic_element_response(eid,u);
            update_demand(demands_[i],r.axial_compression,std::max(std::abs(r.shear_i),std::abs(r.shear_j)));
        }
    }
    ++committed_samples_;
}

std::vector<RCBuildingColumnDemandObservation> RCFrame2DColumnDemandRecorder::observations() const{
    std::vector<RCBuildingColumnDemandObservation> out;out.reserve(bindings_.size());
    for(std::size_t i=0;i<bindings_.size();++i)out.push_back({bindings_[i].component_id,demands_[i]});
    return out;
}

RCFrame3DColumnDemandRecorder::RCFrame3DColumnDemandRecorder(
    const CompiledFrame3D& model,std::vector<RCColumnElementBinding> bindings)
    :model_(&model),bindings_(std::move(bindings)){
    canonicalize_and_validate_bindings(model,bindings_);
    demands_.resize(bindings_.size());
    observe(std::vector<double>(static_cast<std::size_t>(model.dof()),0.0),model.initial_nonlinear_state());
}

void RCFrame3DColumnDemandRecorder::observe(const std::vector<double>& u,const std::vector<double>& committed_state){
    if(static_cast<int>(committed_state.size())!=model_->nonlinear_state_size()) throw std::invalid_argument("3D native recorder committed-state size");
    for(std::size_t i=0;i<bindings_.size();++i){
        for(int eid:bindings_[i].elastic_element_ids){
            const auto r=model_->elastic_element_response(eid,u);
            const double vi=std::hypot(r.shear_y_i,r.shear_z_i);
            const double vj=std::hypot(r.shear_y_j,r.shear_z_j);
            update_demand(demands_[i],r.axial_compression,std::max(vi,vj));
        }
    }
    ++committed_samples_;
}

std::vector<RCBuildingColumnDemandObservation> RCFrame3DColumnDemandRecorder::observations() const{
    std::vector<RCBuildingColumnDemandObservation> out;out.reserve(bindings_.size());
    for(std::size_t i=0;i<bindings_.size();++i)out.push_back({bindings_[i].component_id,demands_[i]});
    return out;
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_compiled_frame2d_two_pass_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    RCNativeCompiledFrame2D& frame,const RCNativeFrameAnalysisSettings& settings){
    canonicalize_and_validate_bindings(frame.model,frame.bindings);
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        apply_asce_material_field(frame.model,frame.bindings,models);
        return analyze_compiled_frame2d(frame.model,frame.bindings,settings,models);
    };
    return RCBuildingASCE41Iteration::run_two_pass_asce41_23_aci369_1_22(columns,resolver,runner);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_compiled_frame2d_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    RCNativeCompiledFrame2D& frame,const RCNativeFrameAnalysisSettings& settings,
    const RCColumnAxialIterationOptions& iteration_options){
    canonicalize_and_validate_bindings(frame.model,frame.bindings);
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        apply_asce_material_field(frame.model,frame.bindings,models);
        return analyze_compiled_frame2d(frame.model,frame.bindings,settings,models);
    };
    return RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(columns,resolver,runner,iteration_options);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_compiled_frame3d_two_pass_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    RCNativeCompiledFrame3D& frame,const RCNativeFrameAnalysisSettings& settings){
    canonicalize_and_validate_bindings(frame.model,frame.bindings);
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        apply_asce_material_field(frame.model,frame.bindings,models);
        return analyze_compiled_frame3d(frame.model,frame.bindings,settings,models);
    };
    return RCBuildingASCE41Iteration::run_two_pass_asce41_23_aci369_1_22(columns,resolver,runner);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_compiled_frame3d_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    RCNativeCompiledFrame3D& frame,const RCNativeFrameAnalysisSettings& settings,
    const RCColumnAxialIterationOptions& iteration_options){
    canonicalize_and_validate_bindings(frame.model,frame.bindings);
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        apply_asce_material_field(frame.model,frame.bindings,models);
        return analyze_compiled_frame3d(frame.model,frame.bindings,settings,models);
    };
    return RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(columns,resolver,runner,iteration_options);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_frame2d_two_pass_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    const RCNativeFrame2DFactory& factory,const RCNativeFrameAnalysisSettings& settings){
    if(!factory) throw std::invalid_argument("native 2D RC model factory is empty");
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t iteration){
        return analyze_native_frame2d(factory(models,iteration),settings,models);
    };
    return RCBuildingASCE41Iteration::run_two_pass_asce41_23_aci369_1_22(columns,resolver,runner);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_frame2d_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    const RCNativeFrame2DFactory& factory,const RCNativeFrameAnalysisSettings& settings,
    const RCColumnAxialIterationOptions& iteration_options){
    if(!factory) throw std::invalid_argument("native 2D RC model factory is empty");
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t iteration){
        return analyze_native_frame2d(factory(models,iteration),settings,models);
    };
    return RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(columns,resolver,runner,iteration_options);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_frame3d_two_pass_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    const RCNativeFrame3DFactory& factory,const RCNativeFrameAnalysisSettings& settings){
    if(!factory) throw std::invalid_argument("native 3D RC model factory is empty");
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t iteration){
        return analyze_native_frame3d(factory(models,iteration),settings,models);
    };
    return RCBuildingASCE41Iteration::run_two_pass_asce41_23_aci369_1_22(columns,resolver,runner);
}

RCBuildingASCE41IterationResult RCBuildingASCE41NativeFrame::run_frame3d_asce41_23_aci369_1_22(
    const std::vector<RCBuildingColumnDefinition>& columns,const RCColumnRulesResolver& resolver,
    const RCNativeFrame3DFactory& factory,const RCNativeFrameAnalysisSettings& settings,
    const RCColumnAxialIterationOptions& iteration_options){
    if(!factory) throw std::invalid_argument("native 3D RC model factory is empty");
    RCBuildingResponseRunner runner=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t iteration){
        return analyze_native_frame3d(factory(models,iteration),settings,models);
    };
    return RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(columns,resolver,runner,iteration_options);
}

} // namespace quake
