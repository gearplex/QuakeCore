#include "quake/model_comparison.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace quake {
namespace {

std::vector<double> drift_ratios(const std::vector<double>& floor,
                                 const std::vector<double>& elevations){
    if(floor.size()!=elevations.size()) throw std::invalid_argument("comparison story geometry size mismatch");
    std::vector<double> out(floor.size());double lower_u=0.0,lower_z=0.0;
    for(std::size_t i=0;i<floor.size();++i){
        const double dz=elevations[i]-lower_z;
        if(dz<=0.0) throw std::invalid_argument("comparison story elevations must increase");
        out[i]=(floor[i]-lower_u)/dz;lower_u=floor[i];lower_z=elevations[i];
    }
    return out;
}

void update_first(double& dst,bool flag,double time){if(flag&&dst<0.0)dst=time;}



ModelComparisonColumnSample column_sample(const CompiledFrame2D& frame,const RCColumnElementBinding& b,
                                          std::size_t step,double time,const std::vector<double>& u){
    ModelComparisonColumnSample out;out.step=step;out.time=time;
    for(int eid:b.elastic_element_ids){
        const auto r=frame.elastic_element_response(eid,u);
        out.max_compression_kip=std::max(out.max_compression_kip,std::max(0.0,r.axial_compression));
        out.max_tension_kip=std::max(out.max_tension_kip,std::max(0.0,-r.axial_compression));
        out.max_abs_shear_kip=std::max(out.max_abs_shear_kip,std::max(std::abs(r.shear_i),std::abs(r.shear_j)));
    }
    return out;
}

ModelComparisonColumnSample column_sample(const CompiledFrame3D& frame,const RCColumnElementBinding& b,
                                          std::size_t step,double time,const std::vector<double>& u){
    ModelComparisonColumnSample out;out.step=step;out.time=time;
    for(int eid:b.elastic_element_ids){
        const auto r=frame.elastic_element_response(eid,u);
        out.max_compression_kip=std::max(out.max_compression_kip,std::max(0.0,r.axial_compression));
        out.max_tension_kip=std::max(out.max_tension_kip,std::max(0.0,-r.axial_compression));
        out.max_abs_shear_kip=std::max(out.max_abs_shear_kip,std::max(std::hypot(r.shear_y_i,r.shear_z_i),std::hypot(r.shear_y_j,r.shear_z_j)));
    }
    return out;
}

struct PeakState {
    bool valid{};
    std::size_t step{};
    double time{};
    std::vector<double> u;
    std::vector<double> state;
};



template<class Frame>
void compute_story_tangent(const Frame& frame,const SparseMatrixCSC& kt,
                           ModelComparisonTangentSnapshot& snap){
    const int n=frame.dof();const int ns=static_cast<int>(frame.story_elevations().size());
    if(ns<=0)return;
    std::vector<std::vector<double>> c(static_cast<std::size_t>(ns),std::vector<double>(static_cast<std::size_t>(n),0.0));
    std::vector<double> e(static_cast<std::size_t>(n),0.0);
    for(int k=0;k<n;++k){
        std::fill(e.begin(),e.end(),0.0);e[static_cast<std::size_t>(k)]=1.0;
        const auto y=frame.story_response_values(e);if(static_cast<int>(y.size())!=ns)return;
        for(int i=0;i<ns;++i)c[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]=y[static_cast<std::size_t>(i)];
    }
    try{
        std::vector<double> flex(static_cast<std::size_t>(ns*ns),0.0);
        for(int j=0;j<ns;++j){
            const auto x=superlu_solve_once(kt,c[static_cast<std::size_t>(j)]);
            for(int i=0;i<ns;++i){double v=0.0;for(int k=0;k<n;++k)v+=c[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]*x[static_cast<std::size_t>(k)];flex[static_cast<std::size_t>(i*ns+j)]=v;}
        }
        std::vector<Triplet> ft;ft.reserve(static_cast<std::size_t>(ns*ns));
        for(int i=0;i<ns;++i)for(int j=0;j<ns;++j)ft.push_back({i,j,flex[static_cast<std::size_t>(i*ns+j)]});
        const auto F=SparseMatrixCSC::from_triplets(ns,ns,ft,0.0);
        std::vector<double> kfloor(static_cast<std::size_t>(ns*ns),0.0);
        for(int j=0;j<ns;++j){std::vector<double> rhs(static_cast<std::size_t>(ns),0.0);rhs[static_cast<std::size_t>(j)]=1.0;const auto x=superlu_solve_once(F,rhs);for(int i=0;i<ns;++i)kfloor[static_cast<std::size_t>(i*ns+j)]=x[static_cast<std::size_t>(i)];}
        // q_floor = T * delta_story, T(i,j)=1 for j<=i.
        std::vector<double> kd(static_cast<std::size_t>(ns*ns),0.0);
        for(int a=0;a<ns;++a)for(int b=0;b<ns;++b){double v=0.0;for(int i=a;i<ns;++i)for(int j=b;j<ns;++j)v+=kfloor[static_cast<std::size_t>(i*ns+j)];kd[static_cast<std::size_t>(a*ns+b)]=v;}
        snap.story_tangent_available=true;snap.floor_tangent_matrix=std::move(kfloor);snap.interstory_tangent_matrix=std::move(kd);snap.interstory_tangent_diagonal.resize(static_cast<std::size_t>(ns));
        for(int i=0;i<ns;++i)snap.interstory_tangent_diagonal[static_cast<std::size_t>(i)]=snap.interstory_tangent_matrix[static_cast<std::size_t>(i*ns+i)];
    }catch(const std::exception&){
        snap.story_tangent_available=false;snap.floor_tangent_matrix.clear();snap.interstory_tangent_matrix.clear();snap.interstory_tangent_diagonal.clear();
    }
}

template<class Frame>
class ComparisonRecorder {
public:
    ComparisonRecorder(const Frame& frame,std::vector<ModelComparisonComponent> components,
                       ModelComparisonOptions options)
        :frame_(&frame),components_(std::move(components)),options_(std::move(options)){
        const auto elevations=frame_->story_elevations();
        peak_story_.assign(elevations.size(),0.0);peak_story_time_.assign(elevations.size(),0.0);
        peak_states_.resize(elevations.size());
        results_.reserve(components_.size());
        std::set<int> ids;
        for(const auto& c:components_){
            if(c.label.empty()) throw std::invalid_argument("comparison component label is empty");
            if(!ids.insert(c.component_id).second) throw std::invalid_argument("duplicate comparison component id");
            if(frame_->nonlinear_component_index(c.component_id)<0) throw std::invalid_argument("comparison references unknown nonlinear component");
            ModelComparisonComponentResult r;r.label=c.label;r.component_id=c.component_id;r.minimum_tangent=std::numeric_limits<double>::infinity();results_.push_back(std::move(r));
        }
        std::set<std::string> column_ids;
        for(const auto& b:options_.column_force_bindings){
            if(b.component_id.empty()||!column_ids.insert(b.component_id).second) throw std::invalid_argument("comparison column bindings require unique component_id");
            if(b.elastic_element_ids.empty()) throw std::invalid_argument("comparison column binding requires elastic elements");
            for(int eid:b.elastic_element_ids) if(frame_->elastic_element_index(eid)<0) throw std::invalid_argument("comparison column binding references unknown elastic element");
            ModelComparisonColumnResult cr;cr.component_id=b.component_id;column_results_.push_back(std::move(cr));
        }
        for(auto step:options_.additional_tangent_mode_steps) requested_steps_.insert(step);
    }

    void observe_envelope(std::size_t step,double time,const std::vector<double>& u,const std::vector<double>& state){
        const auto floor=frame_->story_response_values(u);const auto dr=drift_ratios(floor,frame_->story_elevations());
        for(std::size_t i=0;i<dr.size();++i){
            const double a=std::abs(dr[i]);if(a>peak_story_[i]){peak_story_[i]=a;peak_story_time_[i]=time;peak_states_[i]={true,step,time,u,state};}
        }
        const double roof=frame_->response_value(u);if(std::abs(roof)>peak_roof_){peak_roof_=std::abs(roof);peak_roof_time_=time;}
        for(std::size_t i=0;i<components_.size();++i)update_component(results_[i],frame_->nonlinear_component_snapshot(components_[i].component_id,u,state),time);
        for(std::size_t i=0;i<options_.column_force_bindings.size();++i)update_column(column_results_[i],column_sample(*frame_,options_.column_force_bindings[i],step,time,u),false);
    }

    void observe_output(std::size_t step,double time,const std::vector<double>& u,const std::vector<double>& state){
        observe_envelope(step,time,u,state);
        if(options_.record_story_history){
            const auto floor=frame_->story_response_values(u);story_history_.push_back({step,time,floor,drift_ratios(floor,frame_->story_elevations())});
        }
        if(options_.record_component_history){
            for(std::size_t i=0;i<components_.size();++i){
                const auto s=frame_->nonlinear_component_snapshot(components_[i].component_id,u,state);
                results_[i].history.push_back({step,time,s.deformation,s.force,s.tangent,s.diagnostics});
            }
        }
        for(std::size_t i=0;i<options_.column_force_bindings.size();++i)
            column_results_[i].history.push_back(column_sample(*frame_,options_.column_force_bindings[i],step,time,u));
        if(requested_steps_.count(step)) requested_states_[step]={true,step,time,u,state};
        final_u_=u;final_state_=state;
    }

    void finalize(ModelComparisonVariantResult& out){
        out.story_history=std::move(story_history_);out.peak_abs_story_drift_ratio=peak_story_;out.time_at_peak_abs_story_drift=peak_story_time_;
        out.peak_abs_roof_response=peak_roof_;out.time_at_peak_abs_roof_response=peak_roof_time_;out.components=std::move(results_);out.columns=std::move(column_results_);
        if(!final_u_.empty()){
            const auto floor=frame_->story_response_values(final_u_);out.residual_story_drift_ratio=drift_ratios(floor,frame_->story_elevations());out.residual_roof_response=frame_->response_value(final_u_);
        }
        if(options_.tangent_modes_at_story_peaks){
            for(std::size_t i=0;i<peak_states_.size();++i)if(peak_states_[i].valid){
                const auto key=std::make_pair(peak_states_[i].step,peak_states_[i].time);
                auto it=mode_states_.find(key);
                if(it==mode_states_.end())mode_states_.emplace(key,std::make_pair(std::string("story_")+std::to_string(i+1)+"_peak",peak_states_[i]));
                else it->second.first+=",story_"+std::to_string(i+1)+"_peak";
            }
        }
        for(const auto& [step,ps]:requested_states_)if(ps.valid){
            const auto key=std::make_pair(ps.step,ps.time);auto it=mode_states_.find(key);
            if(it==mode_states_.end())mode_states_.emplace(key,std::make_pair(std::string("requested_step_")+std::to_string(step),ps));
            else it->second.first+=",requested_step_"+std::to_string(step);
        }
        for(const auto& [key,item]:mode_states_){
            const auto& ps=item.second;std::vector<double> force,tangents,trial;
            frame_->internal_force_and_tangent(ps.u,ps.state,force,tangents,trial);
            auto kt=frame_->effective_state_tangent_matrix_with_state(ps.u,tangents,ps.state,0.0,0.0);
            auto modes=modal_analysis_from_stiffness(*frame_,kt,options_.tangent_mode_count);
            ModelComparisonTangentSnapshot snap;snap.step=ps.step;snap.time=ps.time;snap.trigger=item.first;
            for(auto& m:modes){
                ModelComparisonTangentMode mm;mm.period=m.period;mm.eigenvalue=m.eigenvalue;mm.reduced_shape=m.shape;mm.story_shape=frame_->story_response_values(m.shape);
                if(!mm.story_shape.empty()&&std::abs(mm.story_shape.back())>1e-14){const double d=mm.story_shape.back();for(double& v:mm.story_shape)v/=d;}
                snap.modes.push_back(std::move(mm));
            }
            compute_story_tangent(*frame_,kt,snap);
            out.tangent_snapshots.push_back(std::move(snap));
        }
    }

private:
    static void update_component(ModelComparisonComponentResult& r,const NonlinearComponentSnapshot& s,double time){
        const double aq=std::abs(s.deformation);if(aq>r.max_abs_deformation){r.max_abs_deformation=aq;r.time_at_max_abs_deformation=time;}
        const double af=std::abs(s.force);if(af>r.max_abs_force){r.max_abs_force=af;r.time_at_max_abs_force=time;}
        if(s.tangent<r.minimum_tangent){r.minimum_tangent=s.tangent;r.time_at_minimum_tangent=time;}
        update_first(r.first_io_time,s.diagnostics.at_or_beyond_io,time);update_first(r.first_ls_time,s.diagnostics.at_or_beyond_ls,time);
        update_first(r.first_cp_time,s.diagnostics.at_or_beyond_cp,time);update_first(r.first_lateral_loss_time,s.diagnostics.lateral_resistance_lost,time);
        update_first(r.first_failure_time,s.diagnostics.failed,time);
    }
    static void update_column(ModelComparisonColumnResult& r,const ModelComparisonColumnSample& s,bool){
        if(s.max_compression_kip>r.max_compression_kip){r.max_compression_kip=s.max_compression_kip;r.time_at_max_compression=s.time;}
        if(s.max_tension_kip>r.max_tension_kip){r.max_tension_kip=s.max_tension_kip;r.time_at_max_tension=s.time;}
        if(s.max_abs_shear_kip>r.max_abs_shear_kip){r.max_abs_shear_kip=s.max_abs_shear_kip;r.time_at_max_abs_shear=s.time;}
    }

    const Frame* frame_{};std::vector<ModelComparisonComponent> components_;ModelComparisonOptions options_;
    std::vector<ModelComparisonComponentResult> results_;std::vector<ModelComparisonColumnResult> column_results_;std::vector<ModelComparisonStorySample> story_history_;
    std::vector<double> peak_story_,peak_story_time_;std::vector<PeakState> peak_states_;
    double peak_roof_{},peak_roof_time_{};std::set<std::size_t> requested_steps_;std::map<std::size_t,PeakState> requested_states_;
    std::map<std::pair<std::size_t,double>,std::pair<std::string,PeakState>> mode_states_;
    std::vector<double> final_u_,final_state_;
};

template<class Frame>
std::vector<std::pair<int,NonlinearMaterial>> capture_baseline_for_variants(
    const Frame& frame,const std::vector<ModelComparisonVariant>& variants){
    std::set<int> ids;for(const auto& v:variants)for(const auto& [id,m]:v.replacements){(void)m;ids.insert(id);}
    std::vector<std::pair<int,NonlinearMaterial>> out;out.reserve(ids.size());
    for(int id:ids){const int idx=frame.nonlinear_component_index(id);if(idx<0)throw std::invalid_argument("comparison variant references unknown nonlinear component");out.push_back({id,frame.materials()[static_cast<std::size_t>(idx)]});}
    return out;
}

template<class Frame>
ModelComparisonResult run_comparison(Frame& frame,const std::vector<ModelComparisonVariant>& variants,
                                     const std::vector<ModelComparisonComponent>& monitored,
                                     const RCNativeFrameAnalysisSettings& settings,const ModelComparisonOptions& options){
    if(variants.empty())throw std::invalid_argument("comparison requires at least one variant");
    if(settings.dt<=0.0||settings.ground_accel.empty())throw std::invalid_argument("comparison analysis settings invalid");
    std::set<std::string> names;for(const auto& v:variants){if(v.name.empty()||!names.insert(v.name).second)throw std::invalid_argument("comparison variant names must be unique/nonempty");}
    const auto baseline=capture_baseline_for_variants(frame,variants);
    PreparedRobustNewmark prepared(frame,settings.dt,settings.strategy,settings.newmark_options.max_subdivisions);
    ModelComparisonResult out;out.variants.reserve(variants.size());
    for(const auto& v:variants){
        if(!baseline.empty())frame.replace_nonlinear_material_field(baseline);
        if(!v.replacements.empty())frame.replace_nonlinear_material_field(v.replacements);
        ModelComparisonVariantResult vr;vr.name=v.name;vr.provenance=v.provenance;
        ComparisonRecorder<Frame> recorder(frame,monitored,options);
        auto ro=settings.newmark_options;const auto chained_state=ro.accepted_state_observer;const auto chained_sub=ro.accepted_substep_state_observer;
        ro.accepted_state_observer=[&](std::size_t step,double time,double ag,const std::vector<double>& u,const std::vector<double>& vel,const std::vector<double>& acc,const std::vector<double>& state){
            recorder.observe_output(step,time,u,state);if(chained_state)chained_state(step,time,ag,u,vel,acc,state);
        };
        if(options.include_accepted_substeps_in_envelopes||chained_sub)ro.accepted_substep_state_observer=[&](std::size_t step,std::size_t sub,int depth,double time,double ag,const std::vector<double>& u,const std::vector<double>& vel,const std::vector<double>& acc,const std::vector<double>& state){
            if(options.include_accepted_substeps_in_envelopes) recorder.observe_envelope(step,time,u,state);
            if(chained_sub) chained_sub(step,sub,depth,time,ag,u,vel,acc,state);
        };
        vr.analysis=prepared.run(settings.ground_accel,ro);vr.reused_prepared_solver=prepared.last_run_reused_preparation();vr.preparation_count_after_run=prepared.preparation_count();recorder.finalize(vr);out.variants.push_back(std::move(vr));
    }
    if(!baseline.empty())frame.replace_nonlinear_material_field(baseline);
    out.prepared_solver_preparations=prepared.preparation_count();out.prepared_solver_setup_factorizations=prepared.setup_factorizations();return out;
}

} // namespace

ModelComparisonVariant make_rc_column_model_comparison_variant(
    std::string name,std::string provenance,const std::vector<RCColumnElementBinding>& bindings,
    const std::vector<RCBuildingColumnModel>& models){
    std::map<std::string,const RCBuildingColumnModel*> by_id;
    for(const auto& m:models){
        if(m.component_id.empty()) throw std::invalid_argument("RC comparison model id is empty");
        if(!by_id.emplace(m.component_id,&m).second) throw std::invalid_argument("duplicate RC comparison model id");
    }
    ModelComparisonVariant out;out.name=std::move(name);out.provenance=std::move(provenance);std::set<int> hinge_ids;std::set<std::string> bound_models;
    for(const auto& b:bindings){
        if(b.component_id.empty()) throw std::invalid_argument("RC comparison binding id is empty");
        auto it=by_id.find(b.component_id);
        if(it==by_id.end()) throw std::invalid_argument("RC comparison binding missing model: "+b.component_id);
        bound_models.insert(b.component_id);
        for(int hid:b.hinge_component_ids){
            if(!hinge_ids.insert(hid).second) throw std::invalid_argument("RC comparison hinge assigned more than once");
            out.replacements.push_back({hid,NonlinearMaterial(ASCE41HingeMaterial(it->second->model.hinge))});
        }
    }
    for(const auto& [id,model]:by_id){
        (void)model;
        if(!bound_models.count(id)) throw std::invalid_argument("RC comparison model has no physical binding: "+id);
    }
    if(out.replacements.empty()) throw std::invalid_argument("RC comparison variant contains no hinge replacements");
    return out;
}

ModelComparisonResult RCFrameModelComparison::run_frame2d(CompiledFrame2D& frame,const std::vector<ModelComparisonVariant>& variants,
    const std::vector<ModelComparisonComponent>& monitored,const RCNativeFrameAnalysisSettings& settings,const ModelComparisonOptions& options){return run_comparison(frame,variants,monitored,settings,options);}

ModelComparisonResult RCFrameModelComparison::run_frame3d(CompiledFrame3D& frame,const std::vector<ModelComparisonVariant>& variants,
    const std::vector<ModelComparisonComponent>& monitored,const RCNativeFrameAnalysisSettings& settings,const ModelComparisonOptions& options){return run_comparison(frame,variants,monitored,settings,options);}

} // namespace quake
