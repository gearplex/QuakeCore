#include "quake/ida.hpp"
#include "quake/analysis_input.hpp"
#include <memory>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace quake {
namespace {

double pga(const std::vector<double>& a){double m=0.0;for(double v:a)m=std::max(m,std::abs(v));return m;}
std::vector<double> scaled(const std::vector<double>& a,double s){std::vector<double> out(a.size());for(std::size_t i=0;i<a.size();++i)out[i]=s*a[i];return out;}

IDARunResult one_run(const NonlinearDynamicModel& model,const GroundMotionRecord& rec,std::size_t ri,double scale,const IDAOptions& o,PreparedRobustNewmark* prepared){
    auto gm=scaled(rec.acceleration,scale);auto analysis_options=o.analysis;analysis_options.return_numerical_failure=true;auto a=prepared?prepared->run(gm,analysis_options):run_newmark_robust(model,gm,rec.dt,o.strategy,analysis_options);IDARunResult r;r.record_index=ri;r.record_name=rec.name;r.scale_factor=scale;r.pga=pga(rec.acceleration)*scale;r.termination=a.termination;r.collapse_mechanism=a.collapse_mechanism;r.termination_step=a.termination_step;r.termination_time=a.termination_time;r.max_story_drift_ratio=a.max_story_drift_ratio;r.max_roof_abs=a.stats.max_roof_abs;r.newton_iterations=a.stats.newton_iterations;r.global_factorizations=a.stats.global_factorizations;r.direct_fallbacks=a.stats.direct_fallbacks;r.elapsed_seconds=a.stats.elapsed_seconds;r.termination_reason=a.termination_reason;return r;
}

struct RecordWork { IDARecordSummary summary; std::vector<IDARunResult> runs; };
RecordWork run_record(const NonlinearDynamicModel& model,const GroundMotionRecord& rec,std::size_t ri,const IDAOptions& o){
    validate_transient_input(rec.acceleration,rec.dt,o.analysis.tolerance,o.analysis.max_iterations);
    RecordWork w;w.summary.record_index=ri;w.summary.record_name=rec.name;w.summary.base_pga=pga(rec.acceleration);
    std::vector<double> scales=o.scale_factors;
    if(scales.empty() || std::any_of(scales.begin(),scales.end(),[](double v){return !std::isfinite(v)||v<=0.0;}))
        throw std::invalid_argument("IDA requires finite positive scale factors");
    std::sort(scales.begin(),scales.end());scales.erase(std::unique(scales.begin(),scales.end()),scales.end());
    std::unique_ptr<PreparedRobustNewmark> prepared;
    if(o.reuse_preparation)prepared=std::make_unique<PreparedRobustNewmark>(model,rec.dt,o.strategy,o.analysis.max_subdivisions);
    auto run=[&](double scale){return one_run(model,rec,ri,scale,o,prepared.get());};
    double low=0.0,high=std::numeric_limits<double>::quiet_NaN();CollapseMechanism high_cause=CollapseMechanism::None;
    for(double scale:scales){
        auto r=run(scale);w.runs.push_back(r);
        if(r.termination==AnalysisTermination::NumericalFailure||r.termination==AnalysisTermination::InitialInstability){
            w.summary.numerical_failure=true;if(r.termination==AnalysisTermination::InitialInstability)break;continue;
        }
        if(r.termination==AnalysisTermination::PhysicalCollapse){
            if(!std::isfinite(high)){high=scale;high_cause=r.collapse_mechanism;}
            if(o.stop_after_first_collapse)break;
        }else if(!std::isfinite(high))low=scale;
        else w.summary.nonmonotonic_response=true;
    }
    if(std::isfinite(high)&&low>0.0&&o.collapse_refinement_steps>0){
        for(int k=0;k<o.collapse_refinement_steps;++k){
            const double mid=std::exp(0.5*(std::log(low)+std::log(high)));
            if(!(low<mid&&mid<high))break;
            auto r=run(mid);w.runs.push_back(r);
            if(r.termination==AnalysisTermination::PhysicalCollapse){high=mid;high_cause=r.collapse_mechanism;}
            else if(r.termination==AnalysisTermination::Completed)low=mid;
            else {w.summary.numerical_failure=true;break;}
        }
    }
    if(std::isfinite(high))for(const auto& r:w.runs)
        if(r.scale_factor>low&&r.scale_factor<high&&
           (r.termination==AnalysisTermination::NumericalFailure||r.termination==AnalysisTermination::InitialInstability))
            w.summary.bracket_has_numerical_gap=true;
    w.summary.last_noncollapse_scale=low;
    if(std::isfinite(high)){w.summary.collapsed=true;w.summary.first_collapse_scale=high;w.summary.collapse_pga=high*w.summary.base_pga;w.summary.collapse_mechanism=high_cause;}
    else {w.summary.right_censored=!w.summary.numerical_failure;}
    std::sort(w.runs.begin(),w.runs.end(),[](const IDARunResult&a,const IDARunResult&b){return a.scale_factor<b.scale_factor;});return w;
}

} // namespace

IDASuiteResult run_ida_suite(const NonlinearDynamicModel& model,const std::vector<GroundMotionRecord>& records,const IDAOptions& options){
    if(options.workers<=0||options.collapse_refinement_steps<0)
        throw std::invalid_argument("invalid IDA options");
    IDASuiteResult out;
    out.records.resize(records.size());
    if(records.empty()){out.workers=0;return out;}
    const int workers=std::max(1,std::min(options.workers,static_cast<int>(records.size())));
    out.workers=workers;
    std::vector<RecordWork> work(records.size());std::atomic<std::size_t> next{0};std::exception_ptr failure;std::mutex fm;const auto start=std::chrono::steady_clock::now();std::vector<std::thread> threads;threads.reserve(static_cast<std::size_t>(workers));
    for(int t=0;t<workers;++t)threads.emplace_back([&]{try{while(true){const std::size_t i=next.fetch_add(1);if(i>=records.size())break;work[i]=run_record(model,records[i],i,options);}}catch(...){std::lock_guard<std::mutex> lock(fm);if(!failure)failure=std::current_exception();next.store(records.size());}});
    for(auto&t:threads) t.join();
    if(failure) std::rethrow_exception(failure);
    out.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::vector<double> logcaps;
    for(std::size_t i=0;i<work.size();++i){out.records[i]=work[i].summary;out.runs.insert(out.runs.end(),work[i].runs.begin(),work[i].runs.end());if(work[i].summary.collapsed){++out.physical_collapses;if(work[i].summary.collapse_pga>0.0)logcaps.push_back(std::log(work[i].summary.collapse_pga));}if(work[i].summary.numerical_failure)++out.numerical_failures;if(work[i].summary.right_censored)++out.right_censored_records;}
    std::sort(out.runs.begin(),out.runs.end(),[](const IDARunResult&a,const IDARunResult&b){return a.record_index<b.record_index||(a.record_index==b.record_index&&a.scale_factor<b.scale_factor);});
    if(!logcaps.empty()){const double mean=std::accumulate(logcaps.begin(),logcaps.end(),0.0)/logcaps.size();out.collapse_pga_median=std::exp(mean);if(logcaps.size()>1){double ss=0.0;for(double x:logcaps)ss+=(x-mean)*(x-mean);out.collapse_pga_log_stddev=std::sqrt(ss/(logcaps.size()-1));}else out.collapse_pga_log_stddev=0.0;}
    return out;
}

} // namespace quake
