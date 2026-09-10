#define main quakecore_phase8_embedded_main
#include "run_validation_phase8.cpp"
#undef main
#include "quake/modal_damping.hpp"
#include "json_helpers.hpp"
#include <memory>
int main(int argc,char**argv){try{
 if(argc<5){std::cerr<<"Usage: ucb_phase9i motion.csv result.json case duration_s [full|same_pattern|woodbury] [exact_modal_tangent=1] [dt_refinement=1] [kinematic|previous_displacement]\n";return 2;}
 std::string variant=argv[3];if(variant!="epp_rayleigh"&&variant!="mroz_rayleigh"&&variant!="hybrid"&&variant!="modal3")throw std::invalid_argument("unknown Berkeley case");
 const double duration=std::stod(argv[4]);if(!std::isfinite(duration)||duration<=0)throw std::invalid_argument("invalid duration");auto motion=read_motion(argv[1]);
 std::size_t n=std::min(motion.ag.size(),static_cast<std::size_t>(std::floor(duration/motion.dt+1e-8)));motion.ag.resize(n);
 int refinement=argc>7?std::stoi(argv[7]):1;if(refinement<1||refinement>32)throw std::invalid_argument("dt refinement must be 1..32");
 if(refinement>1){std::vector<double> fine;double prev=0;for(double next:motion.ag){for(int j=1;j<=refinement;++j)fine.push_back(prev+(next-prev)*j/refinement);prev=next;}motion.ag=std::move(fine);motion.dt/=refinement;}
 const double w=2*M_PI/.48,alpha=(variant=="mroz_rayleigh"||variant=="epp_rayleigh")?2*.025*w:0,beta=variant=="modal3"?0:2*.005/w;
 auto raw=build_frame_phase8(.934616,gravity_column_preloads(),true,true,false,true,false,3,alpha,beta,true,true,-.005,.20,.55,-1,true,variant!="epp_rayleigh");
 auto& frame=raw.model;std::unique_ptr<FixedModalDampingModel> modal;bool exact=argc>6?std::stoi(argv[6])!=0:true;
 if(variant=="hybrid"||variant=="modal3")modal=std::make_unique<FixedModalDampingModel>(frame,variant=="hybrid"?.025:.030,3,exact);
 const NonlinearDynamicModel& model=modal?static_cast<const NonlinearDynamicModel&>(*modal):static_cast<const NonlinearDynamicModel&>(frame);
 RobustNewmarkOptions o;o.tolerance=1e-7;o.max_iterations=50;o.max_backtracks=10;o.nonmonotone_line_search_factor=1.05;o.max_subdivisions=6;o.return_numerical_failure=true;
 const std::string guess=argc>8?argv[8]:"previous_displacement";if(guess!="kinematic"&&guess!="previous_displacement")throw std::invalid_argument("invalid initial guess");o.kinematic_initial_guess=guess=="kinematic";
 std::vector<double> peak(3,0),peak_time(3,0),at_s2;double b1_rotation=0,b1_axial=0;Json history=Json::array();
 auto drift=[&](const auto& u){auto d=frame.story_response_values(u);double lo=0;for(double& x:d){double old=x;x=(x-lo)/39.*100;lo=old;}return d;};
 o.accepted_substep_state_observer=[&](std::size_t,std::size_t,int,double time,double,const auto& u,const auto&,const auto&,const auto& state){auto d=drift(u);for(int i=0;i<3;++i)if(std::abs(d[i])>peak[i]){peak[i]=std::abs(d[i]);peak_time[i]=time;if(i==1)at_s2=d;}auto h=frame.pm_interaction_snapshot(3,u,state);b1_rotation=std::max(b1_rotation,std::abs(h.plastic_rotation));b1_axial=std::max(b1_axial,std::abs(h.plastic_axial_deformation));};
 o.accepted_state_observer=[&](std::size_t,double time,double,const auto& u,const auto&,const auto&,const auto& state){auto h=frame.pm_interaction_snapshot(3,u,state);history.push_back({{"time_s",time},{"story_drift_percent",drift(u)},{"b1_rotation",h.plastic_rotation},{"b1_axial",h.plastic_axial_deformation},{"b1_P",h.compression},{"b1_M",h.moment}});};
 auto strategy=parse_strategy(argc>5?argv[5]:"same_pattern");const auto ar=run_newmark_robust(model,motion.ag,motion.dt,strategy,o);Json result=analysis_json(ar);
 result["newton_initial_guess"]=guess;result["final_material_state"]=ar.final_state;result["initial_acceleration_sample_g"]=motion.initial_acceleration_g;result["initial_state_convention"]="u=v=a=0; t=0 acceleration omitted, preserving archived proxy convention";
 result["phase"]="9I";result["case"]=variant;result["validation_status"]="proxy_input_only";result["input_file"]=argv[1];result["dt_s"]=motion.dt;result["requested_duration_s"]=duration;result["available_duration_s"]=motion.ag.size()*motion.dt;result["drift_denominator_in"]=39;result["exact_modal_tangent"]=exact;result["peak_story_drift_percent"]=peak;result["peak_time_s"]=peak_time;result["simultaneous_at_story2_peak_percent"]=at_s2;result["b1_peak_plastic_rotation_rad"]=b1_rotation;result["b1_peak_plastic_axial_in"]=b1_axial;result["history"]=history;
 result["periods_s"]=Json::array();for(const auto& md:modal_analysis(frame,3))result["periods_s"].push_back(md.period);
 write_json(argv[2],result);std::cout<<variant<<' '<<result["termination"]<<' '<<peak[0]<<','<<peak[1]<<','<<peak[2]<<'\n';return ar.termination==AnalysisTermination::Completed?0:3;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
