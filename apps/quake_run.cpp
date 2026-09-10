#include "json_helpers.hpp"
#include "wall_json.hpp"
#include "steel_json.hpp"
#include "soil_json.hpp"
#include "quake/frame2d.hpp"
#include "quake/ida.hpp"
#include "quake/modal.hpp"
#include "quake/static_analysis.hpp"
#include <chrono>
#include <iostream>
#include <map>
#include <set>
using namespace quake;
static void keys(const Json& j,std::initializer_list<std::string> allowed){if(!j.is_object())throw std::invalid_argument("expected JSON object");for(auto it=j.begin();it!=j.end();++it)if(std::find(allowed.begin(),allowed.end(),it.key())==allowed.end())throw std::invalid_argument("unknown field: "+it.key());}
static int checked_id(const Json& j){if(!j.is_number_integer()||j.get<double>()<std::numeric_limits<int>::min()||j.get<double>()>std::numeric_limits<int>::max())throw std::invalid_argument("identifiers must be 32-bit integers");return j.get<int>();}
static Dof2D dof(const std::string& s){if(s=="UX")return Dof2D::UX;if(s=="UY")return Dof2D::UY;if(s=="RZ")return Dof2D::RZ;throw std::invalid_argument("unsupported 2D DOF "+s);}
int main(int argc,char** argv){try{
 if(argc!=3){std::cerr<<"Usage: quake_run job.json result.json\n";return 2;}
 const auto job=read_json(argv[1]);keys(job,{"schema","name","units","model","analysis","records","provenance","gravity"});
 if(job.at("schema")!="quakecore.job.v1")throw std::invalid_argument("unsupported job schema");
 keys(job.at("units"),{"force","length","time"});if(job.at("units").at("time")!="s")throw std::invalid_argument("time unit must be s");
 for(const auto* name:{"force","length"})if(!job.at("units").at(name).is_string()||job.at("units").at(name).get<std::string>().empty())throw std::invalid_argument("unit labels must be nonempty strings");
 const auto& m=job.at("model");keys(m,{"type","nodes","fixities","equal_dofs","members","hinges","rayleigh","response_node","story_nodes","story_cut_members","walls","steel_members","panel_zones","brbs","viscous_dampers","soil_springs","soil_dashpots","pile_lines"});
 if(m.at("type")!="frame2d")throw std::invalid_argument("v1 job runner supports frame2d");
 Frame2DBuilder b;std::map<int,std::pair<double,double>> coords;std::set<int> elements;std::vector<int> hinge_ids,panel_zone_ids,brb_ids,damper_ids,steel_ids,soil_spring_ids,soil_dashpot_ids;std::vector<PileLineInfo> pile_lines;
 for(const auto& n:m.at("nodes")){
  if(n.size()!=6)throw std::invalid_argument("node must be [id,x,y,mx,my,mr]");int id=checked_id(n[0]);if(coords.contains(id))throw std::invalid_argument("duplicate node id");
  double mx=n[3],my=n[4],mr=n[5];if(mx<0||my<0||mr<0)throw std::invalid_argument("negative nodal mass");
  coords[id]={n[1],n[2]};b.add_node(id,n[1],n[2],mx,my,mr);
 }
 if(m.contains("pile_lines")&&!m.at("pile_lines").is_array())throw std::invalid_argument("pile_lines must be an array");
 {std::set<std::string> pile_names;for(const auto& p:m.value("pile_lines",Json::array())){auto info=add_vertical_pile_line(p,b,coords,elements,soil_spring_ids);if(!pile_names.insert(info.id).second)throw std::invalid_argument("duplicate pile id");pile_lines.push_back(std::move(info));}}
 for(const auto& f:m.at("fixities")){if(f.size()!=4)throw std::invalid_argument("fixity must be [id,ux,uy,rz]");b.fix(checked_id(f[0]),f[1].get<bool>(),f[2].get<bool>(),f[3].get<bool>());}
 for(const auto& e:m.value("equal_dofs",Json::array())){if(e.size()!=3)throw std::invalid_argument("equal DOF must be [retained,constrained,dof]");b.equal_dof(checked_id(e[0]),checked_id(e[1]),dof(e[2]));}
 for(const auto& e:m.at("members")){
  if(e.size()!=7&&e.size()!=8)throw std::invalid_argument("member must be [id,i,j,E,A,I,constant_compression,pdelta_transformation?]");int id=checked_id(e[0]);if(!elements.insert(id).second)throw std::invalid_argument("duplicate element id");
  if(e[3].get<double>()<=0||e[4].get<double>()<=0||e[5].get<double>()<=0)throw std::invalid_argument("member E,A,I must be positive");
  b.add_elastic_frame(id,checked_id(e[1]),checked_id(e[2]),e[3],e[4],e[5],e[6],e.size()==8?e[7].get<bool>():false);
 }
 for(const auto& h:m.value("hinges",Json::array())){
  keys(h,{"id","i","j","type","k","fy","hardening_ratio","a","b","f","c","io","ls","cp","hardening_stiffness","provenance"});
  int id=checked_id(h.at("id"));checked_id(h.at("i"));checked_id(h.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate element/hinge id");hinge_ids.push_back(id);
  const auto pi=coords.at(h.at("i").get<int>()),pj=coords.at(h.at("j").get<int>());
  if(std::hypot(pi.first-pj.first,pi.second-pj.second)>1e-10)throw std::invalid_argument("v1 rotational hinges must connect coincident nodes");
  std::string type=h.at("type");
  if(type=="bilinear")b.add_rotational_spring(id,h.at("i"),h.at("j"),h.at("k"),h.at("fy"),h.at("hardening_ratio"));
  else if(type=="asce41_parameterized"){
   if(h.value("provenance",std::string()).empty())throw std::invalid_argument("ASCE-style hinge requires parameter provenance");
   ASCE41HingeParams p;p.Ke=h.at("k");p.posFy=p.negFy=h.at("fy");p.pos_a=p.neg_a=h.at("a");p.pos_b=p.neg_b=h.at("b");p.pos_f=p.neg_f=h.value("f",0.0);p.pos_c=p.neg_c=h.value("c",0.0);
   p.pos_io=p.neg_io=h.value("io",0.0);p.pos_ls=p.neg_ls=h.value("ls",0.0);p.pos_cp=p.neg_cp=h.value("cp",0.0);p.hardening_stiffness=h.value("hardening_stiffness",-1.0);p.hardening_ratio=h.value("hardening_ratio",0.0);p.backbone_shape=ASCE41BackboneShape::StraightCE;
   b.add_asce41_hinge(id,h.at("i"),h.at("j"),p);
  }else throw std::invalid_argument("unsupported hinge type "+type);
 }
 if(m.contains("steel_members")&&!m.at("steel_members").is_array())throw std::invalid_argument("steel_members must be an array");
 for(const auto& e:m.value("steel_members",Json::array())){
  int id=checked_id(e.at("id")),i=checked_id(e.at("i")),j=checked_id(e.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate steel member/element id");
  if(e.at("type")!="concentrated_plasticity")throw std::invalid_argument("unsupported steel member type");
  const std::string role=e.at("role");SteelMember2DRole r;if(role=="beam")r=SteelMember2DRole::Beam;else if(role=="column")r=SteelMember2DRole::Column;else throw std::invalid_argument("steel member role must be beam or column");
  if(e.value("provenance",std::string()).empty())throw std::invalid_argument("steel member requires parameter provenance");
  b.add_steel_member(id,i,j,r,steel_member_properties(e));steel_ids.push_back(id);
 }
 if(m.contains("panel_zones")&&!m.at("panel_zones").is_array())throw std::invalid_argument("panel_zones must be an array");
 for(const auto& p:m.value("panel_zones",Json::array())){
  steel_keys(p,{"id","i","j","material","provenance"});int id=checked_id(p.at("id")),i=checked_id(p.at("i")),j=checked_id(p.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate panel zone/element id");
  if(p.value("provenance",std::string()).empty())throw std::invalid_argument("panel zone requires parameter provenance");
  const auto pi=coords.at(i),pj=coords.at(j);if(std::hypot(pi.first-pj.first,pi.second-pj.second)>1e-10)throw std::invalid_argument("panel zone nodes must be coincident");
  b.add_panel_zone(id,i,j,steel_nonlinear_material(p.at("material")));panel_zone_ids.push_back(id);
 }
 if(m.contains("brbs")&&!m.at("brbs").is_array())throw std::invalid_argument("brbs must be an array");
 for(const auto& d:m.value("brbs",Json::array())){
  steel_keys(d,{"id","i","j","k","fy","hardening_ratio","provenance"});int id=checked_id(d.at("id")),i=checked_id(d.at("i")),j=checked_id(d.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate BRB/element id");
  if(d.value("provenance",std::string()).empty())throw std::invalid_argument("BRB requires parameter provenance");
  b.add_brb(id,i,j,d.at("k"),d.at("fy"),d.at("hardening_ratio"));brb_ids.push_back(id);
 }
 if(m.contains("viscous_dampers")&&!m.at("viscous_dampers").is_array())throw std::invalid_argument("viscous_dampers must be an array");
 for(const auto& d:m.value("viscous_dampers",Json::array())){
  steel_keys(d,{"id","i","j","coefficient","alpha","regularization_velocity","provenance"});int id=checked_id(d.at("id")),i=checked_id(d.at("i")),j=checked_id(d.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate viscous damper/element id");
  if(d.value("provenance",std::string()).empty())throw std::invalid_argument("viscous damper requires parameter provenance");
  b.add_viscous_damper(id,i,j,{d.at("coefficient"),d.value("alpha",1.0),d.value("regularization_velocity",0.0)});damper_ids.push_back(id);
 }
 if(m.contains("soil_springs")&&!m.at("soil_springs").is_array())throw std::invalid_argument("soil_springs must be an array");
 for(const auto& s:m.value("soil_springs",Json::array())){
  soil_keys(s,{"id","i","j","direction","dof","type","k","fy","hardening_ratio","positive_capacity","negative_capacity","ultimate_resistance_per_length","ultimate_interface_stress","ultimate_bearing_stress","tributary_length","pile_perimeter","toe_area","displacement_50","suction_ratio","compression_sign","provenance"});require_soil_provenance(s);
  int id=soil_id(s.at("id")),i=soil_id(s.at("i")),j=soil_id(s.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate soil spring/element id");if(s.contains("direction")==s.contains("dof"))throw std::invalid_argument("soil spring requires exactly one of direction or dof");if(s.contains("dof")){if(s.at("dof")!="RZ"||(s.at("type")!="bilinear"&&s.at("type")!="asymmetric_elastic_perfectly_plastic"))throw std::invalid_argument("rotational soil spring supports dof RZ and direct material types only");b.add_rotational_soil_spring(id,i,j,soil_material(s));}else{auto d=soil_direction(s.at("direction"));b.add_soil_spring(id,i,j,d.first,d.second,soil_material(s));}soil_spring_ids.push_back(id);
 }
 if(m.contains("soil_dashpots")&&!m.at("soil_dashpots").is_array())throw std::invalid_argument("soil_dashpots must be an array");
 for(const auto& d:m.value("soil_dashpots",Json::array())){
  soil_keys(d,{"id","i","j","direction","dof","coefficient","alpha","regularization_velocity","provenance"});require_soil_provenance(d);int id=soil_id(d.at("id")),i=soil_id(d.at("i")),j=soil_id(d.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate soil dashpot/element id");if(d.contains("direction")==d.contains("dof"))throw std::invalid_argument("soil dashpot requires exactly one of direction or dof");ViscousDamper2DProperties vp{d.at("coefficient"),d.value("alpha",1.0),d.value("regularization_velocity",0.0)};if(d.contains("dof")){if(d.at("dof")!="RZ")throw std::invalid_argument("rotational soil dashpot dof must be RZ");b.add_rotational_soil_dashpot(id,i,j,vp);}else{auto v=soil_direction(d.at("direction"));b.add_directional_viscous_damper(id,i,j,v.first,v.second,vp);}soil_dashpot_ids.push_back(id);
 }
 if(m.contains("walls")&&!m.at("walls").is_array())throw std::invalid_argument("walls must be an array");
 for(const auto& w:m.value("walls",Json::array())){
  int id=checked_id(w.at("id")),i=checked_id(w.at("i")),j=checked_id(w.at("j"));if(!elements.insert(id).second)throw std::invalid_argument("duplicate wall/element id");
  std::string t=w.at("type");if(t=="mvlem")b.add_mvlem(id,i,j,mvlem_properties(w));else if(t=="sfi_mvlem")b.add_sfi_mvlem(id,i,j,sfi_mvlem_properties(w));else throw std::invalid_argument("unsupported wall type");
 }
 auto damping=m.value("rayleigh",std::vector<double>{0,0});if(damping.size()!=2||damping[0]<0||damping[1]<0)throw std::invalid_argument("Rayleigh requires nonnegative [alphaM,betaKinitial]");b.set_rayleigh(damping[0],damping[1]);b.set_response_node(checked_id(m.at("response_node")));
 for(const auto& n:m.at("story_nodes"))checked_id(n);auto floors=m.at("story_nodes").get<std::vector<int>>();if(floors.empty())throw std::invalid_argument("story_nodes must be nonempty");double lower_y=0;for(int node:floors){double y=coords.at(node).second;if(y<=lower_y)throw std::invalid_argument("story nodes must have strictly increasing positive elevations");lower_y=y;}b.set_story_nodes(floors);auto model=b.compile();
 std::vector<std::vector<int>> cuts=m.value("story_cut_members",std::vector<std::vector<int>>{});
 for(const auto& cut:cuts)for(int id:cut){bool found=false;for(const auto& e:m.at("members"))if(e[0].get<int>()==id){auto pi=coords.at(e[1].get<int>()),pj=coords.at(e[2].get<int>());if(std::abs(pi.first-pj.first)>1e-10||pj.second<=pi.second)throw std::invalid_argument("story cuts require upward vertical members");found=true;}for(const auto& e:m.value("steel_members",Json::array()))if(e.at("id").get<int>()==id){auto pi=coords.at(e.at("i").get<int>()),pj=coords.at(e.at("j").get<int>());if(e.at("role")!="column"||std::abs(pi.first-pj.first)>1e-10||pj.second<=pi.second)throw std::invalid_argument("steel story cuts require upward vertical column members");found=true;}if(!found)throw std::invalid_argument("unknown story cut member");}
 auto shears=[&](const auto& u,const auto& state){std::vector<double> v;for(const auto& cut:cuts){double sum=0;for(int id:cut){if(model.elastic_element_index(id)>=0)sum-=model.elastic_element_response(id,u).shear_i;else sum-=model.steel_member_response(id,u,state).force[0];}v.push_back(sum);}return v;};
 const auto& control=job.at("analysis");for(const auto* name:{"max_iterations","max_subdivisions","workers","refinements","history_stride"})if(control.contains(name))checked_id(control.at(name));keys(control,{"type","strategy","tolerance","relative_tolerance","initial_guess","max_iterations","max_subdivisions","line_search","history_stride","output_mode","scales","refinements","workers","drift_limit","stop_after_first_collapse","check_initial_stability"});
 RobustNewmarkOptions o;o.tolerance=control.value("tolerance",1e-8);o.relative_force_tolerance=control.value("relative_tolerance",true);o.max_iterations=control.value("max_iterations",35);o.max_subdivisions=control.value("max_subdivisions",5);o.line_search=control.value("line_search",true);o.return_numerical_failure=true;o.collapse.max_story_drift_ratio=control.value("drift_limit",0.0);o.collapse.check_initial_stability=control.value("check_initial_stability",true);if(o.collapse.max_story_drift_ratio<0)throw std::invalid_argument("negative drift_limit");
 const std::string guess=control.value("initial_guess",std::string("kinematic"));if(guess!="kinematic"&&guess!="previous_displacement")throw std::invalid_argument("unknown initial_guess");o.kinematic_initial_guess=guess=="kinematic";
 auto strategy=parse_strategy(control.value("strategy",std::string("woodbury")));
 std::vector<double> gravity_load(static_cast<std::size_t>(model.dof()),0.0);StaticAnalysisResult gravity_result;bool gravity_performed=false;
 if(job.contains("gravity")){
  const auto& g=job.at("gravity");keys(g,{"loads","steps","tolerance","max_iterations","provenance"});
  if(g.value("provenance",std::string()).empty())throw std::invalid_argument("gravity analysis requires provenance");
  for(const auto& load:g.at("loads")){
   if(!load.is_array()||load.size()!=4)throw std::invalid_argument("gravity load must be [node,fx,fy,mz]");
   const int node=checked_id(load[0]);coords.at(node);
   for(int j=0;j<3;++j){const int rd=model.reduced_dof(node,static_cast<Dof2D>(j));if(rd>=0)gravity_load[static_cast<std::size_t>(rd)]+=load[static_cast<std::size_t>(j+1)].get<double>();}
  }
  const int steps=g.value("steps",100),maxit=g.value("max_iterations",100);if(steps<1||maxit<1)throw std::invalid_argument("invalid gravity step controls");
  gravity_result=solve_static_load(model,gravity_load,steps,g.value("tolerance",1e-8),maxit);
  if(!gravity_result.converged)throw std::runtime_error("gravity analysis failed to converge");
  gravity_performed=true;o.initial_displacement=gravity_result.displacement;o.initial_committed_state=gravity_result.committed_state;o.constant_load=gravity_load;
 }
 std::vector<GroundMotionRecord> records;std::set<std::string> record_names;
 for(const auto& r:job.at("records")){
  keys(r,{"name","dt","acceleration","sample_convention","provenance"});if(r.at("sample_convention")!="step_end")throw std::invalid_argument("acceleration samples must be a(dt), a(2dt), ...; the initial state assumes a(0)=0");
  GroundMotionRecord rec{r.at("name"),r.at("acceleration").get<std::vector<double>>(),r.at("dt")};if(!record_names.insert(rec.name).second)throw std::invalid_argument("duplicate record name");records.push_back(std::move(rec));
 }
 if(records.empty())throw std::invalid_argument("at least one record is required");
 Json out={{"schema","quakecore.result.v1"},{"name",job.value("name",std::string())},{"units",job.at("units")},{"engine","QuakeCore Phase 9L soil springs and piles research"},{"code_compliance","not_assessed"},{"model",m},{"dof",model.dof()},{"provenance",job.value("provenance",Json::object())}};
 const auto start=std::chrono::steady_clock::now();out["modes"]=Json::array();std::vector<ModeShape> modes;
 if(gravity_performed){std::vector<double> force,tangents,trial;model.internal_force_and_tangent(gravity_result.displacement,gravity_result.committed_state,force,tangents,trial);auto kt=model.effective_state_tangent_matrix_with_state(gravity_result.displacement,tangents,gravity_result.committed_state,0.0,0.0);modes=modal_analysis_from_stiffness(model,kt,std::min(3,static_cast<int>(floors.size())));}else modes=modal_analysis(model,std::min(3,static_cast<int>(floors.size())));
 for(const auto& mode:modes)out["modes"].push_back({{"period_s",mode.period},{"floor_shape",model.story_response_values(mode.shape)}});
 if(gravity_performed)out["gravity"]={{"converged",true},{"steps",gravity_result.load_steps_completed},{"newton_iterations",gravity_result.newton_iterations},{"residual_inf_norm",gravity_result.residual_inf_norm},{"reduced_load",gravity_load},{"displacement",gravity_result.displacement}};
 if(control.at("type")=="nrha"){
  const std::string output_mode=control.value("output_mode",std::string("full"));if(output_mode!="full"&&output_mode!="summary"&&output_mode!="minimal")throw std::invalid_argument("unknown output_mode");out["output_mode"]=output_mode;
  int stride=control.value("history_stride",1);if(stride<1)throw std::invalid_argument("history_stride must be positive");out["runs"]=Json::array();
  for(const auto& rec:records){
   Json history=Json::array();std::vector<double> peak_shear(cuts.size(),0.0);std::vector<double> peak(floors.size(),0.0),peak_acc(peak),peak_time(peak);Json hinges=Json::object(),walls=Json::object(),steel_members=Json::object(),panel_zones=Json::object(),brbs=Json::object(),dampers=Json::object(),soil_springs=Json::object(),soil_dashpots=Json::object(),pile_peaks=Json::object();
   auto drift=[&](const auto& u){auto f=model.story_response_values(u);double lo=0,ly=0;for(std::size_t i=0;i<f.size();++i){double x=f[i],y=coords.at(floors[i]).second;f[i]=(x-lo)/(y-ly);lo=x;ly=y;}return f;};
   o.accepted_substep_state_observer=[&](std::size_t,std::size_t,int,double time,double ag,const auto& u,const auto& vel,const auto& a,const auto& s){auto d=drift(u),acc=model.story_response_values(a),shear=shears(u,s);for(std::size_t i=0;i<shear.size();++i)peak_shear[i]=std::max(peak_shear[i],std::abs(shear[i]));for(std::size_t i=0;i<d.size();++i){if(std::abs(d[i])>peak[i]){peak[i]=std::abs(d[i]);peak_time[i]=time;}peak_acc[i]=std::max(peak_acc[i],std::abs(acc[i]+ag));}
    for(int id:model.wall_ids()){auto r=model.wall_response(id,u,s);auto& w=walls[std::to_string(id)];if(w.is_null())w=Json::object();w["peak_abs_shear"]=std::max(w.value("peak_abs_shear",0.0),std::abs(r.force[3]));w["peak_abs_axial_force"]=std::max(w.value("peak_abs_axial_force",0.0),std::abs(r.force[4]));w["peak_abs_bottom_moment"]=std::max(w.value("peak_abs_bottom_moment",0.0),std::abs(r.force[2]));w["peak_abs_curvature"]=std::max(w.value("peak_abs_curvature",0.0),std::abs(r.curvature));w["peak_abs_shear_deformation"]=std::max(w.value("peak_abs_shear_deformation",0.0),std::abs(r.shear_deformation));}
    for(int id:hinge_ids){auto h=model.nonlinear_component_snapshot(id,u,s);auto& v=hinges[std::to_string(id)];if(v.is_null())v=Json::object();v["peak_abs_rotation"]=std::max(v.value("peak_abs_rotation",0.0),std::abs(h.deformation));v["peak_abs_moment"]=std::max(v.value("peak_abs_moment",0.0),std::abs(h.force));}
    for(int id:steel_ids){auto x=model.steel_member_response(id,u,s);auto& v=steel_members[std::to_string(id)];if(v.is_null())v=Json::object();v["role"]=model.steel_member_role(id)==SteelMember2DRole::Beam?"beam":"column";v["peak_abs_axial_force"]=std::max(v.value("peak_abs_axial_force",0.0),std::abs(x.axial_force));for(int k=0;k<2;++k){const std::string n=std::to_string(k+1);v["peak_abs_hinge_rotation_"+n]=std::max(v.value("peak_abs_hinge_rotation_"+n,0.0),std::abs(x.hinge_rotation[k]));v["peak_abs_end_moment_"+n]=std::max(v.value("peak_abs_end_moment_"+n,0.0),std::abs(x.end_moment[k]));}}
    for(int id:panel_zone_ids){auto x=model.nonlinear_component_snapshot(id,u,s);auto& v=panel_zones[std::to_string(id)];if(v.is_null())v=Json::object();v["peak_abs_rotation"]=std::max(v.value("peak_abs_rotation",0.0),std::abs(x.deformation));v["peak_abs_moment"]=std::max(v.value("peak_abs_moment",0.0),std::abs(x.force));}
    for(int id:brb_ids){auto x=model.nonlinear_component_snapshot(id,u,s);auto& v=brbs[std::to_string(id)];if(v.is_null())v=Json::object();v["peak_abs_axial_deformation"]=std::max(v.value("peak_abs_axial_deformation",0.0),std::abs(x.deformation));v["peak_abs_axial_force"]=std::max(v.value("peak_abs_axial_force",0.0),std::abs(x.force));}
    for(int id:damper_ids){auto x=model.viscous_damper_response(id,vel);auto& v=dampers[std::to_string(id)];if(v.is_null())v=Json::object();v["peak_abs_deformation_rate"]=std::max(v.value("peak_abs_deformation_rate",0.0),std::abs(x.deformation_rate));v["peak_abs_force"]=std::max(v.value("peak_abs_force",0.0),std::abs(x.force));}
    for(int id:soil_spring_ids){auto x=model.nonlinear_component_snapshot(id,u,s);auto& v=soil_springs[std::to_string(id)];if(v.is_null())v=Json::object();v["peak_abs_deformation"]=std::max(v.value("peak_abs_deformation",0.0),std::abs(x.deformation));v["peak_abs_force"]=std::max(v.value("peak_abs_force",0.0),std::abs(x.force));}
    for(int id:soil_dashpot_ids){auto x=model.viscous_damper_response(id,vel);auto& v=soil_dashpots[std::to_string(id)];if(v.is_null())v=Json::object();v["peak_abs_deformation_rate"]=std::max(v.value("peak_abs_deformation_rate",0.0),std::abs(x.deformation_rate));v["peak_abs_force"]=std::max(v.value("peak_abs_force",0.0),std::abs(x.force));}
    for(const auto& p:pile_lines)accumulate_pile_peak_json(pile_peaks[p.id],pile_line_response_json(p,model,u,s));
   };
   o.accepted_state_observer=[&](std::size_t step,double time,double ag,const auto& u,const auto& vel,const auto& a,const auto& s){if(step%stride==0||step+1==rec.acceleration.size()){auto aa=model.story_response_values(a);for(auto& v:aa)v+=ag;Json wall_history=Json::object(),steel_history=Json::object(),panel_history=Json::object(),brb_history=Json::object(),damper_history=Json::object(),soil_history=Json::object(),soil_dashpot_history=Json::object(),pile_history=Json::object();for(int id:model.wall_ids())wall_history[std::to_string(id)]=wall_response_json(model.wall_response(id,u,s));for(int id:steel_ids)steel_history[std::to_string(id)]=steel_member_response_json(model.steel_member_response(id,u,s));for(int id:panel_zone_ids){auto x=model.nonlinear_component_snapshot(id,u,s);panel_history[std::to_string(id)]={{"rotation",x.deformation},{"moment",x.force},{"tangent",x.tangent}};}for(int id:brb_ids){auto x=model.nonlinear_component_snapshot(id,u,s);brb_history[std::to_string(id)]={{"axial_deformation",x.deformation},{"axial_force",x.force},{"tangent",x.tangent}};}for(int id:damper_ids)damper_history[std::to_string(id)]=viscous_damper_response_json(model.viscous_damper_response(id,vel));for(int id:soil_spring_ids){auto x=model.nonlinear_component_snapshot(id,u,s);soil_history[std::to_string(id)]={{"deformation",x.deformation},{"force",x.force},{"tangent",x.tangent}};}for(int id:soil_dashpot_ids)soil_dashpot_history[std::to_string(id)]=viscous_damper_response_json(model.viscous_damper_response(id,vel));for(const auto& p:pile_lines)pile_history[p.id]=pile_line_response_json(p,model,u,s);history.push_back({{"time_s",time},{"floor_displacement",model.story_response_values(u)},{"story_drift_ratio",drift(u)},{"floor_absolute_acceleration",aa},{"story_restoring_shear",shears(u,s)},{"walls",wall_history},{"steel_members",steel_history},{"panel_zones",panel_history},{"brbs",brb_history},{"viscous_dampers",damper_history},{"soil_springs",soil_history},{"soil_dashpots",soil_dashpot_history},{"pile_lines",pile_history}});}};
   if(output_mode!="full")o.accepted_state_observer={};if(output_mode=="minimal")o.accepted_substep_state_observer={};
   auto ar=run_newmark_robust(model,rec.acceleration,rec.dt,strategy,o);auto r=analysis_json(ar);r["record"]=rec.name;r["peak_story_drift_ratio"]=peak;r["peak_time_s"]=peak_time;r["peak_floor_absolute_acceleration"]=peak_acc;r["peak_story_restoring_shear"]=peak_shear;r["hinges"]=hinges;r["walls"]=walls;r["steel_members"]=steel_members;r["panel_zones"]=panel_zones;r["brbs"]=brbs;r["viscous_dampers"]=dampers;r["soil_springs"]=soil_springs;r["soil_dashpots"]=soil_dashpots;r["pile_lines"]=pile_peaks;r["history"]=history;if(output_mode=="minimal")for(auto key:{"peak_story_drift_ratio","peak_time_s","peak_floor_absolute_acceleration","peak_story_restoring_shear","hinges","walls","steel_members","panel_zones","brbs","viscous_dampers","soil_springs","soil_dashpots","pile_lines","history"})r.erase(key);out["runs"].push_back(r);
  }
 }else if(control.at("type")=="ida"){
  IDAOptions io;io.analysis=o;io.strategy=strategy;io.scale_factors=control.at("scales").get<std::vector<double>>();io.workers=control.value("workers",1);io.collapse_refinement_steps=control.value("refinements",3);io.stop_after_first_collapse=control.value("stop_after_first_collapse",true);
  if(o.collapse.max_story_drift_ratio<=0)throw std::invalid_argument("v1 IDA requires a documented positive drift_limit");
  const auto r=run_ida_suite(model,records,io);out["runs"]=Json::array();out["records"]=Json::array();
  for(const auto& a:r.runs)out["runs"].push_back({{"record",a.record_name},{"scale",a.scale_factor},{"pga",a.pga},{"termination",termination_name(a.termination)},{"reason",a.termination_reason},{"peak_story_drift_ratio",a.max_story_drift_ratio},{"peak_roof_displacement",a.max_roof_abs},{"elapsed_seconds",a.elapsed_seconds}});
  for(const auto& a:r.records)out["records"].push_back({{"record",a.record_name},{"lower_scale",a.last_noncollapse_scale},{"upper_scale",a.first_collapse_scale},{"right_censored",a.right_censored},{"numerical_failure",a.numerical_failure},{"bracket_has_numerical_gap",a.bracket_has_numerical_gap},{"nonmonotonic_response",a.nonmonotonic_response}});
  out["intensity_measure"]="PGA in supplied length/time^2 units; scale factor is not Sa(T1)";
 }else throw std::invalid_argument("analysis type must be nrha or ida");
 if(model.wall_count())out["wall_scope"]={{"kinematics","two_node_in_plane_small_displacement"},{"sfi_panel_law","explicit_input; fixed_angle_rc is not FSAM"},{"wall_geometric_stiffness",false},{"wall_acceptance_criteria","not_assessed"},{"solver","coupled wall tangent uses direct sparse factorization"}};
 if(model.steel_member_count()||!panel_zone_ids.empty()||!brb_ids.empty()||!damper_ids.empty())out["steel_scope"]={{"steel_members","two-node condensed concentrated-plasticity beams/columns"},{"panel_zones","zero-length rotational joint-shear spring; user supplies backbone"},{"brbs","small-displacement axial bilinear element; fatigue and fracture not modeled"},{"viscous_dampers","memoryless axial power-law dashpot with consistent velocity tangent"},{"member_geometric_stiffness","optional constant compression preload only"},{"acceptance_criteria","not_assessed"},{"solver","coupled member and velocity-dependent tangents use same-pattern sparse refactorization"}};
 if(!soil_spring_ids.empty()||!soil_dashpot_ids.empty()||!pile_lines.empty())out["soil_scope"]={{"soil_springs","directional translational or relative-rotational zero-length springs with explicit user-supplied backbone provenance"},{"py_tz_qz_adapters","bilinear tributary-force conversions; not OpenSees Simple1 material implementations"},{"pile_lines","vertical small-displacement elastic beam-column meshes with coincident fixed soil nodes"},{"radiation_damping","optional local translational or rotational power-law dashpots; no free-field interaction"},{"excluded","liquefaction, pore-pressure generation, pile-group shadowing, kinematic interaction, and SSI-compatible input motion"},{"acceptance_criteria","not_assessed"}};
 out["gravity_analysis_performed"]=gravity_performed;out["geometric_formulation"]=gravity_performed?"gravity_state_transfer_with_optional_member_pdelta":"small_displacement_with_optional_constant_preload_or_member_pdelta";out["analysis_wall_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();write_json(argv[2],out);
 bool failed=false;for(const auto& r:out["runs"])if(r["termination"]=="numerical_failure"||r["termination"]=="initial_instability")failed=true;
 std::cout<<argv[2]<<'\n';return failed?3:0;
}catch(const std::exception& e){std::cerr<<"QuakeCore job error: "<<e.what()<<'\n';return 2;}}
