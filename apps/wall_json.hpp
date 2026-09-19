#pragma once
#include "json_helpers.hpp"
#include "quake/wall2d.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
namespace quake {
inline void wall_keys(const Json& j,std::initializer_list<std::string> allowed){
    if(!j.is_object())throw std::invalid_argument("wall input must be an object");
    for(auto it=j.begin();it!=j.end();++it)if(std::find(allowed.begin(),allowed.end(),it.key())==allowed.end())throw std::invalid_argument("unknown wall field: "+it.key());
}
inline std::array<double,4> wall_four(const Json& j,const char* name){
    const auto& a=j.at(name);if(!a.is_array()||a.size()!=4)throw std::invalid_argument(std::string(name)+" must contain four numbers");
    std::array<double,4> out{};for(std::size_t i=0;i<4;++i){if(!a.at(i).is_number())throw std::invalid_argument(std::string(name)+" must contain four numbers");out[i]=a.at(i).get<double>();}return out;
}
inline std::array<Pinching4Point,4> wall_pinching_points(const Json& j,const char* name){
    const auto& a=j.at(name);if(!a.is_array()||a.size()!=4)throw std::invalid_argument(std::string(name)+" must contain four [deformation,force] points");
    std::array<Pinching4Point,4> out{};for(std::size_t i=0;i<4;++i){const auto& q=a.at(i);if(!q.is_array()||q.size()!=2||!q.at(0).is_number()||!q.at(1).is_number())throw std::invalid_argument(std::string(name)+" must contain four [deformation,force] points");out[i]={q.at(0).get<double>(),q.at(1).get<double>()};}return out;
}
inline WallUniaxial wall_uniaxial(const Json& j){
    std::string t=j.at("type");
    if(t=="elastic"){wall_keys(j,{"type","E"});return WallUniaxial::elastic(j.at("E"));}
    if(t=="steel_bilinear"){wall_keys(j,{"type","E","fy","b"});return WallUniaxial::steel(j.at("E"),j.at("fy"),j.at("b"));}
    if(t=="concrete01"){wall_keys(j,{"type","fc","epsc","fcu","epsu"});return WallUniaxial::concrete01(j.at("fc"),j.at("epsc"),j.at("fcu"),j.at("epsu"));}
    if(t=="concrete_cm"){
        wall_keys(j,{"type","fc","epsc","Ec","rc","xcrn","ft","et","rt","xcrp","gap_close"});
        ConcreteCMParameters p{j.at("fc"),j.at("epsc"),j.at("Ec"),j.at("rc"),j.at("xcrn"),j.at("ft"),j.at("et"),j.at("rt"),j.at("xcrp"),j.at("gap_close").get<bool>()};
        return WallUniaxial::concrete_cm(p);
    }
    if(t=="pinching4"){
        wall_keys(j,{"type","positive","negative","r_disp_positive","r_force_positive","u_force_positive","r_disp_negative","r_force_negative","u_force_negative","gamma_k","gamma_k_limit","gamma_d","gamma_d_limit","gamma_f","gamma_f_limit","gamma_e","damage_mode","admitted_reversal_count"});
        Pinching4CyclicParameters p{};p.envelope.positive=wall_pinching_points(j,"positive");p.envelope.negative=wall_pinching_points(j,"negative");
        p.r_disp_positive=j.at("r_disp_positive");p.r_force_positive=j.at("r_force_positive");p.u_force_positive=j.at("u_force_positive");
        p.r_disp_negative=j.at("r_disp_negative");p.r_force_negative=j.at("r_force_negative");p.u_force_negative=j.at("u_force_negative");
        p.gamma_k=wall_four(j,"gamma_k");p.gamma_k_limit=j.at("gamma_k_limit");p.gamma_d=wall_four(j,"gamma_d");p.gamma_d_limit=j.at("gamma_d_limit");p.gamma_f=wall_four(j,"gamma_f");p.gamma_f_limit=j.at("gamma_f_limit");p.gamma_e=j.at("gamma_e");
        const std::string damage=j.at("damage_mode");if(damage=="energy")p.damage_mode=Pinching4DamageMode::Energy;else if(damage=="cycle")p.damage_mode=Pinching4DamageMode::Cycle;else throw std::invalid_argument("pinching4 damage_mode must be energy or cycle");
        const auto& count=j.at("admitted_reversal_count");if(!count.is_number_integer()||count.get<long long>()<2)throw std::invalid_argument("pinching4 admitted_reversal_count must be an integer >=2");p.admitted_reversal_count=count.get<unsigned>();
        return WallUniaxial::pinching4(p);
    }
    if(t=="minmax"){
        wall_keys(j,{"type","min","max","material"});return WallUniaxial::minmax(wall_uniaxial(j.at("material")),j.at("min"),j.at("max"));
    }
    if(t=="parallel"){
        wall_keys(j,{"type","materials"});const auto& a=j.at("materials");if(!a.is_array()||a.empty())throw std::invalid_argument("parallel materials must be a nonempty array");std::vector<WallUniaxial> materials;for(const auto& child:a)materials.push_back(wall_uniaxial(child));return WallUniaxial::parallel(std::move(materials));
    }
    throw std::invalid_argument("unsupported wall uniaxial law: "+t);
}
inline WallPanel wall_panel(const Json& j){
    std::string t=j.at("type");const double rad=std::acos(-1.0)/180;
    if(t=="elastic_plane_stress"){wall_keys(j,{"type","E","nu"});return WallPanel(j.at("E"),j.at("nu"));}
    if(t=="fixed_angle_rc"){
        wall_keys(j,{"type","background_E","background_nu","angle_deg","concrete","steel_x","steel_y","rho_x","rho_y"});
        double rx=j.at("rho_x"),ry=j.at("rho_y");if(!std::isfinite(rx)||!std::isfinite(ry)||rx<=0||ry<=0||rx>=1||ry>=1)throw std::invalid_argument("fixed-angle RC requires rho_x,rho_y in (0,1)");
        double a=j.at("angle_deg").get<double>()*rad;auto c=wall_uniaxial(j.at("concrete"));
        return WallPanel(j.at("background_E"),j.at("background_nu"),{{a,1,c},{a+std::acos(-1.0)/2,1,c},{0,rx,wall_uniaxial(j.at("steel_x"))},{std::acos(-1.0)/2,ry,wall_uniaxial(j.at("steel_y"))}});
    }
    if(t=="layered_plane_stress"){
        wall_keys(j,{"type","background_E","background_nu","layers"});std::vector<WallPanelLayer> layers;
        if(!j.at("layers").is_array())throw std::invalid_argument("panel layers must be an array");
        for(const auto& l:j.at("layers")){wall_keys(l,{"angle_deg","weight","material"});layers.push_back({l.at("angle_deg").get<double>()*rad,l.at("weight"),wall_uniaxial(l.at("material"))});}
        return WallPanel(j.at("background_E"),j.at("background_nu"),std::move(layers));
    }
    throw std::invalid_argument("unsupported wall panel law (FSAM is not implemented): "+t);
}
inline MVLEMProperties mvlem_properties(const Json& j){
    wall_keys(j,{"id","type","i","j","c","density","fibers","shear"});MVLEMProperties p;p.c=j.value("c",.4);p.density=j.value("density",0.0);
    if(!j.at("fibers").is_array())throw std::invalid_argument("MVLEM fibers must be an array");
    for(const auto& f:j.at("fibers")){wall_keys(f,{"width","thickness","rho","concrete","steel"});p.fibers.push_back({f.at("width"),f.at("thickness"),f.at("rho"),wall_uniaxial(f.at("concrete")),wall_uniaxial(f.at("steel"))});}
    const auto& sh=j.at("shear");std::string t=sh.at("type");
    if(t=="bilinear"){wall_keys(sh,{"type","stiffness","yield_force","hardening_ratio"});p.shear=WallUniaxial::steel(sh.at("stiffness"),sh.at("yield_force"),sh.at("hardening_ratio"));}
    else if(t=="elastic"){wall_keys(sh,{"type","stiffness"});p.shear=WallUniaxial::elastic(sh.at("stiffness"));}
    else p.shear=wall_uniaxial(sh);return p;
}
inline SFIMVLEMProperties sfi_mvlem_properties(const Json& j){
    wall_keys(j,{"id","type","i","j","c","density","panels","local_max_iterations","local_relative_tolerance"});SFIMVLEMProperties p;p.c=j.value("c",.4);p.density=j.value("density",0.0);
    if(j.contains("local_max_iterations")){const auto& n=j.at("local_max_iterations");if(!n.is_number_integer()||n.get<double>()<1||n.get<double>()>1000)throw std::invalid_argument("invalid local_max_iterations");p.local_max_iterations=n.get<int>();}
    p.local_relative_tolerance=j.value("local_relative_tolerance",1e-10);
    if(!j.at("panels").is_array())throw std::invalid_argument("SFI panels must be an array");
    for(const auto& f:j.at("panels")){wall_keys(f,{"width","thickness","material"});p.panels.push_back({f.at("width"),f.at("thickness"),wall_panel(f.at("material"))});}return p;
}
inline Json wall_response_json(const Wall2DResponse& r){
    return {{"nodal_force",r.force},{"curvature",r.curvature},{"shear_deformation",r.shear_deformation},{"fiber_strain",r.fiber_strain},{"concrete_stress",r.concrete_stress},{"steel_stress",r.steel_stress},{"panel_strain",r.panel_strain},{"panel_stress",r.panel_stress},{"transverse_stress_residual",r.maximum_transverse_stress_residual},{"local_iterations",r.local_iterations}};
}
} // namespace quake
