#pragma once
#include "json_helpers.hpp"
#include "quake/wall2d.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
namespace quake {
inline void wall_keys(const Json& j,std::initializer_list<std::string> allowed){
    if(!j.is_object())throw std::invalid_argument("wall input must be an object");
    for(auto it=j.begin();it!=j.end();++it)if(std::find(allowed.begin(),allowed.end(),it.key())==allowed.end())throw std::invalid_argument("unknown wall field: "+it.key());
}
inline WallUniaxial wall_uniaxial(const Json& j){
    std::string t=j.at("type");
    if(t=="elastic"){wall_keys(j,{"type","E"});return WallUniaxial::elastic(j.at("E"));}
    if(t=="steel_bilinear"){wall_keys(j,{"type","E","fy","b"});return WallUniaxial::steel(j.at("E"),j.at("fy"),j.at("b"));}
    if(t=="concrete01"){wall_keys(j,{"type","fc","epsc","fcu","epsu"});return WallUniaxial::concrete01(j.at("fc"),j.at("epsc"),j.at("fcu"),j.at("epsu"));}
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
    if(t=="elastic"){wall_keys(sh,{"type","stiffness"});p.shear=WallUniaxial::elastic(sh.at("stiffness"));}
    else if(t=="bilinear"){wall_keys(sh,{"type","stiffness","yield_force","hardening_ratio"});p.shear=WallUniaxial::steel(sh.at("stiffness"),sh.at("yield_force"),sh.at("hardening_ratio"));}
    else throw std::invalid_argument("unsupported MVLEM shear law");return p;
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
