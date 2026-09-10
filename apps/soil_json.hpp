#pragma once

#include "json_helpers.hpp"
#include "quake/frame2d.hpp"
#include "quake/soil2d.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

struct PileLineInfo {
    std::string id;
    std::vector<double> depths;
    std::vector<int> node_ids;
    std::vector<int> member_ids;
    std::vector<int> lateral_spring_ids;
    std::vector<int> shaft_spring_ids;
    int toe_spring_id{-1};
};

inline void soil_keys(const Json& j,std::initializer_list<std::string> allowed){
    if(!j.is_object())throw std::invalid_argument("expected soil JSON object");
    for(auto it=j.begin();it!=j.end();++it)
        if(std::find(allowed.begin(),allowed.end(),it.key())==allowed.end())
            throw std::invalid_argument("unknown soil field: "+it.key());
}
inline int soil_id(const Json& j){
    if(!j.is_number_integer()||j.get<double>()<std::numeric_limits<int>::min()||j.get<double>()>std::numeric_limits<int>::max())
        throw std::invalid_argument("soil identifiers must be 32-bit integers");
    return j.get<int>();
}
inline std::pair<double,double> soil_direction(const Json& j){
    if(!j.is_array()||j.size()!=2)throw std::invalid_argument("soil direction must be [x,y]");
    const double x=j[0],y=j[1],L=std::hypot(x,y);
    if(!std::isfinite(L)||L<=0.0)throw std::invalid_argument("soil direction must be finite and nonzero");
    return {x/L,y/L};
}
inline void require_soil_provenance(const Json& j){
    if(j.value("provenance",std::string()).empty())throw std::invalid_argument("soil component requires parameter provenance");
}
inline quake::NonlinearMaterial soil_material(const Json& j){
    using namespace quake;
    const std::string type=j.at("type");
    if(type=="bilinear"){
        soil_keys(j,{"id","i","j","direction","dof","type","k","fy","hardening_ratio","provenance"});
        return NonlinearMaterial(BilinearSpring(j.at("k"),j.at("fy"),j.value("hardening_ratio",0.0)));
    }
    if(type=="asymmetric_elastic_perfectly_plastic"){
        soil_keys(j,{"id","i","j","direction","dof","type","k","positive_capacity","negative_capacity","provenance"});
        return NonlinearMaterial(AsymmetricElasticPerfectlyPlasticSpring(j.at("k"),j.at("positive_capacity"),j.at("negative_capacity")));
    }
    if(type=="py_bilinear"){
        soil_keys(j,{"id","i","j","direction","type","ultimate_resistance_per_length","tributary_length","displacement_50","hardening_ratio","provenance"});
        if(j.value("hardening_ratio",0.0)!=0.0)throw std::invalid_argument("p-y bilinear adapter requires zero hardening to preserve the ultimate force cap");
        auto x=py_bilinear_backbone(j.at("ultimate_resistance_per_length"),j.at("tributary_length"),j.at("displacement_50"));
        return NonlinearMaterial(BilinearSpring(x.initial_stiffness,x.nodal_capacity,0.0));
    }
    if(type=="tz_bilinear"){
        soil_keys(j,{"id","i","j","direction","type","ultimate_interface_stress","pile_perimeter","tributary_length","displacement_50","hardening_ratio","provenance"});
        if(j.value("hardening_ratio",0.0)!=0.0)throw std::invalid_argument("t-z bilinear adapter requires zero hardening to preserve the ultimate force cap");
        auto x=tz_bilinear_backbone(j.at("ultimate_interface_stress"),j.at("pile_perimeter"),j.at("tributary_length"),j.at("displacement_50"));
        return NonlinearMaterial(BilinearSpring(x.initial_stiffness,x.nodal_capacity,0.0));
    }
    if(type=="qz_bilinear"){
        soil_keys(j,{"id","i","j","direction","type","ultimate_bearing_stress","toe_area","displacement_50","suction_ratio","compression_sign","provenance"});
        auto x=qz_bilinear_backbone(j.at("ultimate_bearing_stress"),j.at("toe_area"),j.at("displacement_50"));
        const double suction=j.value("suction_ratio",0.0);if(!std::isfinite(suction)||suction<0.0||suction>1.0)
            throw std::invalid_argument("q-z suction_ratio must be between zero and one");
        const int sign=j.value("compression_sign",-1);if(sign!=1&&sign!=-1)throw std::invalid_argument("q-z compression_sign must be -1 or 1");
        const double pos=sign>0?x.nodal_capacity:suction*x.nodal_capacity;
        const double neg=sign<0?x.nodal_capacity:suction*x.nodal_capacity;
        return NonlinearMaterial(AsymmetricElasticPerfectlyPlasticSpring(x.initial_stiffness,pos,neg));
    }
    throw std::invalid_argument("unsupported soil spring type "+type);
}

inline double soil_profile_value(const Json& j,const char* key,std::size_t i,std::size_t n){
    const auto& v=j.at(key);if(v.is_number())return v.get<double>();
    if(!v.is_array()||v.size()!=n)throw std::invalid_argument(std::string(key)+" must be a scalar or one value per pile node");
    return v[i].get<double>();
}
inline std::vector<int> soil_id_array(const Json& j,const char* key,std::size_t n){
    const auto& a=j.at(key);if(!a.is_array()||a.size()!=n)throw std::invalid_argument(std::string(key)+" has wrong length");
    std::vector<int> out;out.reserve(n);for(const auto& x:a)out.push_back(soil_id(x));return out;
}
inline std::string pile_label(const Json& id){
    if(id.is_string()&&!id.get<std::string>().empty())return id.get<std::string>();
    if(id.is_number_integer())return std::to_string(soil_id(id));
    throw std::invalid_argument("pile id must be a nonempty string or integer");
}

inline PileLineInfo add_vertical_pile_line(const Json& p,quake::Frame2DBuilder& b,
    std::map<int,std::pair<double,double>>& coords,std::set<int>& elements,std::vector<int>& soil_spring_ids){
    using namespace quake;
    soil_keys(p,{"id","head_node","embedded_length","segments","pile_node_ids","soil_node_ids","member_ids","E","A","I","constant_compression","lateral","shaft","toe","provenance"});
    require_soil_provenance(p);PileLineInfo info;info.id=pile_label(p.at("id"));
    const int head=soil_id(p.at("head_node"));const auto hit=coords.find(head);if(hit==coords.end())throw std::invalid_argument("unknown pile head node");
    const int segments=p.at("segments");if(segments<1)throw std::invalid_argument("pile segments must be positive");
    const double length=p.at("embedded_length");if(!std::isfinite(length)||length<=0.0)throw std::invalid_argument("pile embedded_length must be positive");
    const std::size_t nn=static_cast<std::size_t>(segments+1);info.node_ids=soil_id_array(p,"pile_node_ids",nn);
    if(info.node_ids.front()!=head)throw std::invalid_argument("first pile_node_ids entry must equal head_node");
    auto anchors=soil_id_array(p,"soil_node_ids",nn);info.member_ids=soil_id_array(p,"member_ids",segments);
    info.depths.resize(nn);const double dz=length/segments;
    for(std::size_t i=0;i<nn;++i)info.depths[i]=i*dz;
    const auto trib=nodal_tributary_lengths(info.depths);
    for(std::size_t i=1;i<nn;++i){int id=info.node_ids[i];if(coords.contains(id))throw std::invalid_argument("duplicate generated pile node id");coords[id]={hit->second.first,hit->second.second-info.depths[i]};b.add_node(id,coords[id].first,coords[id].second);}
    for(std::size_t i=0;i<nn;++i){int id=anchors[i];if(coords.contains(id))throw std::invalid_argument("duplicate generated soil node id");coords[id]={hit->second.first,hit->second.second-info.depths[i]};b.add_node(id,coords[id].first,coords[id].second);b.fix(id);}
    const double E=p.at("E"),A=p.at("A"),I=p.at("I"),P=p.value("constant_compression",0.0);
    for(int i=0;i<segments;++i){int id=info.member_ids[static_cast<std::size_t>(i)];if(!elements.insert(id).second)throw std::invalid_argument("duplicate pile member/element id");b.add_elastic_frame(id,info.node_ids[static_cast<std::size_t>(i)],info.node_ids[static_cast<std::size_t>(i+1)],E,A,I,P);}
    if(p.contains("lateral")){
        const auto& x=p.at("lateral");soil_keys(x,{"spring_ids","ultimate_resistance_per_length","displacement_50","hardening_ratio","provenance"});require_soil_provenance(x);
        if(x.value("hardening_ratio",0.0)!=0.0)throw std::invalid_argument("pile p-y bilinear adapter requires zero hardening to preserve the ultimate force cap");info.lateral_spring_ids=soil_id_array(x,"spring_ids",nn);
        for(std::size_t i=0;i<nn;++i){const int id=info.lateral_spring_ids[i];if(!elements.insert(id).second)throw std::invalid_argument("duplicate pile soil spring/element id");
            auto bb=py_bilinear_backbone(soil_profile_value(x,"ultimate_resistance_per_length",i,nn),trib[i],soil_profile_value(x,"displacement_50",i,nn));
            b.add_soil_spring(id,anchors[i],info.node_ids[i],1,0,NonlinearMaterial(BilinearSpring(bb.initial_stiffness,bb.nodal_capacity,0.0)));soil_spring_ids.push_back(id);}
    }
    if(p.contains("shaft")){
        const auto& x=p.at("shaft");soil_keys(x,{"spring_ids","ultimate_interface_stress","pile_perimeter","displacement_50","hardening_ratio","provenance"});require_soil_provenance(x);
        if(x.value("hardening_ratio",0.0)!=0.0)throw std::invalid_argument("pile t-z bilinear adapter requires zero hardening to preserve the ultimate force cap");info.shaft_spring_ids=soil_id_array(x,"spring_ids",nn);
        const double perimeter=x.at("pile_perimeter");
        for(std::size_t i=0;i<nn;++i){const int id=info.shaft_spring_ids[i];if(!elements.insert(id).second)throw std::invalid_argument("duplicate pile soil spring/element id");
            auto bb=tz_bilinear_backbone(soil_profile_value(x,"ultimate_interface_stress",i,nn),perimeter,trib[i],soil_profile_value(x,"displacement_50",i,nn));
            b.add_soil_spring(id,anchors[i],info.node_ids[i],0,1,NonlinearMaterial(BilinearSpring(bb.initial_stiffness,bb.nodal_capacity,0.0)));soil_spring_ids.push_back(id);}
    }
    if(p.contains("toe")){
        const auto& x=p.at("toe");soil_keys(x,{"spring_id","ultimate_bearing_stress","toe_area","displacement_50","suction_ratio","compression_sign","provenance"});require_soil_provenance(x);
        info.toe_spring_id=soil_id(x.at("spring_id"));if(!elements.insert(info.toe_spring_id).second)throw std::invalid_argument("duplicate pile toe spring/element id");
        auto bb=qz_bilinear_backbone(x.at("ultimate_bearing_stress"),x.at("toe_area"),x.at("displacement_50"));const double suction=x.value("suction_ratio",0.0);
        if(!std::isfinite(suction)||suction<0.0||suction>1.0)throw std::invalid_argument("pile toe suction_ratio must be between zero and one");
        const int sign=x.value("compression_sign",-1);if(sign!=1&&sign!=-1)throw std::invalid_argument("pile toe compression_sign must be -1 or 1");
        const double pos=sign>0?bb.nodal_capacity:suction*bb.nodal_capacity,neg=sign<0?bb.nodal_capacity:suction*bb.nodal_capacity;
        b.add_soil_spring(info.toe_spring_id,anchors.back(),info.node_ids.back(),0,1,NonlinearMaterial(AsymmetricElasticPerfectlyPlasticSpring(bb.initial_stiffness,pos,neg)));soil_spring_ids.push_back(info.toe_spring_id);
    }
    return info;
}

inline double pile_dof_value(const quake::CompiledFrame2D& model,int node,quake::Dof2D dof,
                             const std::vector<double>& x){
    const int r=model.reduced_dof(node,dof);return r<0?0.0:x.at(static_cast<std::size_t>(r));
}
inline Json pile_line_response_json(const PileLineInfo& p,const quake::CompiledFrame2D& model,
                                    const std::vector<double>& u,const std::vector<double>& state){
    Json j={{"depth",p.depths},{"node_ids",p.node_ids},{"member_ids",p.member_ids}};
    j["lateral_displacement"]=Json::array();j["vertical_displacement"]=Json::array();
    for(int id:p.node_ids){j["lateral_displacement"].push_back(pile_dof_value(model,id,quake::Dof2D::UX,u));j["vertical_displacement"].push_back(pile_dof_value(model,id,quake::Dof2D::UY,u));}
    auto springs=[&](const std::vector<int>& ids,const char* deformation,const char* force){j[deformation]=Json::array();j[force]=Json::array();for(int id:ids){auto x=model.nonlinear_component_snapshot(id,u,state);j[deformation].push_back(x.deformation);j[force].push_back(x.force);}};
    springs(p.lateral_spring_ids,"lateral_spring_deformation","lateral_soil_force");
    springs(p.shaft_spring_ids,"shaft_spring_deformation","shaft_soil_force");
    j["member_moment_i"]=Json::array();j["member_moment_j"]=Json::array();j["member_shear_i"]=Json::array();j["member_shear_j"]=Json::array();
    for(int id:p.member_ids){auto x=model.elastic_element_response(id,u);j["member_moment_i"].push_back(x.moment_i);j["member_moment_j"].push_back(x.moment_j);j["member_shear_i"].push_back(x.shear_i);j["member_shear_j"].push_back(x.shear_j);}
    if(p.toe_spring_id>=0){auto x=model.nonlinear_component_snapshot(p.toe_spring_id,u,state);j["toe_spring_deformation"]=x.deformation;j["toe_soil_force"]=x.force;}
    return j;
}
inline void accumulate_pile_peak_json(Json& peak,const Json& current){
    if(peak.is_null())peak={{"depth",current.at("depth")},{"node_ids",current.at("node_ids")},{"member_ids",current.at("member_ids")}};
    for(const auto* key:{"lateral_displacement","vertical_displacement","lateral_spring_deformation","lateral_soil_force","shaft_spring_deformation","shaft_soil_force","member_moment_i","member_moment_j","member_shear_i","member_shear_j"}){
        if(!current.contains(key))continue;const std::string out=std::string("peak_abs_")+key;
        if(!peak.contains(out))peak[out]=std::vector<double>(current.at(key).size(),0.0);
        for(std::size_t i=0;i<current.at(key).size();++i)peak[out][i]=std::max(peak[out][i].get<double>(),std::abs(current.at(key)[i].get<double>()));
    }
    for(const auto* key:{"toe_spring_deformation","toe_soil_force"})if(current.contains(key)){
        const std::string out=std::string("peak_abs_")+key;peak[out]=std::max(peak.value(out,0.0),std::abs(current.at(key).get<double>()));
    }
}
