#pragma once

#include "json_helpers.hpp"
#include "quake/steel2d.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

inline void steel_keys(const Json& j,std::initializer_list<std::string> allowed){
    if(!j.is_object())throw std::invalid_argument("expected steel component object");
    for(auto it=j.begin();it!=j.end();++it)
        if(std::find(allowed.begin(),allowed.end(),it.key())==allowed.end())
            throw std::invalid_argument("unknown steel component field: "+it.key());
}

inline quake::NonlinearMaterial steel_nonlinear_material(const Json& j){
    steel_keys(j,{"type","k","fy","hardening_ratio","hardening_stiffness","a","b","f","c","io","ls","cp",
                  "up","upc","uu","fcap_fy","fres_fy","lambda_s","lambda_c","lambda_a","lambda_k",
                  "cs","cc","ca","ck","provenance"});
    const std::string type=j.at("type");
    if(type=="bilinear")
        return quake::NonlinearMaterial(quake::BilinearSpring(j.at("k"),j.at("fy"),j.at("hardening_ratio")));
    if(j.value("provenance",std::string()).empty())
        throw std::invalid_argument("degrading steel component material requires nonempty provenance");
    if(type=="asce41_parameterized"){
        quake::ASCE41HingeParams p;p.Ke=j.at("k");p.posFy=p.negFy=j.at("fy");
        p.pos_a=p.neg_a=j.at("a");p.pos_b=p.neg_b=j.at("b");p.pos_f=p.neg_f=j.value("f",0.0);
        p.pos_c=p.neg_c=j.value("c",0.0);p.pos_io=p.neg_io=j.value("io",0.0);
        p.pos_ls=p.neg_ls=j.value("ls",0.0);p.pos_cp=p.neg_cp=j.value("cp",0.0);
        p.hardening_stiffness=j.value("hardening_stiffness",-1.0);
        p.hardening_ratio=j.value("hardening_ratio",0.0);
        p.backbone_shape=quake::ASCE41BackboneShape::StraightCE;
        return quake::NonlinearMaterial(quake::ASCE41HingeMaterial(p));
    }
    if(type=="imk_peak_oriented"){
        quake::IMKPeakOrientedParams p;p.Ke=j.at("k");p.posFy=p.negFy=j.at("fy");
        p.posUp=p.negUp=j.at("up");p.posUpc=p.negUpc=j.at("upc");p.posUu=p.negUu=j.at("uu");
        p.posFcapFy=p.negFcapFy=j.at("fcap_fy");p.posFresFy=p.negFresFy=j.at("fres_fy");
        p.lambdaS=j.value("lambda_s",1e30);p.lambdaC=j.value("lambda_c",1e30);
        p.lambdaA=j.value("lambda_a",1e30);p.lambdaK=j.value("lambda_k",1e30);
        p.cS=j.value("cs",1.0);p.cC=j.value("cc",1.0);p.cA=j.value("ca",1.0);p.cK=j.value("ck",1.0);
        return quake::NonlinearMaterial(quake::IMKPeakOrientedMaterial(p));
    }
    throw std::invalid_argument("unsupported steel component material type "+type);
}

inline quake::SteelMember2DProperties steel_member_properties(const Json& j){
    steel_keys(j,{"id","type","i","j","role","E","A","I","constant_compression","hinge_i","hinge_j",
                  "local_max_iterations","local_relative_tolerance","provenance"});
    quake::SteelMember2DProperties p;p.E=j.at("E");p.A=j.at("A");p.I=j.at("I");
    p.axial_compression=j.value("constant_compression",0.0);
    p.hinge_i=steel_nonlinear_material(j.at("hinge_i"));p.hinge_j=steel_nonlinear_material(j.at("hinge_j"));
    p.local_max_iterations=j.value("local_max_iterations",50);
    p.local_relative_tolerance=j.value("local_relative_tolerance",1e-11);
    return p;
}

inline Json steel_member_response_json(const quake::SteelMember2DResponse& r){
    return {{"nodal_force",r.force},{"axial_deformation",r.axial_deformation},{"axial_force",r.axial_force},
            {"chord_rotation",r.chord_rotation},{"end_rotation",r.end_rotation},{"hinge_rotation",r.hinge_rotation},
            {"end_moment",r.end_moment},{"local_iterations",r.local_iterations}};
}

inline Json viscous_damper_response_json(const quake::ViscousDamper2DTrial& r){
    return {{"deformation_rate",r.deformation_rate},{"force",r.force},{"velocity_tangent",r.tangent}};
}
