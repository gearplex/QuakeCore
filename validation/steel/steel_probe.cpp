#include "json_helpers.hpp"
#include "quake/steel2d.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using namespace quake;

static std::vector<double> protocol(){
    std::vector<double> out;double previous=0.0;
    for(double target:{.006,-.012,.02,-.03,.04,-.05,.06,-.04,.02,0.0}){
        for(int i=1;i<=74;++i)out.push_back(previous+(target-previous)*i/74.0);
        previous=target;
    }
    return out;
}

int main(int argc,char** argv){try{
    if(argc!=2){std::cerr<<"Usage: steel_probe output.json\n";return 2;}
    const auto q=protocol();Json out={{"schema","quakecore.phase9k.steel-probe.v1"},{"panel_zone",Json::array()},{"brb",Json::array()},{"steel_member",Json::array()},{"viscous_damper",Json::array()}};
    NonlinearMaterial panel(BilinearSpring(5000.0,50.0,.02)),brb(BilinearSpring(1000.0,20.0,.01));
    auto ps=std::vector<double>(static_cast<std::size_t>(panel.state_size()),0.0),bs=std::vector<double>(static_cast<std::size_t>(brb.state_size()),0.0);panel.initialize_state(ps.data());brb.initialize_state(bs.data());
    for(double x:q){
        std::vector<double> pt(ps.size()),bt(bs.size());auto pr=panel.trial(x,ps.data(),pt.data());auto br=brb.trial(x,bs.data(),bt.data());
        out["panel_zone"].push_back({{"deformation",x},{"force",pr.force},{"tangent",pr.tangent}});
        out["brb"].push_back({{"deformation",x},{"force",br.force},{"tangent",br.tangent}});
        ps=std::move(pt);bs=std::move(bt);
    }
    SteelMember2DProperties p;p.E=2e8;p.A=.02;p.I=8e-4;
    p.hinge_i=NonlinearMaterial(BilinearSpring(5e6,300,.02));p.hinge_j=NonlinearMaterial(BilinearSpring(5e6,300,.02));
    SteelMember2D member(0,0,6,0,p);auto state=member.initial_state();
    for(double x:q){
        auto r=member.trial({0,0,x,0,0,0},state.data());
        out["steel_member"].push_back({{"rotation",x},{"end_moment",r.end_moment},{"hinge_rotation",r.hinge_rotation},{"nodal_force",r.force}});
        state=std::move(r.state);
    }
    ViscousDamper2D damper({125.0,.5,1e-12});
    for(double v:{-1.0,-.5,-.1,-.01,-.001,.001,.01,.1,.5,1.0}){
        auto r=damper.trial(v);out["viscous_damper"].push_back({{"velocity",v},{"force",r.force},{"tangent",r.tangent}});
    }
    write_json(argv[1],out);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
