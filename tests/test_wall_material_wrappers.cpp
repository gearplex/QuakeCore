#include "quake/wall_material.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace quake;
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double x,double y,double tol,const char* why){check(std::isfinite(x)&&std::abs(x-y)<=tol,why);}

Pinching4CyclicParameters yori_steel(){
    Pinching4CyclicParameters p{};
    p.envelope={{{{0.00204828,59.4},{0.00705517,66.66},{0.0523366,75.9},{0.0654208,33.0}}},
                {{{-0.00204828,-59.4},{-0.00705517,-66.66},{-0.0523366,-75.9},{-0.0654208,-33.0}}}}};
    p.r_disp_positive=p.r_disp_negative=0.6;p.r_force_positive=p.r_force_negative=0.99;p.u_force_positive=p.u_force_negative=0.4;
    p.gamma_d={0.1,0.0,0.0,0.0};p.gamma_d_limit=2.0;p.gamma_e=10000.0;p.damage_mode=Pinching4DamageMode::Energy;p.admitted_reversal_count=10;
    return p;
}
std::vector<double> steel_protocol_to_240(){
    std::vector<double> v{0.0};
    auto segment=[&](double a,double b){for(int i=1;i<=20;++i)v.push_back(a+(b-a)*i/20.0);};
    segment(0,.001);segment(.001,.004);segment(.004,0);segment(0,-.004);segment(-.004,0);segment(0,.012);segment(.012,-.012);
    segment(-.012,.035);segment(.035,-.035);segment(-.035,.060);segment(.060,-.060);segment(-.060,.080);
    return v;
}
int main(){try{
    auto par=WallUniaxial::parallel({WallUniaxial::elastic(100),WallUniaxial::elastic(20)});
    check(par.state_size()==12,"parallel state size");near(par.initial_tangent(),120,1e-12,"parallel initial tangent");
    std::vector<double> s(par.state_size()),t(par.state_size());par.initialize(s.data());auto pr=par.trial(.01,s.data(),t.data());near(pr.stress,1.2,1e-12,"parallel stress");near(pr.tangent,120,1e-12,"parallel tangent");
    auto mm=WallUniaxial::minmax(WallUniaxial::elastic(100),-.02,.02);check(mm.state_size()==7,"minmax state size");std::vector<double> ms(mm.state_size()),mt(mm.state_size());mm.initialize(ms.data());
    auto edge=mm.trial(.02,ms.data(),mt.data());near(edge.stress,0,1e-12,"MinMax max limit fails");near(edge.tangent,1e-6,1e-15,"MinMax failed tangent");check(mt.back()==1,"MinMax limit failure flag");
    auto failed=mm.trial(.021,ms.data(),mt.data());near(failed.stress,0,1e-12,"MinMax fail stress");near(failed.tangent,1e-6,1e-15,"MinMax fail tangent");check(mt.back()==1,"MinMax trial failure flag");
    std::vector<double> recovered(mm.state_size());auto rr=mm.trial(.01,ms.data(),recovered.data());near(rr.stress,1,1e-12,"MinMax rollback stress");near(rr.tangent,100,1e-12,"MinMax rollback tangent");check(recovered.back()==0,"MinMax rollback flag");
    auto committed_failed=mt;std::vector<double> after(mm.state_size());auto ar=mm.trial(.01,committed_failed.data(),after.data());near(ar.stress,0,1e-12,"MinMax permanent fail stress");near(ar.tangent,1e-6,1e-15,"MinMax permanent fail tangent");check(after.back()==1,"MinMax permanent flag");
    auto nested=WallUniaxial::parallel({WallUniaxial::minmax(WallUniaxial::elastic(100),-.02,.02),WallUniaxial::elastic(18)});check(nested.state_size()==13,"nested state size");std::vector<double> ns(nested.state_size()),nt(nested.state_size());nested.initialize(ns.data());
    auto n1=nested.trial(.01,ns.data(),nt.data());near(n1.stress,1.18,1e-12,"nested intact stress");near(n1.tangent,118,1e-12,"nested intact tangent");
    auto nf=nested.trial(.021,ns.data(),nt.data());near(nf.stress,.378,1e-12,"nested fail stress");near(nf.tangent,18.000001,1e-12,"nested fail tangent");auto ncomm=nt;auto n2=nested.trial(.01,ncomm.data(),nt.data());near(n2.stress,.18,1e-12,"nested permanent branch stress");near(n2.tangent,18.000001,1e-12,"nested permanent branch tangent");
    WallPanel panel(0,.2,{{0,1,nested},{1.5707963267948966,1,WallUniaxial::elastic(5)}});check(panel.state_size()==19,"panel dynamic state size");auto ps=panel.initial_state();auto pv=panel.trial({.01,0,0},ps.data());near(pv.stress[0],1.18,1e-12,"panel wrapper stress");

    // Gate 4 integration seam: Pinching4 is now a wall material and therefore
    // composes through the already verified MinMax + Parallel lifecycle. The
    // exact received MinMax limits are not in the checked-in source input; a
    // protocol-scoped limit inside the oracle's (0.073, 0.080] failure bracket
    // is sufficient to verify the sampled step-239 -> step-240 transition.
    auto raw=WallUniaxial::pinching4(yori_steel());
    check(raw.state_size()==17,"Pinching4 wall state size");near(raw.initial_tangent(),28999.941414259767,1e-10,"Pinching4 wall initial tangent");
    auto wrapped=WallUniaxial::parallel({WallUniaxial::minmax(raw,-.075,.075),WallUniaxial::elastic(.01)});
    std::vector<double> ws(wrapped.state_size()),wt(wrapped.state_size());wrapped.initialize(ws.data());
    auto protocol=steel_protocol_to_240();
    for(std::size_t step=0;step<protocol.size();++step){
        auto r=wrapped.trial(protocol[step],ws.data(),wt.data());
        if(step==238){near(r.stress,33.00066002921643,1e-10,"wrapped oracle step238 stress");near(r.tangent,.01005044272310952,1e-12,"wrapped oracle step238 tangent");}
        if(step==239){near(r.stress,33.00073038231548,1e-10,"wrapped oracle step239 stress");near(r.tangent,.01005044272310952,1e-12,"wrapped oracle step239 tangent");}
        if(step==240){near(r.stress,.0008,1e-12,"wrapped oracle step240 stress");near(r.tangent,.010289999414142598,1e-12,"wrapped oracle step240 tangent");}
        ws=wt;
    }
    auto post=wrapped.trial(.072,ws.data(),wt.data());near(post.stress,.00072,1e-12,"wrapped permanent failure stress");near(post.tangent,.010289999414142598,1e-12,"wrapped permanent failure tangent");
    std::cout<<"wrapper state, Pinching4 bridge, rollback, and failure checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
