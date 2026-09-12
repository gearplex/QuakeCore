#include "quake/wall_material.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace quake;
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double x,double y,double tol,const char* why){check(std::isfinite(x)&&std::abs(x-y)<=tol,why);}
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
    std::cout<<"wrapper state, rollback, and failure checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
