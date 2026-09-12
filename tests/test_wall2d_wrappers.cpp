#include "quake/wall2d.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace quake;
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double a,double b,double tol,const char* why){if(!std::isfinite(a)||std::abs(a-b)>tol)throw std::runtime_error(why);}
int main(){try{
    MVLEMProperties p;
    p.shear=WallUniaxial::parallel({WallUniaxial::minmax(WallUniaxial::elastic(10),-.02,.02),WallUniaxial::elastic(2)});
    for(int i=0;i<2;++i)p.fibers.push_back({.5,.2,.1,WallUniaxial::elastic(100),WallUniaxial::parallel({WallUniaxial::elastic(200),WallUniaxial::elastic(20)})});
    Wall2D w(3,p);
    const int fiber_state=6+12;
    const int shear_state=7+6;
    check(w.state_size()==2*fiber_state+shear_state,"dynamic wall state size");
    auto s=w.initial_state();check((int)s.size()==w.state_size(),"initial state size");
    auto r=w.trial({0,0,0,.003,-.001,0},s.data());check((int)r.state.size()==w.state_size(),"trial state size");
    auto old=s;auto fail=w.trial({0,0,0,.07,-.001,0},s.data());(void)fail;check(s==old,"wall trial mutated committed state");
    auto a=w.trial({0,0,0,.003,-.001,0},s.data());auto b=w.trial({0,0,0,.003,-.001,0},old.data());near(a.force[3],b.force[3],1e-12,"wall rollback mismatch");
    std::cout<<"wall dynamic wrapper state checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
