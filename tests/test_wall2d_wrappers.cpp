#include "quake/wall2d.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace quake;
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double a,double b,double tol,const char* why){if(!std::isfinite(a)||std::abs(a-b)>tol)throw std::runtime_error(why);}

ConcreteCMParameters unconfined(){return {-6.5,-.002,4595.486916530173,7.0,1.030,.0604669,2.0*.0604669/4595.486916530173,1.2,10000.0,true};}
ConcreteCMParameters confined(){return {-8.01435,-.00432976,5102.805419570689,7.0,1.015,.0671422,2.0*.0671422/5102.805419570689,1.2,10000.0,true};}
Pinching4CyclicParameters steel(){Pinching4CyclicParameters p{};p.envelope.positive={{{.00204828,59.4},{.00705517,66.66},{.0523366,75.9},{.0654208,33.0}}};p.envelope.negative={{{-.00204828,-59.4},{-.00705517,-66.66},{-.0523366,-75.9},{-.0654208,-33.0}}};p.r_disp_positive=p.r_disp_negative=.6;p.r_force_positive=p.r_force_negative=.99;p.u_force_positive=p.u_force_negative=.4;p.gamma_d={.1,0,0,0};p.gamma_d_limit=2;p.gamma_e=10000;p.damage_mode=Pinching4DamageMode::Energy;p.admitted_reversal_count=10;return p;}
Pinching4CyclicParameters shear(){Pinching4CyclicParameters p{};p.envelope.positive={{{.00010787,171.319},{.288,285.531},{.72,299.808},{1.44,57.1062}}};p.envelope.negative={{{-.00010787,-171.319},{-.288,-285.531},{-.72,-299.808},{-1.44,-57.1062}}};p.r_disp_positive=p.r_disp_negative=.6;p.r_force_positive=p.r_force_negative=.99;p.u_force_positive=p.u_force_negative=.4;p.gamma_d={.1,0,0,0};p.gamma_d_limit=2;p.gamma_e=10000;p.damage_mode=Pinching4DamageMode::Energy;p.admitted_reversal_count=6;return p;}
WallUniaxial received_steel(){return WallUniaxial::parallel({WallUniaxial::minmax(WallUniaxial::pinching4(steel()),-.075,.075),WallUniaxial::elastic(.01)});}
MVLEMProperties story1(){
    MVLEMProperties p;p.c=.4;p.density=0;p.shear=WallUniaxial::pinching4(shear());
    auto cs=WallUniaxial::concrete_cm(confined()),us=WallUniaxial::concrete_cm(unconfined());auto ss=received_steel();
    for(int i=0;i<3;++i)p.fibers.push_back({4.08,12.0,.0333474,cs,ss});
    for(int i=0;i<3;++i)p.fibers.push_back({15.84,12.0,.0025641,us,ss});
    for(int i=0;i<3;++i)p.fibers.push_back({4.08,12.0,.0333474,cs,ss});
    return p;
}
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

    // Story-1 Gate 4 MVLEM admission. Geometry and reinforcement ratios are
    // preserved from the checked-in archetype 10001 matched-law job. The
    // +/-0.075 MinMax limits remain protocol-scoped placeholders because the
    // exact received m[30]/m[31] values are not retained in the repository.
    Wall2D yori(72.0,story1());auto ys=yori.initial_state();
    check(yori.state_size()==9*(28+24)+17,"story1 received-material MVLEM state size");
    // Pure axial virgin trial at eps=-0.0004. Expected resultant aggregates
    // frozen OpenSees ConcreteCM stresses with the symmetric virgin steel law.
    auto axial=yori.trial({0,0,0,0,-.0288,0},ys.data());
    near(axial.force[4],-1632.2038326115676,1e-9,"story1 MVLEM axial oracle resultant");
    near(axial.force[1],1632.2038326115676,1e-9,"story1 MVLEM axial equilibrium");
    near(axial.tangent[28],52795.2153542112,1e-8,"story1 MVLEM axial tangent");
    // Pure shear virgin trial at 0.0001 in shear displacement matches the
    // frozen OpenSees shear Pinching4 point exactly at protocol step 20.
    auto shear_trial=yori.trial({0,0,0,.0001,0,0},ys.data());
    near(shear_trial.force[3],158.8198757763975,1e-10,"story1 MVLEM shear oracle resultant");
    near(shear_trial.force[0],-158.8198757763975,1e-10,"story1 MVLEM shear equilibrium");
    check(ys==yori.initial_state(),"story1 MVLEM trial mutated committed state");
    std::cout<<"wall dynamic wrapper and story1 Gate4 MVLEM material-stack checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
