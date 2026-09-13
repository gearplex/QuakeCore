#include "quake/pinching4_cycle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace quake;

namespace {
void check(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void near(double x,double y,double tol,const char* why){
    check(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tol*std::max({1.0,std::abs(x),std::abs(y)}),why);
}
template<class F>void rejects(F f){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,"expected Pinching4 rejection did not occur");}

Pinching4EnvelopeParameters envelope(){
    std::array<Pinching4Point,4> p{{
        {0.00204828,59.4},{0.00705517,66.66},{0.0523366,75.9},{0.0654208,33.0}}};
    std::array<Pinching4Point,4> n{{
        {-0.00204828,-59.4},{-0.00705517,-66.66},{-0.0523366,-75.9},{-0.0654208,-33.0}}};
    return {p,n};
}

Pinching4 material(){
    Pinching4CyclicParameters p{};
    p.envelope=envelope();
    p.r_disp_positive=p.r_disp_negative=0.6;
    p.r_force_positive=p.r_force_negative=0.99;
    p.u_force_positive=p.u_force_negative=0.4;
    p.gamma_d={0.1,0.0,0.0,0.0};
    p.gamma_d_limit=2.0;
    p.gamma_e=10000.0;
    p.damage_mode=Pinching4DamageMode::Energy;
    p.admitted_reversal_count=3;
    return Pinching4(p);
}
}

int main(){try{
    auto m=material();
    auto state=m.initial_state();
    std::vector<double> protocol{0.0};
    for(int i=1;i<=20;++i)protocol.push_back(0.001*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.001+(0.004-0.001)*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.004*(1.0-i/20.0));
    for(int i=1;i<=20;++i)protocol.push_back(-0.004*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(-0.004*(1.0-i/20.0));
    for(int i=1;i<=20;++i)protocol.push_back(0.012*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.012*(1.0-i/10.0));

    struct Point{int step;double force;double tangent;double work;Pinching4StateKind kind;};
    const Point points[]={
        {120,67.66902796576875,204.05715985559664,1.1074272065762027,Pinching4StateKind::PositiveEnvelope},
        {121,32.86909826865704,28999.941414259767,1.0471043308355472,Pinching4StateKind::PositiveToNegative},
        {122,-1.9308314284546704,28999.941414259767,1.0285413707314257,Pinching4StateKind::PositiveToNegative},
        {123,-30.980860695107893,2826.180958549829,1.0482883860055632,Pinching4StateKind::PositiveToNegative},
        {132,-61.50361504744605,2826.180958549829,1.5477045550153545,Pinching4StateKind::PositiveToNegative},
        {133,-62.52449846910956,356.87499126204204,1.6221214231252878,Pinching4StateKind::PositiveToNegative},
        {134,-63.38999922107336,1450.0018973854026,1.6976701217393977,Pinching4StateKind::NegativeEnvelope},
        {136,-66.68955359846188,204.05715985559664,1.8538738552286418,Pinching4StateKind::NegativeEnvelope},
        {140,-67.66902796576875,204.05715985559664,2.1763344509827953,Pinching4StateKind::NegativeEnvelope}};

    double work=0.0,previous_force=0.0,previous_deformation=0.0;
    for(int step=0;step<=140;++step){
        const double deformation=protocol[static_cast<std::size_t>(step)];
        if(step==121){
            const auto a=m.trial(deformation,state),b=m.trial(deformation,state);
            near(a.response.force,b.response.force,1e-13,"Pinching4 step140 trial purity");
            check(state.kind==Pinching4StateKind::PositiveEnvelope,"Pinching4 step140 trial mutated committed state");
        }
        const auto trial=m.trial(deformation,state);
        if(step>0)work+=0.5*(previous_force+trial.response.force)*(deformation-previous_deformation);
        for(const auto& point:points)if(point.step==step){
            near(trial.response.force,point.force,1e-11,"Pinching4 step140 force oracle");
            near(trial.response.tangent,point.tangent,1e-11,"Pinching4 step140 tangent oracle");
            near(work,point.work,1e-10,"Pinching4 step140 work oracle");
            check(trial.state.kind==point.kind,"Pinching4 step140 state mismatch");
        }
        previous_force=trial.response.force;
        previous_deformation=deformation;
        state=trial.state;
    }
    check(state.reversal_count==3,"Pinching4 step140 reversal history mismatch");
    rejects([&]{(void)m.trial(-0.00965,state);});
    std::cout<<"Pinching4 frozen steel_raw parity through step 140 passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
