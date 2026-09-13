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

Pinching4 material(){
    Pinching4CyclicParameters p{};
    p.envelope.positive={{{0.00010787,171.319},{0.288,285.531},{0.72,299.808},{1.2,138.0068}}};
    p.envelope.negative={{{-0.00010787,-171.319},{-0.288,-285.531},{-0.72,-299.808},{-1.2,-138.0068}}};
    p.r_disp_positive=p.r_disp_negative=0.6;
    p.r_force_positive=p.r_force_negative=0.99;
    p.u_force_positive=p.u_force_negative=0.4;
    p.gamma_d={0.1,0.0,0.0,0.0};
    p.gamma_d_limit=2.0;
    p.gamma_e=10000.0;
    p.damage_mode=Pinching4DamageMode::Energy;
    p.admitted_reversal_count=6;
    return Pinching4(p);
}
}

int main(){try{
    auto m=material();
    auto state=m.initial_state();
    std::vector<double> protocol{0.0};
    auto segment=[&](double a,double b){for(int i=1;i<=20;++i)protocol.push_back(a+(b-a)*i/20.0);};
    segment(0.0,0.0001);segment(0.0001,0.010);segment(0.010,0.10);segment(0.10,0.30);
    segment(0.30,0.0);segment(0.0,-0.30);segment(-0.30,0.72);segment(0.72,-0.72);
    segment(-0.72,1.20);segment(1.20,-1.20);segment(-1.20,0.0);

    struct Point{int step;double force;double tangent;double work;Pinching4StateKind kind;};
    const Point points[]={
        {20,158.8198757763975,1588198.7577639748,0.007940993788819876,Pinching4StateKind::PositiveEnvelope},
        {21,171.51225325621093,396.7180346333192,0.08969819572439047,Pinching4StateKind::PositiveEnvelope},
        {79,285.5970972222222,33.04861111111108,66.33617503860908,Pinching4StateKind::PositiveEnvelope},
        {80,285.92758333333336,33.04861111111108,69.19379844138686,Pinching4StateKind::PositiveEnvelope},
        {81,-122.36672033808712,165.7246583540572,67.96709196892252,Pinching4StateKind::PositiveToNegative},
        {100,-169.5982479689934,165.7246583540572,109.5720999526815,Pinching4StateKind::PositiveToNegative},
        {101,-177.22697654510387,396.7180346333192,112.17328913653722,Pinching4StateKind::NegativeEnvelope},
        {120,-285.92758333333336,33.04861111111108,178.7619583853205,Pinching4StateKind::NegativeEnvelope},
        {121,136.655717846209,329.74079061250774,174.95552581539883,Pinching4StateKind::NegativeToPositive},
        {130,284.3106867424242,21.736291035353766,272.321310683455,Pinching4StateKind::NegativeToPositive},
        {132,286.5277884280303,21.736291035353766,301.4340729171482,Pinching4StateKind::NegativeToPositive},
        {133,288.00964583333337,33.04861111111108,316.084777490813,Pinching4StateKind::PositiveEnvelope},
        {140,299.808,33.04861111111108,421.010227272063,Pinching4StateKind::PositiveEnvelope},
        {141,-132.75232177048485,178.83868732803202,414.99622285580045,Pinching4StateKind::PositiveToNegative},
        {153,-284.4411044886363,21.736291035353766,596.3436088805149,Pinching4StateKind::PositiveToNegative},
        {154,-286.0061174431818,21.736291035353766,616.8797088700604,Pinching4StateKind::PositiveToNegative},
        {155,-287.9105,33.04861111111108,637.5407070980149,Pinching4StateKind::NegativeEnvelope},
        {160,-299.808,33.04861111111108,743.3300370980149,Pinching4StateKind::NegativeEnvelope},
        {161,132.1699586728815,127.92256470807718,735.2834111143133,Pinching4StateKind::NegativeToPositive},
        {173,273.24167149999994,8.697532196969682,972.1243132402508,Pinching4StateKind::NegativeToPositive},
        {175,274.91159768181814,8.697532196969682,1024.7470270817053,Pinching4StateKind::NegativeToPositive},
        {176,267.44775999999996,-337.0858333333333,1050.7802762504325,Pinching4StateKind::PositiveEnvelope},
        {180,138.0068,-337.0858333333333,1128.6275517704325,Pinching4StateKind::PositiveEnvelope},
        {181,-130.85924352066,91.25720053969904,1128.198698381672,Pinching4StateKind::PositiveToNegative},
        {194,-272.82418995454543,8.697532196969682,1443.3571027061248,Pinching4StateKind::PositiveToNegative},
        {196,-274.91159768181814,8.697532196969682,1509.0853972224884,Pinching4StateKind::PositiveToNegative},
        {197,-259.3576999999999,-337.0858333333333,1541.1415550833974,Pinching4StateKind::NegativeEnvelope},
        {200,-138.0068,-337.0858333333333,1612.6671650833973,Pinching4StateKind::NegativeEnvelope},
        {201,25.059881508471527,37.01917882584288,1609.2787575286516,Pinching4StateKind::NegativeToPositive},
        {220,67.26174536993241,37.01917882584288,1661.9020848493417,Pinching4StateKind::NegativeToPositive}};

    double work=0.0,previous_force=0.0,previous_deformation=0.0;
    for(int step=0;step<=220;++step){
        const double deformation=protocol[static_cast<std::size_t>(step)];
        if(step==81||step==121||step==141||step==161||step==181||step==201){
            const auto a=m.trial(deformation,state),b=m.trial(deformation,state);
            near(a.response.force,b.response.force,1e-13,"shear Pinching4 trial purity");
        }
        const auto trial=m.trial(deformation,state);
        if(step>0)work+=0.5*(previous_force+trial.response.force)*(deformation-previous_deformation);
        for(const auto& point:points)if(point.step==step){
            near(trial.response.force,point.force,1e-11,"shear Pinching4 force oracle");
            near(trial.response.tangent,point.tangent,1e-11,"shear Pinching4 tangent oracle");
            near(work,point.work,1e-10,"shear Pinching4 work oracle");
            check(trial.state.kind==point.kind,"shear Pinching4 state mismatch");
        }
        previous_force=trial.response.force;
        previous_deformation=deformation;
        state=trial.state;
    }
    check(state.reversal_count==6,"shear Pinching4 reversal history mismatch");
    check(state.kind==Pinching4StateKind::NegativeToPositive,"shear Pinching4 terminal state mismatch");
    std::cout<<"Pinching4 full frozen shear oracle parity through step 220 passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
