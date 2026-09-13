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

Pinching4EnvelopeParameters envelope(){
    std::array<Pinching4Point,4> p{{
        {0.00204828,59.4},{0.00705517,66.66},{0.0523366,75.9},{0.0654208,33.0}}};
    std::array<Pinching4Point,4> n{{
        {-0.00204828,-59.4},{-0.00705517,-66.66},{-0.0523366,-75.9},{-0.0654208,-33.0}}};
    return {p,n};
}

Pinching4 steel_material(){
    Pinching4CyclicParameters p{};
    p.envelope=envelope();
    p.r_disp_positive=p.r_disp_negative=0.6;
    p.r_force_positive=p.r_force_negative=0.99;
    p.u_force_positive=p.u_force_negative=0.4;
    p.gamma_d={0.1,0.0,0.0,0.0};
    p.gamma_d_limit=2.0;
    p.gamma_e=10000.0;
    p.damage_mode=Pinching4DamageMode::Energy;
    p.admitted_reversal_count=10;
    return Pinching4(p);
}

Pinching4 shear_material(){
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
    auto m=steel_material();
    auto state=m.initial_state();
    std::vector<double> protocol{0.0};
    for(int i=1;i<=20;++i)protocol.push_back(0.001*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.001+(0.004-0.001)*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.004*(1.0-i/20.0));
    for(int i=1;i<=20;++i)protocol.push_back(-0.004*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(-0.004*(1.0-i/20.0));
    for(int i=1;i<=20;++i)protocol.push_back(0.012*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.012*(1.0-i/10.0));
    for(int i=1;i<=20;++i)protocol.push_back(-0.012+(0.035+0.012)*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.035*(1.0-i/10.0));
    for(int i=1;i<=20;++i)protocol.push_back(-0.035+(0.06+0.035)*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.06*(1.0-i/10.0));
    for(int i=1;i<=20;++i)protocol.push_back(-0.06+(0.08+0.06)*i/20.0);
    for(int i=1;i<=20;++i)protocol.push_back(0.08*(1.0-i/10.0));
    for(int i=1;i<=20;++i)protocol.push_back(-0.08*(1.0-i/20.0));

    struct Point{int step;double force;double tangent;double work;Pinching4StateKind kind;};
    const Point points[]={
        {120,67.66902796576875,204.05715985559664,1.1074272065762027,Pinching4StateKind::PositiveEnvelope},
        {121,32.86909826865704,28999.941414259767,1.0471043308355472,Pinching4StateKind::PositiveToNegative},
        {134,-63.38999922107336,1450.0018973854026,1.6976701217393977,Pinching4StateKind::NegativeEnvelope},
        {140,-67.66902796576875,204.05715985559664,2.1763344509827953,Pinching4StateKind::NegativeEnvelope},
        {141,0.48083435774169914,28999.941414259767,2.0973883234933637,Pinching4StateKind::NegativeToPositive},
        {150,67.69523439064866,128.6248040863558,3.142247203416292,Pinching4StateKind::NegativeToPositive},
        {151,68.0465337115016,204.05715985559664,3.301743780936319,Pinching4StateKind::PositiveEnvelope},
        {160,72.36234264244747,204.05715985559664,4.78656764837933,Pinching4StateKind::PositiveEnvelope},
        {161,-29.137452307461775,28999.94141425976,4.710924090293105,Pinching4StateKind::PositiveToNegative},
        {173,-67.56660958656231,128.6248040863558,6.804603349626328,Pinching4StateKind::PositiveToNegative},
        {174,-68.07714228547994,204.05715985559664,7.041979915402402,Pinching4StateKind::NegativeEnvelope},
        {180,-72.36234264244747,204.05715985559664,8.51659450714564,Pinching4StateKind::NegativeEnvelope},
        {181,31.28951312606798,769.5644638960573,8.41904653704424,Pinching4StateKind::NegativeToPositive},
        {193,72.51897817158633,47.45230045580652,11.446446447012438,Pinching4StateKind::NegativeToPositive},
        {195,72.96977502591649,47.45230045580652,12.137518024700576,Pinching4StateKind::NegativeToPositive},
        {196,73.58668560158105,204.05715985559664,12.485589618690883,Pinching4StateKind::PositiveEnvelope},
        {200,50.77352226349338,-3278.763699729444,13.808982081319918,Pinching4StateKind::PositiveEnvelope},
        {201,-32.03430469329189,522.8465963895039,13.752764428609314,Pinching4StateKind::PositiveToNegative},
        {214,-72.38848434533287,47.45230045580652,17.84065597630393,Pinching4StateKind::PositiveToNegative},
        {216,-72.95791195080254,47.45230045580652,18.71273435408074,Pinching4StateKind::PositiveToNegative},
        {217,-73.79074276143665,204.05715985559664,19.15298031821746,Pinching4StateKind::NegativeEnvelope},
        {220,-50.773522263493376,-3278.763699729444,20.39944025438689,Pinching4StateKind::NegativeEnvelope},
        {221,14.158369125362219,199.90962098879086,20.27128721840343,Pinching4StateKind::NegativeToPositive},
        {235,32.737500028984016,12.500000011066907,22.615231350840624,Pinching4StateKind::NegativeToPositive},
        {237,32.912500029138954,12.500000011066907,23.074781351247488,Pinching4StateKind::NegativeToPositive},
        {238,33.00000002921642,5.0442723109520276e-05,23.30547510145173,Pinching4StateKind::PositiveEnvelope},
        {240,33.000000735414545,5.0442723109520276e-05,23.767475106804145,Pinching4StateKind::PositiveEnvelope},
        {241,-14.257075645284722,164.9903618605693,23.692503406443628,Pinching4StateKind::PositiveToNegative},
        {255,-32.675000028928686,12.500000011066907,26.323871443726045,Pinching4StateKind::PositiveToNegative},
        {258,-32.97500002919429,12.500000011066907,27.11167144442352,Pinching4StateKind::PositiveToNegative},
        {259,-33.00000033187276,5.0442723109520276e-05,27.375571445867788,Pinching4StateKind::NegativeEnvelope},
        {260,-33.000000735414545,5.0442723109520276e-05,27.639571450136938,Pinching4StateKind::NegativeEnvelope},
        {261,13.557162754476575,148.39160227554498,27.600685774175062,Pinching4StateKind::NegativeToPositive},
        {270,18.899260436396194,148.39160227554498,28.184901391610765,Pinching4StateKind::NegativeToPositive},
        {280,24.834924527417993,148.39160227554498,29.059585090887055,Pinching4StateKind::NegativeToPositive}};

    double work=0.0,previous_force=0.0,previous_deformation=0.0;
    for(int step=0;step<=280;++step){
        const double deformation=protocol[static_cast<std::size_t>(step)];
        if(step==121||step==141||step==161||step==181||step==201||step==221||step==241||step==261){
            const auto a=m.trial(deformation,state),b=m.trial(deformation,state);
            near(a.response.force,b.response.force,1e-13,"Pinching4 focused trial purity");
            const bool from_negative=step==141||step==181||step==221||step==261;
            const auto expected=from_negative?Pinching4StateKind::NegativeEnvelope:Pinching4StateKind::PositiveEnvelope;
            check(state.kind==expected,"Pinching4 focused trial mutated committed state");
        }
        const auto trial=m.trial(deformation,state);
        if(step>0)work+=0.5*(previous_force+trial.response.force)*(deformation-previous_deformation);
        for(const auto& point:points)if(point.step==step){
            near(trial.response.force,point.force,1e-11,"Pinching4 focused force oracle");
            near(trial.response.tangent,point.tangent,1e-11,"Pinching4 focused tangent oracle");
            near(work,point.work,1e-10,"Pinching4 focused work oracle");
            check(trial.state.kind==point.kind,"Pinching4 focused state mismatch");
        }
        previous_force=trial.response.force;
        previous_deformation=deformation;
        state=trial.state;
    }
    check(state.reversal_count==10,"Pinching4 full steel_raw reversal history mismatch");
    check(state.kind==Pinching4StateKind::NegativeToPositive,"Pinching4 full steel_raw terminal state mismatch");

    auto sm=shear_material();
    auto ss=sm.initial_state();
    std::vector<double> sp{0.0};
    auto segment=[&](double a,double b){for(int i=1;i<=20;++i)sp.push_back(a+(b-a)*i/20.0);};
    segment(0.0,0.0001);segment(0.0001,0.010);segment(0.010,0.10);segment(0.10,0.30);
    segment(0.30,0.0);segment(0.0,-0.30);segment(-0.30,0.72);segment(0.72,-0.72);
    segment(-0.72,1.20);segment(1.20,-1.20);segment(-1.20,0.0);
    const Point shear_points[]={
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
    double sw=0.0,spf=0.0,spd=0.0;
    for(int step=0;step<=220;++step){
        const double deformation=sp[static_cast<std::size_t>(step)];
        if(step==81||step==121||step==141||step==161||step==181||step==201){
            const auto a=sm.trial(deformation,ss),b=sm.trial(deformation,ss);
            near(a.response.force,b.response.force,1e-13,"shear Pinching4 trial purity");
        }
        const auto trial=sm.trial(deformation,ss);
        if(step>0)sw+=0.5*(spf+trial.response.force)*(deformation-spd);
        for(const auto& point:shear_points)if(point.step==step){
            near(trial.response.force,point.force,1e-11,"shear Pinching4 force oracle");
            near(trial.response.tangent,point.tangent,1e-11,"shear Pinching4 tangent oracle");
            near(sw,point.work,1e-10,"shear Pinching4 work oracle");
            check(trial.state.kind==point.kind,"shear Pinching4 state mismatch");
        }
        spf=trial.response.force;spd=deformation;ss=trial.state;
    }
    check(ss.reversal_count==6,"shear Pinching4 reversal history mismatch");
    check(ss.kind==Pinching4StateKind::NegativeToPositive,"shear Pinching4 terminal state mismatch");
    std::cout<<"Pinching4 full frozen steel_raw and shear oracle parity passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
