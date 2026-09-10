#include "quake/newmark.hpp"
#include "quake/ida.hpp"
#include "quake/analysis_input.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/shear_building.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace quake;
static void check(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}

// An exactly elastic SDOF with deliberately injected local integration faults.
// Its state stores the last accepted displacement, so rejected trials can be
// distinguished from legitimate small increments without hidden mutation.
class FaultModel final:public NonlinearDynamicModel {
    SparseMatrixCSC k_=SparseMatrixCSC::from_triplets(1,1,{{0,0,100.0}});
    SparseUpdateBasis b_=SparseUpdateBasis::from_dense(1,1,{1.0});
    std::vector<double> m_{1.0},t_{100.0};
public:
    enum Mode {Normal,NaNForce,IncrementFailure,ProgrammingError,InputGap};
    Mode mode;
    explicit FaultModel(Mode v=Normal):mode(v){}
    int dof()const override{return 1;}
    int nonlinear_count()const override{return 1;}
    int nonlinear_state_size()const override{return 1;}
    const SparseMatrixCSC& K_initial()const override{return k_;}
    const SparseUpdateBasis& nonlinear_basis()const override{return b_;}
    const std::vector<double>& mass()const override{return m_;}
    const std::vector<double>& initial_nonlinear_tangents()const override{return t_;}
    std::vector<double> mass_multiply(const std::vector<double>& v)const override{return v;}
    std::vector<double> damping_multiply(const std::vector<double>&)const override{return {0};}
    SparseMatrixCSC effective_initial_matrix(double a0,double)const override{return SparseMatrixCSC::from_triplets(1,1,{{0,0,100+a0}});}
    SparseMatrixCSC effective_tangent_matrix(const std::vector<double>&,double a0,double a1)const override{return effective_initial_matrix(a0,a1);}
    void evaluate_nonlinear_deformations(const std::vector<double>& u,const std::vector<double>& c,
        std::vector<double>& f,std::vector<double>& t,std::vector<double>& s)const override{
        if(mode==ProgrammingError)throw std::invalid_argument("bad model definition");
        if(mode==IncrementFailure&&std::abs(u[0]-c[0])>0.04)throw ConstitutiveIntegrationError("fixture increment too large");
        f={mode==NaNForce?std::numeric_limits<double>::quiet_NaN():100*u[0]};t=t_;s=u;
    }
    void internal_force_and_tangent(const std::vector<double>& u,const std::vector<double>& c,
        std::vector<double>& f,std::vector<double>& t,std::vector<double>& s)const override{evaluate_nonlinear_deformations(u,c,f,t,s);}
    std::vector<double> base_excitation(double ag)const override{
        if(mode==InputGap&&ag>1.4&&ag<1.6)return {std::numeric_limits<double>::quiet_NaN()};
        return {-ag};
    }
    double response_value(const std::vector<double>& u)const override{return u[0];}
    double max_drift_measure(const std::vector<double>& u)const override{return std::abs(u[0]);}
    double max_drift_ratio(const std::vector<double>& u)const override{return std::abs(u[0]);}
};

int main(){try{
    RobustNewmarkOptions o;o.max_subdivisions=4;o.return_numerical_failure=true;
    FaultModel finite;
    for(auto strategy:{LinearStrategy::FullFactorization,LinearStrategy::SamePatternRefactorization,LinearStrategy::Woodbury}){
        bool threw=false;try{run_newmark_robust(finite,{std::numeric_limits<double>::quiet_NaN()},.1,strategy,o);}catch(const std::invalid_argument&){threw=true;}
        check(threw,"NaN excitation must be rejected");
        FaultModel nan(FaultModel::NaNForce);auto bad=run_newmark_robust(nan,{1},.1,strategy,o);
        check(bad.termination==AnalysisTermination::NumericalFailure&&bad.stats.steps==0&&bad.stats.nonfinite_evaluations>0,"NaN residual must never converge");
        FaultModel recover(FaultModel::IncrementFailure);std::size_t calls=0;double last=0;
        o.accepted_substep_state_observer=[&](std::size_t,std::size_t,int,double time,double,const auto&,const auto&,const auto&,const auto&){check(time>last,"observer exposed a rolled-back substep");last=time;++calls;};
        auto good=run_newmark_robust(recover,{10},.2,strategy,o);
        check(good.termination==AnalysisTermination::Completed&&good.stats.constitutive_integration_failures>0&&calls>1,"local return failure must recover by subdivision");
        o.accepted_substep_state_observer={};
    }
    FaultModel invalid(FaultModel::ProgrammingError);bool threw=false;
    try{run_newmark_robust(invalid,{1},.1,LinearStrategy::Woodbury,o);}catch(const std::invalid_argument&){threw=true;}
    check(threw,"programming errors must propagate");
    IDAOptions ida;ida.scale_factors={.5,1,2,3};ida.stop_after_first_collapse=false;ida.collapse_refinement_steps=0;
    ida.analysis.collapse.max_story_drift_ratio=.003;ida.analysis.max_subdivisions=0;
    auto r=run_ida_suite(finite,{{"fixture",{1},.1}},ida);
    check(r.records[0].first_collapse_scale==2&&r.records[0].last_noncollapse_scale==1,"IDA overwrote first collapse bracket");
    ida.scale_factors={.5,1.5,2};FaultModel gap(FaultModel::InputGap);
    r=run_ida_suite(gap,{{"gap",{1},.1}},ida);
    check(r.records[0].collapsed&&r.records[0].numerical_failure&&r.records[0].bracket_has_numerical_gap,"IDA lost unresolved numerical gap");
    ida.scale_factors={1,-1};threw=false;
    try{run_ida_suite(finite,{{"bad",{1},.1}},ida);}catch(const std::invalid_argument&){threw=true;}
    check(threw,"IDA must reject invalid scales instead of silently dropping them");
    auto a=run_newmark_robust(finite,{1,2,-1,0},.01,LinearStrategy::Woodbury,o);
    check(a.stats.residual_evaluations==a.stats.newton_iterations,"line search evaluation was repeated at next Newton iteration");
    // Singular numerical factors must be discarded, then recover on the same
    // sparsity pattern. A stale factor must never remain callable.
    auto matrix=[](double a,double b,double c,double d){return SparseMatrixCSC::from_triplets(2,2,{{0,0,a},{0,1,b},{1,0,c},{1,1,d}});};
    auto nonsingular=matrix(2,1,1,3),singular=matrix(1,1,1,1);
    SuperLUSamePatternSolver direct(nonsingular);threw=false;
    try{direct.refactor_and_solve(singular,{1,2});}catch(const std::runtime_error&){threw=true;}
    check(threw,"singular refactor must fail");threw=false;
    try{direct.solve_current({1,2});}catch(const std::runtime_error&){threw=true;}
    check(threw,"singular refactor left stale factors active");
    auto x=direct.refactor_and_solve(nonsingular,{1,2});check(std::abs(x[0]-.2)<1e-12&&std::abs(x[1]-.6)<1e-12,"same-pattern recovery produced incorrect solution");
    ShearBuilding sh(2,1,20,{0,1},100,1,.02,0,.002);
    PreparedAdaptiveNewmark adaptive(sh,.01,0);
    auto gm=synthetic_ground_motion(200,.01,15);
    (void)adaptive.run(gm);auto again=adaptive.run(gm);
    auto reference=run_newmark(sh,gm,.01,LinearStrategy::FullFactorization);
    for(std::size_t i=0;i<gm.size();++i)check(std::abs(again.roof_history[i]-reference.roof_history[i])<1e-9,"prepared adaptive retained stale factors across runs");
    ida.scale_factors={.5,1,2};ida.reuse_preparation=true;auto reused=run_ida_suite(finite,{{"reuse",{1},.1}},ida);
    ida.reuse_preparation=false;auto fresh=run_ida_suite(finite,{{"reuse",{1},.1}},ida);
    for(std::size_t i=0;i<fresh.runs.size();++i)check(fresh.runs[i].max_roof_abs==reused.runs[i].max_roof_abs&&fresh.runs[i].termination==reused.runs[i].termination,"IDA preparation changed results");
    std::cout<<"Reliability regressions passed: finite input/state, recoverable integration, rollback, IDA brackets/gaps, and evaluation reuse.\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
