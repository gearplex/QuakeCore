#include "quake/frame2d.hpp"
#include "quake/fsc_shear_damage.hpp"
#include "quake/fsc_shear_spring.hpp"
#include "quake/frame3d.hpp"
#include "quake/corotational3d.hpp"
#include "quake/imk_peak_oriented.hpp"
#include "quake/asce41_hinge.hpp"
#include "quake/rc_column_asce41.hpp"
#include "quake/rc_building_asce41.hpp"
#include "quake/rc_frame_asce41.hpp"
#include "quake/model_comparison.hpp"
#include "quake/stability.hpp"
#include "quake/calibration.hpp"
#include "quake/substructure.hpp"
#include "quake/low_rank_solver.hpp"
#include "quake/generalized_woodbury.hpp"
#include "quake/newmark.hpp"
#include "quake/modal.hpp"
#include "quake/modal_damping.hpp"
#include "quake/reduction.hpp"
#include "quake/shear_building.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/suite.hpp"
#include "quake/ida.hpp"
#include "quake/static_analysis.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

using namespace quake;

static void require(bool cond, const char* msg) {
    if (!cond) throw std::runtime_error(msg);
}

static double rel_err(const std::vector<double>& a, const std::vector<double>& b) {
    double num=0.0, den=0.0;
    for (std::size_t i=0;i<a.size();++i) { num=std::max(num,std::abs(a[i]-b[i])); den=std::max(den,std::abs(a[i])); }
    return num/std::max(1.0,den);
}

static void test_sparse_solve() {
    auto A = SparseMatrixCSC::from_triplets(3,3,{
        {0,0,4},{1,0,-1},{0,1,-1},{1,1,4},{2,1,-1},{1,2,-1},{2,2,3}});
    std::vector<double> xtrue{1,2,3};
    auto b=A.multiply(xtrue);
    auto x=superlu_solve_once(A,b);
    require(rel_err(xtrue,x)<1e-12,"SuperLU solve accuracy");
    SuperLUSamePatternSolver repeated(A);
    auto xr=repeated.solve_current(b);
    require(rel_err(xtrue,xr)<1e-12,"SuperLU same-pattern current solve accuracy");
    auto A2=A; A2.values()[0]+=0.7; A2.values()[3]-=0.2;
    auto b2=A2.multiply(xtrue);
    auto xr2=repeated.refactor_and_solve(A2,b2);
    require(rel_err(xtrue,xr2)<2e-12,"SuperLU same-pattern refactor accuracy");
    require(repeated.factorizations()==2,"SuperLU same-pattern factorization count");
}

static void test_woodbury_equivalence() {
    const int n=8, m=3;
    std::vector<Triplet> t;
    for(int i=0;i<n;++i){ t.push_back({i,i,8.0}); if(i+1<n){t.push_back({i,i+1,-1.0});t.push_back({i+1,i,-1.0});} }
    auto A=SparseMatrixCSC::from_triplets(n,n,t);
    std::vector<double> B(static_cast<std::size_t>(n*m),0.0);
    for(int j=0;j<m;++j){int i=1+2*j;B[static_cast<std::size_t>(j*n+i)]=1;B[static_cast<std::size_t>(j*n+i-1)]=-1;}
    LowRankWoodburySolver w(A,B,m);
    std::vector<double> d{-1.7,0.0,-0.9};
    std::vector<Triplet> tu=t;
    for(int j=0;j<m;++j){
        for(int r=0;r<n;++r) for(int c=0;c<n;++c){
            const double v=d[static_cast<std::size_t>(j)]*B[static_cast<std::size_t>(j*n+r)]*B[static_cast<std::size_t>(j*n+c)];
            if(v!=0.0) tu.push_back({r,c,v});
        }
    }
    auto K=SparseMatrixCSC::from_triplets(n,n,tu,1e-18);
    std::vector<double> rhs{1,2,-1,0.5,3,-2,1,0.25};
    auto xf=superlu_solve_once(K,rhs);
    auto xw=w.solve(rhs,d);
    require(rel_err(xf,xw)<2e-12,"Woodbury exact equivalence");

    LazyLowRankWoodburySolver lazy(A,B,m);
    auto xl=lazy.solve(rhs,d);
    require(rel_err(xf,xl)<2e-12,"Lazy Woodbury exact equivalence");
    require(lazy.cached_update_count()==2,"Lazy Woodbury should cache only active directions");
    std::vector<double> d2{0.0,-0.4,0.0};
    std::vector<Triplet> tu2=t;
    for(int r=0;r<n;++r) for(int c=0;c<n;++c){
        const double v=d2[1]*B[static_cast<std::size_t>(n+r)]*B[static_cast<std::size_t>(n+c)];
        if(v!=0.0) tu2.push_back({r,c,v});
    }
    auto K2=SparseMatrixCSC::from_triplets(n,n,tu2,1e-18);
    auto xf2=superlu_solve_once(K2,rhs);
    auto xl2=lazy.solve(rhs,d2);
    require(rel_err(xf2,xl2)<2e-12,"Lazy Woodbury newly activated column equivalence");
    require(lazy.cached_update_count()==3,"Lazy Woodbury should grow cache on first activation");
}


static void test_generalized_woodbury_coupled_update(){
    const int n=10,r=4;std::vector<Triplet> t;
    for(int i=0;i<n;++i){t.push_back({i,i,8.0});if(i+1<n){t.push_back({i,i+1,-1.0});t.push_back({i+1,i,-1.0});}}
    auto A=SparseMatrixCSC::from_triplets(n,n,t);
    std::vector<std::vector<std::pair<int,double>>> cols(static_cast<std::size_t>(r));
    cols[0]={{1,1.0},{2,-1.0}};cols[1]={{4,1.0},{5,-1.0}};cols[2]={{6,0.7},{7,-0.4}};cols[3]={{2,0.5},{8,0.8}};
    auto U=SparseUpdateBasis::from_columns(n,cols);GeneralizedWoodburySolver w(A,U);
    std::vector<double> C(static_cast<std::size_t>(r*r),0.0);
    C[0*r+0]=-1.2;C[1*r+1]=-0.8;C[2*r+2]=0.4;C[0*r+3]=C[3*r+0]=0.25;C[1*r+2]=C[2*r+1]=-0.15;
    std::vector<Triplet> tt=t;
    for(int a=0;a<r;++a)for(int b=0;b<r;++b){const double cab=C[static_cast<std::size_t>(a*r+b)];if(cab==0)continue;for(const auto&[i,ui]:cols[static_cast<std::size_t>(a)])for(const auto&[j,uj]:cols[static_cast<std::size_t>(b)])tt.push_back({i,j,cab*ui*uj});}
    auto K=SparseMatrixCSC::from_triplets(n,n,tt);std::vector<double> rhs(static_cast<std::size_t>(n));for(int i=0;i<n;++i)rhs[static_cast<std::size_t>(i)]=0.3+0.17*i;
    auto x=w.solve(rhs,C),xr=superlu_solve_once(K,rhs);require(rel_err(xr,x)<2e-12,"generalized Woodbury coupled update equivalence");require(w.last_active_dimension()==4,"generalized Woodbury active dimension");
}

static void test_nrha_equivalence() {
    std::vector<int> nl{0,2,4,6,8,10,12,14,16,18};
    ShearBuilding model(20,1.0,220.0,nl,100.0,0.20,0.02,0.02,0.001);
    auto gm=synthetic_ground_motion(450,0.01,5.0);
    auto full=run_newmark(model,gm,0.01,LinearStrategy::FullFactorization,1e-8,25);
    auto modified=run_newmark(model,gm,0.01,LinearStrategy::ModifiedNewton,1e-8,40);
    auto wood=run_newmark(model,gm,0.01,LinearStrategy::Woodbury,1e-8,25);
    require(rel_err(full.roof_history,wood.roof_history)<2e-9,"NRHA roof histories differ");
    require(rel_err(full.roof_history,modified.roof_history)<2e-8,"Modified-Newton NRHA roof histories differ");
    require(full.stats.global_factorizations>modified.stats.global_factorizations,"Modified Newton should reduce factorizations");
    require(modified.stats.global_factorizations>wood.stats.global_factorizations,"Woodbury should reduce factorizations below modified Newton");
    require(wood.stats.global_factorizations==1,"Woodbury should use one global factorization");
}



static void test_robust_subdivision_recovers_failed_step() {
    std::vector<int> nl; for(int s=0;s<20;s+=2) nl.push_back(s);
    ShearBuilding model(20,1.0,220.0,nl,100.0,0.20,0.01,0.02,0.001);
    auto gm=synthetic_ground_motion(120,0.02,10.0);
    bool failed=false;
    try { (void)run_newmark(model,gm,0.02,LinearStrategy::Woodbury,1e-8,2); }
    catch(const std::runtime_error&) { failed=true; }
    require(failed,"coarse two-iteration reference case should fail without subdivision");
    RobustNewmarkOptions o; o.max_iterations=2; o.max_subdivisions=5; o.max_backtracks=8;
    std::size_t committed_substeps=0; bool saw_subdivision_depth=false;
    o.accepted_substep_state_observer=[&](std::size_t,std::size_t,int depth,double,double,
        const std::vector<double>&,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&){
        ++committed_substeps; saw_subdivision_depth=saw_subdivision_depth||depth>0;
    };
    auto robust=run_newmark_robust(model,gm,0.02,LinearStrategy::Woodbury,o);
    require(robust.stats.failed_steps==0,"robust subdivision should recover reference case");
    require(robust.stats.subdivided_steps>0,"robust reference case should exercise subdivision");
    require(robust.stats.minimum_dt<0.02,"robust reference case should reduce dt");
    require(robust.roof_history.size()==gm.size(),"robust response history should retain original record grid");
    require(committed_substeps>gm.size()&&saw_subdivision_depth,
            "rollback-safe observer should expose committed internal subdivision states");

    RobustNewmarkOptions fail=o;fail.max_iterations=1;fail.max_subdivisions=0;fail.return_numerical_failure=true;
    std::size_t max_observed_output_step=0;bool observed_any=false;
    fail.accepted_substep_state_observer=[&](std::size_t step,std::size_t,int,double,double,
        const std::vector<double>&,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&){
        observed_any=true;max_observed_output_step=std::max(max_observed_output_step,step);
    };
    auto failed_run=run_newmark_robust(model,gm,0.02,LinearStrategy::Woodbury,fail);
    require(failed_run.termination==AnalysisTermination::NumericalFailure,
            "rollback observer fixture should terminate numerically without subdivision");
    if(observed_any) require(max_observed_output_step<failed_run.stats.steps,
            "rolled-back failing output step must never be emitted as a committed substep");
}

static CompiledFrame2D make_cantilever(double axial_compression = 0.0) {
    Frame2DBuilder b;
    b.add_node(1,0.0,0.0);
    b.add_node(2,0.0,4.0);
    b.add_elastic_frame(1,1,2,30000.0,2.0,0.5,axial_compression);
    b.fix(1);
    b.set_response_node(2);
    b.set_story_nodes({2});
    return b.compile();
}

static void test_frame2d_cantilever() {
    auto model=make_cantilever();
    const int ux=model.reduced_dof(2,Dof2D::UX);
    require(ux>=0,"cantilever UX mapping");
    std::vector<double> rhs(static_cast<std::size_t>(model.dof()),0.0);
    rhs[static_cast<std::size_t>(ux)]=1.0;
    auto u=superlu_solve_once(model.K_initial(),rhs);
    const double exact=4.0*4.0*4.0/(3.0*30000.0*0.5);
    require(std::abs(u[static_cast<std::size_t>(ux)]-exact)<1e-11,"2D frame cantilever stiffness");

    auto pdelta=make_cantilever(120.0);
    auto up=superlu_solve_once(pdelta.K_initial(),rhs);
    require(up[static_cast<std::size_t>(ux)]>u[static_cast<std::size_t>(ux)],"compression should soften cantilever tangent");
}

static void test_frame2d_constraints_and_hinge_basis() {
    Frame2DBuilder b;
    b.add_node(1,0,0);
    b.add_node(2,0,3,1.0,0.0,0.0);
    b.add_node(3,5,3,2.0,0.0,0.0);
    b.add_elastic_frame(1,1,2,25000,1.0,0.4);
    b.fix(1);
    b.equal_dof(2,3,Dof2D::UX);
    b.add_rotational_spring(1,1,2,1000.0,5.0,0.02);
    b.set_response_node(2);
    b.set_story_nodes({2});
    auto model=b.compile();
    require(model.reduced_dof(2,Dof2D::UX)==model.reduced_dof(3,Dof2D::UX),"equalDOF compilation");
    require(model.dof()==5,"constraint-eliminated DOF count");
    require(model.nonlinear_count()==1,"compiled hinge count");
    const int rz=model.reduced_dof(2,Dof2D::RZ);
    require(std::abs(model.nonlinear_basis().value(rz,0)-1.0)<1e-15,"hinge generalized deformation basis");
    const int ux=model.reduced_dof(2,Dof2D::UX);
    require(std::abs(model.mass()[static_cast<std::size_t>(ux)]-3.0)<1e-15,"rigid-floor mass accumulation");
}

static CompiledFrame2D make_frame_nrha_model(int stories) {
    Frame2DBuilder b;
    constexpr double h=3.5, bay=6.0;
    constexpr double E=25000.0, Ac=2.0, Ic=0.55, Ab=1.5, Ib=0.35;
    for(int s=0;s<=stories;++s) {
        for(int c=0;c<2;++c) {
            const int id=100*s+c;
            const double mx=s==0?0.0:0.5;
            b.add_node(id,c*bay,s*h,mx,0.0,0.0);
            if(s==0) b.fix(id);
        }
        if(s>0) b.rigid_floor_x(100*s,{100*s+1});
    }
    for(int s=1;s<=stories;++s) {
        const double p=25.0*(stories-s+1);
        for(int c=0;c<2;++c) b.add_elastic_frame(1000+10*s+c,100*(s-1)+c,100*s+c,E,Ac,Ic,p);
        const int dl=10000+10*s, dr=dl+1;
        b.add_node(dl,0.0,s*h); b.add_node(dr,bay,s*h);
        b.equal_dof(100*s,dl,Dof2D::UX); b.equal_dof(100*s,dl,Dof2D::UY);
        b.equal_dof(100*s+1,dr,Dof2D::UX); b.equal_dof(100*s+1,dr,Dof2D::UY);
        b.add_elastic_frame(2000+s,dl,dr,E,Ab,Ib);
        b.add_rotational_spring(3000+2*s,100*s,dl,50000.0,5.0,0.02);
        b.add_rotational_spring(3000+2*s+1,100*s+1,dr,50000.0,5.0,0.02);
    }
    b.set_rayleigh(0.02,0.0005);
    b.set_response_node(100*stories);
    std::vector<int> story_nodes;
    for(int s=1;s<=stories;++s) story_nodes.push_back(100*s);
    b.set_story_nodes(std::move(story_nodes));
    return b.compile();
}

static void test_frame2d_nrha_equivalence() {
    auto model=make_frame_nrha_model(6);
    auto gm=synthetic_ground_motion(300,0.005,4.0);
    auto full=run_newmark(model,gm,0.005,LinearStrategy::FullFactorization,1e-8,30);
    auto wood=run_newmark(model,gm,0.005,LinearStrategy::Woodbury,1e-8,30);
    require(rel_err(full.roof_history,wood.roof_history)<2e-9,"frame NRHA roof histories differ");
    require(full.stats.newton_iterations==wood.stats.newton_iterations,"frame Newton path should match");
    require(wood.stats.global_factorizations==1,"frame Woodbury should use one global factorization");

    PreparedWoodburyNewmark prepared(model,0.005);
    auto cached=prepared.run(gm,1e-8,30);
    require(rel_err(wood.roof_history,cached.roof_history)<2e-10,"prepared Woodbury changes frame response");
    require(cached.stats.global_factorizations==0,"prepared record should not refactor baseline");
    require(prepared.setup_factorizations()==1,"prepared suite should have one setup factorization");
}


static void test_frame3d_cantilever() {
    Frame3DBuilder b;
    b.add_node(1,0,0,0);
    b.add_node(2,0,0,4);
    b.add_elastic_frame(1,1,2,30000.0,12000.0,2.0,0.2,0.5,0.8,0,1,0);
    b.fix(1);
    b.set_response(2,Dof3D::UX);
    b.set_story_nodes({2},Dof3D::UX);
    auto model=b.compile();
    std::vector<double> rhs(static_cast<std::size_t>(model.dof()),0.0);
    const int ux=model.reduced_dof(2,Dof3D::UX);
    const int uy=model.reduced_dof(2,Dof3D::UY);
    rhs[static_cast<std::size_t>(ux)]=1.0;
    auto u=superlu_solve_once(model.K_initial(),rhs);
    const double exact_x=64.0/(3.0*30000.0*0.5);
    require(std::abs(u[static_cast<std::size_t>(ux)]-exact_x)<1e-11,"3D cantilever x-bending stiffness");
    std::fill(rhs.begin(),rhs.end(),0.0);
    rhs[static_cast<std::size_t>(uy)]=1.0;
    u=superlu_solve_once(model.K_initial(),rhs);
    const double exact_y=64.0/(3.0*30000.0*0.8);
    require(std::abs(u[static_cast<std::size_t>(uy)]-exact_y)<1e-11,"3D cantilever y-bending stiffness");
}



static void test_frame3d_vector_hinge_operator() {
    Frame3DBuilder b;
    b.add_node(1,0,0,0); b.add_node(2,0,0,0);
    b.fix(1);
    b.fix_dof(2,Dof3D::UX); b.fix_dof(2,Dof3D::UY); b.fix_dof(2,Dof3D::UZ);
    b.add_rotational_vector_spring(1,1,2,1.0,1.0,0.0,1000.0,5.0,0.02);
    b.set_response(2,Dof3D::RX); b.set_story_nodes({2},Dof3D::RX);
    auto model=b.compile();
    const int rx=model.reduced_dof(2,Dof3D::RX), ry=model.reduced_dof(2,Dof3D::RY), rz=model.reduced_dof(2,Dof3D::RZ);
    require(rx>=0 && ry>=0 && rz>=0,"vector hinge rotational DOFs");
    const double c=1.0/std::sqrt(2.0);
    require(std::abs(model.nonlinear_basis().value(rx,0)-c)<1e-14,"vector hinge RX projection");
    require(std::abs(model.nonlinear_basis().value(ry,0)-c)<1e-14,"vector hinge RY projection");
    require(std::abs(model.nonlinear_basis().value(rz,0))<1e-14,"vector hinge RZ projection");
}

static void test_frame3d_rigid_diaphragm_mpc() {
    Frame3DBuilder b;
    b.add_node(1,0,0,3,1.0,2.0,0,0,0,0);
    b.add_node(2,2,0,3,3.0,4.0,0,0,0,0);
    b.add_node(3,0,3,3,5.0,6.0,0,0,0,0);
    // Retain only the in-plane rigid-diaphragm generalized DOFs.
    for (int id : {1,2,3}) {
        b.fix_dof(id,Dof3D::UZ); b.fix_dof(id,Dof3D::RX); b.fix_dof(id,Dof3D::RY);
    }
    b.rigid_diaphragm_z(1,{2,3});
    b.set_ground_direction(Dof3D::UX);
    b.set_response(3,Dof3D::UX); // dependent response: ux_master - 3*rz_master
    b.set_story_nodes({3},Dof3D::UX);
    auto model=b.compile();
    require(model.dof()==3,"rigid diaphragm should condense to UX/UY/RZ master DOFs");
    const int ux=model.reduced_dof(1,Dof3D::UX), uy=model.reduced_dof(1,Dof3D::UY), rz=model.reduced_dof(1,Dof3D::RZ);
    require(ux>=0 && uy>=0 && rz>=0,"rigid diaphragm master mappings");
    require(model.reduced_dof(3,Dof3D::UX)==-1,"offset slave UX must be a multi-term MPC row");

    std::vector<double> q(static_cast<std::size_t>(model.dof()),0.0);
    q[static_cast<std::size_t>(ux)]=0.10; q[static_cast<std::size_t>(rz)]=0.02;
    require(std::abs(model.response_value(q)-0.04)<1e-14,"rigid diaphragm kinematic coupling");

    std::vector<double> ex(static_cast<std::size_t>(model.dof()),0.0); ex[static_cast<std::size_t>(ux)]=1.0;
    auto mex=model.mass_multiply(ex);
    require(std::abs(mex[static_cast<std::size_t>(ux)]-9.0)<1e-13,"rigid diaphragm generalized Mxx");
    require(std::abs(mex[static_cast<std::size_t>(rz)]+15.0)<1e-13,"rigid diaphragm generalized M-rz coupling");
    std::vector<double> er(static_cast<std::size_t>(model.dof()),0.0); er[static_cast<std::size_t>(rz)]=1.0;
    auto mer=model.mass_multiply(er);
    require(std::abs(mer[static_cast<std::size_t>(rz)]-61.0)<1e-13,"rigid diaphragm generalized rotational mass");
    require(std::abs(mer[static_cast<std::size_t>(uy)]-8.0)<1e-13,"rigid diaphragm generalized My-rz coupling");
    auto p=model.base_excitation(1.0);
    require(std::abs(p[static_cast<std::size_t>(ux)]+9.0)<1e-13,"rigid diaphragm base translation load");
    require(std::abs(p[static_cast<std::size_t>(rz)]-15.0)<1e-13,"rigid diaphragm torsional base load");
}

static int n3(int story,int ix,int iy){return 10000*story+100*ix+iy;}

static CompiledFrame3D make_space_frame_nrha_model(int stories) {
    Frame3DBuilder b;
    constexpr double h=3.5,L=6.0,E=25000.0,G=10000.0;
    for(int s=0;s<=stories;++s){
        for(int ix=0;ix<2;++ix) for(int iy=0;iy<2;++iy){
            const int id=n3(s,ix,iy);
            const double m=s==0?0.0:0.25;
            b.add_node(id,ix*L,iy*L,s*h,m,m,0.05*m,0,0,0);
            if(s==0)b.fix(id);
        }
    }
    int eid=1,sid=1,did=500000;
    for(int s=1;s<=stories;++s){
        const double p=8.0*(stories-s+1);
        for(int ix=0;ix<2;++ix) for(int iy=0;iy<2;++iy)
            b.add_elastic_frame(eid++,n3(s-1,ix,iy),n3(s,ix,iy),E,G,2.0,0.30,0.55,0.55,1,0,0,p);
        // X-direction beams: hinge rotation about global Y.
        for(int iy=0;iy<2;++iy){
            const int jl=n3(s,0,iy), jr=n3(s,1,iy); const int dl=did++,dr=did++;
            b.add_node(dl,0,iy*L,s*h);b.add_node(dr,L,iy*L,s*h);
            for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RX,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}
            b.add_elastic_frame(eid++,dl,dr,E,G,1.5,0.20,0.35,0.35,0,1,0);
            b.add_bilinear_spring(sid++,jl,Dof3D::RY,dl,Dof3D::RY,50000,8.0,0.50);
            b.add_bilinear_spring(sid++,jr,Dof3D::RY,dr,Dof3D::RY,50000,8.0,0.50);
        }
        // Y-direction beams: hinge rotation about global X.
        for(int ix=0;ix<2;++ix){
            const int jl=n3(s,ix,0), jr=n3(s,ix,1); const int dl=did++,dr=did++;
            b.add_node(dl,ix*L,0,s*h);b.add_node(dr,ix*L,L,s*h);
            for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RY,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}
            b.add_elastic_frame(eid++,dl,dr,E,G,1.5,0.20,0.35,0.35,1,0,0);
            b.add_bilinear_spring(sid++,jl,Dof3D::RX,dl,Dof3D::RX,50000,8.0,0.50);
            b.add_bilinear_spring(sid++,jr,Dof3D::RX,dr,Dof3D::RX,50000,8.0,0.50);
        }
    }
    b.set_rayleigh(0.02,0.0005);
    b.set_ground_direction(Dof3D::UX);
    b.set_response(n3(stories,0,0),Dof3D::UX);
    std::vector<int> story_nodes;for(int s=1;s<=stories;++s)story_nodes.push_back(n3(s,0,0));
    b.set_story_nodes(std::move(story_nodes),Dof3D::UX);
    return b.compile();
}

static void test_frame3d_nrha_equivalence() {
    auto model=make_space_frame_nrha_model(3);
    require(model.nonlinear_basis().nnz()==2*model.nonlinear_count(),"3D concentrated hinge basis should remain sparse");
    auto gm=synthetic_ground_motion(240,0.005,1.25);
    auto full=run_newmark(model,gm,0.005,LinearStrategy::FullFactorization,1e-8,35);
    auto same=run_newmark(model,gm,0.005,LinearStrategy::SamePatternRefactorization,1e-8,35);
    auto wood=run_newmark(model,gm,0.005,LinearStrategy::Woodbury,1e-8,35);
    RobustNewmarkOptions ro; ro.max_iterations=35; ro.max_subdivisions=0;
    auto robust=run_newmark_robust(model,gm,0.005,LinearStrategy::Woodbury,ro);
    require(rel_err(full.roof_history,same.roof_history)<3e-9,"3D same-pattern NRHA roof histories differ");
    require(rel_err(full.roof_history,wood.roof_history)<3e-9,"3D frame NRHA roof histories differ");
    require(rel_err(full.roof_history,robust.roof_history)<3e-9,"robust 3D frame NRHA roof histories differ");
    require(full.stats.newton_iterations==same.stats.newton_iterations,"3D same-pattern Newton path should match");
    require(full.stats.newton_iterations==wood.stats.newton_iterations,"3D frame Newton path should match");
    require(wood.stats.global_factorizations==1,"3D frame Woodbury should use one baseline factorization");
}

static void test_modal_analysis_singular_mass() {
    // 1-DOF analytical oscillator: lambda=k/m=4, omega=2, period=pi.
    ShearBuilding sdof(1,2.0,8.0,{},100.0,1.0,0.01,0.0,0.0);
    auto modes=modal_analysis(sdof,1);
    require(modes.size()==1,"SDOF modal mode count");
    require(std::abs(modes[0].eigenvalue-4.0)<1e-10,"SDOF modal eigenvalue");
    require(std::abs(modes[0].omega-2.0)<1e-10,"SDOF modal omega");
    require(std::abs(modes[0].period-3.14159265358979323846)<1e-10,"SDOF modal period");

    // A frame contains massless rotational/internal DOFs. DGGEV should retain
    // the finite physical modes rather than requiring positive-definite M.
    auto frame=make_space_frame_nrha_model(2);
    auto fm=modal_analysis(frame,3);
    require(!fm.empty(),"frame modal analysis with singular mass");
    for(const auto& m:fm) require(std::isfinite(m.period)&&m.period>0.0,"finite positive frame period");
}



static void test_fixed_modal_damping_force_side_wrapper() {
    // Analytical SDOF: m=2, k=8 -> omega=2 rad/s.  With zeta=5%,
    // c_modal = 2*zeta*omega*m = 0.4.  Modal damping remains physically
    // force-side and does not densify the model's sparse effective matrix;
    // the robust Newton solver may nevertheless use its exact rank-1
    // velocity Jacobian through a Woodbury correction.
    ShearBuilding base(1,2.0,8.0,{},100.0,1.0,0.01,0.0,0.0);
    FixedModalDampingModel modal(base,0.05,1);
    require(modal.modal_count()==1,"modal damping SDOF mode count");
    std::vector<double> v{3.0};
    auto f=modal.damping_multiply(v);
    require(f.size()==1 && std::abs(f[0]-1.2)<1e-11,"modal damping SDOF force");
    auto A0=base.effective_initial_matrix(123.0,4.5);
    auto A1=modal.effective_initial_matrix(123.0,4.5);
    require(A0.col_ptr()==A1.col_ptr() && A0.row_ind()==A1.row_ind(),"modal damping effective matrix topology");
    require(rel_err(A0.values(),A1.values())<1e-14,"modal damping must remain outside sparse model matrix");
    require(!modal.has_additional_effective_low_rank_update(),"modal damping defaults to Perform/OpenSees force-side Newton behavior");
    FixedModalDampingModel modal_exact(base,0.05,1,true);
    require(modal_exact.has_additional_effective_low_rank_update(),"modal damping opt-in exact low-rank Newton update");
    const auto& U=modal_exact.additional_effective_low_rank_basis();
    auto C=modal_exact.additional_effective_low_rank_coefficients(123.0,4.5);
    require(U.rows()==1&&U.cols()==1&&C.size()==1,"modal damping low-rank dimensions");
    const double extra=C[0]*U.value(0,0)*U.value(0,0);
    require(std::abs(extra-1.8)<1e-11,"modal damping exact Newmark tangent contribution");
}

static void test_guyan_static_condensation() {
    ShearBuilding full(4,1.0,200.0,{},100.0,1.0,0.01,0.0,0.0);
    ReductionBuildInfo info;
    auto reduced=craig_bampton_reduce(full,{0,3},0,&info);
    require(info.full_dof==4 && info.retained_dof==2 && info.reduced_dof==2,"Guyan reduction dimensions");
    std::vector<double> ff(4,0.0); ff[3]=1.0;
    auto uf=superlu_solve_once(full.K_initial(),ff);
    std::vector<double> fr(2,0.0); fr[1]=1.0;
    auto qr=superlu_solve_once(reduced.K_initial(),fr);
    auto ur=reduced.expand(qr);
    require(std::abs(uf[3]-ur[3])<1e-11,"Guyan retained static displacement");
    require(std::abs(uf[0]-ur[0])<1e-11,"Guyan retained lower displacement");
}

static void test_craig_bampton_modal_enrichment() {
    ShearBuilding full(8,1.0,220.0,{0,7},100.0,10.0,0.1,0.0,0.001);
    auto retained=nonlinear_support_dofs(full);
    auto guyan=craig_bampton_reduce(full,retained,0);
    auto cb=craig_bampton_reduce(full,retained,3);
    auto mf=modal_analysis(full,2);
    auto mg=modal_analysis(guyan,2);
    auto mc=modal_analysis(cb,2);
    require(!mf.empty()&&!mg.empty()&&!mc.empty(),"Craig-Bampton modal results");
    const double eg=std::abs(mg[0].period-mf[0].period)/mf[0].period;
    const double ec=std::abs(mc[0].period-mf[0].period)/mf[0].period;
    require(ec<eg,"fixed-interface modes should improve first-mode period over Guyan-only reduction");
    require(ec<0.03,"Craig-Bampton first-mode period accuracy");
}


static void test_parallel_suite_equivalence() {
    ShearBuilding model(12,1.0,220.0,{0,3,6,9},100.0,0.20,0.02,0.02,0.001);
    std::vector<std::vector<double>> motions;
    for(int r=0;r<4;++r) motions.push_back(synthetic_ground_motion(100,0.01,3.0+0.4*r));
    PreparedWoodburyNewmark prepared(model,0.01);
    std::vector<AnalysisResult> sequential;
    for(const auto& gm:motions) sequential.push_back(prepared.run(gm,1e-8,25));
    auto parallel=run_record_suite_parallel(model,motions,0.01,2,1e-8,25);
    require(parallel.records.size()==sequential.size(),"parallel suite result count");
    require(parallel.workers==2,"parallel suite worker count");
    for(std::size_t i=0;i<sequential.size();++i)
        require(rel_err(sequential[i].roof_history,parallel.records[i].roof_history)<2e-10,"parallel suite changed response");
}


static IMKPeakOrientedParams basic_imk_params(double lambda=1e9) {
    IMKPeakOrientedParams p;
    p.Ke=1000.0;
    p.posUp=p.negUp=0.02;
    p.posUpc=p.negUpc=0.04;
    p.posUu=p.negUu=0.12;
    p.posFy=p.negFy=10.0;
    p.posFcapFy=p.negFcapFy=1.20;
    p.posFresFy=p.negFresFy=0.20;
    p.lambdaS=p.lambdaC=p.lambdaA=p.lambdaK=lambda;
    p.cS=p.cC=p.cA=p.cK=1.0;
    p.Dpos=p.Dneg=1.0;
    return p;
}

static void test_imk_monotonic_backbone() {
    IMKPeakOrientedMaterial imk(basic_imk_params());
    std::vector<double> c(IMKPeakOrientedMaterial::kStateSize),t(IMKPeakOrientedMaterial::kStateSize);
    imk.initialize_state(c.data());
    auto a=imk.trial(0.005,c.data(),t.data());
    require(std::abs(a.force-5.0)<1e-12 && std::abs(a.tangent-1000.0)<1e-12,"IMK elastic branch");
    c=t; auto b=imk.trial(0.020,c.data(),t.data());
    require(std::abs(b.force-11.0)<1e-10 && std::abs(b.tangent-100.0)<1e-10,"IMK post-yield branch");
    c=t; auto d=imk.trial(0.050,c.data(),t.data());
    require(std::abs(d.force-7.0)<1e-9 && std::abs(d.tangent+250.0)<1e-9,"IMK post-capping branch");
    c=t; auto e=imk.trial(0.080,c.data(),t.data());
    require(std::abs(e.force-2.0)<1e-9 && std::abs(e.tangent)<1e-12,"IMK residual branch");
    c=t; auto f=imk.trial(0.121,c.data(),t.data());
    require(std::abs(f.force)<1e-12 && (f.events&IMK_EVENT_FAILURE),"IMK ultimate failure");
}

static void test_imk_reversal_and_deterioration() {
    IMKPeakOrientedMaterial imk(basic_imk_params(0.20));
    std::vector<double> c(IMKPeakOrientedMaterial::kStateSize),t(IMKPeakOrientedMaterial::kStateSize);
    imk.initialize_state(c.data());
    c=t=c;
    auto p1=imk.trial(0.030,c.data(),t.data()); c=t;
    auto p2=imk.trial(0.045,c.data(),t.data()); c=t;
    (void)p1; (void)p2;
    auto rev=imk.trial(0.040,c.data(),t.data());
    require((rev.events&IMK_EVENT_REVERSAL)!=0,"IMK reversal event");
    require((rev.events&IMK_EVENT_DETERIORATION)!=0,"IMK deterioration event");
    require(t[4] < 1000.0,"IMK unloading stiffness deterioration");
    require(t[10] < 1.0 || t[12] < 1.0,"IMK opposite-direction backbone deterioration");
}

static void test_frame3d_mixed_material_state_bank() {
    Frame3DBuilder b;
    b.add_node(1,0,0,0); b.add_node(2,0,0,3,1,1,1); b.fix(1);
    b.add_elastic_frame(1,1,2,25000,10000,2.0,0.3,0.5,0.5,1,0,0);
    b.add_bilinear_spring(1,1,Dof3D::RY,2,Dof3D::RY,5000,5,0.05);
    auto p=basic_imk_params(); p.Ke=4000; p.posFy=p.negFy=4.0;
    b.add_imk_peak_oriented_spring(2,1,Dof3D::RZ,2,Dof3D::RZ,p);
    b.set_response(2,Dof3D::UX); b.set_story_nodes({2},Dof3D::UX);
    auto model=b.compile();
    require(model.nonlinear_count()==2,"mixed material component count");
    require(model.nonlinear_state_size()==2+IMKPeakOrientedMaterial::kStateSize,"mixed material state size");
    auto st=model.initial_nonlinear_state();
    require(std::abs(st[2+4]-p.Ke)<1e-12,"IMK initial unload stiffness in flat state bank");
    std::vector<double> q{0.002,0.005},cf,k,trial;
    model.evaluate_nonlinear_deformations(q,st,cf,k,trial);
    require(cf.size()==2 && k.size()==2 && trial.size()==st.size(),"mixed material evaluation dimensions");
}


static void test_story_block_schur_exactness() {
    auto model=make_space_frame_nrha_model(6);
    const double dt=0.005,a0=1.0/(0.25*dt*dt),a1=0.5/(0.25*dt);
    auto A=model.effective_initial_matrix(a0,a1);
    auto partition=make_story_block_partition(model,2);
    BlockSchurFactor block(A,partition);
    SuperLUFactor full(A);
    std::vector<double> rhs(static_cast<std::size_t>(model.dof()));
    for(int i=0;i<model.dof();++i) rhs[static_cast<std::size_t>(i)]=std::sin(0.031*(i+1));
    auto xf=full.solve(rhs), xb=block.solve(rhs);
    require(rel_err(xf,xb)<2e-10,"story-block Schur solve differs from monolithic factorization");
    require(block.stats().interface_dof<model.dof(),"story-block partition should eliminate some linear DOFs");

    LazyLowRankWoodburySolver wfull(A,model.nonlinear_basis());
    SubstructuredLazyWoodburySolver wblock(A,model.nonlinear_basis(),make_story_block_partition(model,2));
    std::vector<double> dk(static_cast<std::size_t>(model.nonlinear_count()),0.0);
    for(int j=0;j<std::min(8,model.nonlinear_count());++j) dk[static_cast<std::size_t>(j)]=-0.8*model.initial_nonlinear_tangents()[static_cast<std::size_t>(j)];
    auto xw=wfull.solve(rhs,dk), xbw=wblock.solve(rhs,dk);
    require(rel_err(xw,xbw)<3e-10,"substructured Woodbury differs from monolithic Woodbury");
}

static void test_solver_calibration_smoke() {
    auto model=make_space_frame_nrha_model(3);
    auto c=calibrate_solver_crossover(model,0.005,{0,1,4,8},1);
    require(c.points.size()==4,"solver calibration point count");
    require(c.dof==model.dof() && c.nonlinear_count==model.nonlinear_count(),"solver calibration model metadata");
    require(c.recommended_woodbury_rank_limit>=-1 && c.recommended_woodbury_rank_limit<=model.nonlinear_count(),"solver calibration rank limit bounds");
    for(const auto& p:c.points) require(p.woodbury_seconds>=0.0&&p.same_pattern_refactor_seconds>=0.0,"solver calibration timing values");
}


static void test_prepared_adaptive_equivalence() {
    auto model=make_space_frame_nrha_model(3);
    auto gm=synthetic_ground_motion(220,0.005,2.0);
    auto wood=run_newmark(model,gm,0.005,LinearStrategy::Woodbury,1e-8,35);
    PreparedAdaptiveNewmark adaptive(model,0.005,2); // force both paths during yielding
    auto a=adaptive.run(gm,1e-8,35);
    require(rel_err(wood.roof_history,a.roof_history)<5e-9,"adaptive calibrated solver changed NRHA response");
    require(a.stats.adaptive_woodbury_solves+a.stats.adaptive_direct_solves==a.stats.linear_solves,"adaptive solve accounting");
}



static ASCE41HingeParams basic_asce41_params(){
    ASCE41HingeParams p;
    p.Ke=1000.0;p.posFy=p.negFy=10.0;p.hardening_ratio=0.02;
    p.pos_a=p.neg_a=0.02;p.pos_b=p.neg_b=0.08;p.pos_f=p.neg_f=0.11;p.pos_c=p.neg_c=0.20;
    p.pos_drop_span=p.neg_drop_span=0.005;
    p.pos_io=p.neg_io=0.01;p.pos_ls=p.neg_ls=0.025;p.pos_cp=p.neg_cp=0.05;
    return p;
}

static void test_asce41_physical_hardening_stiffness(){
    auto p=basic_asce41_params();
    p.hardening_ratio=0.50; // deliberately very different from the explicit slope
    p.hardening_stiffness=3.0;
    ASCE41HingeMaterial h(p);
    std::vector<double> c(ASCE41HingeMaterial::kStateSize),t(c.size());
    h.initialize_state(c.data());
    const double uy=p.posFy/p.Ke;
    const double q=uy+0.01;
    auto e=h.trial(q,c.data(),t.data());
    require(std::abs(e.tangent-3.0)<1e-12,"ASCE41 explicit physical hardening tangent");
    require(std::abs(e.force-(p.posFy+3.0*0.01))<1e-10,"ASCE41 explicit physical hardening force");
}

static void test_asce41_degrading_backbone_and_limits(){
    auto p=basic_asce41_params();ASCE41HingeMaterial h(p);
    std::vector<double> c(ASCE41HingeMaterial::kStateSize),t(c.size());h.initialize_state(c.data());
    const double uy=p.posFy/p.Ke;
    auto e=h.trial(0.5*uy,c.data(),t.data());require(std::abs(e.tangent-p.Ke)<1e-12,"ASCE41 elastic branch");
    c=t; e=h.trial(uy+p.pos_a,c.data(),t.data());require(e.force>p.posFy,"ASCE41 hardening to C");
    c=t; e=h.trial(uy+p.pos_a+0.5*p.pos_drop_span,c.data(),t.data());require(e.tangent<0.0,"ASCE41 C-D negative degradation tangent");
    c=t; e=h.trial(uy+p.pos_a+p.pos_drop_span+0.01,c.data(),t.data());require(std::abs(e.force-p.pos_c*p.posFy)<1e-8,"ASCE41 residual branch");
    require(h.performance_level(uy+0.012)==ASCE41PerformanceLevel::IO,"ASCE41 IO classification");
    require(h.performance_level(uy+0.03)==ASCE41PerformanceLevel::LS,"ASCE41 LS classification");
    require(h.performance_level(uy+0.05)==ASCE41PerformanceLevel::CP,"ASCE41 CP classification");
    require(h.performance_level(uy+0.06)==ASCE41PerformanceLevel::BeyondCP,"ASCE41 beyond-CP classification");
    c=t;e=h.trial(uy+p.pos_b+1e-4,c.data(),t.data());require((e.events&ASCE41_EVENT_LATERAL_LOSS)!=0&&!(e.events&ASCE41_EVENT_FAILURE)&&std::abs(e.force)<1e-14,"ASCE41 E-point lateral resistance loss should precede F");
    c=t;e=h.trial(uy+p.pos_f+1e-4,c.data(),t.data());require((e.events&ASCE41_EVENT_FAILURE)!=0&&std::abs(e.force)<1e-14,"ASCE41 F-point effective/gravity resistance loss");
}


static void test_asce41_straight_ce_and_explicit_C_strength(){
    auto p=basic_asce41_params();
    p.backbone_shape=ASCE41BackboneShape::StraightCE;
    p.posMc=p.negMc=11.5;
    p.pos_f=p.pos_b; p.neg_f=p.neg_b;
    ASCE41HingeMaterial h(p);
    std::vector<double> c(ASCE41HingeMaterial::kStateSize),t(c.size());h.initialize_state(c.data());
    const double uy=p.posFy/p.Ke;
    auto at_c=h.trial(uy+p.pos_a,c.data(),t.data());
    require(std::abs(at_c.force-p.posMc)<1e-10,"ASCE41 explicit C strength");
    const double qp=0.5*(p.pos_a+p.pos_b);
    auto mid=h.trial(uy+qp,c.data(),t.data());
    const double expected=p.posMc*(p.pos_b-qp)/(p.pos_b-p.pos_a);
    require(std::abs(mid.force-expected)<1e-10,"ASCE41 straight C-E interpolation");
    require(mid.tangent<0.0,"ASCE41 straight C-E negative tangent");
    auto at_e=h.trial(uy+p.pos_b+1e-8,c.data(),t.data());
    require((at_e.events&ASCE41_EVENT_LATERAL_LOSS)!=0 && (at_e.events&ASCE41_EVENT_FAILURE)!=0,
            "ASCE41 straight benchmark E should coincide with F when configured");
}

static void test_rc_column_asce41_provider_separates_code_backbone_from_hysteresis(){
    RCColumnResolvedParameters r;
    r.numerical_hinge_Ke=100000.0;
    r.posMy=r.negMy=100.0; r.posMc=r.negMc=105.0;
    r.pos_a=r.neg_a=0.01; r.pos_b=r.neg_b=0.04;
    r.pos_c=r.neg_c=0.20;
    r.pos_io=r.neg_io=0.004; r.pos_ls=r.neg_ls=0.012; r.pos_cp=r.neg_cp=0.022;
    auto nist=RCColumnASCE41Provider::nist_asce41_17_benchmark(r,"unit-test resolved benchmark values");
    require(nist.hinge.backbone_shape==ASCE41BackboneShape::StraightCE,"NIST provider straight C-E topology");
    require(nist.hinge.pos_f==nist.hinge.pos_b&&nist.hinge.neg_f==nist.hinge.neg_b,"NIST provider E=F semantics");
    require(!nist.code_coefficients_embedded,"provider must not imply embedded copyrighted code tables");

    r.pos_f=r.neg_f=0.06; r.pos_drop_span=r.neg_drop_span=0.005; r.pos_e_drop_span=r.neg_e_drop_span=0.005;
    auto prod=RCColumnASCE41Provider::asce41_23_aci369_1_22(
        r,ASCE41BackboneShape::ResearchExtendedCDE,"authorized ACI 369.1-22 resolved parameters");
    require(prod.hinge.backbone_shape==ASCE41BackboneShape::ResearchExtendedCDE,"production provider preserves resolved topology");
    require(prod.hinge.pos_f==0.06,"production provider preserves separately resolved F semantics");
}


static RCColumnResolvedParameters synthetic_rc_column_params(double compression_kip){
    RCColumnResolvedParameters r;
    r.numerical_hinge_Ke=100000.0;
    // Synthetic demand dependence for architecture tests only: larger
    // compression reduces flexural strength and deformation capacity.
    const double scale=std::max(0.60,1.0-0.001*compression_kip);
    r.posMy=r.negMy=120.0*scale;
    r.posMc=r.negMc=126.0*scale;
    r.pos_a=r.neg_a=0.030*scale;
    r.pos_b=r.neg_b=0.080*scale;
    r.pos_c=r.neg_c=0.20;
    r.pos_io=r.neg_io=0.008*scale;
    r.pos_ls=r.neg_ls=0.025*scale;
    r.pos_cp=r.neg_cp=0.050*scale;
    return r;
}

static void test_rc_column_section_rules_boundary_and_audit(){
    RCColumnSectionInput sec;
    sec.component_id="C1";
    sec.width_in=12.0; sec.depth_in=12.0; sec.clear_length_in=120.0;
    sec.expected_fc_ksi=5.0; sec.expected_fy_long_ksi=70.0;
    sec.longitudinal_bar_count=8; sec.longitudinal_bar_area_in2=0.11;
    sec.gravity_axial_compression_kip=80.0;

    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput& s,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(d.max_compression_kip);
        rr.demand_used=d;
        rr.provenance="authorized-rules-test-fixture";
        rr.audit.push_back({"pos_a",rr.resolved.pos_a,"fixture/table-row-1",
                            s.conforming_transverse_reinforcement?"conforming":"nonconforming"});
        return rr;
    });
    RCColumnDemandState d; d.max_compression_kip=100.0;
    auto rr=resolver.resolve(sec,d);
    auto spec=RCColumnASCE41Provider::asce41_23_aci369_1_22(rr,ASCE41BackboneShape::StraightCE);
    require(spec.provenance=="authorized-rules-test-fixture","RC rules provenance retained");
    require(spec.audit.size()==1&&spec.audit[0].parameter=="pos_a","RC parameter audit retained");
    require(std::abs(spec.demand_used.max_compression_kip-100.0)<1e-12,"RC demand provenance retained");
    require(!spec.code_coefficients_embedded,"RC provider must not claim embedded code coefficients");
}

static void test_rc_column_axial_iteration_converges_and_regenerates(){
    RCColumnSectionInput sec;
    sec.component_id="B1";
    sec.width_in=12.0; sec.depth_in=12.0; sec.clear_length_in=120.0;
    sec.expected_fc_ksi=5.0; sec.expected_fy_long_ksi=70.0;
    sec.longitudinal_bar_count=8; sec.longitudinal_bar_area_in2=0.049;
    sec.gravity_axial_compression_kip=80.0;
    sec.initial_shear_demand_kip=20.0;

    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(d.max_compression_kip);
        rr.demand_used=d;
        rr.provenance="synthetic axial-sensitive rules fixture";
        rr.audit.push_back({"posMy",rr.resolved.posMy,"fixture","max compression controls"});
        return rr;
    });

    // First analysis discovers 200 kip max compression. Reanalysis with the
    // regenerated hinge sees the same envelope, so the fixed point is reached.
    RCColumnResponseRunner runner=[](const RCColumnModelSpec&,std::size_t){
        RCColumnDemandState d;
        d.max_compression_kip=200.0;
        d.max_tension_kip=10.0;
        d.max_abs_shear_kip=30.0;
        return d;
    };
    RCColumnAxialIterationOptions opt;
    opt.max_iterations=5; opt.min_iterations=2;
    opt.demand_relative_tolerance=1e-12;
    opt.parameter_relative_tolerance=1e-12;
    auto result=RCColumnAxialIteration::run_asce41_23_aci369_1_22(
        sec,resolver,runner,ASCE41BackboneShape::StraightCE,opt);
    require(result.converged,"RC axial iteration should converge");
    require(result.steps.size()==2,"RC axial iteration should require regeneration/reanalysis");
    require(result.steps[0].parameter_relative_change>0.0,"RC axial iteration should detect parameter change");
    require(result.steps[1].parameter_relative_change<1e-12,"RC axial iteration parameters should stabilize");
    require(std::abs(result.final_model.demand_used.max_compression_kip-200.0)<1e-12,
            "RC final model should be generated at converged max compression");
    const auto initial=synthetic_rc_column_params(80.0);
    const auto converged=synthetic_rc_column_params(200.0);
    require(result.final_model.hinge.posFy<initial.posMy,"RC max compression should regenerate lower synthetic strength");
    require(std::abs(result.final_model.hinge.posFy-converged.posMy)<1e-12,"RC final strength should match converged resolver state");
}


static void test_rc_column_two_pass_convenience(){
    RCColumnSectionInput sec;
    sec.component_id="A1";
    sec.width_in=12.0; sec.depth_in=12.0; sec.clear_length_in=120.0;
    sec.expected_fc_ksi=5.0; sec.expected_fy_long_ksi=70.0;
    sec.longitudinal_bar_count=8; sec.longitudinal_bar_area_in2=0.049;
    sec.gravity_axial_compression_kip=60.0;
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(d.max_compression_kip);
        rr.provenance="two-pass fixture";
        return rr;
    });
    RCColumnResponseRunner runner=[](const RCColumnModelSpec&,std::size_t){
        RCColumnDemandState d; d.max_compression_kip=180.0; return d;
    };
    auto result=RCColumnAxialIteration::run_two_pass_asce41_23_aci369_1_22(
        sec,resolver,runner,ASCE41BackboneShape::StraightCE);
    require(result.steps.size()==2,"RC two-pass convenience must execute exactly two analyses");
    require(std::abs(result.steps[1].demand_input.max_compression_kip-180.0)<1e-12,
            "RC two-pass second hinge must be regenerated at first-pass max compression");
}

static void test_rc_column_axial_iteration_under_relaxation(){
    RCColumnSectionInput sec;
    sec.component_id="C2";
    sec.width_in=12.0; sec.depth_in=12.0; sec.clear_length_in=120.0;
    sec.expected_fc_ksi=5.0; sec.expected_fy_long_ksi=64.0;
    sec.longitudinal_bar_count=8; sec.longitudinal_bar_area_in2=0.11;
    sec.gravity_axial_compression_kip=100.0;

    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(d.max_compression_kip);
        rr.provenance="synthetic relaxation fixture";
        return rr;
    });
    RCColumnResponseRunner runner=[](const RCColumnModelSpec&,std::size_t){
        RCColumnDemandState d; d.max_compression_kip=300.0; return d;
    };
    RCColumnAxialIterationOptions opt;
    opt.max_iterations=3; opt.min_iterations=2; opt.relaxation=0.5;
    opt.demand_relative_tolerance=0.0; opt.parameter_relative_tolerance=0.0;
    auto result=RCColumnAxialIteration::run_asce41_23_aci369_1_22(
        sec,resolver,runner,ASCE41BackboneShape::StraightCE,opt);
    require(!result.converged&&result.steps.size()==3,"RC under-relaxed fixed point may remain unconverged at iteration cap");
    require(std::abs(result.steps[1].demand_input.max_compression_kip-200.0)<1e-12,
            "RC under-relaxation should update halfway toward observed demand");
}


static RCColumnSectionInput synthetic_building_column(
    const std::string& id,double gravity,double bar_area){
    RCColumnSectionInput sec;
    sec.component_id=id;
    sec.width_in=12.0; sec.depth_in=12.0; sec.clear_length_in=120.0;
    sec.expected_fc_ksi=5.0; sec.expected_fy_long_ksi=70.0;
    sec.longitudinal_bar_count=8; sec.longitudinal_bar_area_in2=bar_area;
    sec.gravity_axial_compression_kip=gravity;
    sec.initial_shear_demand_kip=15.0;
    return sec;
}

static void test_rc_building_iteration_synchronous_and_order_invariant(){
    RCBuildingColumnDefinition b1{synthetic_building_column("B1",80.0,0.049),
                                  ASCE41BackboneShape::StraightCE};
    RCBuildingColumnDefinition c1{synthetic_building_column("C1",100.0,0.11),
                                  ASCE41BackboneShape::StraightCE};

    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput& s,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(d.max_compression_kip);
        // Give the components distinct strengths so cross-coupling is visible.
        if(s.component_id=="C1"){
            rr.resolved.posMy*=1.25; rr.resolved.negMy*=1.25;
            rr.resolved.posMc*=1.25; rr.resolved.negMc*=1.25;
        }
        rr.provenance="whole-building synchronous fixture";
        rr.audit.push_back({"posMy",rr.resolved.posMy,"fixture",s.component_id});
        return rr;
    });

    RCBuildingResponseRunner runner=[](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        require(models.size()==2,"whole-building runner should receive all coordinated models together");
        // Coordinator canonicalizes by component_id, so this response couples
        // each column to the OTHER column's current model. A sequential updater
        // would produce a different second-column demand.
        require(models[0].component_id=="B1"&&models[1].component_id=="C1",
                "whole-building models should be in canonical deterministic order");
        RCBuildingAnalysisObservation o;
        o.status="synthetic global analysis complete";
        const double b_obs=50.0+models[1].model.hinge.posFy;
        const double c_obs=70.0+models[0].model.hinge.posFy;
        // Deliberately return observations in reverse order to verify id-based matching.
        o.column_demands.push_back({"C1",{c_obs,5.0,25.0}});
        o.column_demands.push_back({"B1",{b_obs,4.0,24.0}});
        return o;
    };

    RCColumnAxialIterationOptions opt;
    opt.max_iterations=2; opt.min_iterations=2; opt.relaxation=1.0;
    opt.demand_relative_tolerance=0.0; opt.parameter_relative_tolerance=0.0;
    auto r1=RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(
        {c1,b1},resolver,runner,opt);
    auto r2=RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22(
        {b1,c1},resolver,runner,opt);

    require(r1.steps.size()==2&&r2.steps.size()==2,"whole-building coordinator should execute requested iterations");
    require(r1.steps[0].models[0].component_id=="B1"&&r1.steps[0].models[1].component_id=="C1",
            "whole-building coordinator should canonicalize input order");
    require(r1.final_models.size()==2&&r2.final_models.size()==2,"whole-building final model field size");
    for(std::size_t i=0;i<2;++i){
        require(r1.final_models[i].component_id==r2.final_models[i].component_id,
                "whole-building result ids should be order invariant");
        require(std::abs(r1.final_models[i].model.hinge.posFy-r2.final_models[i].model.hinge.posFy)<1e-12,
                "whole-building synchronous result should be input-order invariant");
    }

    const double c_initial_my=synthetic_rc_column_params(100.0).posMy*1.25;
    const double b_initial_my=synthetic_rc_column_params(80.0).posMy;
    const double expected_b_target=50.0+c_initial_my;
    const double expected_c_target=70.0+b_initial_my;
    const auto& metrics=r1.steps[1].component_metrics;
    require(metrics.size()==2&&metrics[0].component_id=="B1"&&metrics[1].component_id=="C1",
            "whole-building metrics should retain canonical component mapping");
    require(std::abs(metrics[0].demand_input.max_compression_kip-expected_b_target)<1e-12,
            "B1 second-pass demand must use C1 model from the same prior global analysis");
    require(std::abs(metrics[1].demand_input.max_compression_kip-expected_c_target)<1e-12,
            "C1 second-pass demand must use B1 model from the same prior global analysis");
}

static void test_rc_building_iteration_converges_on_global_worst_column(){
    RCBuildingColumnDefinition a{synthetic_building_column("A1",60.0,0.049),
                                 ASCE41BackboneShape::StraightCE};
    RCBuildingColumnDefinition d{synthetic_building_column("D1",90.0,0.11),
                                 ASCE41BackboneShape::StraightCE};
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& demand){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(demand.max_compression_kip);
        rr.provenance="whole-building convergence fixture";
        return rr;
    });
    RCBuildingResponseRunner runner=[](const std::vector<RCBuildingColumnModel>&,std::size_t){
        RCBuildingAnalysisObservation o;
        o.column_demands={{"A1",{180.0,0.0,20.0}}, {"D1",{240.0,0.0,35.0}}};
        return o;
    };
    RCColumnAxialIterationOptions opt;
    opt.max_iterations=5; opt.min_iterations=2;
    opt.demand_relative_tolerance=1e-12; opt.parameter_relative_tolerance=1e-12;
    auto r=RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22({d,a},resolver,runner,opt);
    require(r.converged&&r.termination==RCBuildingASCE41Termination::Converged,
            "whole-building parameter field should converge");
    require(r.steps.size()==2,"constant observed building envelope should converge after regeneration/reanalysis");
    require(r.steps[0].global_parameter_relative_change>0.0,
            "whole-building first pass should detect at least one changed hinge");
    require(r.steps[1].global_parameter_relative_change<1e-12,
            "whole-building regenerated parameter field should stabilize");
    require(r.final_models[0].component_id=="A1"&&r.final_models[1].component_id=="D1",
            "whole-building final models should remain canonical");
    require(std::abs(r.final_models[0].model.demand_used.max_compression_kip-180.0)<1e-12,
            "A1 final hinge should use its converged demand");
    require(std::abs(r.final_models[1].model.demand_used.max_compression_kip-240.0)<1e-12,
            "D1 final hinge should use its converged demand");
}

static void test_rc_building_two_pass_and_failure_semantics(){
    RCBuildingColumnDefinition b{synthetic_building_column("B2",75.0,0.049),
                                 ASCE41BackboneShape::StraightCE};
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& demand){
        RCColumnRuleResolution rr;
        rr.resolved=synthetic_rc_column_params(demand.max_compression_kip);
        rr.provenance="building two-pass fixture";
        return rr;
    });
    RCBuildingResponseRunner ok=[](const std::vector<RCBuildingColumnModel>&,std::size_t){
        RCBuildingAnalysisObservation o;
        o.column_demands={{"B2",{190.0,2.0,31.0}}};
        return o;
    };
    auto two=RCBuildingASCE41Iteration::run_two_pass_asce41_23_aci369_1_22({b},resolver,ok);
    require(two.steps.size()==2,"whole-building two-pass helper must execute exactly two global analyses");
    require(std::abs(two.steps[1].component_metrics[0].demand_input.max_compression_kip-190.0)<1e-12,
            "whole-building two-pass helper should regenerate all hinges at pass-one demand");

    RCBuildingResponseRunner failed=[](const std::vector<RCBuildingColumnModel>&,std::size_t){
        RCBuildingAnalysisObservation o;
        o.analysis_succeeded=false;
        o.status="synthetic NRHA nonconvergence";
        return o;
    };
    auto f=RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22({b},resolver,failed);
    require(!f.converged&&f.termination==RCBuildingASCE41Termination::AnalysisFailure,
            "whole-building analysis failure must not be mislabeled parameter nonconvergence");
    require(f.termination_detail=="synthetic NRHA nonconvergence",
            "whole-building analysis failure detail should be preserved");
    require(f.final_models.size()==1,"failure result should retain the model field that was analyzed");

    RCBuildingResponseRunner collapsed=[](const std::vector<RCBuildingColumnModel>&,std::size_t){
        RCBuildingAnalysisObservation o;o.analysis_succeeded=false;o.physical_collapse=true;o.status="synthetic physical collapse";return o;
    };
    auto pc=RCBuildingASCE41Iteration::run_asce41_23_aci369_1_22({b},resolver,collapsed);
    require(pc.termination==RCBuildingASCE41Termination::PhysicalCollapse,
            "whole-building physical collapse must remain distinct from numerical analysis failure");
}


static void test_native_element_force_recovery_2d_and_3d(){
    {
        Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,0,4);b.fix(1);
        b.add_elastic_frame(11,1,2,30000.0,2.0,0.5,20.0);b.set_response_node(2);b.set_story_nodes({2});
        auto m=b.compile();std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0);
        const int uy=m.reduced_dof(2,Dof2D::UY),ux=m.reduced_dof(2,Dof2D::UX);
        u[static_cast<std::size_t>(uy)]=-0.01;u[static_cast<std::size_t>(ux)]=0.001;
        const auto r=m.elastic_element_response(11,u);
        require(std::abs(r.axial_compression-170.0)<1e-10,"2D native axial force recovery must include preload plus elastic increment");
        require(std::max(std::abs(r.shear_i),std::abs(r.shear_j))>0.0,"2D native shear recovery");
    }
    {
        Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,4);b.fix(1);
        b.add_elastic_frame(21,1,2,30000.0,12000.0,2.0,0.3,0.5,0.5,0,1,0,20.0);
        b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);
        auto m=b.compile();std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0);
        const int uz=m.reduced_dof(2,Dof3D::UZ),ux=m.reduced_dof(2,Dof3D::UX);
        u[static_cast<std::size_t>(uz)]=-0.01;u[static_cast<std::size_t>(ux)]=0.001;
        const auto r=m.elastic_element_response(21,u);
        require(std::abs(r.axial_compression-170.0)<1e-10,"3D native axial force recovery must include preload plus elastic increment");
        const double vmax=std::max(std::hypot(r.shear_y_i,r.shear_z_i),std::hypot(r.shear_y_j,r.shear_z_j));
        require(vmax>0.0,"3D native vector shear recovery");
    }
}


static void test_native_frame2d_building_iteration_uses_asce_material_bank(){
    RCBuildingColumnDefinition c{synthetic_building_column("C2D",20.0,0.11),ASCE41BackboneShape::StraightCE};
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;rr.resolved=synthetic_rc_column_params(d.max_compression_kip);rr.provenance="native 2D fixture";return rr;
    });
    RCNativeFrame2DFactory factory=[](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,120,0,1.0,0,0);b.fix(1);
        b.add_elastic_frame(101,1,2,3000.0,36.0,108.0,20.0);
        b.add_asce41_hinge(201,1,2,models.at(0).model.hinge);
        b.set_response_node(2);b.set_story_nodes({2});
        return RCNativeFrame2DBuild{b.compile(),{{"C2D",{101},{201}}}};
    };
    RCNativeFrameAnalysisSettings settings;settings.dt=0.01;settings.strategy=LinearStrategy::SamePatternRefactorization;
    settings.ground_accel={0.0,5.0,-5.0,0.0,3.0,-3.0,0.0};settings.newmark_options.max_iterations=20;settings.newmark_options.max_subdivisions=3;
    auto r=RCBuildingASCE41NativeFrame::run_frame2d_two_pass_asce41_23_aci369_1_22({c},resolver,factory,settings);
    require(r.steps.size()==2,"native 2D ASCE material-bank workflow should execute two passes");
    require(r.steps[0].component_metrics[0].demand_observed.max_compression_kip>20.0,
            "native 2D recorder should capture axial compression from committed element response");
}

static void test_native_frame3d_building_iteration_records_committed_demands(){
    RCBuildingColumnDefinition c{synthetic_building_column("C1",20.0,0.11),ASCE41BackboneShape::StraightCE};
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;rr.resolved=synthetic_rc_column_params(d.max_compression_kip);rr.provenance="native frame fixture";return rr;
    });
    std::size_t factory_calls=0;
    RCNativeFrame3DFactory factory=[&](const std::vector<RCBuildingColumnModel>& models,std::size_t){
        ++factory_calls;require(models.size()==1&&models[0].component_id=="C1","native factory receives coordinated model field");
        Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,120,0,0,1.0);b.fix(1);
        b.add_elastic_frame(100,1,2,3000.0,1200.0,36.0,20.0,108.0,108.0,0,1,0,20.0);
        b.add_asce41_hinge(200,1,Dof3D::RY,2,Dof3D::RY,models[0].model.hinge);
        b.set_ground_direction(Dof3D::UZ);b.set_response(2,Dof3D::UZ);b.set_story_nodes({2},Dof3D::UZ);
        RCNativeFrame3DBuild out{b.compile(),{{"C1",{100},{200}}}};return out;
    };
    RCNativeFrameAnalysisSettings settings;settings.dt=0.01;settings.strategy=LinearStrategy::SamePatternRefactorization;
    settings.ground_accel={0.0,5.0,-5.0,0.0,3.0,-3.0,0.0};settings.newmark_options.max_iterations=20;settings.newmark_options.max_subdivisions=3;
    auto r=RCBuildingASCE41NativeFrame::run_frame3d_two_pass_asce41_23_aci369_1_22({c},resolver,factory,settings);
    require(factory_calls==2&&r.steps.size()==2,"native frame two-pass workflow should rebuild and analyze twice");
    require(r.steps[0].component_metrics[0].demand_observed.max_compression_kip>20.0,
            "native committed recorder should capture dynamic compression above gravity preload");
    require(r.steps[1].component_metrics[0].demand_input.max_compression_kip>20.0,
            "native second pass should regenerate the hinge from recorded peak compression");
}

static void test_material_field_replacement_preserves_compiled_topology(){
    auto p0=basic_asce41_params();
    Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,120,0,1.0,0,0);b.fix(1);
    b.add_elastic_frame(101,1,2,3000.0,36.0,108.0,20.0);b.add_asce41_hinge(201,1,2,p0);
    b.set_response_node(2);b.set_story_nodes({2});auto m=b.compile();
    const auto Klin=m.K_linear().values();const auto bcp=m.nonlinear_basis().col_ptr();
    const auto bri=m.nonlinear_basis().row_ind();const auto bval=m.nonlinear_basis().values();
    const int nd=m.dof(),ns=m.nonlinear_state_size();const double oldk=m.initial_nonlinear_tangents().at(0);
    auto p1=p0;p1.Ke=1.7*p0.Ke;p1.posFy*=0.85;p1.negFy*=0.85;
    m.replace_nonlinear_material(201,NonlinearMaterial{ASCE41HingeMaterial(p1)});
    require(m.dof()==nd&&m.nonlinear_state_size()==ns,"material-field replacement must preserve DOF/state topology");
    require(m.K_linear().values()==Klin,"material-field replacement must preserve elastic topology");
    require(m.nonlinear_basis().col_ptr()==bcp&&m.nonlinear_basis().row_ind()==bri&&m.nonlinear_basis().values()==bval,
            "material-field replacement must preserve nonlinear update basis");
    require(std::abs(m.initial_nonlinear_tangents().at(0)-p1.Ke)<1e-12&&std::abs(oldk-p1.Ke)>1e-9,
            "material-field replacement must update physical/numerical initial tangent");
    const auto* hp=m.materials().at(0).asce41_params();require(hp&&std::abs(hp->posFy-p1.posFy)<1e-12,
            "material-field replacement must install the regenerated ASCE hinge");

    bool rejected=false;
    try{m.replace_nonlinear_material(201,NonlinearMaterial{BilinearSpring(100.0,10.0,0.01)});}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"immutable compiled state layout must reject material kinds with different state size");
}

static void test_frame3d_material_field_matches_fresh_compile(){
    auto p0=basic_asce41_params(),p1=p0;p1.Ke*=1.35;p1.posFy*=0.82;p1.negFy*=0.82;p1.posMc=p1.posFy*1.05;p1.negMc=p1.negFy*1.05;
    auto build=[&](const ASCE41HingeParams& hp){
        Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,120,0,0,1.0);b.fix(1);
        b.add_elastic_frame(100,1,2,3000.0,1200.0,36.0,20.0,108.0,108.0,0,1,0,20.0);
        b.add_asce41_hinge(200,1,Dof3D::RY,2,Dof3D::RY,hp);
        b.set_ground_direction(Dof3D::UZ);b.set_response(2,Dof3D::UZ);b.set_story_nodes({2},Dof3D::UZ);return b.compile();
    };
    auto reused=build(p0),fresh=build(p1);const auto basis_before=reused.nonlinear_basis().values();
    reused.replace_nonlinear_material(200,NonlinearMaterial{ASCE41HingeMaterial(p1)});
    require(reused.nonlinear_basis().values()==basis_before,"3D material replacement must preserve nonlinear basis");
    require(reused.K_initial().col_ptr()==fresh.K_initial().col_ptr()&&reused.K_initial().row_ind()==fresh.K_initial().row_ind(),
            "3D material replacement and fresh compile must have identical Kinitial pattern");
    require(rel_err(reused.K_initial().values(),fresh.K_initial().values())<1e-14,
            "3D in-place material field must reproduce fresh-compile Kinitial exactly");
    auto sr=reused.initial_nonlinear_state(),sf=fresh.initial_nonlinear_state();
    std::vector<double> u(static_cast<std::size_t>(reused.dof()),0.0);const int ry=reused.reduced_dof(2,Dof3D::RY);u[static_cast<std::size_t>(ry)]=0.015;
    std::vector<double> fr,ff,tr,tf,xr,xf;reused.internal_force_and_tangent(u,sr,fr,tr,xr);fresh.internal_force_and_tangent(u,sf,ff,tf,xf);
    require(rel_err(fr,ff)<1e-14&&rel_err(tr,tf)<1e-14,"3D reused material field response must match fresh compile");
}

static void test_compiled_topology_asce_iteration_reuses_one_frame(){
    RCBuildingColumnDefinition c{synthetic_building_column("CT",20.0,0.11),ASCE41BackboneShape::StraightCE};
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;rr.resolved=synthetic_rc_column_params(d.max_compression_kip);rr.provenance="compiled topology fixture";return rr;
    });
    auto initial=RCColumnASCE41Provider::asce41_23_aci369_1_22(
        synthetic_rc_column_params(20.0),ASCE41BackboneShape::StraightCE,"compiled topology seed");
    Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,120,0,1.0,0,0);b.fix(1);
    b.add_elastic_frame(101,1,2,3000.0,36.0,108.0,20.0);b.add_asce41_hinge(201,1,2,initial.hinge);
    b.set_response_node(2);b.set_story_nodes({2});
    RCNativeCompiledFrame2D frame{b.compile(),{{"CT",{101},{201}}}};
    const auto basis_before=frame.model.nonlinear_basis().values();const int dof_before=frame.model.dof();
    RCNativeFrameAnalysisSettings settings;settings.dt=0.01;settings.strategy=LinearStrategy::SamePatternRefactorization;
    settings.ground_accel={0.0,5.0,-5.0,0.0,3.0,-3.0,0.0};settings.newmark_options.max_iterations=20;settings.newmark_options.max_subdivisions=3;
    auto r=RCBuildingASCE41NativeFrame::run_compiled_frame2d_two_pass_asce41_23_aci369_1_22({c},resolver,frame,settings);
    require(r.steps.size()==2,"compiled-topology ASCE workflow must execute two global analyses");
    require(frame.model.dof()==dof_before&&frame.model.nonlinear_basis().values()==basis_before,
            "compiled-topology ASCE workflow must retain topology/basis across passes");
    const int hi=frame.model.nonlinear_component_index(201);const auto* hp=frame.model.materials().at(static_cast<std::size_t>(hi)).asce41_params();
    require(hp!=nullptr&&std::abs(hp->posFy-r.final_models.at(0).model.hinge.posFy)<1e-12,
            "compiled frame must retain the final regenerated material field");
}

static void test_compiled_topology_frame3d_two_pass(){
    RCBuildingColumnDefinition c{synthetic_building_column("C3T",20.0,0.11),ASCE41BackboneShape::StraightCE};
    CallbackRCColumnRulesResolver resolver([](const RCColumnSectionInput&,const RCColumnDemandState& d){
        RCColumnRuleResolution rr;rr.resolved=synthetic_rc_column_params(d.max_compression_kip);rr.provenance="compiled 3D topology fixture";return rr;
    });
    auto initial=RCColumnASCE41Provider::asce41_23_aci369_1_22(synthetic_rc_column_params(20.0),ASCE41BackboneShape::StraightCE,"compiled 3D seed");
    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,120,0,0,1.0);b.fix(1);
    b.add_elastic_frame(100,1,2,3000.0,1200.0,36.0,20.0,108.0,108.0,0,1,0,20.0);
    b.add_asce41_hinge(200,1,Dof3D::RY,2,Dof3D::RY,initial.hinge);
    b.set_ground_direction(Dof3D::UZ);b.set_response(2,Dof3D::UZ);b.set_story_nodes({2},Dof3D::UZ);
    RCNativeCompiledFrame3D frame{b.compile(),{{"C3T",{100},{200}}}};const auto pattern=frame.model.K_initial().col_ptr();
    RCNativeFrameAnalysisSettings settings;settings.dt=0.01;settings.strategy=LinearStrategy::SamePatternRefactorization;
    settings.ground_accel={0.0,5.0,-5.0,0.0,3.0,-3.0,0.0};settings.newmark_options.max_iterations=20;settings.newmark_options.max_subdivisions=3;
    auto r=RCBuildingASCE41NativeFrame::run_compiled_frame3d_two_pass_asce41_23_aci369_1_22({c},resolver,frame,settings);
    require(r.steps.size()==2&&frame.model.K_initial().col_ptr()==pattern,"compiled 3D ASCE two-pass must reuse fixed sparse topology");
    require(r.steps[0].component_metrics[0].demand_observed.max_compression_kip>20.0,"compiled 3D native demand recorder must remain active");
}

static void test_corotational_section_force_recovery_is_objective(){
    CorotationalFrame3DProperties p;
    p.xi=0;p.yi=0;p.zi=0;p.xj=3;p.yj=0;p.zj=0;
    p.E=30000;p.G=12000;p.A=2.0;p.J=0.3;p.Iy=0.5;p.Iz=0.7;p.reference={0,1,0};p.axial_compression=50.0;
    const double th=0.63,c=std::cos(th),ss=std::sin(th);
    std::array<double,12> rigid{};rigid[6]=3*c-3;rigid[7]=3*ss;rigid[5]=th;rigid[11]=th;
    auto rr=corotational3d_section_response(p,rigid);
    require(std::abs(rr.axial_compression-50.0)<2e-8,"objective corotational recorder must retain preload under rigid rotation");
    require(std::max({std::abs(rr.shear_y_i),std::abs(rr.shear_z_i),std::abs(rr.moment_y_i),std::abs(rr.moment_z_i)})<1e-7,
            "rigid-body rotation must not generate corotational section shear/moment");

    const double ext=0.001;
    std::array<double,12> ur=rigid;ur[6]=(3+ext)*c-3;ur[7]=(3+ext)*ss;
    std::array<double,12> ua{};ua[6]=ext;
    const auto ar=corotational3d_section_response(p,ur),aa=corotational3d_section_response(p,ua);
    const double expected=50.0-p.E*p.A/3.0*ext;
    require(std::abs(ar.axial_compression-expected)<2e-7&&std::abs(aa.axial_compression-expected)<2e-9,
            "corotational physical axial force must be objective under superposed rigid rotation");
    require(std::abs(ar.axial_compression-aa.axial_compression)<2e-7,"rotated/unrotated section demand parity");

    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,3,0,0);b.fix(1);
    b.add_elastic_frame(77,1,2,p.E,p.G,p.A,p.J,p.Iy,p.Iz,0,1,0,p.axial_compression);
    b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);b.set_corotational(true);auto m=b.compile();
    std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0);const int ux=m.reduced_dof(2,Dof3D::UX);u[static_cast<std::size_t>(ux)]=ext;
    const auto er=m.elastic_element_response(77,u);
    require(std::abs(er.axial_compression-expected)<2e-7,"compiled corotational frame must expose native objective section demand");
}


static PMInteractionHingeParams basic_pm_return_params(){
    PMInteractionHingeParams p;
    p.hinge=basic_asce41_params();
    p.hinge.Ke=1000.0;p.hinge.posFy=p.hinge.negFy=10.0;
    p.hinge.posMc=p.hinge.negMc=0.0;
    p.hinge.pos_a=p.hinge.neg_a=0.02;
    p.hinge.pos_b=p.hinge.neg_b=0.06;
    p.hinge.pos_f=p.hinge.neg_f=0.09;
    p.hinge.backbone_shape=ASCE41BackboneShape::StraightCE;
    p.axial_preload=40.0;p.axial_stiffness=2000.0;
    p.axial_force_points={0.0,40.0,80.0,120.0};
    p.moment_capacity_points={12.0,10.0,7.0,4.0};
    p.return_tolerance=1e-11;p.max_return_iterations=40;
    return p;
}

static void test_pm_interaction_return_mapping_and_consistent_tangent(){
    auto p=basic_pm_return_params();p.surface_evolution=PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic;p.enforce_deformation_capacity=false;PMInteractionHinge2D h(p);
    std::vector<double> c(PMInteractionHinge2D::kStateSize),t(c.size());h.initialize_state(c.data());
    auto e=h.trial(0.0,0.004,c.data(),t.data());
    require(e.converged&&(e.events&PM_EVENT_ELASTIC)!=0,"P-M return hinge elastic predictor");
    require(std::abs(e.tangent[0]-p.axial_stiffness)<1e-12&&std::abs(e.tangent[3]-p.hinge.Ke)<1e-12,"P-M elastic tangent");

    auto r=h.trial(-0.001,0.015,c.data(),t.data());
    require(r.converged&&(r.events&PM_EVENT_YIELD)!=0,"P-M return mapping must converge into yield surface");
    require(r.positive_plastic_rotation>0.0,"P-M return mapping plastic rotation");
    require(std::abs(r.plastic_axial_deformation)>1e-9,"associative P-M yielding must generate axial plastic deformation");
    
    const double hd=2e-7,hr=2e-8;
    auto eval=[&](double da,double th){std::vector<double> st(c.size());return h.trial(da,th,c.data(),st.data());};
    const auto ap=eval(-0.001+hd,0.015),am=eval(-0.001-hd,0.015);
    const auto rp=eval(-0.001,0.015+hr),rm=eval(-0.001,0.015-hr);
    const double k00=(ap.axial_force-am.axial_force)/(2*hd),k10=(ap.moment-am.moment)/(2*hd);
    const double k01=(rp.axial_force-rm.axial_force)/(2*hr),k11=(rp.moment-rm.moment)/(2*hr);
        require(std::abs(k00-r.tangent[0])<3e-5*std::max(1.0,std::abs(k00)),"P-M local KNN finite-difference parity");
    require(std::abs(k01-r.tangent[1])<3e-5*std::max(1.0,std::abs(k01)),"P-M local KNM finite-difference parity");
    require(std::abs(k10-r.tangent[2])<3e-5*std::max(1.0,std::abs(k10)),"P-M local KMN finite-difference parity");
    require(std::abs(k11-r.tangent[3])<3e-5*std::max(1.0,std::abs(k11)),"P-M local KMM finite-difference parity");
}


static void test_perform_concrete_pm_surface_return_mapping(){
    auto p=basic_pm_return_params();
    p.surface_shape=PMInteractionSurfaceShape::PerformConcrete;
    p.surface_evolution=PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic;
    p.enforce_deformation_capacity=false;
    p.p_balance=45.0;p.py_tension=-30.0;p.py_compression=140.0;p.my_balance=125.0;
    p.alpha_tension=1.50;p.alpha_compression=1.70;p.beta_pm=1.10;
    p.axial_preload=20.0;p.axial_stiffness=5000.0;p.hinge.Ke=5000.0;
    PMInteractionHinge2D h(p);std::vector<double> c(PMInteractionHinge2D::kStateSize),t(c.size());h.initialize_state(c.data());
    const auto r=h.trial(-0.001,0.03,c.data(),t.data());
    require(r.converged&&(r.events&PM_EVENT_YIELD)!=0,"PERFORM concrete P-M surface return mapping");
    require(std::abs(r.plastic_axial_deformation)>1e-8,"PERFORM concrete normal flow must generate axial plastic deformation");
    const double hd=1e-7,hr=1e-8;
    auto eval=[&](double da,double th){std::vector<double> st(c.size());return h.trial(da,th,c.data(),st.data());};
    const auto ap=eval(-0.001+hd,0.03),am=eval(-0.001-hd,0.03),rp=eval(-0.001,0.03+hr),rm=eval(-0.001,0.03-hr);
    const double k00=(ap.axial_force-am.axial_force)/(2*hd),k10=(ap.moment-am.moment)/(2*hd);
    const double k01=(rp.axial_force-rm.axial_force)/(2*hr),k11=(rp.moment-rm.moment)/(2*hr);
    require(std::abs(k00-r.tangent[0])<2e-4*std::max(1.0,std::abs(k00)),"PERFORM surface KNN finite-difference parity");
    require(std::abs(k01-r.tangent[1])<2e-4*std::max(1.0,std::abs(k01)),"PERFORM surface KNM finite-difference parity");
    require(std::abs(k10-r.tangent[2])<2e-4*std::max(1.0,std::abs(k10)),"PERFORM surface KMN finite-difference parity");
    require(std::abs(k11-r.tangent[3])<2e-4*std::max(1.0,std::abs(k11)),"PERFORM surface KMM finite-difference parity");
}


static void test_perform_mroz_two_surface_translation(){
    auto p=basic_pm_return_params();
    p.surface_shape=PMInteractionSurfaceShape::PerformConcrete;
    p.surface_evolution=PMInteractionSurfaceEvolution::MrozTwoSurface;
    p.enforce_deformation_capacity=false;
    p.p_balance=45.0;p.py_tension=-30.0;p.py_compression=140.0;p.my_balance=25.0;
    p.alpha_tension=1.50;p.alpha_compression=1.70;p.beta_pm=1.10;
    p.axial_preload=20.0;p.axial_stiffness=5000.0;p.hinge.Ke=5000.0;
    p.hinge.posFy=p.hinge.negFy=20.0;p.hinge.hardening_stiffness=200.0;
    p.hinge.pos_a=p.hinge.neg_a=0.02;p.mroz_outer_scale=1.20;
    PMInteractionHinge2D h(p);std::vector<double> c(PMInteractionHinge2D::kStateSize),t(c.size());h.initialize_state(c.data());
    auto r=h.trial(-0.0005,0.010,c.data(),t.data());
    require(r.converged&&(r.events&PM_EVENT_YIELD)!=0,"Mroz P-M first plastic step");
    require(std::hypot(t[9],t[10])>1e-8,"Mroz inner Y surface must translate after yielding");
    std::vector<double> t2(c.size());auto r2=h.trial(-0.0005,0.020,t.data(),t2.data());
    require(r2.converged,"Mroz P-M second plastic step");
    require(std::abs(r2.moment)>=std::abs(r.moment)-1e-5,"Mroz hardening should not lose monotonic strength before U surface");
    const double hd=1e-7,hr=1e-8;auto eval=[&](double da,double th){std::vector<double> st(c.size());return h.trial(da,th,c.data(),st.data());};
    auto ap=eval(-0.0005+hd,0.010),am=eval(-0.0005-hd,0.010),rp=eval(-0.0005,0.010+hr),rm=eval(-0.0005,0.010-hr);
    require(ap.converged&&am.converged&&rp.converged&&rm.converged,"Mroz tangent perturbation convergence");
    const double k00=(ap.axial_force-am.axial_force)/(2*hd),k10=(ap.moment-am.moment)/(2*hd),k01=(rp.axial_force-rm.axial_force)/(2*hr),k11=(rp.moment-rm.moment)/(2*hr);
    require(std::abs(k00-r.tangent[0])<5e-4*std::max(1.0,std::abs(k00)),"Mroz KNN finite-difference parity");
    require(std::abs(k01-r.tangent[1])<5e-4*std::max(1.0,std::abs(k01)),"Mroz KNM finite-difference parity");
    require(std::abs(k10-r.tangent[2])<5e-4*std::max(1.0,std::abs(k10)),"Mroz KMN finite-difference parity");
    require(std::abs(k11-r.tangent[3])<5e-4*std::max(1.0,std::abs(k11)),"Mroz KMM finite-difference parity");
}

static void test_frame3d_true_pm_compound_consistent_jacobian(){
    Frame3DBuilder b;
    b.add_node(1,0,0,0);b.add_node(2,0,0,0);b.add_node(3,0,0,3,1,0,0);
    b.fix(1);b.fix_dof(2,Dof3D::UY);b.fix_dof(2,Dof3D::RX);b.fix_dof(2,Dof3D::RZ);
    b.fix_dof(3,Dof3D::UY);b.fix_dof(3,Dof3D::RX);b.fix_dof(3,Dof3D::RZ);
    b.add_elastic_frame(1,2,3,30000,12000,2.0,0.3,0.5,0.5,1,0,0,40.0);
    auto p=basic_pm_return_params();p.axial_stiffness=2.0e5;p.hinge.Ke=1.0e5;
    b.add_linear_pm_interaction_hinge(10,{{1,Dof3D::RY,-1.0},{2,Dof3D::RY,1.0}},
                                           {{1,Dof3D::UZ,-1.0},{2,Dof3D::UZ,1.0}},p);
    b.set_response(3,Dof3D::UX);b.set_story_nodes({3},Dof3D::UX);auto m=b.compile();
    require(m.has_state_dependent_global_tangent(),"true P-M hinge must use state-aware global tangent");
    require(!m.has_generalized_state_update(),"true P-M block must stay on exact/direct state tangent path");
    auto state=m.initial_nonlinear_state();std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0),v(u.size(),0.0);
    const int uz=m.reduced_dof(2,Dof3D::UZ),ry=m.reduced_dof(2,Dof3D::RY);
    require(uz>=0&&ry>=0,"true P-M reduced hinge DOFs");u[static_cast<std::size_t>(uz)]=-0.0001;u[static_cast<std::size_t>(ry)]=0.00015;
    v[static_cast<std::size_t>(uz)]=0.37;v[static_cast<std::size_t>(ry)]=-0.52;
    std::vector<double> f,t,trial;m.internal_force_and_tangent(u,state,f,t,trial);auto K=m.effective_state_tangent_matrix_with_state(u,t,state,0,0);auto kv=K.multiply(v);
    const auto ps=m.pm_interaction_snapshot(10,u,state);
    const double h=1e-7;auto up=u,um=u;for(std::size_t i=0;i<u.size();++i){up[i]+=h*v[i];um[i]-=h*v[i];}
    std::vector<double> fp,fm,tp,tm,sp,sm;m.internal_force_and_tangent(up,state,fp,tp,sp);m.internal_force_and_tangent(um,state,fm,tm,sm);
    std::vector<double> fd(kv.size());for(std::size_t i=0;i<fd.size();++i)fd[i]=(fp[i]-fm[i])/(2*h);
    require(rel_err(kv,fd)<2e-5,"true P-M assembled tangent must match global internal-force Jacobian");
    std::vector<double> ez(u.size(),0.0),er(u.size(),0.0);ez[static_cast<std::size_t>(uz)]=1;er[static_cast<std::size_t>(ry)]=1;
    auto cz=K.multiply(ez),cr=K.multiply(er);
        require(std::abs(cz[static_cast<std::size_t>(ry)])>1e-5&&std::abs(cr[static_cast<std::size_t>(uz)])>1e-5,"true P-M tangent must contain reciprocal axial-flexural coupling");
    require(std::abs(cz[static_cast<std::size_t>(ry)]-cr[static_cast<std::size_t>(uz)])<2e-5*std::max(1.0,std::abs(cz[static_cast<std::size_t>(ry)])),"true P-M global tangent reciprocity");
}

static void test_frame3d_axial_coupled_asce41_consistent_jacobian(){
    Frame3DBuilder b;
    b.add_node(1,0,0,0);
    b.add_node(2,0,0,3,1,1,1);
    b.fix(1);
    b.add_elastic_frame(1,1,2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,40.0);

    AxialCoupledASCE41Params cp;
    cp.hinge=basic_asce41_params();
    cp.axial_preload=40.0;
    cp.axial_stiffness=200.0;
    cp.axial_force_points={0.0,50.0,100.0};
    cp.moment_capacity_points={10.0,7.5,5.0};
    b.add_linear_axial_coupled_asce41_hinge(
        10,
        {{2,Dof3D::RY,1.0}},
        {{2,Dof3D::UZ,1.0}},
        cp);
    b.set_response(2,Dof3D::UX);
    b.set_story_nodes({2},Dof3D::UX);
    auto m=b.compile();
    require(m.has_state_dependent_global_tangent(),"coupled P-M hinge must request state-aware global tangent");
    require(!m.has_generalized_state_update(),"coupled P-M hinge must not use scalar/generalized Woodbury until nonsymmetric block support is enabled");

    auto state=m.initial_nonlinear_state();
    std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0),v(u.size(),0.0);
    const int ry=m.reduced_dof(2,Dof3D::RY),uz=m.reduced_dof(2,Dof3D::UZ);
    require(ry>=0&&uz>=0,"coupled P-M reduced DOFs");
    // Compression increases because axial_delta is negative.  Rotation is beyond yield
    // so moment depends on the axial-dependent capacity, exercising dM/dP.
    u[static_cast<std::size_t>(ry)]=0.015;
    u[static_cast<std::size_t>(uz)]=-0.05;
    v[static_cast<std::size_t>(ry)]=0.63;
    v[static_cast<std::size_t>(uz)]=-0.41;

    std::vector<double> f,tang,trial;
    m.internal_force_and_tangent(u,state,f,tang,trial);
    auto kt=m.effective_state_tangent_matrix_with_state(u,tang,state,0.0,0.0);
    auto kv=kt.multiply(v);

    const double h=2e-7;
    auto up=u,um=u;
    for(std::size_t i=0;i<u.size();++i){up[i]+=h*v[i];um[i]-=h*v[i];}
    std::vector<double> fp,fm,tp,tm,sp,sm;
    m.internal_force_and_tangent(up,state,fp,tp,sp);
    m.internal_force_and_tangent(um,state,fm,tm,sm);
    std::vector<double> fd(kv.size());
    for(std::size_t i=0;i<fd.size();++i)fd[i]=(fp[i]-fm[i])/(2*h);
    require(rel_err(kv,fd)<2e-6,"coupled P-M assembled Newton tangent must match internal-force Jacobian");

    // P -> M coupling is directional: axial deformation changes hinge moment, while
    // the intermediate model does not invent a reciprocal hinge axial force.
    std::vector<double> e_uz(u.size(),0.0),e_ry(u.size(),0.0);
    e_uz[static_cast<std::size_t>(uz)]=1.0;
    e_ry[static_cast<std::size_t>(ry)]=1.0;
    const auto col_uz=kt.multiply(e_uz),col_ry=kt.multiply(e_ry);
    require(std::abs(col_uz[static_cast<std::size_t>(ry)])>1e-6,"coupled P-M tangent must contain dM/d axial-deformation cross term");
    require(std::abs(col_ry[static_cast<std::size_t>(uz)])<1e-9,"intermediate P-M model must not invent reciprocal axial hinge force");
}

static void test_frame3d_updated_pdelta_state_tangent(){
    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,3,1,1,1);b.fix(1);
    b.add_elastic_frame(1,1,2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,20.0);
    b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);b.set_updated_pdelta(true);
    auto m=b.compile();require(m.has_state_dependent_global_tangent(),"updated PDelta mode flag");
    auto u=std::vector<double>(static_cast<std::size_t>(m.dof()),0.0);const int uz=m.reduced_dof(2,Dof3D::UZ),ux=m.reduced_dof(2,Dof3D::UX);
    require(uz>=0&&ux>=0,"updated PDelta reduced DOFs");u[static_cast<std::size_t>(uz)]=-0.002;
    std::vector<double> tang=m.initial_nonlinear_tangents();
    auto k0=m.effective_tangent_matrix(tang,0,0);auto k1=m.effective_state_tangent_matrix(u,tang,0,0);
    require(k1.diagonal()[static_cast<std::size_t>(ux)]<k0.diagonal()[static_cast<std::size_t>(ux)],"axial shortening should increase PDelta softening");
}


static void test_asce41_cyclic_path_continuity_and_persistent_E_loss(){
    auto p=basic_asce41_params();p.lambda_strength=0.25;p.lambda_unloading=0.25;
    ASCE41HingeMaterial h(p);std::vector<double> c(ASCE41HingeMaterial::kStateSize),t(c.size());h.initialize_state(c.data());
    // Fine cyclic protocol through yield, capping, reversal, and reloading.
    std::vector<double> targets{0.0,0.045,-0.035,0.060,-0.055,0.075,-0.070};
    double q=0.0,last_f=0.0;const double dq=1e-4;double max_jump=0.0;
    for(double target:targets){
        const double sg=target>=q?1.0:-1.0;
        while((sg>0&&q<target-0.5*dq)||(sg<0&&q>target+0.5*dq)){
            double qn=q+sg*dq;if((sg>0&&qn>target)||(sg<0&&qn<target))qn=target;
            auto e=h.trial(qn,c.data(),t.data());
            max_jump=std::max(max_jump,std::abs(e.force-last_f));
            require(std::isfinite(e.force)&&std::isfinite(e.tangent),"ASCE41 cyclic path finite response");
            c=t;q=qn;last_f=e.force;
        }
    }
    require(max_jump<1.0,"ASCE41 explicit branch state should avoid cyclic force jumps");

    // E-loss is persistent lateral loss, while F remains a later effective/gravity failure.
    h.initialize_state(c.data());const double uy=p.posFy/p.Ke;
    auto e=h.trial(uy+p.pos_b+1e-5,c.data(),t.data());require((e.events&ASCE41_EVENT_LATERAL_LOSS)!=0,"ASCE41 E reached");c=t;
    e=h.trial(0.0,c.data(),t.data());require((e.events&ASCE41_EVENT_LATERAL_LOSS)!=0&&std::abs(e.force)<1e-14,"ASCE41 E lateral loss should persist after reversal");c=t;
    e=h.trial(-(uy+p.neg_f+1e-4),c.data(),t.data());require((e.events&ASCE41_EVENT_FAILURE)!=0,"ASCE41 F failure after persistent E loss");
}


static void test_corotational3d_objectivity_and_small_response(){
    CorotationalFrame3DProperties p;
    p.xi=0;p.yi=0;p.zi=0;p.xj=3;p.yj=0;p.zj=0;
    p.E=30000;p.G=12000;p.A=2.0;p.J=0.3;p.Iy=0.5;p.Iz=0.7;p.reference={0,1,0};
    std::array<double,12> u{};
    const double th=0.65,c=std::cos(th),ss=std::sin(th);
    u[6]=3*c-3;u[7]=3*ss;u[5]=th;u[11]=th;
    auto q=corotational3d_basic_deformation(p,u);
    double qmax=0;for(double x:q)qmax=std::max(qmax,std::abs(x));
    require(qmax<2e-9,"corotational 3D basic deformation must be objective under rigid rotation");
    auto f=corotational3d_internal_force(p,u);double fmax=0;for(double x:f)fmax=std::max(fmax,std::abs(x));
    require(fmax<2e-3,"corotational 3D rigid rotation should have negligible internal force");

    u={};u[6]=0.003;
    auto r=corotational3d_response(p,u);
    const double expected=p.E*p.A/3.0*0.003;
    require(std::abs(r.force[6]-expected)/expected<2e-4,"corotational 3D axial response");

    std::array<double,12> z{};auto rz=corotational3d_response(p,z);
    auto klin=frame3d_global_stiffness(p.xi,p.yi,p.zi,p.xj,p.yj,p.zj,p.E,p.G,p.A,p.J,p.Iy,p.Iz,p.reference,0.0);
    double num=0,den=0;for(int i=0;i<144;++i){num=std::max(num,std::abs(rz.tangent[i]-klin[i]));den=std::max(den,std::abs(klin[i]));}
    require(num/std::max(1.0,den)<2e-3,"corotational initial tangent should recover linear 3D frame stiffness");

    // Arbitrary-axis rigid-body rotation plus translation. This is a stronger
    // objectivity check than a planar spin of an axis-aligned member.
    CorotationalFrame3DProperties pa=p;
    pa.xi=0.3;pa.yi=-0.2;pa.zi=0.1;pa.xj=2.2;pa.yj=1.1;pa.zj=3.7;pa.reference={0.2,1.0,-0.15};
    std::array<double,3> axis{0.35,-0.42,0.837};
    const double an=std::sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
    for(double& x:axis)x/=an;
    const double ang=0.82,ca=std::cos(ang),sa=std::sin(ang),omc=1.0-ca;
    auto rotate=[&](std::array<double,3> x){
        const double d=axis[0]*x[0]+axis[1]*x[1]+axis[2]*x[2];
        std::array<double,3> cr{axis[1]*x[2]-axis[2]*x[1],axis[2]*x[0]-axis[0]*x[2],axis[0]*x[1]-axis[1]*x[0]};
        return std::array<double,3>{ca*x[0]+sa*cr[0]+omc*d*axis[0],ca*x[1]+sa*cr[1]+omc*d*axis[1],ca*x[2]+sa*cr[2]+omc*d*axis[2]};
    };
    const std::array<double,3> tr{0.7,-1.1,0.45};
    auto xi=rotate({pa.xi,pa.yi,pa.zi}),xj=rotate({pa.xj,pa.yj,pa.zj});
    std::array<double,12> ua{};
    ua[0]=xi[0]+tr[0]-pa.xi;ua[1]=xi[1]+tr[1]-pa.yi;ua[2]=xi[2]+tr[2]-pa.zi;
    ua[6]=xj[0]+tr[0]-pa.xj;ua[7]=xj[1]+tr[1]-pa.yj;ua[8]=xj[2]+tr[2]-pa.zj;
    ua[3]=ua[9]=ang*axis[0];ua[4]=ua[10]=ang*axis[1];ua[5]=ua[11]=ang*axis[2];
    auto qa=corotational3d_basic_deformation(pa,ua);qmax=0.0;for(double x:qa)qmax=std::max(qmax,std::abs(x));
    require(qmax<2e-8,"corotational 3D must be objective under arbitrary rigid-body motion");
}

static CompiledFrame3D corotational_cantilever_column(int elements,double axial_compression){
    Frame3DBuilder b;const double L=4.0;
    for(int i=0;i<=elements;++i)b.add_node(i+1,0,0,L*i/elements);
    b.fix(1);
    for(int i=1;i<=elements;++i){const int n=i+1;b.fix_dof(n,Dof3D::UY);b.fix_dof(n,Dof3D::RX);b.fix_dof(n,Dof3D::RZ);}
    for(int i=0;i<elements;++i)b.add_elastic_frame(i+1,i+1,i+2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,axial_compression);
    b.set_response(elements+1,Dof3D::UX);b.set_story_nodes({elements+1},Dof3D::UX);b.set_corotational(true);
    return b.compile();
}

static void test_corotational3d_euler_stability_convergence(){
    const double E=30000.0,I=0.5,L=4.0;
    const double pcr=M_PI*M_PI*E*I/(4.0*L*L);
    auto stable=corotational_cantilever_column(8,0.98*pcr);
    auto unstable=corotational_cantilever_column(8,1.02*pcr);
    require(assess_positive_definiteness(stable.K_initial(),1e-9,256).positive_definite,
            "8-element corotational column should remain stable below Euler load");
    require(!assess_positive_definiteness(unstable.K_initial(),1e-9,256).positive_definite,
            "8-element corotational column should lose stability above Euler load");
}

static void test_updated_pdelta_generalized_woodbury_exactness(){
    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,3,1,1,1);b.fix(1);
    b.add_elastic_frame(1,1,2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,35.0);
    b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);b.set_updated_pdelta(true);auto m=b.compile();
    require(m.has_generalized_state_update(),"updated PDelta should expose generalized low-rank update");
    std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0);const int ux=m.reduced_dof(2,Dof3D::UX),uy=m.reduced_dof(2,Dof3D::UY),uz=m.reduced_dof(2,Dof3D::UZ);
    u[static_cast<std::size_t>(ux)]=0.05;u[static_cast<std::size_t>(uy)]=-0.025;u[static_cast<std::size_t>(uz)]=-0.004;
    auto tang=m.initial_nonlinear_tangents();auto A=m.effective_initial_matrix(0.0,0.0);auto K=m.effective_state_tangent_matrix(u,tang,0.0,0.0);auto C=m.generalized_state_update_coefficients(u,tang);
    GeneralizedWoodburySolver w(A,m.generalized_state_update_basis());std::vector<double> rhs(static_cast<std::size_t>(m.dof()));for(int i=0;i<m.dof();++i)rhs[static_cast<std::size_t>(i)]=0.2+0.13*i;
    auto x=w.solve(rhs,C),xr=superlu_solve_once(K,rhs);require(rel_err(xr,x)<2e-8,"localized geometric generalized Woodbury must match direct updated-PDelta tangent");
}

static void test_frame3d_updated_pdelta_consistent_jacobian(){
    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,3,1,1,1);b.fix(1);
    b.add_elastic_frame(1,1,2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,20.0);
    b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);b.set_updated_pdelta(true);
    auto m=b.compile();auto state=m.initial_nonlinear_state();
    std::vector<double> u(static_cast<std::size_t>(m.dof()),0.0),v(u.size(),0.0);
    const int ux=m.reduced_dof(2,Dof3D::UX),uz=m.reduced_dof(2,Dof3D::UZ);
    require(ux>=0&&uz>=0,"updated PDelta Jacobian reduced DOFs");
    u[static_cast<std::size_t>(ux)]=0.04;u[static_cast<std::size_t>(uz)]=-0.003;
    v[static_cast<std::size_t>(ux)]=0.7;v[static_cast<std::size_t>(uz)]=-0.4;
    std::vector<double> f,tang,trial;m.internal_force_and_tangent(u,state,f,tang,trial);
    auto kt=m.effective_state_tangent_matrix(u,tang,0.0,0.0);auto kv=kt.multiply(v);
    const double h=1e-7;auto up=u,um=u;for(std::size_t i=0;i<u.size();++i){up[i]+=h*v[i];um[i]-=h*v[i];}
    std::vector<double> fp,fm,tp,tm,sp,sm;m.internal_force_and_tangent(up,state,fp,tp,sp);m.internal_force_and_tangent(um,state,fm,tm,sm);
    std::vector<double> fd(kv.size());for(std::size_t i=0;i<fd.size();++i)fd[i]=(fp[i]-fm[i])/(2*h);
    require(rel_err(kv,fd)<2e-7,"updated PDelta Newton tangent must match internal-force Jacobian");
}
static void test_collapse_driver_distinguishes_physical_termination(){
    auto model=make_space_frame_nrha_model(3);auto gm=synthetic_ground_motion(180,0.005,3.0);
    RobustNewmarkOptions o;o.max_iterations=35;o.max_subdivisions=2;o.collapse.max_story_drift_ratio=1e-5;o.return_numerical_failure=true;
    auto r=run_newmark_robust(model,gm,0.005,LinearStrategy::Woodbury,o);
    require(r.termination==AnalysisTermination::PhysicalCollapse,"collapse driver should report physical collapse separately");
    require(r.collapse_mechanism==CollapseMechanism::DriftLimit,"collapse driver should classify drift-triggered collapse");
    require(r.termination_step<gm.size(),"physical collapse termination step");
    require(r.max_story_drift_ratio>=o.collapse.max_story_drift_ratio,"collapse drift ratio recorded");
}


static void test_robust_energy_ledger_and_initial_instability(){
    auto model=make_space_frame_nrha_model(3);auto gm=synthetic_ground_motion(220,0.005,0.5);
    RobustNewmarkOptions o;o.max_iterations=35;o.max_subdivisions=2;o.return_numerical_failure=true;
    auto r=run_newmark_robust(model,gm,0.005,LinearStrategy::Woodbury,o);
    require(r.termination==AnalysisTermination::Completed,"energy-ledger reference analysis should complete");
    require(std::isfinite(r.stats.energy_balance_relative_error)&&r.stats.energy_balance_relative_error<0.02,"robust energy ledger should close within reference tolerance");

    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,3,1,1,1);b.fix(1);
    b.add_elastic_frame(1,1,2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,25000.0);
    b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);auto unstable=b.compile();
    RobustNewmarkOptions uo;uo.collapse.tangent_check_interval=1;uo.collapse.min_tangent_ratio=0.02;uo.return_numerical_failure=true;
    auto ur=run_newmark_robust(unstable,std::vector<double>(2,0.0),0.005,LinearStrategy::Woodbury,uo);
    require(ur.termination==AnalysisTermination::InitialInstability,"initially indefinite model must be distinguished from earthquake-induced collapse");
    require(ur.termination_time==0.0,"initial instability should terminate at time zero");
}


static void test_ida_suite_parallel_and_censoring(){
    auto model=make_space_frame_nrha_model(2);GroundMotionRecord a{"A",synthetic_ground_motion(100,0.005,0.35),0.005},b{"B",synthetic_ground_motion(100,0.005,0.45),0.005};
    IDAOptions o;o.scale_factors={0.25,0.5};o.workers=2;o.collapse_refinement_steps=0;o.analysis.max_iterations=35;o.analysis.max_subdivisions=1;o.analysis.return_numerical_failure=true;
    auto r=run_ida_suite(model,{a,b},o);require(r.records.size()==2,"IDA should return one summary per record");require(r.runs.size()>=2&&r.workers==2,"IDA parallel suite run count/workers");
    for(const auto& rr:r.runs){require(rr.scale_factor>0&&std::isfinite(rr.pga),"IDA run scale/PGA");require(rr.record_index<2,"IDA record index");}
    require(r.physical_collapses+r.right_censored_records<=2,"IDA summary accounting");
}

static void test_sparse_initial_stability_certification(){
    auto spd=SparseMatrixCSC::from_triplets(3,3,{
        {0,0,4.0},{1,1,3.0},{2,2,2.0},{0,1,-1.0},{1,0,-1.0},{1,2,-0.5},{2,1,-0.5}});
    auto a=assess_positive_definiteness(spd,1e-12,16);
    require(a.status==InitialStabilityStatus::PositiveDefinite&&a.positive_definite,"sparse SPD certification should accept positive definite matrix");
    require(a.minimum_pivot_ratio>0.0,"sparse SPD certification should report positive pivot margin");
    require(a.dense_check_performed&&a.dense_negative_eigenvalues==0&&a.dense_minimum_eigenvalue>0.0,"dense validation should agree with sparse SPD certification");

    auto indef=SparseMatrixCSC::from_triplets(3,3,{
        {0,0,2.0},{1,1,-10.0},{2,2,1.0}});
    a=assess_positive_definiteness(indef,1e-12,16);
    require(a.status==InitialStabilityStatus::NotPositiveDefinite&&!a.positive_definite,"sparse SPD certification must reject indefinite matrix even when negative mode is not nearest zero");
    require(a.dense_check_performed&&a.dense_negative_eigenvalues==1&&a.dense_minimum_eigenvalue<0.0,"dense inertia cross-check should identify negative mode");

    auto singular=SparseMatrixCSC::from_triplets(2,2,{{0,0,1.0},{1,1,0.0}});
    a=assess_positive_definiteness(singular,1e-12,16);
    require(a.status==InitialStabilityStatus::NotPositiveDefinite,"sparse SPD certification must reject singular tangent");
}

static void test_tangent_stability_estimator(){
    auto model=make_space_frame_nrha_model(2);auto u=std::vector<double>(static_cast<std::size_t>(model.dof()),0.0);auto s=model.initial_nonlinear_state();
    auto e=estimate_tangent_stability(model,u,s,5);require(e.factorization_ok,"initial tangent stability factorization");require(e.rayleigh_eigenvalue>0.0,"initial tangent should be stable positive");
}


static void test_fsc_shear_spring_law(){
    quake::FSCShearSpringParams p;
    p.Ke=100.0;
    p.post_failure_stiffness=-5.0;
    p.residual_strength_ratio=0.20;
    p.nominal_shear_limit_kip=-1.0;
    p.cyclic_strength_coefficient=0.05;
    quake::FSCShearSpringLaw law(p);
    std::array<double,quake::FSCShearSpringLaw::kStateSize> c{},s{},sp{},sm{},srev1{},srev2{};
    law.initialize_state(c.data());

    const double lim=quake::FSCShearSpringLaw::plastic_rotation_limit(p,25.0,8.0);
    auto e=law.trial(0.04,0.5*lim,25.0,c.data(),s.data());
    require(!e.initiated,"FSC shear spring initiated below remote limit");
    require(std::abs(e.material.force-4.0)<1e-12 && std::abs(e.material.tangent-100.0)<1e-12,
            "FSC pre-failure elastic response");

    auto init=law.trial(0.08,1.05*lim,25.0,c.data(),s.data());
    require(init.initiated && init.material.diagnostics.lateral_resistance_lost,
            "FSC shear initiation/E diagnostic");
    c=s;

    const double q=0.10,h=1e-7;
    auto post=law.trial(q,1.05*lim,25.0,c.data(),s.data());
    auto pp=law.trial(q+h,1.05*lim,25.0,c.data(),sp.data());
    auto pm=law.trial(q-h,1.05*lim,25.0,c.data(),sm.data());
    const double fd=(pp.material.force-pm.material.force)/(2.0*h);
    require(std::abs(fd-post.material.tangent)<1e-6,
            "FSC post-failure tangent must match finite differences");
    require(post.material.tangent<0.0,"FSC post-failure softening branch expected");

    c=s;
    auto rev1=law.trial(-0.10,1.05*lim,25.0,c.data(),srev1.data());
    auto rev2=law.trial(-0.10,1.05*lim,25.0,c.data(),srev2.data());
    require(std::abs(rev1.retained_strength_ratio-rev2.retained_strength_ratio)<1e-14,
            "FSC repeated Newton trials must not double-count cyclic damage");
    require(rev1.retained_strength_ratio<post.retained_strength_ratio,
            "FSC opposite-sign excursion should reduce retained strength");
}

static void test_fsc_shear_damage_diagnostic(){
    quake::FSCShearDamageParams p;
    p.confined_concrete_area_in2=18.0;
    p.nominal_shear_limit_kip=-1.0;
    quake::FSCShearDamageDiagnostic d(p);
    const double lim=quake::FSCShearDamageDiagnostic::plastic_rotation_limit_rad(p,25.0,8.0);
    if(!(lim>0.0 && lim<0.032)) throw std::runtime_error("FSC rotation limit out of range");
    d.update(0,0.0,0.5*lim,0.0,25.0,5.0);
    if(d.state().initiated) throw std::runtime_error("FSC diagnostic initiated too early");
    const double lim6=quake::FSCShearDamageDiagnostic::plastic_rotation_limit_rad(p,25.0,6.0);
    d.update(1,0.01,1.01*lim6,0.0,25.0,6.0);
    if(!d.state().initiated || d.state().initiation_cause!=1) throw std::runtime_error("FSC rotation initiation failed");
    const double before=d.state().retained_strength_ratio;
    d.update(2,0.02,1.01*lim6,0.0,25.0,-6.0);
    if(d.state().opposite_sign_crossings!=1 || !(d.state().retained_strength_ratio<before))
        throw std::runtime_error("FSC cyclic strength decrement failed");
    quake::FSCShearDamageParams p2=p; p2.confined_concrete_area_in2=24.0;
    if(!(quake::FSCShearDamageDiagnostic::cyclic_strength_coefficient(p2) <
         quake::FSCShearDamageDiagnostic::cyclic_strength_coefficient(p)))
        throw std::runtime_error("FSC confinement sensitivity has wrong sign");
}


static CompiledFrame2D make_comparison_frame(const ASCE41HingeParams& hp){
    Frame2DBuilder b;
    b.add_node(1,0,0);b.add_node(2,0,120,1.0,0,0);b.fix(1);
    b.add_elastic_frame(901,1,2,3000.0,36.0,108.0,20.0);
    b.add_asce41_hinge(902,1,2,hp);
    b.set_response_node(2);b.set_story_nodes({2});
    return b.compile();
}

static void test_prepared_robust_reuse_and_material_invalidation(){
    auto hp=basic_asce41_params();hp.Ke=5000.0;hp.posFy=hp.negFy=50.0;
    auto model=make_comparison_frame(hp);
    RobustNewmarkOptions opt;opt.max_iterations=25;opt.max_subdivisions=3;opt.return_numerical_failure=true;
    const auto gm=synthetic_ground_motion(40,0.01,1.0);
    PreparedRobustNewmark prepared(model,0.01,LinearStrategy::SamePatternRefactorization,3);
    auto a=prepared.run(gm,opt);
    require(!prepared.last_run_reused_preparation()&&prepared.preparation_count()==1,
            "first prepared robust run must create preparation");
    auto b=prepared.run(gm,opt);
    require(prepared.last_run_reused_preparation()&&prepared.preparation_count()==1,
            "unchanged prepared robust model should reuse preparation");
    require(rel_err(a.roof_history,b.roof_history)<1e-13,"prepared robust reuse changes response");

    auto strength_only=hp;strength_only.posFy*=0.8;strength_only.negFy*=0.8;
    model.replace_nonlinear_material(902,NonlinearMaterial{ASCE41HingeMaterial(strength_only)});
    (void)prepared.run(gm,opt);
    require(prepared.last_run_reused_preparation()&&prepared.preparation_count()==1,
            "strength-only material change with unchanged Ke should reuse effective-initial preparation");

    auto stiffness_change=strength_only;stiffness_change.Ke*=1.25;
    model.replace_nonlinear_material(902,NonlinearMaterial{ASCE41HingeMaterial(stiffness_change)});
    (void)prepared.run(gm,opt);
    require(!prepared.last_run_reused_preparation()&&prepared.preparation_count()==2,
            "Ke change must invalidate prepared robust solver contexts");
}

static void test_model_comparison_harness_is_aligned_and_order_invariant(){
    auto hp=basic_asce41_params();hp.Ke=5000.0;hp.posFy=hp.negFy=50.0;
    auto model=make_comparison_frame(hp);
    auto same=hp;
    auto changed=hp;changed.Ke*=1.15;changed.posFy*=0.9;changed.negFy*=0.9;
    std::vector<ModelComparisonVariant> variants{
        {"research","captured baseline",{}},
        {"same","identical field",{{902,NonlinearMaterial{ASCE41HingeMaterial(same)}}}},
        {"changed","stiffness/strength perturbation",{{902,NonlinearMaterial{ASCE41HingeMaterial(changed)}}}}
    };
    RCNativeFrameAnalysisSettings settings;settings.dt=0.01;settings.strategy=LinearStrategy::SamePatternRefactorization;
    settings.ground_accel=synthetic_ground_motion(50,settings.dt,1.2);settings.newmark_options.max_iterations=25;settings.newmark_options.max_subdivisions=3;settings.newmark_options.return_numerical_failure=true;
    ModelComparisonOptions options;options.tangent_modes_at_story_peaks=true;options.tangent_mode_count=1;options.column_force_bindings={{"Ccmp",{901},{902}}};
    auto r=RCFrameModelComparison::run_frame2d(model,variants,{{"base_hinge",902}},settings,options);
    require(r.variants.size()==3,"comparison harness variant count");
    require(r.variants[0].story_history.size()==settings.ground_accel.size(),"comparison story histories must align to output steps");
    require(r.variants[0].components.size()==1&&!r.variants[0].components[0].history.empty(),"comparison component history recorder");
    require(r.variants[0].columns.size()==1&&r.variants[0].columns[0].history.size()==settings.ground_accel.size(),"comparison native column P/V history recorder");
    require(!r.variants[0].tangent_snapshots.empty(),"comparison must retain tangent modes at story peak");
    require(r.variants[0].tangent_snapshots[0].story_tangent_available&&!r.variants[0].tangent_snapshots[0].interstory_tangent_diagonal.empty(),"comparison story tangent condensation");
    require(rel_err(r.variants[0].analysis.roof_history,r.variants[1].analysis.roof_history)<1e-13,
            "identical material variant should reproduce baseline response");
    require(r.variants[1].reused_prepared_solver,"identical-Ke comparison variant should reuse prepared solver");
    require(r.prepared_solver_preparations==2,"changed Ke comparison variant should trigger one safe preparation refresh");

    auto model2=make_comparison_frame(hp);
    std::reverse(variants.begin(),variants.end());
    auto q=RCFrameModelComparison::run_frame2d(model2,variants,{{"base_hinge",902}},settings,options);
    auto peak_by_name=[](const ModelComparisonResult& x){std::map<std::string,double> m;for(const auto& v:x.variants)m[v.name]=v.peak_abs_roof_response;return m;};
    const auto p1=peak_by_name(r),p2=peak_by_name(q);
    require(p1.size()==p2.size(),"comparison order-invariance name coverage");
    for(const auto& [name,val]:p1)require(std::abs(val-p2.at(name))<1e-11,"comparison physical result depends on variant execution order");
}


static void test_rc_model_comparison_variant_requires_exact_binding_coverage(){
    auto hp=basic_asce41_params();
    RCBuildingColumnModel a;a.component_id="A";a.model.hinge=hp;
    RCBuildingColumnModel extra;extra.component_id="EXTRA";extra.model.hinge=hp;
    std::vector<RCColumnElementBinding> bindings{{"A",{901},{902}}};
    const auto ok=make_rc_column_model_comparison_variant("ok","coverage test",bindings,{a});
    require(ok.replacements.size()==1&&ok.replacements[0].first==902,"RC comparison exact binding adapter basic mapping");
    bool threw=false;
    try{(void)make_rc_column_model_comparison_variant("bad","coverage test",bindings,{a,extra});}
    catch(const std::invalid_argument&){threw=true;}
    require(threw,"RC comparison adapter must reject unbound extra model");
}

static void test_model_comparison_harness_frame3d_smoke(){
    auto hp=basic_asce41_params();hp.Ke=5000.0;hp.posFy=hp.negFy=50.0;
    Frame3DBuilder b;b.add_node(1,0,0,0);b.add_node(2,0,0,120,1.0,0,0);b.fix(1);
    b.add_elastic_frame(911,1,2,3000.0,1200.0,36.0,20.0,108.0,108.0,0,1,0,20.0);
    b.add_asce41_hinge(912,1,Dof3D::RY,2,Dof3D::RY,hp);
    b.set_ground_direction(Dof3D::UX);b.set_response(2,Dof3D::UX);b.set_story_nodes({2},Dof3D::UX);
    auto model=b.compile();
    std::vector<ModelComparisonVariant> variants{{"research3d","baseline",{}},{"same3d","same",{{912,NonlinearMaterial{ASCE41HingeMaterial(hp)}}}}};
    RCNativeFrameAnalysisSettings settings;settings.dt=0.01;settings.strategy=LinearStrategy::SamePatternRefactorization;
    settings.ground_accel=synthetic_ground_motion(30,settings.dt,0.8);settings.newmark_options.max_iterations=25;settings.newmark_options.max_subdivisions=2;settings.newmark_options.return_numerical_failure=true;
    ModelComparisonOptions options;options.tangent_mode_count=1;options.column_force_bindings={{"C3cmp",{911},{912}}};
    auto r=RCFrameModelComparison::run_frame3d(model,variants,{{"hinge3d",912}},settings,options);
    require(r.variants.size()==2&&r.variants[0].story_history.size()==settings.ground_accel.size(),"3D comparison aligned story history");
    require(rel_err(r.variants[0].analysis.roof_history,r.variants[1].analysis.roof_history)<1e-13,"3D comparison identical variant parity");
    require(!r.variants[0].tangent_snapshots.empty(),"3D comparison tangent modal snapshot");
    require(r.variants[0].columns.size()==1&&!r.variants[0].columns[0].history.empty(),"3D comparison native axial/shear history");
}

static void test_gravity_state_transfer_remains_in_equilibrium(){
    Frame2DBuilder b;b.add_node(1,0,0);b.add_node(2,0,10,1,1,0);b.fix(1);
    b.add_elastic_frame(1,1,2,1000,10,10);b.set_response_node(2);b.set_story_nodes({2});
    auto model=b.compile();std::vector<double> load(static_cast<std::size_t>(model.dof()),0.0);
    const int uy=model.reduced_dof(2,Dof2D::UY);load[static_cast<std::size_t>(uy)]=-10.0;
    auto gravity=solve_static_load(model,load,10,1e-12,10);
    require(gravity.converged,"gravity solve convergence");
    require(std::abs(gravity.displacement[static_cast<std::size_t>(uy)]+.01)<1e-12,"gravity axial displacement");
    RobustNewmarkOptions o;o.constant_load=load;o.initial_displacement=gravity.displacement;
    o.initial_committed_state=gravity.committed_state;o.collapse.check_initial_stability=false;
    o.line_search=false;o.max_subdivisions=0;
    auto a=run_newmark_robust(model,std::vector<double>(10,0.0),.01,LinearStrategy::SamePatternRefactorization,o);
    require(a.termination==AnalysisTermination::Completed,"gravity-state transient completion");
    require(std::abs(a.final_displacement[static_cast<std::size_t>(uy)]+.01)<1e-11,"gravity equilibrium retained in transient");
}

int main(){
    try{
        test_sparse_solve();
        test_woodbury_equivalence();
        test_generalized_woodbury_coupled_update();
        test_nrha_equivalence();
        test_robust_subdivision_recovers_failed_step();
        test_frame2d_cantilever();
        test_frame2d_constraints_and_hinge_basis();
        test_frame2d_nrha_equivalence();
        test_frame3d_cantilever();
        test_frame3d_vector_hinge_operator();
        test_frame3d_rigid_diaphragm_mpc();
        test_frame3d_nrha_equivalence();
        test_modal_analysis_singular_mass();
        test_fixed_modal_damping_force_side_wrapper();
        test_guyan_static_condensation();
        test_craig_bampton_modal_enrichment();
        test_parallel_suite_equivalence();
        test_imk_monotonic_backbone();
        test_imk_reversal_and_deterioration();
        test_frame3d_mixed_material_state_bank();
        test_story_block_schur_exactness();
        test_solver_calibration_smoke();
        test_prepared_adaptive_equivalence();
        test_fsc_shear_damage_diagnostic();
        test_fsc_shear_spring_law();
        test_asce41_physical_hardening_stiffness();
        test_asce41_degrading_backbone_and_limits();
        test_asce41_straight_ce_and_explicit_C_strength();
        test_rc_column_asce41_provider_separates_code_backbone_from_hysteresis();
        test_rc_column_section_rules_boundary_and_audit();
        test_rc_column_axial_iteration_converges_and_regenerates();
        test_rc_column_two_pass_convenience();
        test_rc_column_axial_iteration_under_relaxation();
        test_rc_building_iteration_synchronous_and_order_invariant();
        test_rc_building_iteration_converges_on_global_worst_column();
        test_rc_building_two_pass_and_failure_semantics();
        test_native_element_force_recovery_2d_and_3d();
        test_native_frame2d_building_iteration_uses_asce_material_bank();
        test_native_frame3d_building_iteration_records_committed_demands();
        test_material_field_replacement_preserves_compiled_topology();
        test_frame3d_material_field_matches_fresh_compile();
        test_compiled_topology_asce_iteration_reuses_one_frame();
        test_compiled_topology_frame3d_two_pass();
        test_corotational_section_force_recovery_is_objective();
        test_prepared_robust_reuse_and_material_invalidation();
        test_model_comparison_harness_is_aligned_and_order_invariant();
        test_rc_model_comparison_variant_requires_exact_binding_coverage();
        test_model_comparison_harness_frame3d_smoke();
        test_gravity_state_transfer_remains_in_equilibrium();
        test_asce41_cyclic_path_continuity_and_persistent_E_loss();
        test_pm_interaction_return_mapping_and_consistent_tangent();
        test_perform_concrete_pm_surface_return_mapping();
        test_perform_mroz_two_surface_translation();
        test_frame3d_true_pm_compound_consistent_jacobian();
        test_frame3d_axial_coupled_asce41_consistent_jacobian();
        test_frame3d_updated_pdelta_state_tangent();
        test_corotational3d_objectivity_and_small_response();
        test_corotational3d_euler_stability_convergence();
        test_frame3d_updated_pdelta_consistent_jacobian();
        test_updated_pdelta_generalized_woodbury_exactness();
        test_collapse_driver_distinguishes_physical_termination();
        test_robust_energy_ledger_and_initial_instability();
        test_ida_suite_parallel_and_censoring();
        test_sparse_initial_stability_certification();
        test_tangent_stability_estimator();
        std::cout<<"All quake-core tests passed.\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"TEST FAILURE: "<<e.what()<<"\n";return 1;}
}
