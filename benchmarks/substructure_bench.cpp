#include "quake/substructure.hpp"
#include "quake/frame3d.hpp"
#include "quake/newmark.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;
static int nid(int s,int ix,int iy){return 10000*s+100*ix+iy;}
static CompiledFrame3D make_model(int stories){Frame3DBuilder b;constexpr double h=3.5,L=6,E=25000,G=10000;int eid=1,sid=1,did=500000;for(int s=0;s<=stories;++s)for(int ix=0;ix<2;++ix)for(int iy=0;iy<2;++iy){double m=s?0.25:0;b.add_node(nid(s,ix,iy),ix*L,iy*L,s*h,m,m,.05*m);if(!s)b.fix(nid(s,ix,iy));}for(int s=1;s<=stories;++s){for(int ix=0;ix<2;++ix)for(int iy=0;iy<2;++iy)b.add_elastic_frame(eid++,nid(s-1,ix,iy),nid(s,ix,iy),E,G,2,.3,.55,.55,1,0,0,8*(stories-s+1));for(int iy=0;iy<2;++iy){int jl=nid(s,0,iy),jr=nid(s,1,iy),dl=did++,dr=did++;b.add_node(dl,0,iy*L,s*h);b.add_node(dr,L,iy*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RX,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,.2,.35,.35,0,1,0);b.add_bilinear_spring(sid++,jl,Dof3D::RY,dl,Dof3D::RY,50000,8,.5);b.add_bilinear_spring(sid++,jr,Dof3D::RY,dr,Dof3D::RY,50000,8,.5);}for(int ix=0;ix<2;++ix){int jl=nid(s,ix,0),jr=nid(s,ix,1),dl=did++,dr=did++;b.add_node(dl,ix*L,0,s*h);b.add_node(dr,ix*L,L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RY,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,.2,.35,.35,1,0,0);b.add_bilinear_spring(sid++,jl,Dof3D::RX,dl,Dof3D::RX,50000,8,.5);b.add_bilinear_spring(sid++,jr,Dof3D::RX,dr,Dof3D::RX,50000,8,.5);}}b.set_rayleigh(.02,.0005);b.set_ground_direction(Dof3D::UX);b.set_response(nid(stories,0,0),Dof3D::UX);return b.compile();}
template<class F>double tm(F&&f){auto a=std::chrono::steady_clock::now();f();return std::chrono::duration<double>(std::chrono::steady_clock::now()-a).count();}
int main(int argc,char**argv){
    int stories=argc>1?std::stoi(argv[1]):20,bs=argc>2?std::stoi(argv[2]):4,reps=argc>3?std::stoi(argv[3]):100,rank=argc>4?std::stoi(argv[4]):40;
    auto m=make_model(stories);double dt=.005,a0=1.0/(.25*dt*dt),a1=.5/(.25*dt);auto A=m.effective_initial_matrix(a0,a1);auto p=make_story_block_partition(m,bs);
    SuperLUFactor full(A);BlockSchurFactor block(A,p);std::vector<double>rhs(static_cast<std::size_t>(m.dof()));for(int i=0;i<m.dof();++i)rhs[i]=std::sin(.013*(i+1));
    auto xf=full.solve(rhs),xb=block.solve(rhs);double err=0,den=0;for(std::size_t i=0;i<xf.size();++i){err=std::max(err,std::abs(xf[i]-xb[i]));den=std::max(den,std::abs(xf[i]));}
    volatile double sink=0;double tf=tm([&]{for(int i=0;i<reps;++i){auto x=full.solve(rhs);sink+=x[i%x.size()];}});double tb=tm([&]{for(int i=0;i<reps;++i){auto x=block.solve(rhs);sink+=x[(i+3)%x.size()];}});(void)sink;auto st=block.stats();
    rank=std::max(0,std::min(rank,m.nonlinear_count()));std::vector<double>dk(static_cast<std::size_t>(m.nonlinear_count()),0.0);for(int j=0;j<rank;++j)dk[static_cast<std::size_t>(j)]=-0.8*m.initial_nonlinear_tangents()[static_cast<std::size_t>(j)];
    LazyLowRankWoodburySolver mono(A,m.nonlinear_basis());SubstructuredLazyWoodburySolver sub(A,m.nonlinear_basis(),make_story_block_partition(m,bs));
    auto xm=mono.solve(rhs,dk),xs=sub.solve(rhs,dk);double ew=0,dw=0;for(std::size_t i=0;i<xm.size();++i){ew=std::max(ew,std::abs(xm[i]-xs[i]));dw=std::max(dw,std::abs(xm[i]));}
    sink=0;double twm=tm([&]{for(int i=0;i<reps;++i){auto x=mono.solve(rhs,dk);sink+=x[i%x.size()];}});double tws=tm([&]{for(int i=0;i<reps;++i){auto x=sub.solve(rhs,dk);sink+=x[(i+5)%x.size()];}});(void)sink;
    std::cout<<std::fixed<<std::setprecision(6)<<"dof="<<m.dof()<<" block_size="<<bs<<" interface="<<st.interface_dof<<" interior="<<st.interior_dof<<" blocks="<<st.blocks<<" schur_nnz="<<st.schur_nnz<<" setup="<<st.setup_seconds<<"\n"<<"relative_solve_error="<<(den?err/den:err)<<" full_solve_seconds="<<tf<<" block_solve_seconds="<<tb<<" block_speedup="<<tf/tb<<"x\n"<<"woodbury_rank="<<rank<<" substructured_relative_error="<<(dw?ew/dw:ew)<<" monolithic_wood_seconds="<<twm<<" substructured_wood_seconds="<<tws<<" substructured_wood_speedup="<<twm/tws<<"x\n";
}
