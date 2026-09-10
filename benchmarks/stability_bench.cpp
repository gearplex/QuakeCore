#include "quake/frame3d.hpp"
#include "quake/stability.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;
static int nid(int s,int ix,int iy){return 10000*s+100*ix+iy;}
static CompiledFrame3D make_model(int stories,int nx=2,int ny=2){
    Frame3DBuilder b;constexpr double h=3.5,L=6.0,E=25000.0,G=10000.0;int eid=1,sid=1,did=900000;const int nc=(nx+1)*(ny+1);
    for(int s=0;s<=stories;++s)for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy){double m=s?1.0/nc:0.0;b.add_node(nid(s,ix,iy),ix*L,iy*L,s*h,m,m,.05*m);if(s==0)b.fix(nid(s,ix,iy));}
    for(int s=1;s<=stories;++s){
        const double p=2.0*(stories-s+1); // deliberately stable reference preload
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy)b.add_elastic_frame(eid++,nid(s-1,ix,iy),nid(s,ix,iy),E,G,2,.3,.55,.55,1,0,0,p);
        for(int ix=0;ix<nx;++ix)for(int iy=0;iy<=ny;++iy){int jl=nid(s,ix,iy),jr=nid(s,ix+1,iy),dl=did++,dr=did++;b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,(ix+1)*L,iy*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RX,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,.2,.35,.35,0,1,0);b.add_bilinear_spring(sid++,jl,Dof3D::RY,dl,Dof3D::RY,50000,8,.5);b.add_bilinear_spring(sid++,jr,Dof3D::RY,dr,Dof3D::RY,50000,8,.5);}
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<ny;++iy){int jl=nid(s,ix,iy),jr=nid(s,ix,iy+1),dl=did++,dr=did++;b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,ix*L,(iy+1)*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RY,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,.2,.35,.35,1,0,0);b.add_bilinear_spring(sid++,jl,Dof3D::RX,dl,Dof3D::RX,50000,8,.5);b.add_bilinear_spring(sid++,jr,Dof3D::RX,dr,Dof3D::RX,50000,8,.5);}
    }
    b.set_response(nid(stories,0,0),Dof3D::UX);std::vector<int> sn;for(int s=1;s<=stories;++s)sn.push_back(nid(s,0,0));b.set_story_nodes(std::move(sn),Dof3D::UX);return b.compile();
}
int main(int argc,char**argv){std::vector<int> stories{10,20,40};if(argc>1)stories={std::stoi(argv[1])};std::cout<<std::fixed<<std::setprecision(6);for(int ns:stories){auto m=make_model(ns);auto u=std::vector<double>(static_cast<std::size_t>(m.dof()),0.0);auto t=m.initial_nonlinear_tangents();auto t0=std::chrono::steady_clock::now();auto a=assess_initial_stability_from_tangents(m,u,t,1e-12,0);double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();std::cout<<"stories="<<ns<<" dof="<<m.dof()<<" K_nnz="<<m.K_initial().nnz()<<" status="<<(a.positive_definite?"SPD":"NOT_SPD")<<" min_pivot_ratio="<<a.minimum_pivot_ratio<<" factor_nnz="<<a.factor_nonzeros<<" seconds="<<sec<<" symmetry_error="<<a.symmetry_relative_error<<"\n";}}
