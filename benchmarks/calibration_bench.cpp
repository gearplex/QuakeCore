#include "quake/calibration.hpp"
#include "quake/frame3d.hpp"
#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;
static int nid(int s,int ix,int iy){return 10000*s+100*ix+iy;}
static CompiledFrame3D make_model(int stories,int nx,int ny){
    Frame3DBuilder b; constexpr double h=3.5,L=6,E=25000,G=10000; int eid=1,sid=1,did=700000; const int nc=(nx+1)*(ny+1);
    for(int s=0;s<=stories;++s)for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy){double m=s?1.0/nc:0;b.add_node(nid(s,ix,iy),ix*L,iy*L,s*h,m,m,.05*m);if(!s)b.fix(nid(s,ix,iy));}
    for(int s=1;s<=stories;++s){double p=8.0*(stories-s+1);for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy)b.add_elastic_frame(eid++,nid(s-1,ix,iy),nid(s,ix,iy),E,G,2,.3,.55,.55,1,0,0,p);
        for(int ix=0;ix<nx;++ix)for(int iy=0;iy<=ny;++iy){int jl=nid(s,ix,iy),jr=nid(s,ix+1,iy),dl=did++,dr=did++;b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,(ix+1)*L,iy*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RX,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,.2,.35,.35,0,1,0);b.add_bilinear_spring(sid++,jl,Dof3D::RY,dl,Dof3D::RY,50000,8,.5);b.add_bilinear_spring(sid++,jr,Dof3D::RY,dr,Dof3D::RY,50000,8,.5);}
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<ny;++iy){int jl=nid(s,ix,iy),jr=nid(s,ix,iy+1),dl=did++,dr=did++;b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,ix*L,(iy+1)*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RY,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,.2,.35,.35,1,0,0);b.add_bilinear_spring(sid++,jl,Dof3D::RX,dl,Dof3D::RX,50000,8,.5);b.add_bilinear_spring(sid++,jr,Dof3D::RX,dr,Dof3D::RX,50000,8,.5);}}
    b.set_rayleigh(.02,.0005);b.set_ground_direction(Dof3D::UX);b.set_response(nid(stories,0,0),Dof3D::UX);return b.compile();
}
int main(int argc,char**argv){int stories=argc>1?std::stoi(argv[1]):10,nx=argc>2?std::stoi(argv[2]):2,ny=argc>3?std::stoi(argv[3]):2,reps=argc>4?std::stoi(argv[4]):5;auto m=make_model(stories,nx,ny);auto c=calibrate_solver_crossover(m,.005,{},reps);std::cout<<std::fixed<<std::setprecision(6)<<"dof="<<c.dof<<" nonlinear="<<c.nonlinear_count<<" nnz="<<c.matrix_nnz<<" stories="<<stories<<" bays="<<nx<<"x"<<ny<<"\n";for(auto&p:c.points)std::cout<<"rank="<<p.active_rank<<" wood="<<p.woodbury_seconds<<" direct="<<p.same_pattern_refactor_seconds<<" speedup="<<p.speedup<<"x\n";std::cout<<"recommended_woodbury_rank_limit="<<c.recommended_woodbury_rank_limit<<" fraction_of_dof="<<c.recommended_rank_fraction<<"\n";}
