#include "quake/frame3d.hpp"
#include "quake/modal.hpp"
#include "quake/newmark.hpp"
#include "quake/reduction.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;

static int nid(int s,int ix,int iy){return 10000*s+100*ix+iy;}
static CompiledFrame3D make_space_frame(int stories,int nx=2,int ny=2){
    Frame3DBuilder b;constexpr double h=3.5,L=6.0,E=25000.0,G=10000.0;const int ncols=(nx+1)*(ny+1);
    for(int s=0;s<=stories;++s)for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy){const int id=nid(s,ix,iy);const double m=s?1.0/ncols:0.0;b.add_node(id,ix*L,iy*L,s*h,m,m,0.05*m);if(!s)b.fix(id);}
    int eid=1,sid=1,did=500000;
    for(int s=1;s<=stories;++s){const double p=6.0*(stories-s+1);for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<=ny;++iy)b.add_elastic_frame(eid++,nid(s-1,ix,iy),nid(s,ix,iy),E,G,2.0,0.30,0.55,0.55,1,0,0,p);
        for(int ix=0;ix<nx;++ix)for(int iy=0;iy<=ny;++iy){const int jl=nid(s,ix,iy),jr=nid(s,ix+1,iy),dl=did++,dr=did++;b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,(ix+1)*L,iy*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RX,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,0.20,0.35,0.35,0,1,0);b.add_bilinear_spring(sid++,jl,Dof3D::RY,dl,Dof3D::RY,50000,8.0,0.50);b.add_bilinear_spring(sid++,jr,Dof3D::RY,dr,Dof3D::RY,50000,8.0,0.50);}
        for(int ix=0;ix<=nx;++ix)for(int iy=0;iy<ny;++iy){const int jl=nid(s,ix,iy),jr=nid(s,ix,iy+1),dl=did++,dr=did++;b.add_node(dl,ix*L,iy*L,s*h);b.add_node(dr,ix*L,(iy+1)*L,s*h);for(auto d:{Dof3D::UX,Dof3D::UY,Dof3D::UZ,Dof3D::RY,Dof3D::RZ}){b.equal_dof(jl,dl,d);b.equal_dof(jr,dr,d);}b.add_elastic_frame(eid++,dl,dr,E,G,1.5,0.20,0.35,0.35,1,0,0);b.add_bilinear_spring(sid++,jl,Dof3D::RX,dl,Dof3D::RX,50000,8.0,0.50);b.add_bilinear_spring(sid++,jr,Dof3D::RX,dr,Dof3D::RX,50000,8.0,0.50);}}
    b.set_rayleigh(0.02,0.0005);b.set_ground_direction(Dof3D::UX);b.set_response(nid(stories,0,0),Dof3D::UX);std::vector<int> sn;for(int s=1;s<=stories;++s)sn.push_back(nid(s,0,0));b.set_story_nodes(std::move(sn),Dof3D::UX);return b.compile();
}
static double hist_error(const std::vector<double>& a,const std::vector<double>& b){double num=0,den=0;for(std::size_t i=0;i<a.size();++i){num=std::max(num,std::abs(a[i]-b[i]));den=std::max(den,std::abs(a[i]));}return num/std::max(1e-12,den);}
static double active_fraction(const AnalysisStats& s){return s.nonlinear_component_evaluations?double(s.nonlinear_active_component_evaluations)/double(s.nonlinear_component_evaluations):0.0;}
int main(int argc,char**argv){const int stories=argc>1?std::stoi(argv[1]):8;const int modes=argc>2?std::stoi(argv[2]):12;const int steps=argc>3?std::stoi(argv[3]):300;const double amp=argc>4?std::stod(argv[4]):2.0;const double dt=0.005;
    auto full=make_space_frame(stories);auto retained=nonlinear_support_dofs(full);ReductionBuildInfo gi,ci;auto guyan=craig_bampton_reduce(full,retained,0,&gi);auto cb=craig_bampton_reduce(full,retained,modes,&ci);
    auto mf=modal_analysis(full,3),mg=modal_analysis(guyan,3),mc=modal_analysis(cb,3);auto gm=synthetic_ground_motion(steps,dt,amp);
    auto rf=run_newmark(full,gm,dt,LinearStrategy::Woodbury,1e-8,35);auto rg=run_newmark(guyan,gm,dt,LinearStrategy::Woodbury,1e-8,35);auto rc=run_newmark(cb,gm,dt,LinearStrategy::Woodbury,1e-8,35);
    std::cout<<std::fixed<<std::setprecision(6);std::cout<<"full_dof="<<full.dof()<<" retained="<<retained.size()<<" guyan_dof="<<guyan.dof()<<" cb_modes="<<ci.fixed_interface_modes<<" cb_dof="<<cb.dof()<<" nonlinear="<<full.nonlinear_count()<<"\n";
    for(int i=0;i<std::min({3,(int)mf.size(),(int)mg.size(),(int)mc.size()});++i){std::cout<<"mode="<<i+1<<" full_T="<<mf[i].period<<" guyan_T="<<mg[i].period<<" guyan_err="<<std::abs(mg[i].period-mf[i].period)/mf[i].period<<" cb_T="<<mc[i].period<<" cb_err="<<std::abs(mc[i].period-mf[i].period)/mf[i].period<<"\n";}
    std::cout<<"full_seconds="<<rf.stats.elapsed_seconds<<" guyan_seconds="<<rg.stats.elapsed_seconds<<" cb_seconds="<<rc.stats.elapsed_seconds<<"\n";
    std::cout<<"guyan_speedup="<<rf.stats.elapsed_seconds/rg.stats.elapsed_seconds<<" cb_speedup="<<rf.stats.elapsed_seconds/rc.stats.elapsed_seconds<<"\n";
    std::cout<<"guyan_history_err="<<hist_error(rf.roof_history,rg.roof_history)<<" cb_history_err="<<hist_error(rf.roof_history,rc.roof_history)<<"\n";
    std::cout<<"full_peak="<<rf.stats.max_roof_abs<<" guyan_peak="<<rg.stats.max_roof_abs<<" cb_peak="<<rc.stats.max_roof_abs<<"\n";
    std::cout<<"full_active_fraction="<<active_fraction(rf.stats)<<" cb_active_fraction="<<active_fraction(rc.stats)<<"\n";
}
