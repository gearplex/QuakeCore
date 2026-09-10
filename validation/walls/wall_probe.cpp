#include "wall_json.hpp"
#include "quake/dynamic_model.hpp"
#include <iostream>
using namespace quake;
int main(int argc,char** argv){try{
 if(argc!=3)throw std::invalid_argument("wall_probe input.json output.json");auto j=read_json(argv[1]);Json out;
 if(j.contains("material")){
  auto m=wall_uniaxial(j.at("material"));std::vector<double>s(6),t(6);m.initialize(s.data());out["history"]=Json::array();
  for(double e:j.at("strains")){auto r=m.trial(e,s.data(),t.data());out["history"].push_back({e,r.stress,r.tangent});s=t;}
 }else{
  const auto& w=j.at("wall");double h=j.at("height");Wall2D e=w.at("type")=="mvlem"?Wall2D(h,mvlem_properties(w)):Wall2D(h,sfi_mvlem_properties(w));
  auto s=e.initial_state();std::array<double,6>u{};out["initial_tangent"]=e.initial_tangent();out["history"]=Json::array();
  if(j.contains("deformations")){
   for(const auto& d:j.at("deformations")){u=d.get<std::array<double,6>>();auto r=e.trial(u,s.data());auto row=wall_response_json(r);row["u"]=u;row["tangent"]=r.tangent;out["history"].push_back(row);s=r.state;}
  }else{
   for(const auto& target:j.at("protocol")){
    u[3]=target.at(0);double P=target.at(1),M=target.at(2);bool ok=false;Wall2DResponse r;
    for(int it=0;it<150;++it){
     r=e.trial(u,s.data());double f=r.force[4]-P,g=r.force[5]-M,norm=std::max(std::abs(f),std::abs(g));
     if(norm<=1e-8){ok=true;break;}
     double a=r.tangent[28],b=r.tangent[29],c=r.tangent[34],d=r.tangent[35],det=a*d-b*c;
     if(!std::isfinite(det)||std::abs(det)<1e-15*std::max(std::abs(a*d),std::abs(b*c)))break;
     double dv=(-d*f+b*g)/det,dr=(c*f-a*g)/det,alpha=1;bool accepted=false;
     for(int ls=0;ls<22;++ls){auto v=u;v[4]+=alpha*dv;v[5]+=alpha*dr;
      try{auto z=e.trial(v,s.data());if(std::max(std::abs(z.force[4]-P),std::abs(z.force[5]-M))<norm){u=v;accepted=true;break;}}catch(const ConstitutiveIntegrationError&){}alpha*=.5;}
     if(!accepted)break;
    }
    if(!ok){out["completed"]=false;out["failed_target"]=target;write_json(argv[2],out);return 3;}
    auto row=wall_response_json(r);row["u"]=u;row["tangent"]=r.tangent;out["history"].push_back(row);s=r.state;
   }
  }out["final_state"]=s;
 }out["completed"]=true;write_json(argv[2],out);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
