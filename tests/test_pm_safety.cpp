#define main pm_audit_main
#include "../validation/component_audit/pm_audit.cpp"
#undef main
#include <stdexcept>
int main(){try{
 {auto p=params(1);p.alpha_tension=.5;bool rejected=false;
  try{PMInteractionHinge2D h(p);}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("convex return map accepted nonconvex surface");}
 for(auto ev:{PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic,PMInteractionSurfaceEvolution::MrozTwoSurface}){
  auto p=params(1);p.surface_evolution=ev;PMInteractionHinge2D h(p);std::vector<double> c(h.kStateSize),s(c.size());h.initialize_state(c.data());
  for(double da:{-.05,.03}){auto t=h.trial(da,0,c.data(),s.data());if(!t.converged||t.compression<p.py_tension-1e-8||t.compression>p.p_balance+1.2*(p.py_compression-p.p_balance)+1e-8)throw std::runtime_error("pure axial projection violated the surface");if(std::abs(t.axial_force-p.axial_stiffness*(da-t.plastic_axial_deformation))>1e-7)throw std::runtime_error("axial tip force/state inconsistency");}
 }
 PMInteractionHinge2D a(params(1)),b(params(25.4));std::vector<double> ca(a.kStateSize),cb(a.kStateSize),sa(ca.size()),sb(ca.size());a.initialize_state(ca.data());b.initialize_state(cb.data());
 double error=0;
 for(double rot:{.01,.02,.005,-.01,-.03,.01,.04}){
  auto ra=a.trial(-.0005,rot,ca.data(),sa.data()),rb=b.trial(-.0005*25.4,rot,cb.data(),sb.data());
  if(!ra.converged||!rb.converged)throw std::runtime_error("unit conversion test return failed");
  error=std::max(error,std::abs(ra.moment-rb.moment/25.4));
  if(std::abs(ra.plastic_axial_deformation-rb.plastic_axial_deformation/25.4)>1e-9)throw std::runtime_error("Mroz axial flow depends on length units");
  ca=sa;cb=sb;
 }
 if(error>1e-7)throw std::runtime_error("Mroz cyclic moment depends on length units");
 // Penalty stiffness and strength scales representative of the Berkeley
 // columns; include trials outside both axial intercepts and large moments.
 for(double da:{-.004,-.002,0.,.002,.004})for(double rot:{-.04,-.003,-1e-8,1e-8,.003,.04}){
  auto p=params(1);p.surface_evolution=PMInteractionSurfaceEvolution::ElasticPerfectlyPlastic;
  p.return_tolerance=1e-12;p.axial_stiffness=90000;p.hinge.Ke=1400000;p.my_balance=177.4861105332723;p.p_balance=47.22493333333333;
  p.py_tension=-56.32;p.py_compression=165.562;p.alpha_tension=1.524;p.alpha_compression=1.386;
  PMInteractionHinge2D h(p);std::vector<double> c(h.kStateSize),s(c.size());h.initialize_state(c.data());
  auto r=h.trial(da,rot,c.data(),s.data());if(!r.converged)throw std::runtime_error("penalty-scale P-M projection failed");
  const bool comp=r.compression>=p.p_balance;double den=(comp?p.py_compression:p.py_tension)-p.p_balance;
  const double surface=std::pow(std::abs((r.compression-p.p_balance)/den),comp?p.alpha_compression:p.alpha_tension)+std::pow(std::abs(r.moment/p.my_balance),p.beta_pm);
  if(surface>1+1e-8)throw std::runtime_error("P-M projection outside yield surface");
  if(std::abs(r.axial_force-p.axial_stiffness*(da-r.plastic_axial_deformation))>1e-5)throw std::runtime_error("P-M projection violates elastic compatibility");
  if(std::abs(rot)>.001){
   const double eps[2]{1e-7,1e-6};
   for(int j=0;j<2;++j){std::vector<double> sp(c.size()),sm(c.size());
    auto rp=h.trial(da+(j==0?eps[j]:0),rot+(j==1?eps[j]:0),c.data(),sp.data());
    auto rm=h.trial(da-(j==0?eps[j]:0),rot-(j==1?eps[j]:0),c.data(),sm.data());
    if(!rp.converged||!rm.converged)throw std::runtime_error("P-M tangent reference failed");
    double fd[2]{(rp.axial_force-rm.axial_force)/(2*eps[j]),(rp.moment-rm.moment)/(2*eps[j])};
    for(int i=0;i<2;++i)if(std::abs(fd[i]-r.tangent[2*i+j])>1e-3*std::max(1.0,std::abs(fd[i]))){std::cerr<<"da="<<da<<" theta="<<rot<<" i="<<i<<" j="<<j<<" fd="<<fd[i]<<" tangent="<<r.tangent[2*i+j]<<" P="<<r.compression<<" M="<<r.moment<<"\n";throw std::runtime_error("P-M tangent disagrees with force finite difference");}
   }
  }
 }
 std::cout<<"P-M domain and unit-invariance tests passed; max moment difference="<<error<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
