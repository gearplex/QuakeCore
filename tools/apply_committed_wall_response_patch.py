from pathlib import Path

# Response extraction must not advance hysteretic history.  The transient
# solver already passes the committed state corresponding to the supplied
# displacement; calling trial() again at the same strain can create a second
# reversal/transition in non-idempotent cyclic laws such as ConcreteCM.

wm_h = Path("include/quake/wall_material.hpp")
text = wm_h.read_text()
old = "    Result trial(double strain, const double* committed, double* trial_state) const;\n"
new = old + "    Result committed_response(double strain, const double* committed) const;\n"
if text.count(old) != 1:
    raise SystemExit("wall_material.hpp trial declaration anchor not found exactly once")
wm_h.write_text(text.replace(old, new, 1))

wm_cpp = Path("src/wall_material.cpp")
text = wm_cpp.read_text()
anchor = "WallPanel::WallPanel(double E,double nu,std::vector<WallPanelLayer> layers):layers_(std::move(layers)) {\n"
impl = r'''WallUniaxial::Result WallUniaxial::committed_response(double e,const double* s) const {
    if(!std::isfinite(e))throw std::invalid_argument("nonfinite wall material strain");
    const int nstate=state_size();for(int i=0;i<nstate;++i)if(!std::isfinite(s[i]))throw std::invalid_argument("nonfinite wall material state");
    if(kind_==Kind::ConcreteCM){const auto c=decode_concrete_cm_state(s);return {c.stress,c.tangent};}
    if(kind_==Kind::Pinching4){const auto c=decode_pinching4_state(s);return {c.force,c.tangent};}
    if(kind_==Kind::MinMax){
        const auto& child=children_->front();const int n=child.state_size();
        if(s[n]!=0.0)return {0,1.0e-8*child.initial_tangent()};
        return child.committed_response(e,s);
    }
    if(kind_==Kind::Parallel){
        double stress=0,tangent=0;int o=0;
        for(const auto& child:*children_){auto r=child.committed_response(e,s+o);stress+=r.stress;tangent+=r.tangent;o+=child.state_size();}
        return {stress,tangent};
    }
    if(kind_==Kind::Elastic)return {E_*e,E_};
    if(kind_==Kind::SteelBilinear){
        const double stress=E_*(e-s[0]);
        const double xi=stress-s[1];
        const double H=b_==0.0?0.0:b_*E_/(1.0-b_);
        const double tol=1e-12*std::max(1.0,fy_);
        const double tangent=std::abs(xi)>=fy_-tol ? E_*H/(E_+H) : E_;
        return {stress,tangent};
    }
    // Concrete01 explicitly stores the last committed stress/tangent.
    return {s[4],s[5]==0.0 && e==0.0 && s[0]==0.0 ? E_ : s[5]};
}

'''
if text.count(anchor) != 1:
    raise SystemExit("wall_material.cpp WallPanel anchor not found exactly once")
text = text.replace(anchor, impl + anchor, 1)

panel_anchor = "WallPanelResult WallPanel::trial(const std::array<double,3>& e,const double* s) const {\n"
panel_impl = r'''WallPanelResult WallPanel::committed_response(const std::array<double,3>& e,const double* s) const {
    for(double v:e)if(!std::isfinite(v))throw std::invalid_argument("nonfinite panel strain");
    WallPanelResult r;r.tangent=background_;r.state.assign(s,s+state_size());
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)r.stress[a]+=background_[3*a+b]*e[b];
    for(std::size_t i=0;i<layers_.size();++i){
        const auto& n=directions_[i];double q=n[0]*e[0]+n[1]*e[1]+n[2]*e[2];
        auto v=layers_[i].material.committed_response(q,s+offsets_[i]);
        for(int a=0;a<3;++a){r.stress[a]+=layers_[i].weight*v.stress*n[a];for(int b=0;b<3;++b)r.tangent[3*a+b]+=layers_[i].weight*v.tangent*n[a]*n[b];}
    }
    return r;
}

'''
if text.count(panel_anchor) != 1:
    raise SystemExit("wall_material.cpp WallPanel trial anchor not found exactly once")
wm_cpp.write_text(text.replace(panel_anchor, panel_impl + panel_anchor, 1))

wm_h = Path("include/quake/wall_material.hpp")
text = wm_h.read_text()
old = "    WallPanelResult trial(const std::array<double,3>& strain, const double* committed) const;\n"
new = old + "    WallPanelResult committed_response(const std::array<double,3>& strain, const double* committed) const;\n"
if text.count(old) != 1:
    raise SystemExit("WallPanel trial declaration anchor not found exactly once")
wm_h.write_text(text.replace(old, new, 1))

w_h = Path("include/quake/wall2d.hpp")
text = w_h.read_text()
old = "    Wall2DResponse trial(const std::array<double,6>& u,const double* committed) const;\n"
new = old + "    Wall2DResponse committed_response(const std::array<double,6>& u,const double* committed) const;\n"
if text.count(old) != 1:
    raise SystemExit("wall2d.hpp trial declaration anchor not found exactly once")
w_h.write_text(text.replace(old, new, 1))

w_cpp = Path("src/wall2d.cpp")
text = w_cpp.read_text()
anchor = "Wall2DResponse Wall2D::trial(const std::array<double,6>& u,const double* s) const {\n"
impl = r'''Wall2DResponse Wall2D::committed_response(const std::array<double,6>& u,const double* s) const {
    for(double v:u)if(!std::isfinite(v))throw std::invalid_argument("nonfinite wall displacement");
    Wall2DResponse r;r.state.assign(s,s+state_size_);auto g=shear_B();double gamma=dot(g,u);r.shear_deformation=h_*gamma;r.curvature=(u[5]-u[2])/h_;
    if(!is_sfi()){
        const auto& p=std::get<MVLEMProperties>(properties_);
        for(std::size_t i=0;i<p.fibers.size();++i){
            const auto& f=p.fibers[i];auto y=axial_B(i);double e=dot(y,u);int o=offsets_[i];int co=f.concrete.state_size();
            auto a=f.concrete.committed_response(e,s+o),b=f.steel.committed_response(e,s+o+co);
            double A=f.width*f.thickness,V=h_*A,rho=f.reinforcement_ratio;
            for(int j=0;j<6;++j)r.force[j]+=V*((1-rho)*a.stress+rho*b.stress)*y[j];
            outer(r.tangent,y,y,V*((1-rho)*a.tangent+rho*b.tangent));
            r.fiber_strain.push_back(e);r.concrete_stress.push_back(a.stress);r.steel_stress.push_back(b.stress);
        }
        auto sh=p.shear.committed_response(r.shear_deformation,s+shear_offset_);
        for(int j=0;j<6;++j)r.force[j]+=h_*sh.stress*g[j];
        outer(r.tangent,g,g,h_*h_*sh.tangent);
    }else{
        const auto& p=std::get<SFIMVLEMProperties>(properties_);
        for(std::size_t i=0;i<p.panels.size();++i){
            const auto& f=p.panels[i];auto y=axial_B(i);double ey=dot(y,u);int o=offsets_[i];double ex=s[o];
            if(!std::isfinite(ex))throw std::invalid_argument("nonfinite transverse panel history");
            auto v=f.material.committed_response({ex,ey,gamma},s+o+1);
            const auto D0=f.material.initial_tangent();const auto d=condense(v.tangent,D0[0]);double V=h_*f.width*f.thickness;
            for(int j=0;j<6;++j)r.force[j]+=V*(v.stress[1]*y[j]+v.stress[2]*g[j]);
            panel_stiffness(r.tangent,y,g,d,V);
            r.panel_strain.push_back({ex,ey,gamma});r.panel_stress.push_back(v.stress);
            r.maximum_transverse_stress_residual=std::max(r.maximum_transverse_stress_residual,std::abs(v.stress[0]));
        }
    }
    return r;
}

'''
if text.count(anchor) != 1:
    raise SystemExit("wall2d.cpp trial anchor not found exactly once")
w_cpp.write_text(text.replace(anchor, impl + anchor, 1))

f_cpp = Path("src/frame2d.cpp")
text = f_cpp.read_text()
old = "return w.element.trial(q,s.data()+w.state_offset);"
new = "return w.element.committed_response(q,s.data()+w.state_offset);"
if text.count(old) != 1:
    raise SystemExit("frame2d.cpp wall_response trial anchor not found exactly once")
f_cpp.write_text(text.replace(old, new, 1))
