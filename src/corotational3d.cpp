#include "quake/corotational3d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {
namespace {

using V3=std::array<double,3>;
using M3=std::array<double,9>; // row-major

double dot(const V3&a,const V3&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
V3 cross(const V3&a,const V3&b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double norm(const V3&a){return std::sqrt(dot(a,a));}
V3 normalize(V3 a){const double n=norm(a);if(n<1e-14)throw std::invalid_argument("corotational frame degenerate axis");for(double&v:a)v/=n;return a;}
M3 eye(){return {1,0,0,0,1,0,0,0,1};}
M3 transpose(const M3&a){M3 r{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)r[3*i+j]=a[3*j+i];return r;}
M3 mul(const M3&a,const M3&b){M3 r{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)r[3*i+j]+=a[3*i+k]*b[3*k+j];return r;}
M3 skew(const V3&v){return {0,-v[2],v[1],v[2],0,-v[0],-v[1],v[0],0};}
M3 add_scaled(M3 a,const M3&b,double s){for(int i=0;i<9;++i)a[i]+=s*b[i];return a;}

M3 exp_rot(const V3&w){
    const double th=norm(w); if(th<1e-12){auto K=skew(w);return add_scaled(eye(),K,1.0);}
    V3 a{w[0]/th,w[1]/th,w[2]/th};auto K=skew(a);auto K2=mul(K,K);
    auto R=add_scaled(eye(),K,std::sin(th));return add_scaled(R,K2,1.0-std::cos(th));
}

V3 log_rot(const M3&R){
    double c=0.5*(R[0]+R[4]+R[8]-1.0);c=std::clamp(c,-1.0,1.0);const double th=std::acos(c);
    V3 v{R[7]-R[5],R[2]-R[6],R[3]-R[1]};
    if(th<1e-9){return {0.5*v[0],0.5*v[1],0.5*v[2]};}
    if(M_PI-th<1e-6){
        // Robust axis extraction near pi.
        V3 axis{std::sqrt(std::max(0.0,(R[0]+1.0)*0.5)),std::sqrt(std::max(0.0,(R[4]+1.0)*0.5)),std::sqrt(std::max(0.0,(R[8]+1.0)*0.5))};
        if(R[1]<0) axis[1]=-axis[1];
        if(R[2]<0) axis[2]=-axis[2];
        axis=normalize(axis);
        return {th*axis[0],th*axis[1],th*axis[2]};
    }
    const double s=th/(2.0*std::sin(th));return {s*v[0],s*v[1],s*v[2]};
}

M3 frame_from_axes(const V3&ex,const V3&ey,const V3&ez){
    // Columns are local axes in global coordinates.
    return {ex[0],ey[0],ez[0],ex[1],ey[1],ez[1],ex[2],ey[2],ez[2]};
}

void initial_frame(const CorotationalFrame3DProperties&p,V3&ex,V3&ey,V3&ez,M3&R0,double&L0){
    V3 d{p.xj-p.xi,p.yj-p.yi,p.zj-p.zi};L0=norm(d);if(L0<=0)throw std::invalid_argument("zero-length corotational frame");ex=normalize(d);
    const double proj=dot(p.reference,ex);V3 yr{p.reference[0]-proj*ex[0],p.reference[1]-proj*ex[1],p.reference[2]-proj*ex[2]};ey=normalize(yr);ez=normalize(cross(ex,ey));ey=normalize(cross(ez,ex));R0=frame_from_axes(ex,ey,ez);
}

M3 minimal_rotation(const V3&a0,const V3&a1){
    const V3 a=normalize(a0),b=normalize(a1);const double c=std::clamp(dot(a,b),-1.0,1.0);V3 v=cross(a,b);const double s=norm(v);
    if(s<1e-12){
        if(c>0)return eye();
        V3 trial=std::abs(a[0])<0.8?V3{1,0,0}:V3{0,1,0};V3 axis=normalize(cross(a,trial));return exp_rot({M_PI*axis[0],M_PI*axis[1],M_PI*axis[2]});
    }
    V3 axis{v[0]/s,v[1]/s,v[2]/s};const double th=std::atan2(s,c);return exp_rot({th*axis[0],th*axis[1],th*axis[2]});
}

struct Kinematics{std::array<double,6> q{};double L{};};
Kinematics kinematics(const CorotationalFrame3DProperties&p,const std::array<double,12>&u){
    if(p.E<=0||p.G<=0||p.A<=0||p.J<=0||p.Iy<=0||p.Iz<=0||p.axial_compression<0)throw std::invalid_argument("invalid corotational frame properties");
    V3 ex0,ey0,ez0;M3 R0;double L0;initial_frame(p,ex0,ey0,ez0,R0,L0);
    V3 xi{p.xi+u[0],p.yi+u[1],p.zi+u[2]},xj{p.xj+u[6],p.yj+u[7],p.zj+u[8]};V3 d{xj[0]-xi[0],xj[1]-xi[1],xj[2]-xi[2]};const double L=norm(d);if(L<1e-10*L0)throw std::runtime_error("corotational frame collapsed to zero chord length");const V3 ex=normalize(d);
    const M3 Qalign=minimal_rotation(ex0,ex);M3 Rbar=mul(Qalign,R0);
    V3 thi{u[3],u[4],u[5]},thj{u[9],u[10],u[11]};const M3 Ri=mul(exp_rot(thi),R0),Rj=mul(exp_rot(thj),R0);
    auto rel_i=log_rot(mul(transpose(Rbar),Ri));auto rel_j=log_rot(mul(transpose(Rbar),Rj));
    const double psi=0.5*(rel_i[0]+rel_j[0]);const M3 twist=exp_rot({psi,0,0});M3 Rf=mul(Rbar,twist);
    rel_i=log_rot(mul(transpose(Rf),Ri));rel_j=log_rot(mul(transpose(Rf),Rj));
    Kinematics k;k.L=L;k.q={L-L0,rel_j[0]-rel_i[0],rel_i[1],rel_j[1],rel_i[2],rel_j[2]};return k;
}

double raw_energy(const CorotationalFrame3DProperties&p,const std::array<double,12>&u){
    V3 ex0,ey0,ez0;M3 R0;double L0;initial_frame(p,ex0,ey0,ez0,R0,L0);const auto k=kinematics(p,u);const auto&q=k.q;
    double U=0.5*(p.E*p.A/L0)*q[0]*q[0]-p.axial_compression*q[0]+0.5*(p.G*p.J/L0)*q[1]*q[1];
    const double ky=p.E*p.Iy/L0,kz=p.E*p.Iz/L0;
    U+=0.5*ky*(4*q[2]*q[2]+4*q[3]*q[3]+4*q[2]*q[3]);
    U+=0.5*kz*(4*q[4]*q[4]+4*q[5]*q[5]+4*q[4]*q[5]);
    return U;
}

double scale_for_dof(const CorotationalFrame3DProperties&p,int i){
    V3 d{p.xj-p.xi,p.yj-p.yi,p.zj-p.zi};const double L=norm(d);return (i%6)<3?std::max(1.0,L):1.0;
}

std::array<double,12> energy_gradient(const CorotationalFrame3DProperties&p,const std::array<double,12>&u,double rel){
    std::array<double,12> g{};for(int i=0;i<12;++i){const double h=std::max(1e-9,rel*scale_for_dof(p,i));auto up=u,um=u;up[i]+=h;um[i]-=h;g[i]=(raw_energy(p,up)-raw_energy(p,um))/(2*h);}return g;
}

} // namespace

std::array<double,6> corotational3d_basic_deformation(const CorotationalFrame3DProperties&p,const std::array<double,12>&u){return kinematics(p,u).q;}
double corotational3d_strain_energy(const CorotationalFrame3DProperties&p,const std::array<double,12>&u){const std::array<double,12> z{};return raw_energy(p,u)-raw_energy(p,z);}

CorotationalFrame3DSectionResponse corotational3d_section_response(
    const CorotationalFrame3DProperties&p,const std::array<double,12>&u){
    V3 ex0,ey0,ez0;M3 R0;double L0;initial_frame(p,ex0,ey0,ez0,R0,L0);
    const auto k=kinematics(p,u);const auto&q=k.q;
    CorotationalFrame3DSectionResponse r;r.current_length=k.L;
    // The raw energy contains -P0*deltaL; the solver subtracts the undeformed
    // prestress residual. The physical section force therefore retains P0.
    r.axial_compression=p.axial_compression-(p.E*p.A/L0)*q[0];
    const double T=(p.G*p.J/L0)*q[1];r.torsion_i=-T;r.torsion_j=T;
    const double ky=p.E*p.Iy/L0,kz=p.E*p.Iz/L0;
    r.moment_y_i=ky*(4.0*q[2]+2.0*q[3]);
    r.moment_y_j=ky*(2.0*q[2]+4.0*q[3]);
    r.moment_z_i=kz*(4.0*q[4]+2.0*q[5]);
    r.moment_z_j=kz*(2.0*q[4]+4.0*q[5]);
    // No distributed member loads are present in this element. Local end shear
    // follows exactly from end-moment equilibrium on the current chord.
    r.shear_y_i=(r.moment_z_i+r.moment_z_j)/k.L;r.shear_y_j=-r.shear_y_i;
    r.shear_z_i=-(r.moment_y_i+r.moment_y_j)/k.L;r.shear_z_j=-r.shear_z_i;
    return r;
}

std::array<double,12> corotational3d_internal_force(const CorotationalFrame3DProperties&p,const std::array<double,12>&u,double rel){
    auto f=energy_gradient(p,u,rel);const std::array<double,12> z{};const auto f0=energy_gradient(p,z,rel);for(int i=0;i<12;++i)f[i]-=f0[i];return f;
}

CorotationalFrame3DResponse corotational3d_response(const CorotationalFrame3DProperties&p,const std::array<double,12>&u,double rel){
    CorotationalFrame3DResponse r;r.basic_deformation=corotational3d_basic_deformation(p,u);r.strain_energy=corotational3d_strain_energy(p,u);V3 d{p.xj+u[6]-p.xi-u[0],p.yj+u[7]-p.yi-u[1],p.zj+u[8]-p.zi-u[2]};r.current_length=norm(d);r.force=corotational3d_internal_force(p,u,std::max(2e-8,rel*0.1));
    for(int j=0;j<12;++j){const double h=std::max(2e-8,rel*scale_for_dof(p,j));auto up=u,um=u;up[j]+=h;um[j]-=h;const auto fp=corotational3d_internal_force(p,up,std::max(2e-8,rel*0.05));const auto fm=corotational3d_internal_force(p,um,std::max(2e-8,rel*0.05));for(int i=0;i<12;++i)r.tangent[static_cast<std::size_t>(i*12+j)]=(fp[i]-fm[i])/(2*h);}
    // Numerical differentiation is theoretically symmetric because it comes
    // from an energy potential. Explicit symmetrization suppresses differencing
    // noise and is useful for the stability/inertia checks.
    for(int i=0;i<12;++i)for(int j=i+1;j<12;++j){const double a=0.5*(r.tangent[static_cast<std::size_t>(i*12+j)]+r.tangent[static_cast<std::size_t>(j*12+i)]);r.tangent[static_cast<std::size_t>(i*12+j)]=a;r.tangent[static_cast<std::size_t>(j*12+i)]=a;}
    return r;
}

} // namespace quake
