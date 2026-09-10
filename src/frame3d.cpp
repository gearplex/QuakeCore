#include "quake/frame3d.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <limits>
#include <set>
#include <sstream>
#include <iomanip>
#include <stdexcept>

extern "C" {
void dsyev_(const char* jobz,const char* uplo,const int* n,double* a,const int* lda,
            double* w,double* work,const int* lwork,int* info);
}

namespace quake {
namespace {

int dof_index(int node_slot,Dof3D dof){return 6*node_slot+static_cast<int>(dof);}

double dot3(const std::array<double,3>& a,const std::array<double,3>& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
std::array<double,3> cross3(const std::array<double,3>& a,const std::array<double,3>& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double norm3(const std::array<double,3>& a){return std::sqrt(dot3(a,a));}
std::array<double,3> normalize3(std::array<double,3> a){const double n=norm3(a);if(n<1e-14)throw std::invalid_argument("degenerate 3D local axis");for(double& v:a)v/=n;return a;}

std::array<double,144> transform_12(const std::array<double,144>& K,const std::array<double,9>& R){
    std::array<double,144> T{};
    for(int block=0;block<4;++block){
        const int o=3*block;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)T[static_cast<std::size_t>((o+i)*12+(o+j))]=R[static_cast<std::size_t>(i*3+j)];
    }
    std::array<double,144> KT{},out{};
    for(int i=0;i<12;++i)for(int j=0;j<12;++j)for(int k=0;k<12;++k)
        KT[static_cast<std::size_t>(i*12+j)] += K[static_cast<std::size_t>(i*12+k)]*T[static_cast<std::size_t>(k*12+j)];
    for(int i=0;i<12;++i)for(int j=0;j<12;++j)for(int k=0;k<12;++k)
        out[static_cast<std::size_t>(i*12+j)] += T[static_cast<std::size_t>(k*12+i)]*KT[static_cast<std::size_t>(k*12+j)];
    return out;
}

void add_scaled_bbt3(std::vector<Triplet>& t,const SparseUpdateBasis& basis,int col,double scale){
    if(scale==0.0)return;
    const int b=basis.col_ptr()[static_cast<std::size_t>(col)],e=basis.col_ptr()[static_cast<std::size_t>(col+1)];
    for(int pi=b;pi<e;++pi){const int i=basis.row_ind()[static_cast<std::size_t>(pi)];const double bi=basis.values()[static_cast<std::size_t>(pi)];
        for(int pj=b;pj<e;++pj){const int j=basis.row_ind()[static_cast<std::size_t>(pj)];const double bj=basis.values()[static_cast<std::size_t>(pj)];t.push_back({i,j,scale*bi*bj});}}
}


using ConstraintTerms = std::vector<std::pair<int,double>>;

void accumulate_terms(std::map<int,double>& out, const ConstraintTerms& terms, double scale) {
    for (const auto& [d,c] : terms) out[d] += scale*c;
}

ConstraintTerms compact_terms(const std::map<int,double>& in, double tol=1e-14) {
    ConstraintTerms out;
    for (const auto& [d,c] : in) if (std::abs(c)>tol) out.emplace_back(d,c);
    return out;
}

double eval_terms(const ConstraintTerms& terms, const std::vector<double>& q) {
    double v=0.0; for (const auto& [d,c] : terms) v += c*q[static_cast<std::size_t>(d)]; return v;
}

std::vector<double> monotone_capacity_slopes(const std::vector<double>& x,const std::vector<double>& y){
    if(x.size()<2||x.size()!=y.size())throw std::invalid_argument("coupled P-M capacity table size");
    const std::size_t n=x.size();std::vector<double> d(n-1),h(n-1),m(n,0.0);
    for(std::size_t i=0;i+1<n;++i){if(!(x[i+1]>x[i]))throw std::invalid_argument("coupled P-M axial-force points must increase");h[i]=x[i+1]-x[i];d[i]=(y[i+1]-y[i])/h[i];}
    m[0]=d[0];m[n-1]=d[n-2];
    for(std::size_t i=1;i+1<n;++i){
        if(d[i-1]*d[i]<=0.0)m[i]=0.0;
        else{const double w1=2.0*h[i]+h[i-1],w2=h[i]+2.0*h[i-1];m[i]=(w1+w2)/(w1/d[i-1]+w2/d[i]);}
    }
    if (m[0] * d[0] <= 0.0) m[0] = 0.0;
    if (m[n - 1] * d[n - 2] <= 0.0) m[n - 1] = 0.0;
    return m;
}

struct CapacityEval { double value{},derivative{}; };

// Force on the monotonic ASCE 41 envelope for a supplied parameter set and
// committed cyclic strength scales. This mirrors ASCE41HingeMaterial::envelope
// and is used only to keep a P-dependent reload target on the *current*
// envelope. The branch zero intercept and target deformation remain history
// variables; only the dimensional target force moves continuously with My(P).
double asce41_envelope_force(const ASCE41HingeParams& p,double u,double pos_scale,double neg_scale){
    const bool pos=u>=0.0;
    const double sg=pos?1.0:-1.0;
    const double x=std::abs(u);
    const double fy0=pos?p.posFy:p.negFy;
    const double scale=std::max(1e-6,pos?pos_scale:neg_scale);
    const double a=pos?p.pos_a:p.neg_a;
    const double b=pos?p.pos_b:p.neg_b;
    const double fpar=pos?p.pos_f:p.neg_f;
    const double f=(fpar>0.0)?fpar:b;
    const double c=pos?p.pos_c:p.neg_c;
    const double drop_req=pos?p.pos_drop_span:p.neg_drop_span;
    const double drop=drop_req>0.0?drop_req:0.05*(b-a);
    const double e_req=pos?p.pos_e_drop_span:p.neg_e_drop_span;
    const double e_drop=e_req>0.0?e_req:0.05*(b-a-drop);
    const double fy=std::max(1e-12*fy0,fy0*scale);
    const double uy=fy/p.Ke;
    const double uc=uy+a,ud=uc+drop,ue=uy+b,ue0=ue-e_drop,uf=uy+f;
    const double kh=(p.hardening_stiffness>0.0)?p.hardening_stiffness:p.hardening_ratio*p.Ke;
    const double fc=fy+kh*a;
    const double fr=fy0*c*scale;
    if(x>=uf || x>=ue)return 0.0;
    if(x<=uy)return sg*p.Ke*x;
    if(x<=uc)return sg*(fy+kh*(x-uy));
    if(x<=ud){const double kd=(fr-fc)/drop;return sg*(fc+kd*(x-uc));}
    if(x<=ue0)return sg*fr;
    const double kE=-fr/e_drop;
    return sg*std::max(0.0,fr+kE*(x-ue0));
}

CapacityEval eval_capacity_curve(const std::vector<double>& x,const std::vector<double>& y,const std::vector<double>& m,double q){
    if(x.size()<2||x.size()!=y.size()||m.size()!=x.size())throw std::invalid_argument("coupled P-M capacity curve");
    // Constant continuation outside the tabulated section range prevents an
    // artificial extrapolated strength slope in extreme Newton trials.
    if (q <= x.front()) return {y.front(), 0.0};
    if (q >= x.back()) return {y.back(), 0.0};
    auto it=std::upper_bound(x.begin(),x.end(),q);const std::size_t i=static_cast<std::size_t>(it-x.begin()-1);
    const double h=x[i+1]-x[i],t=(q-x[i])/h,t2=t*t,t3=t2*t;
    const double h00=2*t3-3*t2+1,h10=t3-2*t2+t,h01=-2*t3+3*t2,h11=t3-t2;
    const double v=h00*y[i]+h10*h*m[i]+h01*y[i+1]+h11*h*m[i+1];
    const double dh00=(6*t2-6*t)/h,dh10=3*t2-4*t+1,dh01=(-6*t2+6*t)/h,dh11=3*t2-2*t;
    const double dv=dh00*y[i]+dh10*m[i]+dh01*y[i+1]+dh11*m[i+1];
    return {v,dv};
}

int find_csc_position(const SparseMatrixCSC& a, int row, int col) {
    const int b=a.col_ptr()[static_cast<std::size_t>(col)], e=a.col_ptr()[static_cast<std::size_t>(col+1)];
    auto it=std::lower_bound(a.row_ind().begin()+b,a.row_ind().begin()+e,row);
    if(it==a.row_ind().begin()+e || *it!=row) return -1;
    return static_cast<int>(std::distance(a.row_ind().begin(),it));
}

struct SmallSymEigen { std::vector<double> values; std::vector<double> vectors_col_major; };
SmallSymEigen small_symmetric_eigen(std::vector<double> row_major,int n){
    SmallSymEigen out;if(n<=0)return out;if(static_cast<int>(row_major.size())!=n*n)throw std::invalid_argument("small eigen dimension");
    std::vector<double> a(static_cast<std::size_t>(n*n));for(int i=0;i<n;++i)for(int j=0;j<n;++j)a[static_cast<std::size_t>(j*n+i)]=0.5*(row_major[static_cast<std::size_t>(i*n+j)]+row_major[static_cast<std::size_t>(j*n+i)]);
    out.values.resize(static_cast<std::size_t>(n));const char job='V',uplo='U';const int lda=n;int info=0,lwork=-1;double query=0.0;
    dsyev_(&job,&uplo,&n,a.data(),&lda,out.values.data(),&query,&lwork,&info);if(info!=0)throw std::runtime_error("LAPACK dsyev workspace query failed");
    lwork=std::max(3*n-1,static_cast<int>(std::ceil(query)));std::vector<double> work(static_cast<std::size_t>(lwork));
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)a[static_cast<std::size_t>(j*n+i)]=0.5*(row_major[static_cast<std::size_t>(i*n+j)]+row_major[static_cast<std::size_t>(j*n+i)]);
    dsyev_(&job,&uplo,&n,a.data(),&lda,out.values.data(),work.data(),&lwork,&info);if(info!=0)throw std::runtime_error("LAPACK dsyev failed for geometric block");out.vectors_col_major=std::move(a);return out;
}

} // namespace

std::array<double,144> frame3d_global_stiffness(
    double xi,double yi,double zi,double xj,double yj,double zj,
    double E,double G,double A,double J,double Iy,double Iz,
    const std::array<double,3>& reference,double axial_compression){
    if(E<=0||G<=0||A<=0||J<=0||Iy<=0||Iz<=0||axial_compression<0)throw std::invalid_argument("invalid 3D frame property");
    std::array<double,3> dx{xj-xi,yj-yi,zj-zi};
    const double L=norm3(dx);if(L<=0)throw std::invalid_argument("zero-length 3D frame element");
    auto ex=normalize3(dx);
    const double proj=dot3(reference,ex);
    std::array<double,3> eyraw{reference[0]-proj*ex[0],reference[1]-proj*ex[1],reference[2]-proj*ex[2]};
    auto ey=normalize3(eyraw);
    auto ez=normalize3(cross3(ex,ey));
    // Re-orthogonalize ey to suppress accumulated floating error.
    ey=normalize3(cross3(ez,ex));
    std::array<double,9> R{ex[0],ex[1],ex[2],ey[0],ey[1],ey[2],ez[0],ez[1],ez[2]};

    std::array<double,144> k{};
    auto add=[&](int r,int c,double v){k[static_cast<std::size_t>(r*12+c)]+=v;};
    const double ka=E*A/L, kt=G*J/L;
    add(0,0,ka);add(0,6,-ka);add(6,0,-ka);add(6,6,ka);
    add(3,3,kt);add(3,9,-kt);add(9,3,-kt);add(9,9,kt);

    // Bending about local z: local v / rz.
    const double kz12=12*E*Iz/(L*L*L), kz6=6*E*Iz/(L*L), kz4=4*E*Iz/L, kz2=2*E*Iz/L;
    const int vz[4]{1,5,7,11};
    const double bz[4][4]={{kz12,kz6,-kz12,kz6},{kz6,kz4,-kz6,kz2},{-kz12,-kz6,kz12,-kz6},{kz6,kz2,-kz6,kz4}};
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)add(vz[a],vz[b],bz[a][b]);

    // Bending about local y: local w / ry. Sign follows right-hand rotation convention.
    const double ky12=12*E*Iy/(L*L*L), ky6=6*E*Iy/(L*L), ky4=4*E*Iy/L, ky2=2*E*Iy/L;
    const int wy[4]{2,4,8,10};
    const double by[4][4]={{ky12,-ky6,-ky12,-ky6},{-ky6,ky4,ky6,ky2},{-ky12,ky6,ky12,ky6},{-ky6,ky2,ky6,ky4}};
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)add(wy[a],wy[b],by[a][b]);

    if(axial_compression>0){
        const double q=axial_compression/(30.0*L);
        const double gz[4][4]={{36,3*L,-36,3*L},{3*L,4*L*L,-3*L,-L*L},{-36,-3*L,36,-3*L},{3*L,-L*L,-3*L,4*L*L}};
        const double gy[4][4]={{36,-3*L,-36,-3*L},{-3*L,4*L*L,3*L,-L*L},{-36,3*L,36,3*L},{-3*L,-L*L,3*L,4*L*L}};
        for(int a=0;a<4;++a)for(int b=0;b<4;++b){add(vz[a],vz[b],-q*gz[a][b]);add(wy[a],wy[b],-q*gy[a][b]);}
    }
    return transform_12(k,R);
}


void Frame3DBuilder::add_node(int id,double x,double y,double z,double mx,double my,double mz,double mrx,double mry,double mrz){
    const std::array<double,6> m{mx,my,mz,mrx,mry,mrz};
    for(double v:m) if(v<0) throw std::invalid_argument("negative nodal mass");
    if(std::any_of(nodes_.begin(),nodes_.end(),[&](const Node3D& n){return n.id==id;})) throw std::invalid_argument("duplicate node id");
    nodes_.push_back({id,x,y,z,m});
}

void Frame3DBuilder::add_elastic_frame(int id,int ni,int nj,double E,double G,double A,double J,double Iy,double Iz,double rx,double ry,double rz,double p){
    if(std::any_of(elements_.begin(),elements_.end(),[&](const ElasticFrame3D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 3D elastic element id");
    elements_.push_back({id,ni,nj,E,G,A,J,Iy,Iz,{rx,ry,rz},p});
}
void Frame3DBuilder::add_bilinear_spring(int id,int ni,Dof3D di,int nj,Dof3D dj,double k0,double fy,double ratio){
    add_linear_bilinear_spring(id,{{ni,di,-1.0},{nj,dj,1.0}},k0,fy,ratio);
}
void Frame3DBuilder::add_linear_bilinear_spring(int id,std::vector<MpcTerm3D> terms,double k0,double fy,double ratio){
    if(terms.empty()) throw std::invalid_argument("generalized spring requires deformation terms");
    springs_.push_back({id,std::move(terms),NonlinearMaterial(BilinearSpring(k0,fy,ratio))});
}

void Frame3DBuilder::add_linear_imk_peak_oriented_spring(int id,std::vector<MpcTerm3D> terms,IMKPeakOrientedParams params){
    if(terms.empty()) throw std::invalid_argument("generalized IMK spring requires deformation terms");
    springs_.push_back({id,std::move(terms),NonlinearMaterial(IMKPeakOrientedMaterial(params))});
}
void Frame3DBuilder::add_imk_peak_oriented_spring(int id,int ni,Dof3D di,int nj,Dof3D dj,IMKPeakOrientedParams params){
    add_linear_imk_peak_oriented_spring(id,{{ni,di,-1.0},{nj,dj,1.0}},params);
}
void Frame3DBuilder::add_linear_asce41_hinge(int id,std::vector<MpcTerm3D> terms,ASCE41HingeParams params){
    if(terms.empty()) throw std::invalid_argument("generalized ASCE41 hinge requires deformation terms");
    springs_.push_back({id,std::move(terms),NonlinearMaterial(ASCE41HingeMaterial(params))});
}
void Frame3DBuilder::add_linear_axial_coupled_asce41_hinge(
    int id,std::vector<MpcTerm3D> rotation_terms,std::vector<MpcTerm3D> axial_terms,
    AxialCoupledASCE41Params params){
    if(rotation_terms.empty()||axial_terms.empty())throw std::invalid_argument("coupled ASCE41 hinge requires rotation and axial terms");
    if(params.axial_stiffness<=0.0||params.axial_force_points.size()<2||
       params.axial_force_points.size()!=params.moment_capacity_points.size())
        throw std::invalid_argument("invalid coupled ASCE41 P-M parameters");
    for(double m:params.moment_capacity_points)if(!(m>0.0))throw std::invalid_argument("coupled ASCE41 moment capacities must be positive");
    // The base material provides initial state layout/stiffness.  Its fixed Fy
    // is replaced by the continuously evaluated capacity during global trials.
    const int spring_index=static_cast<int>(springs_.size());
    const double seed=params.moment_capacity_points.front();
    params.hinge.posFy=params.hinge.negFy=seed;
    springs_.push_back({id,std::move(rotation_terms),NonlinearMaterial(ASCE41HingeMaterial(params.hinge))});
    coupled_pm_.push_back({spring_index,std::move(axial_terms),std::move(params)});
}
void Frame3DBuilder::add_linear_pm_interaction_hinge(
    int id,std::vector<MpcTerm3D> rotation_terms,std::vector<MpcTerm3D> axial_terms,
    PMInteractionHingeParams params){
    if(rotation_terms.empty()||axial_terms.empty())throw std::invalid_argument("P-M interaction hinge requires rotation and axial terms");
    PMInteractionHinge2D law(params);
    const int spring_index=static_cast<int>(springs_.size());
    // Placeholder owns the ordinary rotational basis and initial Kr. The PM
    // compound law owns all constitutive state and both generalized forces.
    springs_.push_back({id,std::move(rotation_terms),NonlinearMaterial(BilinearSpring(params.hinge.Ke,1.0e30,0.0))});
    pm_return_.push_back({spring_index,std::move(axial_terms),std::move(params)});
}

void Frame3DBuilder::add_linear_fsc_shear_spring(
    int id,std::vector<MpcTerm3D> shear_terms,
    std::vector<MpcTerm3D> bottom_rotation_terms,
    std::vector<MpcTerm3D> top_rotation_terms,
    std::vector<MpcTerm3D> axial_terms,
    double axial_preload,double axial_stiffness,FSCShearSpringParams params){
    if(shear_terms.empty()||bottom_rotation_terms.empty()||top_rotation_terms.empty()||axial_terms.empty())
        throw std::invalid_argument("FSC shear spring requires shear, end-rotation, and axial terms");
    if(!(axial_stiffness>0.0))throw std::invalid_argument("FSC shear spring axial stiffness must be positive");
    // Keep the generic material bank topology intact. The placeholder state is
    // never used for force evaluation; the compiled FSC law owns the extra
    // state appended to this component. A huge yield force makes the fallback
    // evaluator harmless for diagnostic calls that do not supply global u.
    const int spring_index=static_cast<int>(springs_.size());
    springs_.push_back({id,std::move(shear_terms),NonlinearMaterial(BilinearSpring(params.Ke,1.0e30,0.0))});
    fsc_shear_.push_back({spring_index,std::move(bottom_rotation_terms),std::move(top_rotation_terms),
                          std::move(axial_terms),axial_preload,axial_stiffness,std::move(params)});
}
void Frame3DBuilder::add_asce41_hinge(int id,int ni,Dof3D di,int nj,Dof3D dj,ASCE41HingeParams params){
    add_linear_asce41_hinge(id,{{ni,di,-1.0},{nj,dj,1.0}},params);
}
void Frame3DBuilder::add_rotational_vector_spring(int id,int ni,int nj,double ax,double ay,double az,double k0,double fy,double ratio){
    const double norm=std::sqrt(ax*ax+ay*ay+az*az); if(norm<1e-14) throw std::invalid_argument("rotational spring axis is zero");
    ax/=norm; ay/=norm; az/=norm;
    add_linear_bilinear_spring(id,{{ni,Dof3D::RX,-ax},{ni,Dof3D::RY,-ay},{ni,Dof3D::RZ,-az},
                                   {nj,Dof3D::RX, ax},{nj,Dof3D::RY, ay},{nj,Dof3D::RZ, az}},k0,fy,ratio);
}
void Frame3DBuilder::fix_dof(int node,Dof3D dof){fixed_.push_back({node,dof});}
void Frame3DBuilder::fix(int node,bool ux,bool uy,bool uz,bool rx,bool ry,bool rz){
    const bool f[6]{ux,uy,uz,rx,ry,rz};
    for(int d=0;d<6;++d) if(f[d]) fix_dof(node,static_cast<Dof3D>(d));
}
void Frame3DBuilder::equal_dof(int master,int slave,Dof3D dof){equal_.push_back({master,slave,dof});}
void Frame3DBuilder::linear_constraint(int slave_node,Dof3D slave_dof,std::vector<MpcTerm3D> masters){
    if(masters.empty()) throw std::invalid_argument("linear MPC must have at least one master term");
    linear_.push_back({slave_node,slave_dof,std::move(masters)});
}
void Frame3DBuilder::rigid_diaphragm_z(int master_node,const std::vector<int>& slave_nodes){
    auto find_node=[&](int id)->const Node3D&{
        auto it=std::find_if(nodes_.begin(),nodes_.end(),[&](const Node3D& n){return n.id==id;});
        if(it==nodes_.end()) throw std::invalid_argument("unknown rigid-diaphragm node id "+std::to_string(id));
        return *it;
    };
    const auto& m=find_node(master_node);
    for(int sid:slave_nodes){
        if(sid==master_node) continue;
        const auto& s=find_node(sid);
        const double dx=s.x-m.x, dy=s.y-m.y;
        linear_constraint(sid,Dof3D::UX,{{master_node,Dof3D::UX,1.0},{master_node,Dof3D::RZ,-dy}});
        linear_constraint(sid,Dof3D::UY,{{master_node,Dof3D::UY,1.0},{master_node,Dof3D::RZ, dx}});
        linear_constraint(sid,Dof3D::RZ,{{master_node,Dof3D::RZ,1.0}});
    }
}
void Frame3DBuilder::set_rayleigh(double a,double b){if(a<0||b<0)throw std::invalid_argument("negative Rayleigh coefficient");alpha_m_=a;beta_k_=b;}
void Frame3DBuilder::set_ground_direction(Dof3D d){if(static_cast<int>(d)>2)throw std::invalid_argument("ground direction must be translational");ground_dof_=d;}
void Frame3DBuilder::set_response(int node,Dof3D dof){response_node_id_=node;response_dof_=dof;}
void Frame3DBuilder::set_story_nodes(std::vector<int> ids,Dof3D dof){story_node_ids_=std::move(ids);story_dof_=dof;}

CompiledFrame3D Frame3DBuilder::compile() const{
    if(nodes_.empty()) throw std::invalid_argument("3D frame has no nodes");
    if(corotational_ && updated_pdelta_) throw std::invalid_argument("corotational and updated_pdelta reference modes are mutually exclusive");
    CompiledFrame3D out;
    out.alpha_m_=alpha_m_; out.beta_k_=beta_k_; out.updated_pdelta_=updated_pdelta_; out.corotational_=corotational_; out.elastic_element_count_=static_cast<int>(elements_.size());
    out.node_ids_.reserve(nodes_.size());
    for(std::size_t i=0;i<nodes_.size();++i){out.node_ids_.push_back(nodes_[i].id);out.node_slot_.emplace(nodes_[i].id,static_cast<int>(i));}
    auto slot=[&](int id){auto it=out.node_slot_.find(id);if(it==out.node_slot_.end())throw std::invalid_argument("unknown 3D node id "+std::to_string(id));return it->second;};
    const int full_n=static_cast<int>(6*nodes_.size());

    // Build raw full-DOF constraints. A dependent DOF is expressed as a linear
    // combination of other full DOFs, then recursively expanded to independent q.
    std::vector<bool> fixed(static_cast<std::size_t>(full_n),false);
    for(const auto& f:fixed_) fixed[static_cast<std::size_t>(dof_index(slot(f.node),f.dof))]=true;
    std::vector<std::optional<std::vector<std::pair<int,double>>>> dep(static_cast<std::size_t>(full_n));
    auto set_dep=[&](int slave,std::vector<std::pair<int,double>> masters){
        if(fixed[static_cast<std::size_t>(slave)]) throw std::invalid_argument("DOF cannot be both fixed and dependent");
        if(dep[static_cast<std::size_t>(slave)].has_value()) throw std::invalid_argument("duplicate MPC/equalDOF slave constraint");
        dep[static_cast<std::size_t>(slave)]=std::move(masters);
    };
    for(const auto& e:equal_){
        const int master=dof_index(slot(e.master),e.dof), slave=dof_index(slot(e.slave),e.dof);
        if(master==slave) continue;
        set_dep(slave,{{master,1.0}});
    }
    for(const auto& c:linear_){
        const int slave=dof_index(slot(c.slave),c.dof);
        std::vector<std::pair<int,double>> masters; masters.reserve(c.masters.size());
        for(const auto& t:c.masters){
            if(std::abs(t.coefficient)<1e-16) continue;
            masters.emplace_back(dof_index(slot(t.node),t.dof),t.coefficient);
        }
        if(masters.empty()) throw std::invalid_argument("linear MPC collapsed to zero terms");
        set_dep(slave,std::move(masters));
    }

    std::vector<int> independent(static_cast<std::size_t>(full_n),-1);
    int next=0;
    for(int f=0;f<full_n;++f) if(!fixed[static_cast<std::size_t>(f)] && !dep[static_cast<std::size_t>(f)].has_value()) independent[static_cast<std::size_t>(f)]=next++;
    out.n_=next;
    if(out.n_==0) throw std::invalid_argument("3D frame has no active DOFs");
    out.full_to_terms_.resize(static_cast<std::size_t>(full_n));
    out.full_to_reduced_.assign(static_cast<std::size_t>(full_n),-1);
    std::vector<unsigned char> state(static_cast<std::size_t>(full_n),0); // 0 unseen,1 visiting,2 complete
    std::function<const ConstraintTerms&(int)> expand=[&](int f)->const ConstraintTerms&{
        auto& st=state[static_cast<std::size_t>(f)];
        if(st==2) return out.full_to_terms_[static_cast<std::size_t>(f)];
        if(st==1) throw std::invalid_argument("cyclic multi-point constraint detected");
        st=1;
        std::map<int,double> coeff;
        if(fixed[static_cast<std::size_t>(f)]) {
            // zero row in T
        } else if(independent[static_cast<std::size_t>(f)]>=0) {
            coeff[independent[static_cast<std::size_t>(f)]]=1.0;
        } else {
            for(const auto& [master,c] : *dep[static_cast<std::size_t>(f)]) accumulate_terms(coeff,expand(master),c);
        }
        out.full_to_terms_[static_cast<std::size_t>(f)]=compact_terms(coeff);
        const auto& terms=out.full_to_terms_[static_cast<std::size_t>(f)];
        if(terms.size()==1 && std::abs(terms[0].second-1.0)<1e-14) out.full_to_reduced_[static_cast<std::size_t>(f)]=terms[0].first;
        st=2;
        return out.full_to_terms_[static_cast<std::size_t>(f)];
    };
    for(int f=0;f<full_n;++f) (void)expand(f);

    // Representative elevation of each generalized DOF. Rigid-diaphragm and
    // same-floor MPC rows remain on one elevation; cross-story generalized
    // coordinates are marked NaN and retained as substructure interfaces.
    out.reduced_dof_z_.assign(static_cast<std::size_t>(out.n_), std::numeric_limits<double>::quiet_NaN());
    std::vector<bool> z_seen(static_cast<std::size_t>(out.n_),false), z_mixed(static_cast<std::size_t>(out.n_),false);
    for(std::size_t ni=0;ni<nodes_.size();++ni){
        const double z=nodes_[ni].z;
        for(int d=0;d<6;++d){
            for(const auto& [r,c]:out.full_to_terms_[6*ni+static_cast<std::size_t>(d)]){
                if(std::abs(c)<1e-14) continue;
                if(!z_seen[static_cast<std::size_t>(r)]){out.reduced_dof_z_[static_cast<std::size_t>(r)]=z;z_seen[static_cast<std::size_t>(r)]=true;}
                else if(std::abs(out.reduced_dof_z_[static_cast<std::size_t>(r)]-z)>1e-9) z_mixed[static_cast<std::size_t>(r)]=true;
            }
        }
    }
    for(int r=0;r<out.n_;++r) if(z_mixed[static_cast<std::size_t>(r)]) out.reduced_dof_z_[static_cast<std::size_t>(r)]=std::numeric_limits<double>::quiet_NaN();

    // Condensed generalized mass Mq = T^T M T and influence vector T^T M r.
    std::vector<Triplet> mt;
    out.base_mass_.assign(static_cast<std::size_t>(out.n_),0.0);
    for(std::size_t ni=0;ni<nodes_.size();++ni){
        for(int d=0;d<6;++d){
            const double mass=nodes_[ni].mass[static_cast<std::size_t>(d)]; if(mass==0.0) continue;
            const auto& terms=out.full_to_terms_[6*ni+static_cast<std::size_t>(d)];
            for(const auto& [ra,ca]:terms) for(const auto& [rb,cb]:terms) mt.push_back({ra,rb,mass*ca*cb});
            if(d==static_cast<int>(ground_dof_)) for(const auto& [r,c]:terms) out.base_mass_[static_cast<std::size_t>(r)] += mass*c;
        }
    }
    // Preserve generalized diagonals even if zero so model diagnostics remain stable.
    for(int i=0;i<out.n_;++i) mt.push_back({i,i,0.0});
    out.mass_matrix_=SparseMatrixCSC::from_triplets(out.n_,out.n_,mt,-1.0);
    out.mass_=out.mass_matrix_.diagonal();

    // Assemble Kq = T^T Kfull T element by element.
    std::vector<Triplet> linear;
    out.element_response_data_.reserve(elements_.size());
    for(const auto& e:elements_){
        const int si=slot(e.node_i),sj=slot(e.node_j);
        const auto& a=nodes_[static_cast<std::size_t>(si)]; const auto& b=nodes_[static_cast<std::size_t>(sj)];
        CorotationalFrame3DProperties cp{a.x,a.y,a.z,b.x,b.y,b.z,e.E,e.G,e.A,e.J,e.Iy,e.Iz,e.reference,e.axial_compression};
        std::array<double,144> kg{},kel{},kunit{};
        if(corotational_){
            const std::array<double,12> z{};kg=corotational3d_response(cp,z).tangent;kel=kg;kunit=kg;
        }else{
            kg=frame3d_global_stiffness(a.x,a.y,a.z,b.x,b.y,b.z,e.E,e.G,e.A,e.J,e.Iy,e.Iz,e.reference,e.axial_compression);
            kel=updated_pdelta_ ? frame3d_global_stiffness(a.x,a.y,a.z,b.x,b.y,b.z,e.E,e.G,e.A,e.J,e.Iy,e.Iz,e.reference,0.0) : kg;
            kunit=updated_pdelta_ ? frame3d_global_stiffness(a.x,a.y,a.z,b.x,b.y,b.z,e.E,e.G,e.A,e.J,e.Iy,e.Iz,e.reference,1.0) : kg;
        }
        int fd[12]; for(int d=0;d<6;++d){fd[d]=6*si+d;fd[6+d]=6*sj+d;}
        CompiledFrame3D::ElementResponseData rd; rd.properties=e;
        rd.node_i_xyz={a.x,a.y,a.z}; rd.node_j_xyz={b.x,b.y,b.z};
        for(int d=0;d<12;++d) rd.terms[static_cast<std::size_t>(d)]=out.full_to_terms_[static_cast<std::size_t>(fd[d])];
        if(out.element_index_by_id_.count(e.id)) throw std::invalid_argument("duplicate 3D elastic element id");
        out.element_index_by_id_[e.id]=static_cast<int>(out.element_response_data_.size());
        out.element_response_data_.push_back(std::move(rd));
        if(corotational_){
            CompiledFrame3D::CorotElement ce;ce.properties=cp;
            for(int d=0;d<12;++d)ce.terms[static_cast<std::size_t>(d)]=out.full_to_terms_[static_cast<std::size_t>(fd[d])];
            out.corotational_elements_.push_back(std::move(ce));
        }
        for(int i=0;i<12;++i){
            const auto& ti=out.full_to_terms_[static_cast<std::size_t>(fd[i])]; if(ti.empty()) continue;
            for(int j=0;j<12;++j){
                const double kij=kg[static_cast<std::size_t>(i*12+j)]; if(kij==0.0) continue;
                const auto& tj=out.full_to_terms_[static_cast<std::size_t>(fd[j])];
                for(const auto& [ri,ci]:ti) for(const auto& [rj,cj]:tj) linear.push_back({ri,rj,ci*kij*cj});
            }
        }
        if(updated_pdelta_){
            CompiledFrame3D::GeometricUpdate gu; gu.preload=e.axial_compression;
            const std::array<double,3> dx{b.x-a.x,b.y-a.y,b.z-a.z};const double L=norm3(dx);const auto ex=normalize3(dx);gu.axial_k=e.E*e.A/L;
            std::map<int,double> amap;
            for(int d=0;d<3;++d){accumulate_terms(amap,out.full_to_terms_[static_cast<std::size_t>(6*si+d)],-ex[static_cast<std::size_t>(d)]);accumulate_terms(amap,out.full_to_terms_[static_cast<std::size_t>(6*sj+d)],ex[static_cast<std::size_t>(d)]);}
            gu.axial_terms=compact_terms(amap);
            for(int i=0;i<12;++i){const auto& ti=out.full_to_terms_[static_cast<std::size_t>(fd[i])];if(ti.empty())continue;for(int j=0;j<12;++j){
                const double v=kunit[static_cast<std::size_t>(i*12+j)]-kel[static_cast<std::size_t>(i*12+j)];if(std::abs(v)<1e-18)continue;const auto& tj=out.full_to_terms_[static_cast<std::size_t>(fd[j])];
                for(const auto& [ri,ci]:ti)for(const auto& [rj,cj]:tj)gu.entries.push_back({ri,rj,-1,ci*v*cj});
            }}
            std::vector<int> grows;grows.reserve(gu.entries.size());for(const auto& ge:gu.entries)grows.push_back(ge.row);
            std::sort(grows.begin(),grows.end());grows.erase(std::unique(grows.begin(),grows.end()),grows.end());
            for(int r:grows)for(const auto& [c,ac]:gu.axial_terms)gu.coupling.push_back({r,c,-1,-1,ac});
            out.geometric_updates_.push_back(std::move(gu));
        }
    }
    out.k_linear_=SparseMatrixCSC::from_triplets(out.n_,out.n_,linear,1e-18);

    // Compile nonlinear generalized deformation directions bq = T^T bfull.
    const int m=static_cast<int>(springs_.size());
    std::vector<unsigned char> coupled_state_extra(static_cast<std::size_t>(m),0);
    std::vector<int> pm_return_state_extra(static_cast<std::size_t>(m),0);
    std::vector<int> fsc_state_extra(static_cast<std::size_t>(m),0);
    for(const auto& raw:coupled_pm_){
        if(raw.spring_index<0||raw.spring_index>=m)throw std::runtime_error("coupled P-M spring index");
        if(coupled_state_extra[static_cast<std::size_t>(raw.spring_index)]!=0)throw std::invalid_argument("duplicate coupled P-M spring");
        coupled_state_extra[static_cast<std::size_t>(raw.spring_index)]=1;
    }
    for(const auto& raw:pm_return_){
        if(raw.spring_index<0||raw.spring_index>=m)throw std::runtime_error("P-M return hinge spring index");
        if(coupled_state_extra[static_cast<std::size_t>(raw.spring_index)]!=0||pm_return_state_extra[static_cast<std::size_t>(raw.spring_index)]!=0)
            throw std::invalid_argument("duplicate/incompatible P-M compound spring");
        pm_return_state_extra[static_cast<std::size_t>(raw.spring_index)]=PMInteractionHinge2D::kStateSize;
    }
    for(const auto& raw:fsc_shear_){
        if(raw.spring_index<0||raw.spring_index>=m)throw std::runtime_error("FSC shear spring index");
        if(fsc_state_extra[static_cast<std::size_t>(raw.spring_index)]!=0)throw std::invalid_argument("duplicate FSC shear spring");
        if(coupled_state_extra[static_cast<std::size_t>(raw.spring_index)]!=0||pm_return_state_extra[static_cast<std::size_t>(raw.spring_index)]!=0)throw std::invalid_argument("component cannot be both P-M compound and FSC shear");
        fsc_state_extra[static_cast<std::size_t>(raw.spring_index)]=FSCShearSpringLaw::kStateSize;
    }
    std::vector<std::vector<std::pair<int,double>>> basis_columns(static_cast<std::size_t>(m)); out.materials_.reserve(springs_.size());
    out.material_state_offsets_.resize(static_cast<std::size_t>(m+1),0);
    for(int j=0;j<m;++j){
        const auto& sp=springs_[static_cast<std::size_t>(j)];
        std::map<int,double> bmap;
        for(const auto& term:sp.terms)
            accumulate_terms(bmap,out.full_to_terms_[static_cast<std::size_t>(dof_index(slot(term.node),term.dof))],term.coefficient);
        basis_columns[static_cast<std::size_t>(j)]=compact_terms(bmap);
        if(basis_columns[static_cast<std::size_t>(j)].empty()) throw std::invalid_argument("3D generalized spring has constrained zero deformation");
        out.materials_.push_back(sp.material);
        if(out.component_index_by_id_.count(sp.id))throw std::invalid_argument("duplicate nonlinear spring id");
        out.component_index_by_id_[sp.id]=j;
        out.component_ids_.push_back(sp.id);
        out.initial_tangents_.push_back(sp.material.initial_stiffness());
        out.material_state_offsets_[static_cast<std::size_t>(j+1)] =
            out.material_state_offsets_[static_cast<std::size_t>(j)] + sp.material.state_size()
            + static_cast<int>(coupled_state_extra[static_cast<std::size_t>(j)])
            + pm_return_state_extra[static_cast<std::size_t>(j)]
            + fsc_state_extra[static_cast<std::size_t>(j)];
    }
    out.nonlinear_state_size_ = out.material_state_offsets_.empty() ? 0 : out.material_state_offsets_.back();
    out.basis_=SparseUpdateBasis::from_columns(out.n_,basis_columns,1e-14);

    // Compile auxiliary axial coordinates for the small subset of rotational
    // hinges whose current trial strength depends on member compression.
    out.coupled_pm_by_spring_.assign(static_cast<std::size_t>(m),-1);
    for(const auto& raw:coupled_pm_){
        if(raw.spring_index<0||raw.spring_index>=m)throw std::runtime_error("coupled P-M spring index");
        CompiledFrame3D::CoupledPMHinge cp;cp.spring_index=raw.spring_index;cp.axial_preload=raw.params.axial_preload;
        cp.axial_stiffness=raw.params.axial_stiffness;cp.base_params=raw.params.hinge;
        cp.axial_force_points=raw.params.axial_force_points;cp.moment_capacity_points=raw.params.moment_capacity_points;
        cp.capacity_slopes=monotone_capacity_slopes(cp.axial_force_points,cp.moment_capacity_points);
        std::map<int,double> amap;
        for(const auto& term:raw.axial_terms)
            accumulate_terms(amap,out.full_to_terms_[static_cast<std::size_t>(dof_index(slot(term.node),term.dof))],term.coefficient);
        cp.axial_terms=compact_terms(amap);
        if(cp.axial_terms.empty())throw std::invalid_argument("coupled P-M hinge has constrained zero axial deformation");
        out.coupled_pm_by_spring_[static_cast<std::size_t>(raw.spring_index)]=static_cast<int>(out.coupled_pm_hinges_.size());
        out.coupled_pm_hinges_.push_back(std::move(cp));
    }

    // Compile Phase-9F associative two-force P-M hinges. Their axial terms
    // are a true hinge deformation coordinate and therefore contribute force
    // and stiffness, unlike the legacy member-P monitoring wrapper above.
    out.pm_return_by_spring_.assign(static_cast<std::size_t>(m),-1);
    for(const auto& raw:pm_return_){
        CompiledFrame3D::PMReturnHinge ph;ph.spring_index=raw.spring_index;ph.law.emplace(raw.params);
        std::map<int,double> amap;
        for(const auto& term:raw.axial_terms)
            accumulate_terms(amap,out.full_to_terms_[static_cast<std::size_t>(dof_index(slot(term.node),term.dof))],term.coefficient);
        ph.axial_terms=compact_terms(amap);
        if(ph.axial_terms.empty())throw std::invalid_argument("P-M return hinge has constrained zero axial deformation");
        ph.state_offset=out.material_state_offsets_[static_cast<std::size_t>(raw.spring_index)]
                       + out.materials_[static_cast<std::size_t>(raw.spring_index)].state_size();
        out.pm_return_by_spring_[static_cast<std::size_t>(raw.spring_index)]=static_cast<int>(out.pm_return_hinges_.size());
        out.pm_return_hinges_.push_back(std::move(ph));
    }

    // Compile remote end-rotation and axial coordinates for flexure-shear-
    // critical column springs. Their force direction remains the ordinary
    // scalar shear basis column; the remote coordinates only govern the
    // committed failure state and therefore do not add a global tangent block.
    out.fsc_shear_by_spring_.assign(static_cast<std::size_t>(m),-1);
    for(const auto& raw:fsc_shear_){
        CompiledFrame3D::FSCShearSpring fs{raw.spring_index,-1,raw.axial_preload,raw.axial_stiffness,
                                             FSCShearSpringLaw(raw.params),{},{},{}};
        auto compile_aux=[&](const std::vector<MpcTerm3D>& src){
            std::map<int,double> amap;
            for(const auto& term:src)
                accumulate_terms(amap,out.full_to_terms_[static_cast<std::size_t>(dof_index(slot(term.node),term.dof))],term.coefficient);
            return compact_terms(amap);
        };
        fs.bottom_rotation_terms=compile_aux(raw.bottom_rotation_terms);
        fs.top_rotation_terms=compile_aux(raw.top_rotation_terms);
        fs.axial_terms=compile_aux(raw.axial_terms);
        if(fs.bottom_rotation_terms.empty()||fs.top_rotation_terms.empty()||fs.axial_terms.empty())
            throw std::invalid_argument("FSC shear spring has constrained zero auxiliary coordinate");
        fs.state_offset=out.material_state_offsets_[static_cast<std::size_t>(raw.spring_index)]
                       + out.materials_[static_cast<std::size_t>(raw.spring_index)].state_size();
        out.fsc_shear_by_spring_[static_cast<std::size_t>(raw.spring_index)]=static_cast<int>(out.fsc_shear_springs_.size());
        out.fsc_shear_springs_.push_back(std::move(fs));
    }

    // Compile an exact fixed sparse basis for localized state-dependent
    // geometric tangent changes. Hinge columns come first. Each geometric
    // block adds the nonzero eigenvectors of its fixed unit geometric matrix
    // plus the axial-deformation direction a. The changing tangent can then be
    // written exactly as U C(u) U^T without changing U during the time loop.
    if(updated_pdelta_){
        std::vector<std::vector<std::pair<int,double>>> gcols=basis_columns;
        for(const auto& gu:out.geometric_updates_){
            std::set<int> support;for(const auto& e:gu.entries){support.insert(e.row);support.insert(e.col);}for(const auto& [r,c]:gu.axial_terms)if(std::abs(c)>1e-14)support.insert(r);
            if(support.empty()) continue;
            std::vector<int> ids(support.begin(),support.end());
            std::map<int,int> loc;
            for(std::size_t i=0;i<ids.size();++i) loc[ids[i]]=static_cast<int>(i);
            const int ns=static_cast<int>(ids.size());std::vector<double> Gs(static_cast<std::size_t>(ns*ns),0.0);
            for(const auto& e:gu.entries)Gs[static_cast<std::size_t>(loc[e.row]*ns+loc[e.col])]+=e.value;
            auto eig=small_symmetric_eigen(Gs,ns);double maxeig=0.0;for(double v:eig.values)maxeig=std::max(maxeig,std::abs(v));
            CompiledFrame3D::GeometricLowRankBlock block;block.start_col=static_cast<int>(gcols.size());block.axial_k=gu.axial_k;
            const double etol=std::max(1e-13,1e-10*maxeig);
            for(int k=0;k<ns;++k){const double lam=eig.values[static_cast<std::size_t>(k)];if(std::abs(lam)<=etol)continue;std::vector<std::pair<int,double>> col;for(int i=0;i<ns;++i){const double v=eig.vectors_col_major[static_cast<std::size_t>(k*ns+i)];if(std::abs(v)>1e-12)col.push_back({ids[static_cast<std::size_t>(i)],v});}gcols.push_back(std::move(col));block.eigenvalues.push_back(lam);++block.eigen_count;}
            double an=0.0;for(const auto& [r,c]:gu.axial_terms)an+=c*c;if(an>1e-24){block.axial_col=static_cast<int>(gcols.size());gcols.push_back(gu.axial_terms);}
            if(block.eigen_count>0 && block.axial_col>=0)out.geometric_low_rank_blocks_.push_back(std::move(block));
        }
        out.generalized_state_update_basis_=SparseUpdateBasis::from_columns(out.n_,gcols,1e-14);
    }
    std::vector<Triplet> initial=linear;
    for(int j=0;j<m;++j) add_scaled_bbt3(initial,out.basis_,j,out.initial_tangents_[static_cast<std::size_t>(j)]);
    // Rayleigh stiffness damping follows the physical/pre-existing model
    // stiffness. The additional axial Ka below is a rigid-plastic numerical
    // penalty coordinate and must not participate in beta*K damping.
    {auto damping=initial;for(int i=0;i<out.n_;++i)damping.push_back({i,i,0.0});out.k_damping_reference_=SparseMatrixCSC::from_triplets(out.n_,out.n_,damping,-1.0);}
    for(const auto& ph:out.pm_return_hinges_){
        const double Ka=ph.law->axial_initial_stiffness();
        for(const auto& [ri,ai]:ph.axial_terms)for(const auto& [rj,aj]:ph.axial_terms)
            initial.push_back({ri,rj,Ka*ai*aj});
    }
    for(int i=0;i<out.n_;++i) initial.push_back({i,i,0.0});
    out.k_initial_=SparseMatrixCSC::from_triplets(out.n_,out.n_,initial,-1.0);
    out.spring_initial_scatter_.resize(static_cast<std::size_t>(m));
    for(int j=0;j<m;++j){
        auto& entries=out.spring_initial_scatter_[static_cast<std::size_t>(j)];
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(j)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(j+1)];
        for(int pc=pb;pc<pe;++pc){
            const int col=out.basis_.row_ind()[static_cast<std::size_t>(pc)];const double bj=out.basis_.values()[static_cast<std::size_t>(pc)];
            for(int pr=pb;pr<pe;++pr){
                const int row=out.basis_.row_ind()[static_cast<std::size_t>(pr)];const double bi=out.basis_.values()[static_cast<std::size_t>(pr)];
                const int q=find_csc_position(out.k_initial_,row,col);if(q<0)throw std::runtime_error("missing 3D initial nonlinear scatter");
                entries.emplace_back(q,bi*bj);
            }
        }
    }
    out.spring_damping_scatter_.resize(static_cast<std::size_t>(m));
    for(int j=0;j<m;++j){
        auto& entries=out.spring_damping_scatter_[static_cast<std::size_t>(j)];
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(j)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(j+1)];
        for(int pc=pb;pc<pe;++pc){const int col=out.basis_.row_ind()[static_cast<std::size_t>(pc)];const double bj=out.basis_.values()[static_cast<std::size_t>(pc)];
            for(int pr=pb;pr<pe;++pr){const int row=out.basis_.row_ind()[static_cast<std::size_t>(pr)];const double bi=out.basis_.values()[static_cast<std::size_t>(pr)];
                const int q=find_csc_position(out.k_damping_reference_,row,col);if(q<0)throw std::runtime_error("missing 3D damping-reference nonlinear scatter");entries.emplace_back(q,bi*bj);}}
    }

    // One fixed CSC pattern for every effective tangent: union(Kinitial,M).
    std::vector<Triplet> pattern;
    for(int c=0;c<out.k_initial_.cols();++c) for(int p=out.k_initial_.col_ptr()[static_cast<std::size_t>(c)];p<out.k_initial_.col_ptr()[static_cast<std::size_t>(c+1)];++p)
        pattern.push_back({out.k_initial_.row_ind()[static_cast<std::size_t>(p)],c,0.0});
    for(int c=0;c<out.mass_matrix_.cols();++c) for(int p=out.mass_matrix_.col_ptr()[static_cast<std::size_t>(c)];p<out.mass_matrix_.col_ptr()[static_cast<std::size_t>(c+1)];++p)
        pattern.push_back({out.mass_matrix_.row_ind()[static_cast<std::size_t>(p)],c,0.0});
    for(const auto& gu:out.geometric_updates_){
        for(const auto& e:gu.entries) pattern.push_back({e.row,e.col,0.0});
        for(const auto& e:gu.coupling){pattern.push_back({e.row,e.col,0.0});pattern.push_back({e.col,e.row,0.0});}
    }
    // P->M coupling contributes b_theta (dM/dP)(-EA/L a^T).
    // It is generally nonsymmetric, so only the actual row/column direction
    // is inserted into the fixed sparse pattern.
    for(const auto& cp:out.coupled_pm_hinges_){
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(cp.spring_index)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(cp.spring_index+1)];
        for(int pr=pb;pr<pe;++pr){const int row=out.basis_.row_ind()[static_cast<std::size_t>(pr)];
            for(const auto& [col,ac]:cp.axial_terms)if(std::abs(ac)>1e-18)pattern.push_back({row,col,0.0});}
    }
    // Full 2x2 tangent support for associative P-M hinges.
    for(const auto& ph:out.pm_return_hinges_){
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(ph.spring_index)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(ph.spring_index+1)];
        for(const auto& [ri,ai]:ph.axial_terms)for(const auto& [rj,aj]:ph.axial_terms)if(std::abs(ai*aj)>1e-18)pattern.push_back({ri,rj,0.0});
        for(int pr=pb;pr<pe;++pr){const int rr=out.basis_.row_ind()[static_cast<std::size_t>(pr)];
            for(const auto& [ra,aa]:ph.axial_terms){pattern.push_back({rr,ra,0.0});pattern.push_back({ra,rr,0.0});}}
    }
    for(const auto& ce:out.corotational_elements_){
        for(int i=0;i<12;++i)for(int j=0;j<12;++j)
            for(const auto& [ri,ci]:ce.terms[static_cast<std::size_t>(i)])
                for(const auto& [rj,cj]:ce.terms[static_cast<std::size_t>(j)])
                    if(std::abs(ci*cj)>1e-18) pattern.push_back({ri,rj,0.0});
    }
    out.system_pattern_=SparseMatrixCSC::from_triplets(out.n_,out.n_,pattern,-1.0);
    out.k_to_system_.resize(static_cast<std::size_t>(out.k_initial_.nnz()));
    for(int c=0;c<out.k_initial_.cols();++c) for(int p=out.k_initial_.col_ptr()[static_cast<std::size_t>(c)];p<out.k_initial_.col_ptr()[static_cast<std::size_t>(c+1)];++p){
        const int q=find_csc_position(out.system_pattern_,out.k_initial_.row_ind()[static_cast<std::size_t>(p)],c); if(q<0)throw std::runtime_error("missing K-to-system scatter"); out.k_to_system_[static_cast<std::size_t>(p)]=q;
    }
    out.kd_to_system_.resize(static_cast<std::size_t>(out.k_damping_reference_.nnz()));
    for(int c=0;c<out.k_damping_reference_.cols();++c) for(int p=out.k_damping_reference_.col_ptr()[static_cast<std::size_t>(c)];p<out.k_damping_reference_.col_ptr()[static_cast<std::size_t>(c+1)];++p){
        const int q=find_csc_position(out.system_pattern_,out.k_damping_reference_.row_ind()[static_cast<std::size_t>(p)],c);if(q<0)throw std::runtime_error("missing damping-K-to-system scatter");out.kd_to_system_[static_cast<std::size_t>(p)]=q;
    }
    out.m_to_system_.resize(static_cast<std::size_t>(out.mass_matrix_.nnz()));
    for(int c=0;c<out.mass_matrix_.cols();++c) for(int p=out.mass_matrix_.col_ptr()[static_cast<std::size_t>(c)];p<out.mass_matrix_.col_ptr()[static_cast<std::size_t>(c+1)];++p){
        const int q=find_csc_position(out.system_pattern_,out.mass_matrix_.row_ind()[static_cast<std::size_t>(p)],c); if(q<0)throw std::runtime_error("missing M-to-system scatter"); out.m_to_system_[static_cast<std::size_t>(p)]=q;
    }
    out.spring_scatter_.resize(static_cast<std::size_t>(m));
    for(int j=0;j<m;++j){
        auto& entries=out.spring_scatter_[static_cast<std::size_t>(j)];
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(j)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(j+1)];
        for(int pc=pb;pc<pe;++pc){const int col=out.basis_.row_ind()[static_cast<std::size_t>(pc)];const double bj=out.basis_.values()[static_cast<std::size_t>(pc)];
            for(int pr=pb;pr<pe;++pr){const int row=out.basis_.row_ind()[static_cast<std::size_t>(pr)];const double bi=out.basis_.values()[static_cast<std::size_t>(pr)];
                const int q=find_csc_position(out.system_pattern_,row,col);if(q<0)throw std::runtime_error("missing 3D nonlinear scatter");entries.emplace_back(q,bi*bj);}}
    }

    for(auto& cp:out.coupled_pm_hinges_){
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(cp.spring_index)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(cp.spring_index+1)];
        for(int pr=pb;pr<pe;++pr){const int row=out.basis_.row_ind()[static_cast<std::size_t>(pr)];const double br=out.basis_.values()[static_cast<std::size_t>(pr)];
            for(const auto& [col,ac]:cp.axial_terms){const int q=find_csc_position(out.system_pattern_,row,col);if(q<0)throw std::runtime_error("missing coupled P-M scatter");cp.cross_scatter.push_back({q,br*ac});}}
    }
    for(auto& ph:out.pm_return_hinges_){
        // Local coordinate 0=axial, 1=rotation. Rot-rot is already handled by
        // the ordinary spring tangent; store axial-axial and both cross terms.
        for(const auto& [ri,ai]:ph.axial_terms)for(const auto& [rj,aj]:ph.axial_terms){
            const int q=find_csc_position(out.system_pattern_,ri,rj);if(q<0)throw std::runtime_error("missing P-M axial scatter");
            ph.block_scatter.push_back({q,0,0,ai*aj});
        }
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(ph.spring_index)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(ph.spring_index+1)];
        for(int pr=pb;pr<pe;++pr){const int rr=out.basis_.row_ind()[static_cast<std::size_t>(pr)];const double br=out.basis_.values()[static_cast<std::size_t>(pr)];
            for(const auto& [ra,aa]:ph.axial_terms){
                int q=find_csc_position(out.system_pattern_,ra,rr);if(q<0)throw std::runtime_error("missing P-M N-theta scatter");ph.block_scatter.push_back({q,0,1,aa*br});
                q=find_csc_position(out.system_pattern_,rr,ra);if(q<0)throw std::runtime_error("missing P-M M-axial scatter");ph.block_scatter.push_back({q,1,0,br*aa});
            }
        }
    }

    for(auto& gu:out.geometric_updates_){
        for(auto& e:gu.entries){e.system_pos=find_csc_position(out.system_pattern_,e.row,e.col);if(e.system_pos<0)throw std::runtime_error("missing geometric scatter");}
        for(auto& e:gu.coupling){
            e.system_pos=find_csc_position(out.system_pattern_,e.row,e.col);
            e.transpose_pos=find_csc_position(out.system_pattern_,e.col,e.row);
            if(e.system_pos<0||e.transpose_pos<0)throw std::runtime_error("missing geometric coupling scatter");
        }
    }
    for(auto& ce:out.corotational_elements_){
        for(int i=0;i<12;++i)for(int j=0;j<12;++j){
            for(const auto& [ri,ci]:ce.terms[static_cast<std::size_t>(i)])for(const auto& [rj,cj]:ce.terms[static_cast<std::size_t>(j)]){
                const double c=ci*cj;if(std::abs(c)<1e-18)continue;const int q=find_csc_position(out.system_pattern_,ri,rj);if(q<0)throw std::runtime_error("missing corotational scatter");ce.scatter.push_back({i,j,q,c});
            }
        }
    }

    const int rid=response_node_id_>=0?response_node_id_:nodes_.back().id;
    out.response_terms_=out.full_to_terms_[static_cast<std::size_t>(dof_index(slot(rid),response_dof_))];
    if(out.response_terms_.empty()) throw std::invalid_argument("3D response DOF fixed");
    auto stories=story_node_ids_; if(stories.empty()) stories.push_back(rid);
    for(int id:stories){auto terms=out.full_to_terms_[static_cast<std::size_t>(dof_index(slot(id),story_dof_))];if(terms.empty())throw std::invalid_argument("3D story response DOF fixed");out.story_terms_.push_back(std::move(terms));out.story_z_.push_back(nodes_[static_cast<std::size_t>(slot(id))].z);}
    return out;
}

int CompiledFrame3D::reduced_dof(int node,Dof3D dof) const{
    auto it=node_slot_.find(node); if(it==node_slot_.end()) return -1;
    return full_to_reduced_[static_cast<std::size_t>(dof_index(it->second,dof))];
}
int CompiledFrame3D::nonlinear_component_index(int spring_id) const{
    auto it=component_index_by_id_.find(spring_id);return it==component_index_by_id_.end()?-1:it->second;
}
int CompiledFrame3D::nonlinear_component_id(int index) const{
    if(index<0||index>=nonlinear_count()) return -1;
    return component_ids_[static_cast<std::size_t>(index)];
}

NonlinearComponentSnapshot CompiledFrame3D::nonlinear_component_snapshot(
    int spring_id,const std::vector<double>& u,const std::vector<double>& committed_state) const {
    if(static_cast<int>(u.size())!=n_||static_cast<int>(committed_state.size())!=nonlinear_state_size_)
        throw std::invalid_argument("3D nonlinear component snapshot state size");
    const int idx=nonlinear_component_index(spring_id);
    if(idx<0) throw std::invalid_argument("unknown 3D nonlinear component id");
    const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(idx)];
    const int ci=coupled_pm_by_spring_.empty()?-1:coupled_pm_by_spring_[static_cast<std::size_t>(idx)];
    const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(idx)];
    if(fi>=0||ci>=0||pi>=0) throw std::invalid_argument("compound 3D nonlinear component snapshot requires dedicated recorder");
    const double q=basis_.column_dot(idx,u);
    const int o=material_state_offsets_[static_cast<std::size_t>(idx)];
    std::vector<double> trial(static_cast<std::size_t>(materials_[static_cast<std::size_t>(idx)].state_size()));
    const auto tr=materials_[static_cast<std::size_t>(idx)].trial(q,committed_state.data()+o,trial.data());
    return {spring_id,idx,q,tr.force,tr.tangent,tr.diagnostics};
}

PMInteractionTrialResult CompiledFrame3D::pm_interaction_snapshot(
    int spring_id,const std::vector<double>& u,const std::vector<double>& committed_state) const {
    if(static_cast<int>(u.size())!=n_||static_cast<int>(committed_state.size())!=nonlinear_state_size_)
        throw std::invalid_argument("3D P-M snapshot state size");
    const int idx=nonlinear_component_index(spring_id);if(idx<0)throw std::invalid_argument("unknown P-M component id");
    const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(idx)];
    if(pi<0)throw std::invalid_argument("component is not return-mapped P-M hinge");
    const int o=material_state_offsets_[static_cast<std::size_t>(idx)];const double q=basis_.column_dot(idx,u);
    std::vector<double> scratch(static_cast<std::size_t>(materials_[static_cast<std::size_t>(idx)].state_size()+PMInteractionHinge2D::kStateSize),0.0);
    return evaluate_pm_return(idx,q,u,committed_state.data()+o,scratch.data());
}

void CompiledFrame3D::replace_nonlinear_material(int spring_id,NonlinearMaterial material){
    replace_nonlinear_material_field({{spring_id,std::move(material)}});
}

void CompiledFrame3D::replace_nonlinear_material_field(
    const std::vector<std::pair<int,NonlinearMaterial>>& replacements){
    if(replacements.empty()) return;
    struct Pending{int index{};int component_id{};NonlinearMaterial material{BilinearSpring(1.0,1.0,0.0)};double old_k{},new_k{};};
    std::vector<Pending> pending;pending.reserve(replacements.size());
    std::set<int> seen;
    for(const auto& [id,material]:replacements){
        const int idx=nonlinear_component_index(id);
        if(idx<0) throw std::invalid_argument("unknown 3D nonlinear component id in material field: "+std::to_string(id));
        if(!seen.insert(idx).second) throw std::invalid_argument("duplicate 3D nonlinear component replacement: "+std::to_string(id));
        if(material.state_size()!=materials_[static_cast<std::size_t>(idx)].state_size())
            throw std::invalid_argument("3D material-field replacement changes nonlinear state layout");
        const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(idx)];
        if(fi>=0) throw std::invalid_argument("FSC compound spring cannot be replaced through the scalar material field");
        const int ci=coupled_pm_by_spring_.empty()?-1:coupled_pm_by_spring_[static_cast<std::size_t>(idx)];
        if(ci>=0) throw std::invalid_argument("coupled P-M compound hinge requires dedicated interaction-field regeneration");
        const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(idx)];
        if(pi>=0) throw std::invalid_argument("return-mapped P-M compound hinge requires dedicated interaction-field regeneration");
        pending.push_back({idx,id,material,initial_tangents_[static_cast<std::size_t>(idx)],material.initial_stiffness()});
    }
    for(auto& r:pending){
        const double dk=r.new_k-r.old_k;
        if(dk!=0.0){
            for(const auto& [pos,c]:spring_initial_scatter_[static_cast<std::size_t>(r.index)])
                k_initial_.values()[static_cast<std::size_t>(pos)]+=dk*c;
            for(const auto& [pos,c]:spring_damping_scatter_[static_cast<std::size_t>(r.index)])
                k_damping_reference_.values()[static_cast<std::size_t>(pos)]+=dk*c;
        }
        materials_[static_cast<std::size_t>(r.index)]=std::move(r.material);
        initial_tangents_[static_cast<std::size_t>(r.index)]=r.new_k;
    }
}

int CompiledFrame3D::elastic_element_index(int element_id) const {
    const auto it=element_index_by_id_.find(element_id);
    return it==element_index_by_id_.end()?-1:it->second;
}

ElasticFrame3DResponse CompiledFrame3D::elastic_element_response(
    int element_id,const std::vector<double>& u) const {
    if(static_cast<int>(u.size())!=n_) throw std::invalid_argument("3D element response vector size");
    const int idx=elastic_element_index(element_id);
    if(idx<0) throw std::invalid_argument("unknown 3D elastic element id");
    const auto& e=element_response_data_[static_cast<std::size_t>(idx)];
    if(corotational_){
        std::array<double,12> eu{};
        for(int d=0;d<12;++d) eu[static_cast<std::size_t>(d)]=eval_terms(e.terms[static_cast<std::size_t>(d)],u);
        CorotationalFrame3DProperties cp;
        cp.xi=e.node_i_xyz[0];cp.yi=e.node_i_xyz[1];cp.zi=e.node_i_xyz[2];
        cp.xj=e.node_j_xyz[0];cp.yj=e.node_j_xyz[1];cp.zj=e.node_j_xyz[2];
        cp.E=e.properties.E;cp.G=e.properties.G;cp.A=e.properties.A;cp.J=e.properties.J;cp.Iy=e.properties.Iy;cp.Iz=e.properties.Iz;
        cp.reference=e.properties.reference;cp.axial_compression=e.properties.axial_compression;
        const auto r=corotational3d_section_response(cp,eu);
        return {element_id,r.axial_compression,r.shear_y_i,r.shear_y_j,r.shear_z_i,r.shear_z_j,
                r.torsion_i,r.torsion_j,r.moment_y_i,r.moment_y_j,r.moment_z_i,r.moment_z_j};
    }
    std::array<double,3> dx{e.node_j_xyz[0]-e.node_i_xyz[0],e.node_j_xyz[1]-e.node_i_xyz[1],e.node_j_xyz[2]-e.node_i_xyz[2]};
    const double L=norm3(dx); const auto ex=normalize3(dx);
    const double proj=dot3(e.properties.reference,ex);
    std::array<double,3> eyraw{e.properties.reference[0]-proj*ex[0],e.properties.reference[1]-proj*ex[1],e.properties.reference[2]-proj*ex[2]};
    auto ey=normalize3(eyraw); auto ez=normalize3(cross3(ex,ey)); ey=normalize3(cross3(ez,ex));
    const std::array<double,9> R{ex[0],ex[1],ex[2],ey[0],ey[1],ey[2],ez[0],ez[1],ez[2]};
    std::array<double,12> ug{},q{};
    for(int d=0;d<12;++d) ug[static_cast<std::size_t>(d)]=eval_terms(e.terms[static_cast<std::size_t>(d)],u);
    for(int block=0;block<4;++block){
        const int o=3*block;
        for(int r=0;r<3;++r)for(int c=0;c<3;++c)q[static_cast<std::size_t>(o+r)]+=R[static_cast<std::size_t>(r*3+c)]*ug[static_cast<std::size_t>(o+c)];
    }
    const double ka=e.properties.E*e.properties.A/L;
    const double p=e.properties.axial_compression-ka*(q[6]-q[0]);
    std::array<double,144> k{};
    auto add=[&](int r,int c,double v){k[static_cast<std::size_t>(r*12+c)]+=v;};
    const double kt=e.properties.G*e.properties.J/L;
    add(0,0,ka);add(0,6,-ka);add(6,0,-ka);add(6,6,ka);
    add(3,3,kt);add(3,9,-kt);add(9,3,-kt);add(9,9,kt);
    const double kz12=12*e.properties.E*e.properties.Iz/(L*L*L),kz6=6*e.properties.E*e.properties.Iz/(L*L),kz4=4*e.properties.E*e.properties.Iz/L,kz2=2*e.properties.E*e.properties.Iz/L;
    const int vz[4]{1,5,7,11}; const double bz[4][4]={{kz12,kz6,-kz12,kz6},{kz6,kz4,-kz6,kz2},{-kz12,-kz6,kz12,-kz6},{kz6,kz2,-kz6,kz4}};
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)add(vz[a],vz[b],bz[a][b]);
    const double ky12=12*e.properties.E*e.properties.Iy/(L*L*L),ky6=6*e.properties.E*e.properties.Iy/(L*L),ky4=4*e.properties.E*e.properties.Iy/L,ky2=2*e.properties.E*e.properties.Iy/L;
    const int wy[4]{2,4,8,10}; const double by[4][4]={{ky12,-ky6,-ky12,-ky6},{-ky6,ky4,ky6,ky2},{-ky12,ky6,ky12,ky6},{-ky6,ky2,ky6,ky4}};
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)add(wy[a],wy[b],by[a][b]);
    if(p>0.0){
        const double fac=p/(30.0*L);
        const double gz[4][4]={{36,3*L,-36,3*L},{3*L,4*L*L,-3*L,-L*L},{-36,-3*L,36,-3*L},{3*L,-L*L,-3*L,4*L*L}};
        const double gy[4][4]={{36,-3*L,-36,-3*L},{-3*L,4*L*L,3*L,-L*L},{-36,3*L,36,3*L},{-3*L,-L*L,3*L,4*L*L}};
        for(int a=0;a<4;++a)for(int b=0;b<4;++b){add(vz[a],vz[b],-fac*gz[a][b]);add(wy[a],wy[b],-fac*gy[a][b]);}
    }
    std::array<double,12> f{};
    for(int i=0;i<12;++i)for(int j=0;j<12;++j)f[static_cast<std::size_t>(i)]+=k[static_cast<std::size_t>(i*12+j)]*q[static_cast<std::size_t>(j)];
    return {element_id,p,f[1],f[7],f[2],f[8],f[3],f[9],f[4],f[10],f[5],f[11]};
}
CompiledFrame3D::FSCShearStateSnapshot CompiledFrame3D::fsc_shear_state(
    int spring_index,const std::vector<double>& state) const{
    FSCShearStateSnapshot out;
    if(spring_index<0||spring_index>=nonlinear_count()||static_cast<int>(state.size())!=nonlinear_state_size())return out;
    const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(spring_index)];
    if(fi<0)return out;
    const auto& fs=fsc_shear_springs_[static_cast<std::size_t>(fi)];
    enum : int { Q=0,F=1,INIT=2,QFAIL=3,VFAIL=4,CUMDEC=5,RETAIN=6,RESID=9 };
    const double* x=state.data()+fs.state_offset;
    out.valid=true;out.initiated=x[INIT]>0.5;out.residual_reached=x[RESID]>0.5;
    out.deformation=x[Q];out.force=x[F];out.failure_deformation=x[QFAIL];out.failure_strength=x[VFAIL];
    out.cumulative_strength_decrement=x[CUMDEC];out.retained_strength_ratio=x[RETAIN];
    return out;
}
std::vector<double> CompiledFrame3D::mass_multiply(const std::vector<double>& a) const{
    if(static_cast<int>(a.size())!=n_) throw std::invalid_argument("3D mass vector size");
    return mass_matrix_.multiply(a);
}
std::vector<double> CompiledFrame3D::damping_multiply(const std::vector<double>& v) const{
    if(static_cast<int>(v.size())!=n_) throw std::invalid_argument("3D damping vector size");
    auto mv=mass_matrix_.multiply(v), kv=k_damping_reference_.multiply(v); std::vector<double> out(static_cast<std::size_t>(n_));
    for(int i=0;i<n_;++i) out[static_cast<std::size_t>(i)]=alpha_m_*mv[static_cast<std::size_t>(i)]+beta_k_*kv[static_cast<std::size_t>(i)];
    return out;
}
SparseMatrixCSC CompiledFrame3D::effective_initial_matrix(double a0,double a1) const{
    auto out=system_pattern_; std::fill(out.values().begin(),out.values().end(),0.0);
    const double sm=a0+a1*alpha_m_;
    for(int p=0;p<k_initial_.nnz();++p) out.values()[static_cast<std::size_t>(k_to_system_[static_cast<std::size_t>(p)])] += k_initial_.values()[static_cast<std::size_t>(p)];
    if(a1*beta_k_!=0.0)for(int p=0;p<k_damping_reference_.nnz();++p) out.values()[static_cast<std::size_t>(kd_to_system_[static_cast<std::size_t>(p)])] += a1*beta_k_*k_damping_reference_.values()[static_cast<std::size_t>(p)];
    for(int p=0;p<mass_matrix_.nnz();++p) out.values()[static_cast<std::size_t>(m_to_system_[static_cast<std::size_t>(p)])] += sm*mass_matrix_.values()[static_cast<std::size_t>(p)];
    return out;
}
SparseMatrixCSC CompiledFrame3D::effective_tangent_matrix(const std::vector<double>& tangents,double a0,double a1) const{
    if(static_cast<int>(tangents.size())!=nonlinear_count()) throw std::invalid_argument("3D tangent vector size");
    auto out=effective_initial_matrix(a0,a1);
    for(int j=0;j<nonlinear_count();++j){const double dk=tangents[static_cast<std::size_t>(j)]-initial_tangents_[static_cast<std::size_t>(j)];if(dk==0)continue;for(const auto& [p,c]:spring_scatter_[static_cast<std::size_t>(j)])out.values()[static_cast<std::size_t>(p)]+=dk*c;}
    return out;
}
std::vector<double> CompiledFrame3D::generalized_state_update_coefficients(const std::vector<double>& u,const std::vector<double>& tangents) const{
    if(!has_generalized_state_update())return {};
    if(static_cast<int>(u.size())!=n_||static_cast<int>(tangents.size())!=nonlinear_count())throw std::invalid_argument("generalized state update dimension mismatch");
    const int r=generalized_state_update_basis_.cols();std::vector<double> C(static_cast<std::size_t>(r*r),0.0);
    // Concentrated nonlinear components remain the diagonal special case.
    for(int j=0;j<nonlinear_count();++j)C[static_cast<std::size_t>(j*r+j)]=tangents[static_cast<std::size_t>(j)]-initial_tangents_[static_cast<std::size_t>(j)];
    // For each geometric block: dp G - k[(G u)a^T + a(G u)^T].
    // G = Q Lambda Q^T and columns of Q plus a are fixed in U.
    for(const auto& b:geometric_low_rank_blocks_){
        const double delta=generalized_state_update_basis_.column_dot(b.axial_col,u);const double dp=-b.axial_k*delta;
        for(int k=0;k<b.eigen_count;++k){const int qcol=b.start_col+k;const double lam=b.eigenvalues[static_cast<std::size_t>(k)];const double eta=generalized_state_update_basis_.column_dot(qcol,u);
            C[static_cast<std::size_t>(qcol*r+qcol)]+=dp*lam;const double cross=-b.axial_k*lam*eta;C[static_cast<std::size_t>(qcol*r+b.axial_col)]+=cross;C[static_cast<std::size_t>(b.axial_col*r+qcol)]+=cross;}
    }
    return C;
}

SparseMatrixCSC CompiledFrame3D::effective_state_tangent_matrix(const std::vector<double>& u,const std::vector<double>& tangents,double a0,double a1) const{
    if(static_cast<int>(u.size())!=n_||static_cast<int>(tangents.size())!=nonlinear_count())throw std::invalid_argument("3D state tangent dimension mismatch");
    if(corotational_){
        auto out=system_pattern_;std::fill(out.values().begin(),out.values().end(),0.0);
        // Current structural tangent plus initial-stiffness Rayleigh damping.
        const double sm=a0+a1*alpha_m_;
        if(a1*beta_k_!=0.0)for(int p=0;p<k_damping_reference_.nnz();++p)out.values()[static_cast<std::size_t>(kd_to_system_[static_cast<std::size_t>(p)])]+=a1*beta_k_*k_damping_reference_.values()[static_cast<std::size_t>(p)];
        for(int p=0;p<mass_matrix_.nnz();++p)out.values()[static_cast<std::size_t>(m_to_system_[static_cast<std::size_t>(p)])]+=sm*mass_matrix_.values()[static_cast<std::size_t>(p)];
        for(const auto& ce:corotational_elements_){
            std::array<double,12> eu{};for(int d=0;d<12;++d)eu[static_cast<std::size_t>(d)]=eval_terms(ce.terms[static_cast<std::size_t>(d)],u);
            const auto er=corotational3d_response(ce.properties,eu);
            for(const auto& e:ce.scatter)out.values()[static_cast<std::size_t>(e.system_pos)]+=e.coefficient*er.tangent[static_cast<std::size_t>(e.local_row*12+e.local_col)];
        }
        for(int j=0;j<nonlinear_count();++j)for(const auto& [p,c]:spring_scatter_[static_cast<std::size_t>(j)])out.values()[static_cast<std::size_t>(p)]+=tangents[static_cast<std::size_t>(j)]*c;
        return out;
    }
    auto out=effective_tangent_matrix(tangents,a0,a1);if(!updated_pdelta_)return out;
    for(const auto& gu:geometric_updates_){
        const double delta=eval_terms(gu.axial_terms,u);const double dp=-gu.axial_k*delta;
        std::vector<std::pair<int,double>> guv;guv.reserve(gu.entries.size());
        std::map<int,double> accum;
        for(const auto& e:gu.entries){
            if(std::abs(dp)>=1e-18)out.values()[static_cast<std::size_t>(e.system_pos)]+=dp*e.value;
            accum[e.row]+=e.value*u[static_cast<std::size_t>(e.col)];
        }
        // f_geo = dp(u) G u, dp = -(EA/L) a^T u.  Its exact Jacobian is
        // dp G -(EA/L)(G u) a^T.  The second term is localized but generally
        // nonsymmetric; omitting it gives an inconsistent Newton tangent.
        for(const auto& e:gu.coupling){
            const auto it=accum.find(e.row);if(it==accum.end()||std::abs(it->second)<1e-18)continue;
            const double v=-gu.axial_k*it->second*e.axial_coefficient;
            out.values()[static_cast<std::size_t>(e.system_pos)]+=v;
            if(e.transpose_pos!=e.system_pos)out.values()[static_cast<std::size_t>(e.transpose_pos)]+=v;
        }
    }
    return out;
}
SparseMatrixCSC CompiledFrame3D::effective_state_tangent_matrix_with_state(
    const std::vector<double>& u,const std::vector<double>& tangents,
    const std::vector<double>& committed_state,double a0,double a1) const{
    if(static_cast<int>(committed_state.size())!=nonlinear_state_size())throw std::invalid_argument("3D committed-state tangent dimension mismatch");
    auto out=effective_state_tangent_matrix(u,tangents,a0,a1);
    if(coupled_pm_hinges_.empty()&&pm_return_hinges_.empty())return out;
    for(const auto& cp:coupled_pm_hinges_){
        const int j=cp.spring_index;const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        const double q=basis_.column_dot(j,u);
        std::vector<double> scratch(static_cast<std::size_t>(materials_[static_cast<std::size_t>(j)].state_size()+1),0.0);
        const auto ev=evaluate_coupled_pm(j,q,u,committed_state.data()+o,scratch.data());
        const double scale=-cp.axial_stiffness*ev.dmoment_dP;
        if(std::abs(scale)<1e-18)continue;
        for(const auto& e:cp.cross_scatter)out.values()[static_cast<std::size_t>(e.system_pos)]+=scale*e.coefficient;
    }
    for(const auto& ph:pm_return_hinges_){
        const int j=ph.spring_index;const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        const double q=basis_.column_dot(j,u);
        std::vector<double> scratch(static_cast<std::size_t>(materials_[static_cast<std::size_t>(j)].state_size()+PMInteractionHinge2D::kStateSize),0.0);
        const auto ev=evaluate_pm_return(j,q,u,committed_state.data()+o,scratch.data());
        if(!ev.converged)throw ConstitutiveIntegrationError("P-M local return mapping failed in tangent evaluation");
        // DEBUG9F
        // std::cerr << "pm block " << j << " n=" << ph.block_scatter.size() << " K=" << ev.tangent[0] << "," << ev.tangent[1] << "," << ev.tangent[2] << "," << ev.tangent[3] << "\n";
        for(const auto& e:ph.block_scatter){
            double kval=ev.tangent[2*e.local_row+e.local_col];
            if(e.local_row==0&&e.local_col==0&&!corotational_)kval-=ph.law->axial_initial_stiffness();
            out.values()[static_cast<std::size_t>(e.system_pos)]+=kval*e.coefficient;
        }
    }
    return out;
}

std::vector<double> CompiledFrame3D::initial_nonlinear_state() const {
    std::vector<double> state(static_cast<std::size_t>(nonlinear_state_size_),0.0);
    for(int j=0;j<nonlinear_count();++j){
        const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        materials_[static_cast<std::size_t>(j)].initialize_state(state.data()+o);
        const int ci=coupled_pm_by_spring_.empty()?-1:coupled_pm_by_spring_[static_cast<std::size_t>(j)];
        if(ci>=0){
            const auto& cp=coupled_pm_hinges_[static_cast<std::size_t>(ci)];
            const auto cap=eval_capacity_curve(cp.axial_force_points,cp.moment_capacity_points,cp.capacity_slopes,cp.axial_preload);
            state[static_cast<std::size_t>(o+materials_[static_cast<std::size_t>(j)].state_size())]=std::max(1e-9,cap.value);
        }
        const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(j)];
        if(pi>=0){
            const auto& ph=pm_return_hinges_[static_cast<std::size_t>(pi)];
            ph.law->initialize_state(state.data()+ph.state_offset);
        }
        const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(j)];
        if(fi>=0){
            const auto& fs=fsc_shear_springs_[static_cast<std::size_t>(fi)];
            fs.law.initialize_state(state.data()+fs.state_offset);
        }
    }
    return state;
}

CompiledFrame3D::CoupledPMEval CompiledFrame3D::evaluate_coupled_pm(
    int spring_index,double rotation,const std::vector<double>& u,
    const double* committed,double* trial_state) const {
    if(spring_index<0||spring_index>=nonlinear_count())throw std::invalid_argument("coupled P-M spring index");
    const int ci=coupled_pm_by_spring_.empty()?-1:coupled_pm_by_spring_[static_cast<std::size_t>(spring_index)];
    if(ci<0||ci>=static_cast<int>(coupled_pm_hinges_.size()))throw std::invalid_argument("spring is not coupled P-M");
    const auto& cp=coupled_pm_hinges_[static_cast<std::size_t>(ci)];
    const double axial_delta=eval_terms(cp.axial_terms,u);
    const double P=cp.axial_preload-cp.axial_stiffness*axial_delta;
    const auto cap=eval_capacity_curve(cp.axial_force_points,cp.moment_capacity_points,cp.capacity_slopes,P);
    const double my_trial=std::max(1e-9,cap.value);

    // ASCE41 dimensional reload targets must remain on the same moving
    // P-dependent strength surface as the envelope. Earlier Phase-8A trials
    // either changed Fy while retaining a stale TARGETF (discontinuous at
    // branch re-entry) or froze Fy for the entire branch (stable but lagged).
    // Re-evaluating only TARGETF on the current envelope preserves the
    // committed branch geometry (UZERO/TARGETQ), makes the force continuous in
    // P, and allows a consistent dM/dP tangent. At a fixed committed state the
    // current P reproduces the previously committed target exactly.
    enum : int { POSS=10,NEGS=11,FAIL=12,LATLOSS=13,BRANCH=14,TARGETQ=16,TARGETF=17 };
    const int base_state_size=materials_[static_cast<std::size_t>(spring_index)].state_size();
    auto eval_at=[&](double my_value,double* out_state){
        auto hp=cp.base_params;hp.posFy=hp.negFy=std::max(1e-9,my_value);
        std::array<double,ASCE41HingeMaterial::kStateSize> cwork{};
        std::copy(committed,committed+base_state_size,cwork.begin());
        if(cwork[BRANCH]>0.5 && cwork[LATLOSS]<=0.5 && cwork[FAIL]<=0.5){
            const double ps=cwork[POSS]>0.0?cwork[POSS]:1.0;
            const double ns=cwork[NEGS]>0.0?cwork[NEGS]:1.0;
            cwork[TARGETF]=asce41_envelope_force(hp,cwork[TARGETQ],ps,ns);
        }
        NonlinearMaterial dyn{ASCE41HingeMaterial(hp)};
        return dyn.trial(rotation,cwork.data(),out_state);
    };
    const auto nominal=eval_at(my_trial,trial_state);
    // Store the capacity associated with this accepted trial. It is diagnostic
    // state and also guarantees that a zero-increment trial begins from the
    // same moving strength surface used at the preceding accepted point.
    trial_state[base_state_size]=my_trial;

    // Partial derivative of moment with respect to the current capacity. The
    // branch-target update above is included in this finite difference, so the
    // global cross tangent represents the actual coupled constitutive map.
    const double h=std::max(1e-6,2e-6*std::abs(my_trial));
    std::array<double,ASCE41HingeMaterial::kStateSize> sp{},sm{};
    const auto plus=eval_at(my_trial+h,sp.data());
    const auto minus=eval_at(std::max(1e-9,my_trial-h),sm.data());
    const double den=(my_trial-h>1e-9)?2.0*h:h;
    const double dM_dMy=(my_trial-h>1e-9)?(plus.force-minus.force)/den:(plus.force-nominal.force)/den;
    CoupledPMEval out;out.material=nominal;out.axial_force=P;out.capacity=my_trial;
    out.dcapacity_dP=cap.derivative;out.dmoment_dP=dM_dMy*cap.derivative;
    return out;
}
PMInteractionTrialResult CompiledFrame3D::evaluate_pm_return(
    int spring_index,double rotation,const std::vector<double>& u,
    const double* committed,double* trial_state) const {
    if(spring_index<0||spring_index>=nonlinear_count())throw std::invalid_argument("P-M return spring index");
    const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(spring_index)];
    if(pi<0||pi>=static_cast<int>(pm_return_hinges_.size()))throw std::invalid_argument("spring is not P-M return hinge");
    const auto& ph=pm_return_hinges_[static_cast<std::size_t>(pi)];
    const int bs=materials_[static_cast<std::size_t>(spring_index)].state_size();
    std::copy(committed,committed+bs,trial_state);
    const double axial_delta=eval_terms(ph.axial_terms,u);
    return ph.law->trial(axial_delta,rotation,committed+bs,trial_state+bs);
}

FSCShearSpringTrial CompiledFrame3D::evaluate_fsc_shear(
    int spring_index,double shear_deformation,const std::vector<double>& u,
    const double* committed,double* trial_state) const {
    if(spring_index<0||spring_index>=nonlinear_count())throw std::invalid_argument("FSC shear spring index");
    const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(spring_index)];
    if(fi<0||fi>=static_cast<int>(fsc_shear_springs_.size()))throw std::invalid_argument("spring is not FSC shear");
    const auto& fs=fsc_shear_springs_[static_cast<std::size_t>(fi)];
    const int base_state_size=materials_[static_cast<std::size_t>(spring_index)].state_size();
    std::copy(committed,committed+base_state_size,trial_state);
    const double qb=eval_terms(fs.bottom_rotation_terms,u);
    const double qt=eval_terms(fs.top_rotation_terms,u);
    const double local_rotation=std::abs(qb)>=std::abs(qt)?qb:qt;
    const double axial_delta=eval_terms(fs.axial_terms,u);
    const double P=fs.axial_preload-fs.axial_stiffness*axial_delta;
    return fs.law.trial(shear_deformation,local_rotation,P,
                        committed+base_state_size,trial_state+base_state_size);
}

void CompiledFrame3D::evaluate_nonlinear_deformations(const std::vector<double>& q,
                                                        const std::vector<double>& committed_state,
                                                        std::vector<double>& component_forces,
                                                        std::vector<double>& tangents,
                                                        std::vector<double>& trial_state) const {
    if (static_cast<int>(q.size()) != nonlinear_count() ||
        static_cast<int>(committed_state.size()) != nonlinear_state_size())
        throw std::invalid_argument("3D nonlinear state dimension mismatch");
    component_forces.resize(static_cast<std::size_t>(nonlinear_count()));
    tangents.resize(static_cast<std::size_t>(nonlinear_count()));
    trial_state.resize(static_cast<std::size_t>(nonlinear_state_size()));
    for (int j=0;j<nonlinear_count();++j) {
        const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        const auto tr=materials_[static_cast<std::size_t>(j)].trial(q[static_cast<std::size_t>(j)],
            committed_state.data()+o, trial_state.data()+o);
        const int ci=coupled_pm_by_spring_.empty()?-1:coupled_pm_by_spring_[static_cast<std::size_t>(j)];
        if(ci>=0){
            const int bs=materials_[static_cast<std::size_t>(j)].state_size();
            trial_state[static_cast<std::size_t>(o+bs)]=committed_state[static_cast<std::size_t>(o+bs)];
        }
        const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(j)];
        if(pi>=0){
            const auto& ph=pm_return_hinges_[static_cast<std::size_t>(pi)];
            std::copy(committed_state.begin()+ph.state_offset,
                      committed_state.begin()+ph.state_offset+PMInteractionHinge2D::kStateSize,
                      trial_state.begin()+ph.state_offset);
        }
        const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(j)];
        if(fi>=0){
            // This reduced-coordinate-only API has no access to the remote end
            // rotations/axial deformation that trigger an FSC spring. Preserve
            // its committed FSC state and report the elastic placeholder. Full
            // state-aware response is available through internal_force... and
            // fsc_shear_state().
            const auto& fs=fsc_shear_springs_[static_cast<std::size_t>(fi)];
            std::copy(committed_state.begin()+fs.state_offset,
                      committed_state.begin()+fs.state_offset+FSCShearSpringLaw::kStateSize,
                      trial_state.begin()+fs.state_offset);
        }
        component_forces[static_cast<std::size_t>(j)]=tr.force;
        tangents[static_cast<std::size_t>(j)]=tr.tangent;
    }
}

void CompiledFrame3D::internal_force_and_tangent(const std::vector<double>& u,const std::vector<double>& committed_state,
                                                  std::vector<double>& force,std::vector<double>& tangents,
                                                  std::vector<double>& trial_state) const {
    internal_force_and_tangent_diagnostics(u,committed_state,force,tangents,trial_state,nullptr);
}

void CompiledFrame3D::internal_force_and_tangent_diagnostics(const std::vector<double>& u,
                                                  const std::vector<double>& committed_state,
                                                  std::vector<double>& force,std::vector<double>& tangents,
                                                  std::vector<double>& trial_state,
                                                  NonlinearEvalDiagnostics* diagnostics) const {
    if(static_cast<int>(u.size())!=n_||static_cast<int>(committed_state.size())!=nonlinear_state_size())
        throw std::invalid_argument("3D state dimension mismatch");
    if(corotational_){
        force.assign(static_cast<std::size_t>(n_),0.0);
        for(const auto& ce:corotational_elements_){
            std::array<double,12> eu{};for(int d=0;d<12;++d)eu[static_cast<std::size_t>(d)]=eval_terms(ce.terms[static_cast<std::size_t>(d)],u);
            const auto ef=corotational3d_internal_force(ce.properties,eu);
            for(int d=0;d<12;++d)for(const auto& [r,c]:ce.terms[static_cast<std::size_t>(d)])force[static_cast<std::size_t>(r)]+=c*ef[static_cast<std::size_t>(d)];
        }
    }else force=k_linear_.multiply(u);
    if(updated_pdelta_){for(const auto& gu:geometric_updates_){
        const double delta=eval_terms(gu.axial_terms,u);const double dp=-gu.axial_k*delta;
        std::map<int,double> guv;double quad=0.0;
        for(const auto& e:gu.entries)guv[e.row]+=e.value*u[static_cast<std::size_t>(e.col)];
        for(const auto& [r,v]:guv){if(std::abs(dp)>=1e-18)force[static_cast<std::size_t>(r)]+=dp*v;quad+=u[static_cast<std::size_t>(r)]*v;}
        // Conservative small-rotation update from the potential
        // 1/2 * dp(u) * u^T G u, with dp=-(EA/L)a^T u.  The axial
        // correction below makes the internal force and Newton tangent
        // energetically consistent and preserves tangent symmetry.
        if(std::abs(quad)>=1e-18)for(const auto& [r,ac]:gu.axial_terms)force[static_cast<std::size_t>(r)]+=-0.5*gu.axial_k*quad*ac;
    }}
    tangents.resize(static_cast<std::size_t>(nonlinear_count()));
    trial_state.resize(static_cast<std::size_t>(nonlinear_state_size()));
    for(int j=0;j<nonlinear_count();++j){
        const double q=basis_.column_dot(j,u);
        const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        MaterialTrialResult tr;
        const int ci=coupled_pm_by_spring_.empty()?-1:coupled_pm_by_spring_[static_cast<std::size_t>(j)];
        const int pi=pm_return_by_spring_.empty()?-1:pm_return_by_spring_[static_cast<std::size_t>(j)];
        const int fi=fsc_shear_by_spring_.empty()?-1:fsc_shear_by_spring_[static_cast<std::size_t>(j)];
        if(ci>=0)tr=evaluate_coupled_pm(j,q,u,committed_state.data()+o,trial_state.data()+o).material;
        else if(pi>=0){
            const auto pev=evaluate_pm_return(j,q,u,committed_state.data()+o,trial_state.data()+o);
            if(!pev.converged){
                const auto& ph=pm_return_hinges_[static_cast<std::size_t>(pi)];
                const double da=eval_terms(ph.axial_terms,u);
                const int bs=materials_[static_cast<std::size_t>(j)].state_size();
                const double* ps=committed_state.data()+o+bs;
                std::ostringstream ss;
                ss<<std::setprecision(17)<<"P-M local return mapping failed component="<<component_ids_[static_cast<std::size_t>(j)]
                  <<" index="<<j<<" da="<<da<<" theta="<<q
                  <<" epA="<<ps[0]<<" epR="<<ps[1]<<" kp="<<ps[2]<<" kn="<<ps[3]
                  <<" lastN="<<ps[4]<<" lastM="<<ps[5];
                const auto& prm=ph.law->params();
                ss<<" P0="<<prm.axial_preload<<" Ka="<<prm.axial_stiffness<<" Kr="<<prm.hinge.Ke<<" PB="<<prm.p_balance<<" PT="<<prm.py_tension<<" PC="<<prm.py_compression<<" My="<<prm.my_balance<<" alphaT="<<prm.alpha_tension<<" alphaC="<<prm.alpha_compression;
                throw ConstitutiveIntegrationError(ss.str());
            }
            tr.force=pev.moment;tr.tangent=pev.tangent[3];
            tr.diagnostics.fast_path=(pev.events&PM_EVENT_ELASTIC)!=0;
            tr.diagnostics.tangent_active=std::abs(pev.tangent[3]-materials_[static_cast<std::size_t>(j)].initial_stiffness())>1e-12;
            tr.diagnostics.transition=(pev.events&(PM_EVENT_YIELD|PM_EVENT_CAP|PM_EVENT_STRENGTH_DROP))!=0;
            tr.diagnostics.failed=(pev.events&PM_EVENT_FAILURE)!=0;
            tr.diagnostics.lateral_resistance_lost=(pev.events&PM_EVENT_LATERAL_LOSS)!=0;
            const auto& ph=pm_return_hinges_[static_cast<std::size_t>(pi)];
            for(const auto& [r,a]:ph.axial_terms)force[static_cast<std::size_t>(r)]+=a*pev.axial_force;
        }
        else if(fi>=0)tr=evaluate_fsc_shear(j,q,u,committed_state.data()+o,trial_state.data()+o).material;
        else tr=materials_[static_cast<std::size_t>(j)].trial(q,committed_state.data()+o,trial_state.data()+o);
        tangents[static_cast<std::size_t>(j)]=tr.tangent;
        basis_.axpy_column(j,tr.force,force);
        if(diagnostics){
            ++diagnostics->component_evaluations;
            diagnostics->active_tangent_evaluations += static_cast<std::size_t>(tr.diagnostics.tangent_active);
            diagnostics->fast_path_evaluations += static_cast<std::size_t>(tr.diagnostics.fast_path);
            diagnostics->full_state_evaluations += static_cast<std::size_t>(!tr.diagnostics.fast_path);
            diagnostics->transition_events += static_cast<std::size_t>(tr.diagnostics.transition);
            diagnostics->reversal_events += static_cast<std::size_t>(tr.diagnostics.reversal);
            diagnostics->deterioration_events += static_cast<std::size_t>(tr.diagnostics.deterioration);
            diagnostics->failure_events += static_cast<std::size_t>(tr.diagnostics.failed);
            diagnostics->io_or_beyond_evaluations += static_cast<std::size_t>(tr.diagnostics.at_or_beyond_io);
            diagnostics->ls_or_beyond_evaluations += static_cast<std::size_t>(tr.diagnostics.at_or_beyond_ls);
            diagnostics->cp_or_beyond_evaluations += static_cast<std::size_t>(tr.diagnostics.at_or_beyond_cp);
            diagnostics->beyond_cp_evaluations += static_cast<std::size_t>(tr.diagnostics.beyond_cp);
            diagnostics->lateral_loss_evaluations += static_cast<std::size_t>(tr.diagnostics.lateral_resistance_lost);
        }
    }
}

std::vector<double> CompiledFrame3D::base_excitation(double ag) const{std::vector<double> p(static_cast<std::size_t>(n_));for(int i=0;i<n_;++i)p[static_cast<std::size_t>(i)]=-base_mass_[static_cast<std::size_t>(i)]*ag;return p;}
double CompiledFrame3D::response_value(const std::vector<double>& u) const{if(static_cast<int>(u.size())!=n_)throw std::invalid_argument("3D response size");return eval_terms(response_terms_,u);}
std::vector<double> CompiledFrame3D::story_response_values(const std::vector<double>& x) const{
    if(static_cast<int>(x.size())!=n_)throw std::invalid_argument("3D story response size");
    std::vector<double> out;out.reserve(story_terms_.size());
    for(const auto& terms:story_terms_)out.push_back(eval_terms(terms,x));
    return out;
}
double CompiledFrame3D::max_drift_measure(const std::vector<double>& u) const{if(static_cast<int>(u.size())!=n_)throw std::invalid_argument("3D response size");double lower=0,out=0;for(const auto& terms:story_terms_){const double x=eval_terms(terms,u);out=std::max(out,std::abs(x-lower));lower=x;}return out;}
double CompiledFrame3D::max_drift_ratio(const std::vector<double>& u) const{if(static_cast<int>(u.size())!=n_)throw std::invalid_argument("3D response size");if(story_terms_.empty()||story_z_.size()!=story_terms_.size())return std::numeric_limits<double>::quiet_NaN();double lower_u=0.0,lower_z=0.0,out=0.0;for(std::size_t i=0;i<story_terms_.size();++i){const double x=eval_terms(story_terms_[i],u);const double dz=story_z_[i]-lower_z;if(dz<=0.0)return std::numeric_limits<double>::quiet_NaN();out=std::max(out,std::abs(x-lower_u)/dz);lower_u=x;lower_z=story_z_[i];}return out;}

} // namespace quake
