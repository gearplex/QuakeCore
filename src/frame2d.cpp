#include "quake/frame2d.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <set>
#include <unordered_map>

namespace quake {
namespace {

class UnionFind {
public:
    explicit UnionFind(int n) : p_(static_cast<std::size_t>(n)), rank_(static_cast<std::size_t>(n),0) {
        std::iota(p_.begin(), p_.end(), 0);
    }
    int find(int x) {
        if (p_[static_cast<std::size_t>(x)] != x)
            p_[static_cast<std::size_t>(x)] = find(p_[static_cast<std::size_t>(x)]);
        return p_[static_cast<std::size_t>(x)];
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        auto& ra = rank_[static_cast<std::size_t>(a)];
        auto& rb = rank_[static_cast<std::size_t>(b)];
        if (ra < rb) std::swap(a,b);
        p_[static_cast<std::size_t>(b)] = a;
        if (ra == rb) ++ra;
    }
private:
    std::vector<int> p_;
    std::vector<int> rank_;
};

int dof_index(int node_slot, Dof2D dof) {
    return 3*node_slot + static_cast<int>(dof);
}

std::array<double,36> transpose_multiply(const std::array<double,36>& T,
                                         const std::array<double,36>& K) {
    // Return T^T K T, row-major storage.
    std::array<double,36> KT{};
    std::array<double,36> out{};
    for (int i=0;i<6;++i)
        for (int j=0;j<6;++j)
            for (int k=0;k<6;++k)
                KT[static_cast<std::size_t>(i*6+j)] += K[static_cast<std::size_t>(i*6+k)] * T[static_cast<std::size_t>(k*6+j)];
    for (int i=0;i<6;++i)
        for (int j=0;j<6;++j)
            for (int k=0;k<6;++k)
                out[static_cast<std::size_t>(i*6+j)] += T[static_cast<std::size_t>(k*6+i)] * KT[static_cast<std::size_t>(k*6+j)];
    return out;
}

void add_scaled_bbt(std::vector<Triplet>& t, const SparseUpdateBasis& basis,
                    int col, double scale) {
    if (scale == 0.0) return;
    const int b=basis.col_ptr()[static_cast<std::size_t>(col)];
    const int e=basis.col_ptr()[static_cast<std::size_t>(col+1)];
    for(int pi=b;pi<e;++pi){
        const int i=basis.row_ind()[static_cast<std::size_t>(pi)];
        const double bi=basis.values()[static_cast<std::size_t>(pi)];
        for(int pj=b;pj<e;++pj){
            const int j=basis.row_ind()[static_cast<std::size_t>(pj)];
            const double bj=basis.values()[static_cast<std::size_t>(pj)];
            t.push_back({i,j,scale*bi*bj});
        }
    }
}

} // namespace

std::array<double,36> frame2d_global_stiffness(double xi, double yi,
                                               double xj, double yj,
                                               double E, double A, double I,
                                               double axial_compression) {
    if (E <= 0.0 || A <= 0.0 || I <= 0.0 || axial_compression < 0.0)
        throw std::invalid_argument("invalid frame property");
    const double dx=xj-xi, dy=yj-yi;
    const double L=std::hypot(dx,dy);
    if (L <= 0.0) throw std::invalid_argument("zero-length frame element");
    const double c=dx/L, s=dy/L;
    const double ea=E*A/L;
    const double e12=12.0*E*I/(L*L*L);
    const double e6=6.0*E*I/(L*L);
    const double e4=4.0*E*I/L;
    const double e2=2.0*E*I/L;

    std::array<double,36> k{};
    auto set=[&](int r,int c_,double v){k[static_cast<std::size_t>(r*6+c_)]=v;};
    set(0,0,ea); set(0,3,-ea); set(3,0,-ea); set(3,3,ea);
    set(1,1,e12); set(1,2,e6); set(1,4,-e12); set(1,5,e6);
    set(2,1,e6); set(2,2,e4); set(2,4,-e6); set(2,5,e2);
    set(4,1,-e12); set(4,2,-e6); set(4,4,e12); set(4,5,-e6);
    set(5,1,e6); set(5,2,e2); set(5,4,-e6); set(5,5,e4);

    if (axial_compression > 0.0) {
        // Standard beam-column geometric stiffness. Compression reduces tangent.
        const double q=axial_compression/(30.0*L);
        const int map[4]{1,2,4,5};
        const double g[4][4] = {
            {36.0, 3.0*L, -36.0, 3.0*L},
            {3.0*L, 4.0*L*L, -3.0*L, -L*L},
            {-36.0, -3.0*L, 36.0, -3.0*L},
            {3.0*L, -L*L, -3.0*L, 4.0*L*L}
        };
        for(int a=0;a<4;++a) for(int b=0;b<4;++b)
            k[static_cast<std::size_t>(map[a]*6+map[b])] -= q*g[a][b];
    }

    std::array<double,36> T{};
    T[0]=c; T[1]=s;
    T[6]=-s; T[7]=c;
    T[14]=1.0;
    T[21]=c; T[22]=s;
    T[27]=-s; T[28]=c;
    T[35]=1.0;
    return transpose_multiply(T,k);
}

void Frame2DBuilder::add_node(int id, double x, double y,
                              double mass_x, double mass_y, double mass_r) {
    if (!std::isfinite(x) || !std::isfinite(y)) throw std::invalid_argument("nonfinite 2D node coordinate");
    if (!std::isfinite(mass_x) || !std::isfinite(mass_y) || !std::isfinite(mass_r) ||
        mass_x < 0.0 || mass_y < 0.0 || mass_r < 0.0)
        throw std::invalid_argument("2D nodal mass must be finite and nonnegative");
    if (std::any_of(nodes_.begin(),nodes_.end(),[&](const Node2D& n){return n.id==id;}))
        throw std::invalid_argument("duplicate node id");
    nodes_.push_back({id,x,y,mass_x,mass_y,mass_r});
}

void Frame2DBuilder::add_elastic_frame(int id, int node_i, int node_j,
                                       double E, double A, double I,
                                       double axial_compression) {
    if(std::any_of(elements_.begin(),elements_.end(),[&](const ElasticFrame2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D elastic element id");
    elements_.push_back({id,node_i,node_j,E,A,I,axial_compression});
}

void Frame2DBuilder::add_rotational_spring(int id, int node_i, int node_j,
                                           double k0, double yield_moment,
                                           double post_yield_ratio) {
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D rotational spring id");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::RotationalHinge,NonlinearMaterial(BilinearSpring(k0,yield_moment,post_yield_ratio))});
}
void Frame2DBuilder::add_asce41_hinge(int id,int node_i,int node_j,ASCE41HingeParams params){
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D rotational spring id");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::RotationalHinge,NonlinearMaterial(ASCE41HingeMaterial(params))});
}
void Frame2DBuilder::add_imk_peak_oriented_hinge(int id,int node_i,int node_j,IMKPeakOrientedParams params){
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D rotational spring id");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::RotationalHinge,NonlinearMaterial(IMKPeakOrientedMaterial(params))});
}

void Frame2DBuilder::add_panel_zone(int id,int node_i,int node_j,NonlinearMaterial material){
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D scalar component id");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::PanelZone,std::move(material)});
}

void Frame2DBuilder::add_brb(int id,int node_i,int node_j,double k0,double yield_force,double post_yield_ratio){
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D scalar component id");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::BRB,
                        NonlinearMaterial(BilinearSpring(k0,yield_force,post_yield_ratio))});
}

void Frame2DBuilder::add_soil_spring(int id,int node_i,int node_j,double dx,double dy,NonlinearMaterial material){
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D scalar component id");
    const double L=std::hypot(dx,dy);if(!std::isfinite(L)||L<=0.0)
        throw std::invalid_argument("soil spring direction must be finite and nonzero");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::SoilSpring,std::move(material),true,dx/L,dy/L,0.0});
}

void Frame2DBuilder::add_rotational_soil_spring(int id,int node_i,int node_j,NonlinearMaterial material){
    if(std::any_of(springs_.begin(),springs_.end(),[&](const ScalarComponent2D& e){return e.id==id;}))
        throw std::invalid_argument("duplicate 2D scalar component id");
    springs_.push_back({id,node_i,node_j,NonlinearComponent2DKind::SoilSpring,std::move(material),true,0.0,0.0,1.0});
}

void Frame2DBuilder::add_steel_member(int id,int i,int j,SteelMember2DRole role,SteelMember2DProperties p){
    steel_members_.push_back({id,i,j,role,std::move(p)});
}

void Frame2DBuilder::add_viscous_damper(int id,int i,int j,ViscousDamper2DProperties p){
    dampers_.push_back({id,i,j,p,false,0.0,0.0,0.0});
}

void Frame2DBuilder::add_directional_viscous_damper(int id,int i,int j,double dx,double dy,ViscousDamper2DProperties p){
    const double L=std::hypot(dx,dy);if(!std::isfinite(L)||L<=0.0)
        throw std::invalid_argument("soil dashpot direction must be finite and nonzero");
    dampers_.push_back({id,i,j,p,true,dx/L,dy/L,0.0});
}

void Frame2DBuilder::add_rotational_soil_dashpot(int id,int i,int j,ViscousDamper2DProperties p){
    dampers_.push_back({id,i,j,p,true,0.0,0.0,1.0});
}

void Frame2DBuilder::add_mvlem(int id,int i,int j,MVLEMProperties p){walls_.push_back({id,i,j,std::move(p)});}
void Frame2DBuilder::add_sfi_mvlem(int id,int i,int j,SFIMVLEMProperties p){walls_.push_back({id,i,j,std::move(p)});}

void Frame2DBuilder::fix(int node_id, bool ux, bool uy, bool rz) {
    fixed_node_ids_.push_back(node_id);
    fixed_.push_back({ux,uy,rz});
}

void Frame2DBuilder::equal_dof(int master_node, int slave_node, Dof2D dof) {
    equal_.push_back({master_node,slave_node,dof});
}

void Frame2DBuilder::rigid_floor_x(int master_node, const std::vector<int>& slave_nodes) {
    for (int s : slave_nodes) if (s != master_node) equal_dof(master_node,s,Dof2D::UX);
}

void Frame2DBuilder::set_rayleigh(double alpha_m, double beta_k) {
    if (alpha_m < 0.0 || beta_k < 0.0) throw std::invalid_argument("negative Rayleigh coefficient");
    alpha_m_=alpha_m; beta_k_=beta_k;
}

void Frame2DBuilder::set_response_node(int node_id) { response_node_id_=node_id; }
void Frame2DBuilder::set_story_nodes(std::vector<int> ids) { story_node_ids_=std::move(ids); }

CompiledFrame2D Frame2DBuilder::compile() const {
    if (nodes_.empty()) throw std::invalid_argument("frame has no nodes");
    CompiledFrame2D out;
    out.alpha_m_=alpha_m_; out.beta_k_=beta_k_;
    out.elastic_element_count_=static_cast<int>(elements_.size());
    out.node_ids_.reserve(nodes_.size());
    for (std::size_t i=0;i<nodes_.size();++i) {
        out.node_ids_.push_back(nodes_[i].id);
        out.node_slot_.emplace(nodes_[i].id,static_cast<int>(i));
    }
    auto slot=[&](int id)->int{
        auto it=out.node_slot_.find(id);
        if(it==out.node_slot_.end()) throw std::invalid_argument("unknown node id "+std::to_string(id));
        return it->second;
    };

    const int full_n=static_cast<int>(3*nodes_.size());
    UnionFind uf(full_n);
    for(const auto& e:equal_) uf.unite(dof_index(slot(e.master),e.dof),dof_index(slot(e.slave),e.dof));

    std::vector<bool> fixed_root(static_cast<std::size_t>(full_n),false);
    for(std::size_t i=0;i<fixed_node_ids_.size();++i){
        const int ns=slot(fixed_node_ids_[i]);
        for(int d=0;d<3;++d) if(fixed_[i][static_cast<std::size_t>(d)])
            fixed_root[static_cast<std::size_t>(uf.find(3*ns+d))]=true;
    }
    // Propagate fixed status after all unions.
    for(int i=0;i<full_n;++i) if(fixed_root[static_cast<std::size_t>(i)]) fixed_root[static_cast<std::size_t>(uf.find(i))]=true;

    std::unordered_map<int,int> root_to_reduced;
    out.full_to_reduced_.assign(static_cast<std::size_t>(full_n),-1);
    for(int i=0;i<full_n;++i){
        const int r=uf.find(i);
        if(fixed_root[static_cast<std::size_t>(r)]) continue;
        auto [it,inserted]=root_to_reduced.emplace(r,static_cast<int>(root_to_reduced.size()));
        (void)inserted;
        out.full_to_reduced_[static_cast<std::size_t>(i)]=it->second;
    }
    out.n_=static_cast<int>(root_to_reduced.size());
    if(out.n_==0) throw std::invalid_argument("frame has no active DOFs");
    out.mass_.assign(static_cast<std::size_t>(out.n_),0.0);
    out.base_mass_x_.assign(static_cast<std::size_t>(out.n_),0.0);
    for(std::size_t ni=0;ni<nodes_.size();++ni){
        const double masses[3]{nodes_[ni].mass_x,nodes_[ni].mass_y,nodes_[ni].mass_r};
        for(int d=0;d<3;++d){
            const int rd=out.full_to_reduced_[3*ni+static_cast<std::size_t>(d)];
            if(rd<0) continue;
            out.mass_[static_cast<std::size_t>(rd)] += masses[d];
            if(d==0) out.base_mass_x_[static_cast<std::size_t>(rd)] += nodes_[ni].mass_x;
        }
    }

    std::vector<Triplet> linear;
    out.element_response_data_.reserve(elements_.size());
    for(const auto& e:elements_){
        const int si=slot(e.node_i), sj=slot(e.node_j);
        const auto& ni=nodes_[static_cast<std::size_t>(si)];
        const auto& nj=nodes_[static_cast<std::size_t>(sj)];
        const auto kg=frame2d_global_stiffness(ni.x,ni.y,nj.x,nj.y,e.E,e.A,e.I,e.axial_compression);
        int fd[6]{3*si,3*si+1,3*si+2,3*sj,3*sj+1,3*sj+2};
        CompiledFrame2D::ElementResponseData rd; rd.properties=e; rd.xi=ni.x; rd.yi=ni.y; rd.xj=nj.x; rd.yj=nj.y;
        for(int d=0;d<6;++d) rd.reduced_dofs[static_cast<std::size_t>(d)]=out.full_to_reduced_[static_cast<std::size_t>(fd[d])];
        if(out.element_index_by_id_.count(e.id)) throw std::invalid_argument("duplicate 2D elastic element id");
        out.element_index_by_id_[e.id]=static_cast<int>(out.element_response_data_.size());
        out.element_response_data_.push_back(std::move(rd));
        for(int a=0;a<6;++a){
            const int ra=out.full_to_reduced_[static_cast<std::size_t>(fd[a])];
            if(ra<0) continue;
            for(int b=0;b<6;++b){
                const int rb=out.full_to_reduced_[static_cast<std::size_t>(fd[b])];
                if(rb<0) continue;
                const double v=kg[static_cast<std::size_t>(a*6+b)];
                if(v!=0.0) linear.push_back({ra,rb,v});
            }
        }
    }
    out.k_linear_=SparseMatrixCSC::from_triplets(out.n_,out.n_,linear,1e-18);

    const int m=static_cast<int>(springs_.size());
    std::vector<std::vector<std::pair<int,double>>> basis_columns(static_cast<std::size_t>(m));
    out.materials_.reserve(springs_.size());
    out.material_state_offsets_.resize(static_cast<std::size_t>(m+1),0);
    for(int j=0;j<m;++j){
        const auto& sp=springs_[static_cast<std::size_t>(j)];
        const int si=slot(sp.node_i),sj=slot(sp.node_j);
        const auto& ni=nodes_[static_cast<std::size_t>(si)];const auto& nj=nodes_[static_cast<std::size_t>(sj)];
        if(sp.kind==NonlinearComponent2DKind::SoilSpring){
            if(std::hypot(nj.x-ni.x,nj.y-ni.y)>1e-10)
                throw std::invalid_argument("soil spring nodes must be coincident");
            const int full[6]{3*si,3*si+1,3*si+2,3*sj,3*sj+1,3*sj+2};
            const double coeff[6]{-sp.direction_x,-sp.direction_y,-sp.direction_rz,sp.direction_x,sp.direction_y,sp.direction_rz};
            for(int a=0;a<6;++a)if(coeff[a]!=0.0){const int r=out.full_to_reduced_[static_cast<std::size_t>(full[a])];if(r>=0)basis_columns[static_cast<std::size_t>(j)].push_back({r,coeff[a]});}
        }else if(sp.kind==NonlinearComponent2DKind::BRB){
            const double dx=nj.x-ni.x,dy=nj.y-ni.y,L=std::hypot(dx,dy);
            if(!std::isfinite(L)||L<=0.0)throw std::invalid_argument("BRB requires distinct finite nodes");
            const double c=dx/L,s=dy/L;const int full[4]{3*si,3*si+1,3*sj,3*sj+1};
            const double coeff[4]{-c,-s,c,s};
            for(int a=0;a<4;++a){const int r=out.full_to_reduced_[static_cast<std::size_t>(full[a])];if(r>=0)basis_columns[static_cast<std::size_t>(j)].push_back({r,coeff[a]});}
        }else{
            if(sp.kind==NonlinearComponent2DKind::PanelZone&&std::hypot(nj.x-ni.x,nj.y-ni.y)>1e-10)
                throw std::invalid_argument("panel zone nodes must be coincident");
            const int ri=out.full_to_reduced_[static_cast<std::size_t>(dof_index(si,Dof2D::RZ))];
            const int rj=out.full_to_reduced_[static_cast<std::size_t>(dof_index(sj,Dof2D::RZ))];
            if(ri>=0) basis_columns[static_cast<std::size_t>(j)].push_back({ri,-1.0});
            if(rj>=0) basis_columns[static_cast<std::size_t>(j)].push_back({rj,1.0});
        }
        if(out.spring_index_by_id_.count(sp.id)) throw std::invalid_argument("duplicate 2D rotational spring id");
        out.spring_index_by_id_[sp.id]=j;
        out.spring_ids_.push_back(sp.id);
        out.component_kinds_.push_back(sp.kind);
        out.materials_.push_back(sp.material);
        out.initial_tangents_.push_back(sp.material.initial_stiffness());
        out.material_state_offsets_[static_cast<std::size_t>(j+1)]=out.material_state_offsets_[static_cast<std::size_t>(j)]+sp.material.state_size();
    }
    out.nonlinear_state_size_=out.material_state_offsets_.back();
    out.basis_=SparseUpdateBasis::from_columns(out.n_,basis_columns,1e-14);
    for(int j=0;j<m;++j) if(out.basis_.col_ptr()[static_cast<std::size_t>(j)]==out.basis_.col_ptr()[static_cast<std::size_t>(j+1)])
        throw std::invalid_argument("nonlinear scalar component has constrained zero deformation");
    std::vector<Triplet> initial=linear;
    for(int j=0;j<m;++j) add_scaled_bbt(initial,out.basis_,j,out.initial_tangents_[static_cast<std::size_t>(j)]);
    // Walls contribute all 6x6 structural entries, including zeros: cracking and
    // coupled panel changes must never alter the prepared sparse pattern.
    std::set<int> block_ids;
    for(const auto& w:walls_){
        if(!block_ids.insert(w.id).second || out.element_index_by_id_.count(w.id) || out.spring_index_by_id_.count(w.id))
            throw std::invalid_argument("duplicate wall/element id");
        int si=slot(w.i),sj=slot(w.j);
        const auto& ni=nodes_[si];const auto& nj=nodes_[sj];double h=nj.y-ni.y;
        if(!std::isfinite(h)||h<=0||std::abs(nj.x-ni.x)>1e-10*h)throw std::invalid_argument("wall nodes must be vertical, ordered bottom to top");
        Wall2D e=std::visit([&](const auto& p){return Wall2D(h,p);},w.properties);
        CompiledFrame2D::CompiledWall cw{w.id,out.nonlinear_state_size_,std::move(e),{},{},{}};
        out.nonlinear_state_size_+=cw.element.state_size();cw.k0=cw.element.initial_tangent();cw.scatter.fill(-1);
        int fd[6]{3*si,3*si+1,3*si+2,3*sj,3*sj+1,3*sj+2};
        for(int a=0;a<6;++a){cw.dofs[a]=out.full_to_reduced_[fd[a]];
            if(a%3!=2 && cw.dofs[a]>=0){double mass=cw.element.total_mass()/2;out.mass_[cw.dofs[a]]+=mass;if(a%3==0)out.base_mass_x_[cw.dofs[a]]+=mass;}}
        for(int a=0;a<6;++a)for(int b=0;b<6;++b)if(cw.dofs[a]>=0&&cw.dofs[b]>=0)initial.push_back({cw.dofs[a],cw.dofs[b],cw.k0[6*a+b]});
        out.walls_.push_back(std::move(cw));
    }
    for(const auto& sm:steel_members_){
        if(!block_ids.insert(sm.id).second||out.element_index_by_id_.count(sm.id)||out.spring_index_by_id_.count(sm.id))
            throw std::invalid_argument("duplicate steel member/element id");
        const int si=slot(sm.i),sj=slot(sm.j);const auto& ni=nodes_[static_cast<std::size_t>(si)];const auto& nj=nodes_[static_cast<std::size_t>(sj)];
        SteelMember2D element(ni.x,ni.y,nj.x,nj.y,sm.properties);
        CompiledFrame2D::CompiledSteelMember cs{sm.id,out.nonlinear_state_size_,sm.role,std::move(element),{},{},{}};
        out.nonlinear_state_size_+=cs.element.state_size();cs.k0=cs.element.initial_tangent();cs.scatter.fill(-1);
        const int fd[6]{3*si,3*si+1,3*si+2,3*sj,3*sj+1,3*sj+2};
        for(int a=0;a<6;++a)cs.dofs[a]=out.full_to_reduced_[static_cast<std::size_t>(fd[a])];
        for(int a=0;a<6;++a)for(int b=0;b<6;++b)if(cs.dofs[a]>=0&&cs.dofs[b]>=0)
            initial.push_back({cs.dofs[a],cs.dofs[b],cs.k0[6*a+b]});
        out.steel_members_.push_back(std::move(cs));
    }
    for(const auto& d:dampers_){
        if(!block_ids.insert(d.id).second||out.element_index_by_id_.count(d.id)||out.spring_index_by_id_.count(d.id))
            throw std::invalid_argument("duplicate viscous damper/element id");
        const int si=slot(d.i),sj=slot(d.j);const auto& ni=nodes_[static_cast<std::size_t>(si)];const auto& nj=nodes_[static_cast<std::size_t>(sj)];
        const double dx=nj.x-ni.x,dy=nj.y-ni.y,L=std::hypot(dx,dy);
        if(!d.explicit_direction&&(!std::isfinite(L)||L<=0.0))throw std::invalid_argument("viscous damper requires distinct finite nodes");
        if(d.explicit_direction&&std::hypot(dx,dy)>1e-10)throw std::invalid_argument("directional soil dashpot nodes must be coincident");
        CompiledFrame2D::CompiledDamper cd{d.id,ViscousDamper2D(d.properties),{},{},{},0.0};
        cd.c0=cd.element.initial_tangent();cd.scatter.fill(-1);const double c=d.explicit_direction?d.direction_x:dx/L,s=d.explicit_direction?d.direction_y:dy/L;
        const double rz=d.explicit_direction?d.direction_rz:0.0;cd.b={-c,-s,-rz,c,s,rz};const int fd[6]{3*si,3*si+1,3*si+2,3*sj,3*sj+1,3*sj+2};
        for(int a=0;a<6;++a)cd.dofs[a]=out.full_to_reduced_[static_cast<std::size_t>(fd[a])];
        // Reserve the full axial dashpot block in the immutable CSC pattern.
        for(int a=0;a<6;++a)for(int b=0;b<6;++b)if(cd.dofs[a]>=0&&cd.dofs[b]>=0&&cd.b[a]!=0.0&&cd.b[b]!=0.0)
            initial.push_back({cd.dofs[a],cd.dofs[b],0.0});
        out.dampers_.push_back(std::move(cd));
    }
    out.k_initial_=SparseMatrixCSC::from_triplets(out.n_,out.n_,initial,-1.0);
    // Preserve an explicit diagonal slot for every active DOF so the compiled
    // effective matrix can add mass terms without rebuilding its sparsity.
    {
        std::vector<Triplet> with_diag;
        with_diag.reserve(static_cast<std::size_t>(out.k_initial_.nnz()+out.n_));
        for(int c=0;c<out.k_initial_.cols();++c)
            for(int p=out.k_initial_.col_ptr()[static_cast<std::size_t>(c)];
                p<out.k_initial_.col_ptr()[static_cast<std::size_t>(c+1)];++p)
                with_diag.push_back({out.k_initial_.row_ind()[static_cast<std::size_t>(p)],c,
                                     out.k_initial_.values()[static_cast<std::size_t>(p)]});
        for(int i=0;i<out.n_;++i) with_diag.push_back({i,i,0.0});
        out.k_initial_=SparseMatrixCSC::from_triplets(out.n_,out.n_,with_diag,-1.0);
    }

    auto find_csc_position=[&](int row,int col)->int {
        const int begin=out.k_initial_.col_ptr()[static_cast<std::size_t>(col)];
        const int end=out.k_initial_.col_ptr()[static_cast<std::size_t>(col+1)];
        auto first=out.k_initial_.row_ind().begin()+begin;
        auto last=out.k_initial_.row_ind().begin()+end;
        auto it=std::lower_bound(first,last,row);
        if(it==last || *it!=row) return -1;
        return static_cast<int>(std::distance(out.k_initial_.row_ind().begin(),it));
    };
    out.diagonal_positions_.resize(static_cast<std::size_t>(out.n_));
    for(int i=0;i<out.n_;++i){
        const int pos=find_csc_position(i,i);
        if(pos<0) throw std::runtime_error("compiled frame has active DOF without stiffness diagonal");
        out.diagonal_positions_[static_cast<std::size_t>(i)]=pos;
    }
    out.spring_scatter_.resize(static_cast<std::size_t>(m));
    for(int j=0;j<m;++j){
        std::vector<std::pair<int,double>> entries;
        const int pb=out.basis_.col_ptr()[static_cast<std::size_t>(j)],pe=out.basis_.col_ptr()[static_cast<std::size_t>(j+1)];
        for(int pc=pb;pc<pe;++pc){
            const int col=out.basis_.row_ind()[static_cast<std::size_t>(pc)];const double bj=out.basis_.values()[static_cast<std::size_t>(pc)];
            for(int pr=pb;pr<pe;++pr){
                const int row=out.basis_.row_ind()[static_cast<std::size_t>(pr)];const double bi=out.basis_.values()[static_cast<std::size_t>(pr)];
                const int pos=find_csc_position(row,col);if(pos<0) throw std::runtime_error("missing compiled nonlinear scatter entry");entries.emplace_back(pos,bi*bj);
            }
        }
        out.spring_scatter_[static_cast<std::size_t>(j)]=std::move(entries);
    }

    for(auto& w:out.walls_)for(int a=0;a<6;++a)for(int b=0;b<6;++b)
        if(w.dofs[a]>=0&&w.dofs[b]>=0){w.scatter[6*a+b]=find_csc_position(w.dofs[a],w.dofs[b]);if(w.scatter[6*a+b]<0)throw std::runtime_error("wall scatter missing");}
    for(auto& e:out.steel_members_)for(int a=0;a<6;++a)for(int b=0;b<6;++b)
        if(e.dofs[a]>=0&&e.dofs[b]>=0){e.scatter[6*a+b]=find_csc_position(e.dofs[a],e.dofs[b]);if(e.scatter[6*a+b]<0)throw std::runtime_error("steel member scatter missing");}
    for(auto& d:out.dampers_)for(int a=0;a<6;++a)for(int b=0;b<6;++b)
        if(d.dofs[a]>=0&&d.dofs[b]>=0&&d.b[a]!=0.0&&d.b[b]!=0.0){d.scatter[6*a+b]=find_csc_position(d.dofs[a],d.dofs[b]);if(d.scatter[6*a+b]<0)throw std::runtime_error("viscous damper scatter missing");}

    int response_id=response_node_id_;
    if(response_id<0) response_id=nodes_.back().id;
    out.response_dof_=out.reduced_dof(response_id,Dof2D::UX);
    if(out.response_dof_<0) throw std::invalid_argument("response node UX is constrained");
    auto stories=story_node_ids_;
    if(stories.empty()) stories.push_back(response_id);
    for(int id:stories){
        const int rd=out.reduced_dof(id,Dof2D::UX);
        if(rd<0) throw std::invalid_argument("story node UX is constrained");
        out.story_dofs_.push_back(rd);
        const auto it=out.node_slot_.find(id);
        if(it==out.node_slot_.end()) throw std::invalid_argument("unknown story node");
        out.story_y_.push_back(nodes_[static_cast<std::size_t>(it->second)].y);
    }
    for(std::size_t i=1;i<out.story_y_.size();++i)
        if(out.story_y_[i]<=out.story_y_[i-1]) throw std::invalid_argument("2D story nodes must have increasing elevation");
    return out;
}

int CompiledFrame2D::reduced_dof(int node_id, Dof2D dof) const {
    auto it=node_slot_.find(node_id);
    if(it==node_slot_.end()) return -1;
    return full_to_reduced_[static_cast<std::size_t>(dof_index(it->second,dof))];
}

int CompiledFrame2D::elastic_element_index(int element_id) const {
    const auto it=element_index_by_id_.find(element_id);
    return it==element_index_by_id_.end()?-1:it->second;
}

int CompiledFrame2D::nonlinear_component_index(int spring_id) const {
    const auto it=spring_index_by_id_.find(spring_id);
    return it==spring_index_by_id_.end()?-1:it->second;
}

int CompiledFrame2D::nonlinear_component_id(int index) const {
    if(index<0||index>=nonlinear_count()) return -1;
    return spring_ids_[static_cast<std::size_t>(index)];
}

NonlinearComponent2DKind CompiledFrame2D::nonlinear_component_kind(int id) const{
    const int i=nonlinear_component_index(id);
    if(i<0)throw std::invalid_argument("unknown 2D nonlinear component id");
    return component_kinds_[static_cast<std::size_t>(i)];
}

NonlinearComponentSnapshot CompiledFrame2D::nonlinear_component_snapshot(
    int spring_id,const std::vector<double>& u,const std::vector<double>& committed_state) const {
    if(static_cast<int>(u.size())!=n_||static_cast<int>(committed_state.size())!=nonlinear_state_size_)
        throw std::invalid_argument("2D nonlinear component snapshot state size");
    const int idx=nonlinear_component_index(spring_id);
    if(idx<0) throw std::invalid_argument("unknown 2D nonlinear component id");
    const double q=basis_.column_dot(idx,u);
    const int o=material_state_offsets_[static_cast<std::size_t>(idx)];
    std::vector<double> trial(static_cast<std::size_t>(materials_[static_cast<std::size_t>(idx)].state_size()));
    const auto tr=materials_[static_cast<std::size_t>(idx)].trial(q,committed_state.data()+o,trial.data());
    return {spring_id,idx,q,tr.force,tr.tangent,tr.diagnostics};
}

void CompiledFrame2D::replace_nonlinear_material(int spring_id,NonlinearMaterial material){
    replace_nonlinear_material_field({{spring_id,std::move(material)}});
}

void CompiledFrame2D::replace_nonlinear_material_field(
    const std::vector<std::pair<int,NonlinearMaterial>>& replacements){
    if(replacements.empty()) return;
    struct Pending{int index{};NonlinearMaterial material{BilinearSpring(1.0,1.0,0.0)};double old_k{},new_k{};};
    std::vector<Pending> pending;pending.reserve(replacements.size());
    std::set<int> seen;
    for(const auto& [id,material]:replacements){
        const int idx=nonlinear_component_index(id);
        if(idx<0) throw std::invalid_argument("unknown 2D nonlinear component id in material field: "+std::to_string(id));
        if(!seen.insert(idx).second) throw std::invalid_argument("duplicate 2D nonlinear component replacement: "+std::to_string(id));
        if(material.state_size()!=materials_[static_cast<std::size_t>(idx)].state_size())
            throw std::invalid_argument("2D material-field replacement changes nonlinear state layout");
        pending.push_back({idx,material,initial_tangents_[static_cast<std::size_t>(idx)],material.initial_stiffness()});
    }
    for(auto& r:pending){
        const double dk=r.new_k-r.old_k;
        if(dk!=0.0) for(const auto& [pos,c]:spring_scatter_[static_cast<std::size_t>(r.index)])
            k_initial_.values()[static_cast<std::size_t>(pos)]+=dk*c;
        materials_[static_cast<std::size_t>(r.index)]=std::move(r.material);
        initial_tangents_[static_cast<std::size_t>(r.index)]=r.new_k;
    }
}

ElasticFrame2DResponse CompiledFrame2D::elastic_element_response(
    int element_id,const std::vector<double>& u) const {
    if(static_cast<int>(u.size())!=n_) throw std::invalid_argument("2D element response vector size");
    const int idx=elastic_element_index(element_id);
    if(idx<0) throw std::invalid_argument("unknown 2D elastic element id");
    const auto& e=element_response_data_[static_cast<std::size_t>(idx)];
    const double dx=e.xj-e.xi,dy=e.yj-e.yi,L=std::hypot(dx,dy);
    const double c=dx/L,sn=dy/L;
    std::array<double,6> ug{};
    for(int d=0;d<6;++d){const int r=e.reduced_dofs[static_cast<std::size_t>(d)]; if(r>=0) ug[static_cast<std::size_t>(d)]=u[static_cast<std::size_t>(r)];}
    std::array<double,6> q{};
    q[0]= c*ug[0]+sn*ug[1]; q[1]=-sn*ug[0]+c*ug[1]; q[2]=ug[2];
    q[3]= c*ug[3]+sn*ug[4]; q[4]=-sn*ug[3]+c*ug[4]; q[5]=ug[5];
    const double ea=e.properties.E*e.properties.A/L;
    const double p=e.properties.axial_compression-ea*(q[3]-q[0]);
    const double EI=e.properties.E*e.properties.I;
    const double k12=12.0*EI/(L*L*L),k6=6.0*EI/(L*L),k4=4.0*EI/L,k2=2.0*EI/L;
    std::array<double,4> qb{q[1],q[2],q[4],q[5]},fb{};
    const double kb[4][4]={{k12,k6,-k12,k6},{k6,k4,-k6,k2},{-k12,-k6,k12,-k6},{k6,k2,-k6,k4}};
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)fb[static_cast<std::size_t>(i)]+=kb[i][j]*qb[static_cast<std::size_t>(j)];
    if(p>0.0){
        const double fac=p/(30.0*L);
        const double g[4][4]={{36.0,3.0*L,-36.0,3.0*L},{3.0*L,4.0*L*L,-3.0*L,-L*L},{-36.0,-3.0*L,36.0,-3.0*L},{3.0*L,-L*L,-3.0*L,4.0*L*L}};
        for(int i=0;i<4;++i)for(int j=0;j<4;++j)fb[static_cast<std::size_t>(i)]-=fac*g[i][j]*qb[static_cast<std::size_t>(j)];
    }
    return {element_id,p,fb[0],fb[2],fb[1],fb[3]};
}

std::vector<double> CompiledFrame2D::mass_multiply(const std::vector<double>& a) const {
    if (static_cast<int>(a.size()) != n_) throw std::invalid_argument("mass vector size");
    std::vector<double> out(static_cast<std::size_t>(n_));
    for (int i=0;i<n_;++i) out[static_cast<std::size_t>(i)] = mass_[static_cast<std::size_t>(i)] * a[static_cast<std::size_t>(i)];
    return out;
}

std::vector<double> CompiledFrame2D::damping_multiply(const std::vector<double>& v) const {
    if(static_cast<int>(v.size())!=n_) throw std::invalid_argument("damping vector size");
    auto kv=k_initial_.multiply(v);
    std::vector<double> out(static_cast<std::size_t>(n_));
    for(int i=0;i<n_;++i) out[static_cast<std::size_t>(i)] = alpha_m_*mass_[static_cast<std::size_t>(i)]*v[static_cast<std::size_t>(i)] + beta_k_*kv[static_cast<std::size_t>(i)];
    for(const auto& d:dampers_){
        double rate=0.0;for(int a=0;a<6;++a)if(d.dofs[a]>=0)rate+=d.b[a]*v[static_cast<std::size_t>(d.dofs[a])];
        const auto r=d.element.trial(rate);
        for(int a=0;a<6;++a)if(d.dofs[a]>=0)out[static_cast<std::size_t>(d.dofs[a])]+=d.b[a]*r.force;
    }
    return out;
}

SparseMatrixCSC CompiledFrame2D::effective_initial_matrix(double a0,double a1) const {
    auto out=k_initial_;
    const double stiffness_scale=1.0+a1*beta_k_;
    for(double& v:out.values()) v*=stiffness_scale;
    for(int i=0;i<n_;++i)
        out.values()[static_cast<std::size_t>(diagonal_positions_[static_cast<std::size_t>(i)])] +=
            (a0+a1*alpha_m_)*mass_[static_cast<std::size_t>(i)];
    for(const auto& d:dampers_)for(int a=0;a<6;++a)for(int b=0;b<6;++b){
        const int pos=d.scatter[6*a+b];if(pos>=0)out.values()[static_cast<std::size_t>(pos)]+=a1*d.c0*d.b[a]*d.b[b];
    }
    return out;
}

SparseMatrixCSC CompiledFrame2D::effective_tangent_matrix(const std::vector<double>& tangents,
                                                          double a0,double a1) const {
    if(static_cast<int>(tangents.size())!=nonlinear_count()) throw std::invalid_argument("tangent vector size");
    auto out=effective_initial_matrix(a0,a1);
    for(int j=0;j<nonlinear_count();++j){
        const double dk=tangents[static_cast<std::size_t>(j)]-initial_tangents_[static_cast<std::size_t>(j)];
        if(dk==0.0) continue;
        for(const auto& [pos,coeff]:spring_scatter_[static_cast<std::size_t>(j)])
            out.values()[static_cast<std::size_t>(pos)] += dk*coeff;
    }
    return out;
}

std::vector<double> CompiledFrame2D::initial_nonlinear_state() const{
    std::vector<double> state(static_cast<std::size_t>(nonlinear_state_size_),0.0);
    for(int j=0;j<nonlinear_count();++j){
        const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        materials_[static_cast<std::size_t>(j)].initialize_state(state.data()+o);
    }
    for(const auto& w:walls_){auto s=w.element.initial_state();std::copy(s.begin(),s.end(),state.begin()+w.state_offset);}
    for(const auto& e:steel_members_){auto s=e.element.initial_state();std::copy(s.begin(),s.end(),state.begin()+e.state_offset);}
    return state;
}

void CompiledFrame2D::evaluate_nonlinear_deformations(const std::vector<double>& q,const std::vector<double>& committed_state,std::vector<double>& component_forces,std::vector<double>& tangents,std::vector<double>& trial_state) const{
    if(!walls_.empty()||!steel_members_.empty())throw std::invalid_argument("scalar deformation-bank reduction does not support coupled wall or steel-member elements");
    if(static_cast<int>(q.size())!=nonlinear_count()||static_cast<int>(committed_state.size())!=nonlinear_state_size())throw std::invalid_argument("2D nonlinear state dimension mismatch");
    component_forces.resize(static_cast<std::size_t>(nonlinear_count()));tangents.resize(static_cast<std::size_t>(nonlinear_count()));trial_state.resize(static_cast<std::size_t>(nonlinear_state_size()));
    for(int j=0;j<nonlinear_count();++j){
        const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        const auto tr=materials_[static_cast<std::size_t>(j)].trial(q[static_cast<std::size_t>(j)],committed_state.data()+o,trial_state.data()+o);
        component_forces[static_cast<std::size_t>(j)]=tr.force;tangents[static_cast<std::size_t>(j)]=tr.tangent;
    }
}
void CompiledFrame2D::internal_force_and_tangent(const std::vector<double>& u,const std::vector<double>& committed_state,std::vector<double>& force,std::vector<double>& tangents,std::vector<double>& trial_state) const{
    internal_force_and_tangent_diagnostics(u,committed_state,force,tangents,trial_state,nullptr);
}
void CompiledFrame2D::internal_force_and_tangent_diagnostics(const std::vector<double>& u,const std::vector<double>& committed_state,std::vector<double>& force,std::vector<double>& tangents,std::vector<double>& trial_state,NonlinearEvalDiagnostics* diagnostics) const{
    if(static_cast<int>(u.size())!=n_||static_cast<int>(committed_state.size())!=nonlinear_state_size())throw std::invalid_argument("2D state dimension mismatch");
    force=k_linear_.multiply(u);tangents.resize(static_cast<std::size_t>(nonlinear_count()));trial_state.resize(static_cast<std::size_t>(nonlinear_state_size()));
    for(const auto& w:walls_){
        std::array<double,6> q{};for(int a=0;a<6;++a)if(w.dofs[a]>=0)q[a]=u[w.dofs[a]];
        auto r=w.element.trial(q,committed_state.data()+w.state_offset);
        for(int a=0;a<6;++a)if(w.dofs[a]>=0)force[w.dofs[a]]+=r.force[a];
        std::copy(r.state.begin(),r.state.end(),trial_state.begin()+w.state_offset);
        if(diagnostics){++diagnostics->component_evaluations;++diagnostics->full_state_evaluations;}
    }
    for(const auto& e:steel_members_){
        std::array<double,6> q{};for(int a=0;a<6;++a)if(e.dofs[a]>=0)q[a]=u[static_cast<std::size_t>(e.dofs[a])];
        auto r=e.element.trial(q,committed_state.data()+e.state_offset);
        for(int a=0;a<6;++a)if(e.dofs[a]>=0)force[static_cast<std::size_t>(e.dofs[a])]+=r.force[a];
        std::copy(r.state.begin(),r.state.end(),trial_state.begin()+e.state_offset);
        if(diagnostics){
            diagnostics->component_evaluations+=2;diagnostics->full_state_evaluations+=2;
            for(const auto& d:r.hinge_diagnostics){
                diagnostics->active_tangent_evaluations+=static_cast<std::size_t>(d.tangent_active);
                diagnostics->fast_path_evaluations+=static_cast<std::size_t>(d.fast_path);
                diagnostics->transition_events+=static_cast<std::size_t>(d.transition);
                diagnostics->reversal_events+=static_cast<std::size_t>(d.reversal);
                diagnostics->deterioration_events+=static_cast<std::size_t>(d.deterioration);
                diagnostics->failure_events+=static_cast<std::size_t>(d.failed);
                diagnostics->io_or_beyond_evaluations+=static_cast<std::size_t>(d.at_or_beyond_io);
                diagnostics->ls_or_beyond_evaluations+=static_cast<std::size_t>(d.at_or_beyond_ls);
                diagnostics->cp_or_beyond_evaluations+=static_cast<std::size_t>(d.at_or_beyond_cp);
                diagnostics->beyond_cp_evaluations+=static_cast<std::size_t>(d.beyond_cp);
                diagnostics->lateral_loss_evaluations+=static_cast<std::size_t>(d.lateral_resistance_lost);
            }
        }
    }
    for(int j=0;j<nonlinear_count();++j){
        const double q=basis_.column_dot(j,u);const int o=material_state_offsets_[static_cast<std::size_t>(j)];
        const auto tr=materials_[static_cast<std::size_t>(j)].trial(q,committed_state.data()+o,trial_state.data()+o);
        tangents[static_cast<std::size_t>(j)]=tr.tangent;basis_.axpy_column(j,tr.force,force);
        if(diagnostics){
            ++diagnostics->component_evaluations;
            diagnostics->active_tangent_evaluations+=static_cast<std::size_t>(tr.diagnostics.tangent_active);
            diagnostics->fast_path_evaluations+=static_cast<std::size_t>(tr.diagnostics.fast_path);
            diagnostics->full_state_evaluations+=static_cast<std::size_t>(!tr.diagnostics.fast_path);
            diagnostics->transition_events+=static_cast<std::size_t>(tr.diagnostics.transition);
            diagnostics->reversal_events+=static_cast<std::size_t>(tr.diagnostics.reversal);
            diagnostics->deterioration_events+=static_cast<std::size_t>(tr.diagnostics.deterioration);
            diagnostics->failure_events+=static_cast<std::size_t>(tr.diagnostics.failed);
            diagnostics->io_or_beyond_evaluations+=static_cast<std::size_t>(tr.diagnostics.at_or_beyond_io);
            diagnostics->ls_or_beyond_evaluations+=static_cast<std::size_t>(tr.diagnostics.at_or_beyond_ls);
            diagnostics->cp_or_beyond_evaluations+=static_cast<std::size_t>(tr.diagnostics.at_or_beyond_cp);
            diagnostics->beyond_cp_evaluations+=static_cast<std::size_t>(tr.diagnostics.beyond_cp);
            diagnostics->lateral_loss_evaluations+=static_cast<std::size_t>(tr.diagnostics.lateral_resistance_lost);
        }
    }
}

std::vector<int> CompiledFrame2D::wall_ids() const {std::vector<int> r;for(const auto& w:walls_)r.push_back(w.id);return r;}
Wall2DResponse CompiledFrame2D::wall_response(int id,const std::vector<double>& u,const std::vector<double>& s) const {
    if(static_cast<int>(u.size())!=n_||static_cast<int>(s.size())!=nonlinear_state_size_)throw std::invalid_argument("wall response vector size");
    for(const auto& w:walls_)if(w.id==id){std::array<double,6> q{};for(int a=0;a<6;++a)if(w.dofs[a]>=0)q[a]=u[w.dofs[a]];return w.element.trial(q,s.data()+w.state_offset);}
    throw std::invalid_argument("unknown wall id");
}

std::vector<int> CompiledFrame2D::steel_member_ids() const{std::vector<int> r;for(const auto& e:steel_members_)r.push_back(e.id);return r;}
SteelMember2DRole CompiledFrame2D::steel_member_role(int id) const{for(const auto& e:steel_members_)if(e.id==id)return e.role;throw std::invalid_argument("unknown steel member id");}
SteelMember2DResponse CompiledFrame2D::steel_member_response(int id,const std::vector<double>& u,const std::vector<double>& s) const{
    if(static_cast<int>(u.size())!=n_||static_cast<int>(s.size())!=nonlinear_state_size_)throw std::invalid_argument("steel member response vector size");
    for(const auto& e:steel_members_)if(e.id==id){std::array<double,6> q{};for(int a=0;a<6;++a)if(e.dofs[a]>=0)q[a]=u[static_cast<std::size_t>(e.dofs[a])];return e.element.trial(q,s.data()+e.state_offset);}
    throw std::invalid_argument("unknown steel member id");
}

std::vector<int> CompiledFrame2D::viscous_damper_ids() const{std::vector<int> r;for(const auto& d:dampers_)r.push_back(d.id);return r;}
ViscousDamper2DTrial CompiledFrame2D::viscous_damper_response(int id,const std::vector<double>& v) const{
    if(static_cast<int>(v.size())!=n_)throw std::invalid_argument("viscous damper response vector size");
    for(const auto& d:dampers_)if(d.id==id){double rate=0.0;for(int a=0;a<6;++a)if(d.dofs[a]>=0)rate+=d.b[a]*v[static_cast<std::size_t>(d.dofs[a])];return d.element.trial(rate);}
    throw std::invalid_argument("unknown viscous damper id");
}
SparseMatrixCSC CompiledFrame2D::effective_state_tangent_matrix(const std::vector<double>& u,const std::vector<double>& t,double a0,double a1) const {
    return effective_state_tangent_matrix_with_state(u,t,initial_nonlinear_state(),a0,a1);
}
SparseMatrixCSC CompiledFrame2D::effective_state_tangent_matrix_with_state(const std::vector<double>& u,const std::vector<double>& t,const std::vector<double>& s,double a0,double a1) const {
    if(static_cast<int>(u.size())!=n_||static_cast<int>(s.size())!=nonlinear_state_size_)throw std::invalid_argument("wall tangent vector size");
    auto K=effective_tangent_matrix(t,a0,a1);
    for(const auto& w:walls_){std::array<double,6> q{};for(int a=0;a<6;++a)if(w.dofs[a]>=0)q[a]=u[w.dofs[a]];auto r=w.element.trial(q,s.data()+w.state_offset);
        for(int a=0;a<36;++a)if(w.scatter[a]>=0)K.values()[w.scatter[a]]+=r.tangent[a]-w.k0[a];}
    for(const auto& e:steel_members_){std::array<double,6> q{};for(int a=0;a<6;++a)if(e.dofs[a]>=0)q[a]=u[static_cast<std::size_t>(e.dofs[a])];auto r=e.element.trial(q,s.data()+e.state_offset);
        for(int a=0;a<36;++a)if(e.scatter[a]>=0)K.values()[static_cast<std::size_t>(e.scatter[a])]+=r.tangent[a]-e.k0[a];}
    return K;
}

SparseMatrixCSC CompiledFrame2D::effective_state_tangent_matrix_with_state_and_velocity(
    const std::vector<double>& u,const std::vector<double>& t,const std::vector<double>& s,
    const std::vector<double>& v,double a0,double a1) const{
    if(static_cast<int>(v.size())!=n_)throw std::invalid_argument("velocity-dependent tangent vector size");
    auto K=effective_state_tangent_matrix_with_state(u,t,s,a0,a1);
    for(const auto& d:dampers_){
        double rate=0.0;for(int a=0;a<6;++a)if(d.dofs[a]>=0)rate+=d.b[a]*v[static_cast<std::size_t>(d.dofs[a])];
        const double dc=d.element.trial(rate).tangent-d.c0;
        if(dc==0.0)continue;
        for(int a=0;a<6;++a)for(int b=0;b<6;++b){const int pos=d.scatter[6*a+b];if(pos>=0)K.values()[static_cast<std::size_t>(pos)]+=a1*dc*d.b[a]*d.b[b];}
    }
    return K;
}

std::vector<double> CompiledFrame2D::base_excitation(double ground_accel) const {
    std::vector<double> p(static_cast<std::size_t>(n_));
    for(int i=0;i<n_;++i) p[static_cast<std::size_t>(i)] = -base_mass_x_[static_cast<std::size_t>(i)]*ground_accel;
    return p;
}

double CompiledFrame2D::response_value(const std::vector<double>& u) const {
    if(static_cast<int>(u.size())!=n_) throw std::invalid_argument("response vector size");
    return u[static_cast<std::size_t>(response_dof_)];
}

std::vector<double> CompiledFrame2D::story_response_values(const std::vector<double>& x) const {
    if(static_cast<int>(x.size())!=n_) throw std::invalid_argument("2D story response size");
    std::vector<double> out;out.reserve(story_dofs_.size());
    for(int rd:story_dofs_)out.push_back(x[static_cast<std::size_t>(rd)]);
    return out;
}

double CompiledFrame2D::max_drift_measure(const std::vector<double>& u) const {
    if(static_cast<int>(u.size())!=n_) throw std::invalid_argument("response vector size");
    double lower=0.0, out=0.0;
    for(int rd:story_dofs_){
        const double ui=u[static_cast<std::size_t>(rd)];
        out=std::max(out,std::abs(ui-lower));
        lower=ui;
    }
    return out;
}

double CompiledFrame2D::max_drift_ratio(const std::vector<double>& u) const {
    if(static_cast<int>(u.size())!=n_) throw std::invalid_argument("response vector size");
    if(story_dofs_.empty()||story_y_.size()!=story_dofs_.size())return std::numeric_limits<double>::quiet_NaN();
    double lower_u=0.0,lower_y=0.0,out=0.0;
    for(std::size_t i=0;i<story_dofs_.size();++i){
        const double ui=u[static_cast<std::size_t>(story_dofs_[i])];
        const double dy=story_y_[i]-lower_y;if(dy<=0.0)return std::numeric_limits<double>::quiet_NaN();
        out=std::max(out,std::abs(ui-lower_u)/dy);lower_u=ui;lower_y=story_y_[i];
    }
    return out;
}

} // namespace quake
