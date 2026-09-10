#include "quake/substructure.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace quake {

BlockPartition make_story_block_partition(const CompiledFrame3D& model,int stories_per_block,bool retain_nonlinear_support){
    if(stories_per_block<1) throw std::invalid_argument("stories_per_block must be positive");
    const int n=model.dof(); const auto& z=model.reduced_dof_elevations();
    if(static_cast<int>(z.size())!=n) throw std::invalid_argument("missing reduced DOF elevations");
    std::vector<double> levels; for(double v:z) if(std::isfinite(v)) levels.push_back(v);
    std::sort(levels.begin(),levels.end());
    levels.erase(std::unique(levels.begin(),levels.end(),[](double a,double b){return std::abs(a-b)<1e-9;}),levels.end());
    if(levels.empty()) throw std::invalid_argument("story partition has no finite elevations");
    std::vector<char> support(static_cast<std::size_t>(n),0), boundary_level(levels.size(),0), interface(static_cast<std::size_t>(n),0);
    if(retain_nonlinear_support) for(int r:model.nonlinear_basis().row_ind()) support[static_cast<std::size_t>(r)]=1;
    for(std::size_t i=static_cast<std::size_t>(stories_per_block-1);i<levels.size();i+=static_cast<std::size_t>(stories_per_block)) boundary_level[i]=1;
    if(!boundary_level.empty()) boundary_level[boundary_level.size()-1]=1;
    auto level_index=[&](double zz){auto it=std::lower_bound(levels.begin(),levels.end(),zz-1e-9);if(it==levels.end()||std::abs(*it-zz)>1e-7)return -1;return static_cast<int>(it-levels.begin());};
    const int block_count=(static_cast<int>(levels.size())+stories_per_block-1)/stories_per_block;
    BlockPartition p; p.interior_blocks.resize(static_cast<std::size_t>(block_count));
    for(int d=0;d<n;++d){
        const int li=std::isfinite(z[static_cast<std::size_t>(d)])?level_index(z[static_cast<std::size_t>(d)]):-1;
        if(li<0 || support[static_cast<std::size_t>(d)] || boundary_level[static_cast<std::size_t>(li)]){interface[static_cast<std::size_t>(d)]=1;p.interface_dofs.push_back(d);}
        else p.interior_blocks[static_cast<std::size_t>(li/stories_per_block)].push_back(d);
    }
    p.interior_blocks.erase(std::remove_if(p.interior_blocks.begin(),p.interior_blocks.end(),[](const auto& b){return b.empty();}),p.interior_blocks.end());
    return p;
}

struct BlockSchurFactor::Impl {
    struct Block {
        std::vector<int> dofs;
        std::vector<int> boundary_global;
        std::vector<int> boundary_interface;
        std::vector<int> global_to_local; // n, -1 except block dofs
        std::unique_ptr<SuperLUFactor> factor;
        std::vector<double> Aib; // ni x nb, column-major
        std::vector<double> Abi; // nb x ni, row-major
        std::vector<double> W;   // Aii^-1 Aib, ni x nb column-major
    };
    int n{};
    std::vector<int> interface;
    std::vector<int> interface_map;
    std::vector<int> dof_block;
    std::vector<Block> blocks;
    std::unique_ptr<SuperLUFactor> schur;
    BlockSchurStats stats;
};

BlockSchurFactor::BlockSchurFactor(const SparseMatrixCSC& A,BlockPartition partition):impl_(std::make_unique<Impl>()){
    const auto t0=std::chrono::steady_clock::now();
    if(A.rows()!=A.cols()) throw std::invalid_argument("BlockSchurFactor requires square matrix");
    auto& q=*impl_; q.n=A.rows(); q.interface=std::move(partition.interface_dofs);
    std::sort(q.interface.begin(),q.interface.end());q.interface.erase(std::unique(q.interface.begin(),q.interface.end()),q.interface.end());
    q.interface_map.assign(static_cast<std::size_t>(q.n),-1);q.dof_block.assign(static_cast<std::size_t>(q.n),-2);
    for(int i=0;i<static_cast<int>(q.interface.size());++i){int d=q.interface[static_cast<std::size_t>(i)];if(d<0||d>=q.n)throw std::out_of_range("interface DOF");q.interface_map[static_cast<std::size_t>(d)]=i;q.dof_block[static_cast<std::size_t>(d)]=-1;}
    q.blocks.resize(partition.interior_blocks.size());
    for(int g=0;g<static_cast<int>(partition.interior_blocks.size());++g){
        auto& b=q.blocks[static_cast<std::size_t>(g)];b.dofs=std::move(partition.interior_blocks[static_cast<std::size_t>(g)]);std::sort(b.dofs.begin(),b.dofs.end());
        b.global_to_local.assign(static_cast<std::size_t>(q.n),-1);
        for(int i=0;i<static_cast<int>(b.dofs.size());++i){int d=b.dofs[static_cast<std::size_t>(i)];if(d<0||d>=q.n||q.dof_block[static_cast<std::size_t>(d)]!=-2)throw std::invalid_argument("duplicate/unassigned block DOF");q.dof_block[static_cast<std::size_t>(d)]=g;b.global_to_local[static_cast<std::size_t>(d)]=i;}
    }
    for(int d=0;d<q.n;++d) if(q.dof_block[static_cast<std::size_t>(d)]==-2){q.interface_map[static_cast<std::size_t>(d)]=static_cast<int>(q.interface.size());q.interface.push_back(d);q.dof_block[static_cast<std::size_t>(d)]=-1;}
    // Rebuild interface map after appending any uncovered DOFs.
    std::fill(q.interface_map.begin(),q.interface_map.end(),-1);for(int i=0;i<static_cast<int>(q.interface.size());++i)q.interface_map[static_cast<std::size_t>(q.interface[static_cast<std::size_t>(i)])]=i;

    // Identify local interface neighborhoods and reject cross-block interior couplings.
    std::vector<std::set<int>> bsets(q.blocks.size());
    for(int c=0;c<q.n;++c)for(int p=A.col_ptr()[static_cast<std::size_t>(c)];p<A.col_ptr()[static_cast<std::size_t>(c+1)];++p){int r=A.row_ind()[static_cast<std::size_t>(p)];int bc=q.dof_block[static_cast<std::size_t>(c)],br=q.dof_block[static_cast<std::size_t>(r)];if(br>=0&&bc>=0&&br!=bc&&std::abs(A.values()[static_cast<std::size_t>(p)])>1e-16)throw std::invalid_argument("story blocks have direct cross-block interior coupling");if(br>=0&&bc==-1)bsets[static_cast<std::size_t>(br)].insert(c);if(br==-1&&bc>=0)bsets[static_cast<std::size_t>(bc)].insert(r);}
    for(std::size_t g=0;g<q.blocks.size();++g){auto& b=q.blocks[g];b.boundary_global.assign(bsets[g].begin(),bsets[g].end());for(int d:b.boundary_global)b.boundary_interface.push_back(q.interface_map[static_cast<std::size_t>(d)]);}

    std::vector<Triplet> schur_t;
    // A_BB
    for(int c=0;c<q.n;++c){int ic=q.interface_map[static_cast<std::size_t>(c)];if(ic<0)continue;for(int p=A.col_ptr()[static_cast<std::size_t>(c)];p<A.col_ptr()[static_cast<std::size_t>(c+1)];++p){int r=A.row_ind()[static_cast<std::size_t>(p)],ir=q.interface_map[static_cast<std::size_t>(r)];if(ir>=0)schur_t.push_back({ir,ic,A.values()[static_cast<std::size_t>(p)]});}}

    for(std::size_t g=0;g<q.blocks.size();++g){auto& b=q.blocks[g];const int ni=static_cast<int>(b.dofs.size()),nb=static_cast<int>(b.boundary_global.size());if(ni==0)continue;
        std::unordered_map<int,int> blocal;for(int j=0;j<nb;++j)blocal[b.boundary_global[static_cast<std::size_t>(j)]]=j;
        std::vector<Triplet> ii; b.Aib.assign(static_cast<std::size_t>(ni*nb),0.0);b.Abi.assign(static_cast<std::size_t>(nb*ni),0.0);
        for(int c=0;c<q.n;++c)for(int p=A.col_ptr()[static_cast<std::size_t>(c)];p<A.col_ptr()[static_cast<std::size_t>(c+1)];++p){int r=A.row_ind()[static_cast<std::size_t>(p)];const double v=A.values()[static_cast<std::size_t>(p)];const int lc=b.global_to_local[static_cast<std::size_t>(c)],lr=b.global_to_local[static_cast<std::size_t>(r)];
            if(lc>=0&&lr>=0)ii.push_back({lr,lc,v});
            else if(lr>=0){auto it=blocal.find(c);if(it!=blocal.end())b.Aib[static_cast<std::size_t>(it->second*ni+lr)]+=v;}
            else if(lc>=0){auto it=blocal.find(r);if(it!=blocal.end())b.Abi[static_cast<std::size_t>(it->second*ni+lc)]+=v;}
        }
        auto Aii=SparseMatrixCSC::from_triplets(ni,ni,ii,1e-18);b.factor=std::make_unique<SuperLUFactor>(Aii);b.W=nb?b.factor->solve_multiple(b.Aib,nb):std::vector<double>{};
        for(int a=0;a<nb;++a)for(int bb=0;bb<nb;++bb){double v=0.0;for(int i=0;i<ni;++i)v+=b.Abi[static_cast<std::size_t>(a*ni+i)]*b.W[static_cast<std::size_t>(bb*ni+i)];if(std::abs(v)>1e-18)schur_t.push_back({b.boundary_interface[static_cast<std::size_t>(a)],b.boundary_interface[static_cast<std::size_t>(bb)],-v});}
    }
    auto S=SparseMatrixCSC::from_triplets(static_cast<int>(q.interface.size()),static_cast<int>(q.interface.size()),schur_t,1e-18);q.schur=std::make_unique<SuperLUFactor>(S);
    q.stats.full_dof=q.n;q.stats.interface_dof=static_cast<int>(q.interface.size());q.stats.interior_dof=q.n-q.stats.interface_dof;q.stats.blocks=static_cast<int>(q.blocks.size());q.stats.schur_nnz=S.nnz();q.stats.setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
}
BlockSchurFactor::~BlockSchurFactor()=default;BlockSchurFactor::BlockSchurFactor(BlockSchurFactor&&) noexcept=default;BlockSchurFactor& BlockSchurFactor::operator=(BlockSchurFactor&&) noexcept=default;
int BlockSchurFactor::size() const{return impl_->n;}const BlockSchurStats& BlockSchurFactor::stats() const{return impl_->stats;}
std::vector<double> BlockSchurFactor::solve(const std::vector<double>& rhs) const{
    const auto& q=*impl_;if(static_cast<int>(rhs.size())!=q.n)throw std::invalid_argument("BlockSchur RHS dimension");std::vector<double> rb(q.interface.size());for(std::size_t i=0;i<q.interface.size();++i)rb[i]=rhs[static_cast<std::size_t>(q.interface[i])];std::vector<std::vector<double>> ys(q.blocks.size());
    for(std::size_t g=0;g<q.blocks.size();++g){const auto& b=q.blocks[g];std::vector<double> ri(b.dofs.size());for(std::size_t i=0;i<b.dofs.size();++i)ri[i]=rhs[static_cast<std::size_t>(b.dofs[i])];ys[g]=b.factor->solve(ri);const int ni=static_cast<int>(b.dofs.size());for(int a=0;a<static_cast<int>(b.boundary_global.size());++a){double v=0;for(int i=0;i<ni;++i)v+=b.Abi[static_cast<std::size_t>(a*ni+i)]*ys[g][static_cast<std::size_t>(i)];rb[static_cast<std::size_t>(b.boundary_interface[static_cast<std::size_t>(a)])]-=v;}}
    auto xb=q.schur->solve(rb);std::vector<double>x(static_cast<std::size_t>(q.n));for(std::size_t i=0;i<q.interface.size();++i)x[static_cast<std::size_t>(q.interface[i])]=xb[i];
    for(std::size_t g=0;g<q.blocks.size();++g){const auto& b=q.blocks[g];const int ni=static_cast<int>(b.dofs.size()),nb=static_cast<int>(b.boundary_global.size());for(int i=0;i<ni;++i){double v=ys[g][static_cast<std::size_t>(i)];for(int j=0;j<nb;++j)v-=b.W[static_cast<std::size_t>(j*ni+i)]*xb[static_cast<std::size_t>(b.boundary_interface[static_cast<std::size_t>(j)])];x[static_cast<std::size_t>(b.dofs[static_cast<std::size_t>(i)])]=v;}}
    return x;
}
std::vector<double> BlockSchurFactor::solve_multiple(const std::vector<double>& rhs,int nrhs) const{if(nrhs<=0||static_cast<int>(rhs.size())!=size()*nrhs)throw std::invalid_argument("BlockSchur multi RHS dimension");std::vector<double>out(rhs.size());for(int c=0;c<nrhs;++c){std::vector<double>b(rhs.begin()+static_cast<std::ptrdiff_t>(c*size()),rhs.begin()+static_cast<std::ptrdiff_t>((c+1)*size()));auto x=solve(b);std::copy(x.begin(),x.end(),out.begin()+static_cast<std::ptrdiff_t>(c*size()));}return out;}

SubstructuredLazyWoodburySolver::SubstructuredLazyWoodburySolver(const SparseMatrixCSC& A,SparseUpdateBasis basis,BlockPartition p):n_(A.rows()),m_(basis.cols()),factor_(A,std::move(p)),basis_(std::move(basis)),index_to_cache_(static_cast<std::size_t>(m_),-1){if(A.rows()!=A.cols()||basis_.rows()!=n_)throw std::invalid_argument("substructured Woodbury dimension");}
void SubstructuredLazyWoodburySolver::cache_update(int idx){if(idx<0||idx>=m_)throw std::out_of_range("update index");if(index_to_cache_[static_cast<std::size_t>(idx)]>=0)return;auto w=factor_.solve(basis_.dense_column(idx));const int old=static_cast<int>(cached_indices_.size()),neu=old+1;std::vector<double>G(static_cast<std::size_t>(neu*neu));for(int i=0;i<old;++i)for(int j=0;j<old;++j)G[static_cast<std::size_t>(i*neu+j)]=G_cache_[static_cast<std::size_t>(i*old+j)];for(int c=0;c<old;++c){int cand=cached_indices_[static_cast<std::size_t>(c)];double a=basis_.column_dot(cand,w),b=0;for(int p=basis_.col_ptr()[static_cast<std::size_t>(idx)];p<basis_.col_ptr()[static_cast<std::size_t>(idx+1)];++p){int r=basis_.row_ind()[static_cast<std::size_t>(p)];b+=basis_.values()[static_cast<std::size_t>(p)]*W_cache_[static_cast<std::size_t>(c*n_+r)];}G[static_cast<std::size_t>(c*neu+old)]=a;G[static_cast<std::size_t>(old*neu+c)]=b;}G[static_cast<std::size_t>(old*neu+old)]=basis_.column_dot(idx,w);W_cache_.insert(W_cache_.end(),w.begin(),w.end());G_cache_=std::move(G);index_to_cache_[static_cast<std::size_t>(idx)]=old;cached_indices_.push_back(idx);}
std::vector<double> SubstructuredLazyWoodburySolver::solve(const std::vector<double>& rhs,const std::vector<double>& dk,double tol){if(static_cast<int>(rhs.size())!=n_||static_cast<int>(dk.size())!=m_)throw std::invalid_argument("substructured Woodbury solve dimension");auto y=factor_.solve(rhs);std::vector<int>a;for(int i=0;i<m_;++i)if(std::abs(dk[static_cast<std::size_t>(i)])>tol){cache_update(i);a.push_back(i);}if(a.empty())return y;const int r=static_cast<int>(a.size()),kc=static_cast<int>(cached_indices_.size());std::vector<double>S(static_cast<std::size_t>(r*r)),q(static_cast<std::size_t>(r));for(int i=0;i<r;++i){int ia=a[static_cast<std::size_t>(i)],ca=index_to_cache_[static_cast<std::size_t>(ia)];double di=dk[static_cast<std::size_t>(ia)];q[static_cast<std::size_t>(i)]=di*basis_.column_dot(ia,y);for(int j=0;j<r;++j){int ib=a[static_cast<std::size_t>(j)],cb=index_to_cache_[static_cast<std::size_t>(ib)];S[static_cast<std::size_t>(i*r+j)]=(i==j?1.0:0.0)+di*G_cache_[static_cast<std::size_t>(ca*kc+cb)];}}auto alpha=dense_solve(std::move(S),std::move(q),r);for(int i=0;i<n_;++i){double corr=0;for(int j=0;j<r;++j){int cj=index_to_cache_[static_cast<std::size_t>(a[static_cast<std::size_t>(j)])];corr+=W_cache_[static_cast<std::size_t>(cj*n_+i)]*alpha[static_cast<std::size_t>(j)];}y[static_cast<std::size_t>(i)]-=corr;}return y;}

} // namespace quake
