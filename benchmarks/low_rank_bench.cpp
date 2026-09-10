#include "quake/low_rank_solver.hpp"
#include "quake/sparse.hpp"
#include "quake/superlu_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
using namespace quake;

static SparseMatrixCSC grid_matrix(int nx,int ny,double shift){
    const int n=nx*ny; std::vector<Triplet> t;
    auto id=[nx](int x,int y){return y*nx+x;};
    for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){
        int i=id(x,y); double diag=shift;
        auto link=[&](int xx,int yy){ if(xx<0||xx>=nx||yy<0||yy>=ny)return; int j=id(xx,yy); t.push_back({i,j,-1.0}); diag+=1.0;};
        link(x-1,y);link(x+1,y);link(x,y-1);link(x,y+1); t.push_back({i,i,diag});
    }
    return SparseMatrixCSC::from_triplets(n,n,t,1e-18);
}

static SparseMatrixCSC updated(const SparseMatrixCSC&A,const std::vector<double>&B,int m,const std::vector<double>&d){
    std::vector<Triplet> t; t.reserve(static_cast<std::size_t>(A.nnz()+4*m));
    for(int c=0;c<A.cols();++c)for(int p=A.col_ptr()[c];p<A.col_ptr()[c+1];++p)t.push_back({A.row_ind()[p],c,A.values()[p]});
    const int n=A.rows();
    for(int j=0;j<m;++j){ if(d[j]==0)continue; std::vector<int> ids;std::vector<double> vals; for(int i=0;i<n;++i){double v=B[static_cast<std::size_t>(j*n+i)];if(v!=0){ids.push_back(i);vals.push_back(v);}} for(std::size_t a=0;a<ids.size();++a)for(std::size_t b=0;b<ids.size();++b)t.push_back({ids[a],ids[b],d[j]*vals[a]*vals[b]}); }
    return SparseMatrixCSC::from_triplets(n,n,t,1e-18);
}

int main(int argc,char**argv){
    int nx=argc>1?std::atoi(argv[1]):45; int ny=argc>2?std::atoi(argv[2]):45; int solves=argc>3?std::atoi(argv[3]):100; int active_rank=argc>4?std::atoi(argv[4]):5;
    const int n=nx*ny; const int m=std::min(200,nx*(ny-1));
    auto A=grid_matrix(nx,ny,0.5);
    std::vector<double>B(static_cast<std::size_t>(n*m),0.0);
    for(int j=0;j<m;++j){ int x=(j*17)%nx; int y=(j*31)%(ny-1); int i0=y*nx+x,i1=(y+1)*nx+x; B[static_cast<std::size_t>(j*n+i1)]=1;B[static_cast<std::size_t>(j*n+i0)]=-1; }
    std::vector<double> rhs(static_cast<std::size_t>(n)); for(int i=0;i<n;++i)rhs[i]=std::sin(0.013*i)+0.2*std::cos(0.071*i);
    std::mt19937 rng(42); std::uniform_int_distribution<int> pick(0,m-1);
    std::vector<std::vector<double>> ds; ds.reserve(solves);
    for(int s=0;s<solves;++s){std::vector<double>d(static_cast<std::size_t>(m),0);for(int r=0;r<active_rank;++r)d[static_cast<std::size_t>(pick(rng))]=-0.35-0.05*(r%3);ds.push_back(std::move(d));}

    auto t0=std::chrono::steady_clock::now();
    std::vector<double> last_full;
    for(const auto&d:ds){auto K=updated(A,B,m,d);last_full=superlu_solve_once(K,rhs);} auto t1=std::chrono::steady_clock::now();
    auto t2=std::chrono::steady_clock::now(); LowRankWoodburySolver w(A,B,m); auto tsetup=std::chrono::steady_clock::now(); std::vector<double> last_w; for(const auto&d:ds)last_w=w.solve(rhs,d); auto t3=std::chrono::steady_clock::now();
    double full=std::chrono::duration<double>(t1-t0).count(); double wood_setup=std::chrono::duration<double>(tsetup-t2).count(); double wood_solve=std::chrono::duration<double>(t3-tsetup).count(); double wood=wood_setup+wood_solve;
    double e=0,den=0;for(int i=0;i<n;++i){e=std::max(e,std::abs(last_full[i]-last_w[i]));den=std::max(den,std::abs(last_full[i]));}
    std::cout<<std::fixed<<std::setprecision(6)
             <<"grid_dof="<<n<<" nnz="<<A.nnz()<<" candidate_updates="<<m<<" active_rank="<<active_rank<<" solves="<<solves<<"\n"
             <<"full_factorization_seconds="<<full<<"\nwoodbury_setup_seconds="<<wood_setup<<"\nwoodbury_repeated_solve_seconds="<<wood_solve<<"\nwoodbury_total_seconds="<<wood<<"\nspeedup_including_setup="<<(full/wood)<<"x\n"
             <<"full_global_factorizations="<<solves<<"\nwoodbury_global_factorizations=1\nrelative_solution_error="<<(e/std::max(1.0,den))<<"\n";
}
