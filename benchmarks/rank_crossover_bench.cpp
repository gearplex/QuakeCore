#include "quake/low_rank_solver.hpp"
#include "quake/sparse.hpp"
#include "quake/superlu_solver.hpp"
#include "quake/update_basis.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace quake;

static double err(const std::vector<double>&a,const std::vector<double>&b){double x=0,d=0;for(std::size_t i=0;i<a.size();++i){x=std::max(x,std::abs(a[i]-b[i]));d=std::max(d,std::abs(a[i]));}return x/std::max(1.0,d);}
int main(int argc,char**argv){const int n=argc>1?std::stoi(argv[1]):3000;const int m=argc>2?std::stoi(argv[2]):300;const int reps=argc>3?std::stoi(argv[3]):20;
    std::vector<Triplet> t;for(int i=0;i<n;++i){t.push_back({i,i,8.0});if(i+1<n){t.push_back({i,i+1,-1.0});t.push_back({i+1,i,-1.0});}}
    auto A=SparseMatrixCSC::from_triplets(n,n,t);std::vector<std::vector<std::pair<int,double>>> cols(static_cast<std::size_t>(m));for(int j=0;j<m;++j){int i=1+(j*(n-2))/std::max(1,m-1);cols[static_cast<std::size_t>(j)]={{i-1,-1.0},{i,1.0}};}auto B=SparseUpdateBasis::from_columns(n,cols);
    LazyLowRankWoodburySolver wood(A,B);std::vector<double> rhs(static_cast<std::size_t>(n));for(int i=0;i<n;++i)rhs[static_cast<std::size_t>(i)]=std::sin(0.01*i)+0.2*std::cos(0.037*i);
    std::cout<<std::fixed<<std::setprecision(6)<<"n="<<n<<" candidates="<<m<<" reps="<<reps<<"\n";
    for(int rank:{1,5,10,20,40,80,120,200}){if(rank>m)continue;std::vector<double> dk(static_cast<std::size_t>(m),0.0);for(int j=0;j<rank;++j)dk[static_cast<std::size_t>(j)]=-2.5;
        std::vector<Triplet> tt=t;for(int j=0;j<rank;++j){const auto& c=cols[static_cast<std::size_t>(j)];for(auto [r,br]:c)for(auto [q,bq]:c)tt.push_back({r,q,dk[static_cast<std::size_t>(j)]*br*bq});}auto K=SparseMatrixCSC::from_triplets(n,n,tt,-1.0);
        SuperLUSamePatternSolver same(A); // pattern is tridiagonal, same as device stamps
        // warm both and verify
        auto xw=wood.solve(rhs,dk);auto xd=superlu_solve_once(K,rhs);double e=err(xd,xw);
        auto ts=std::chrono::steady_clock::now();for(int r=0;r<reps;++r)(void)same.refactor_and_solve(K,rhs);double ds=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();
        ts=std::chrono::steady_clock::now();for(int r=0;r<reps;++r)(void)wood.solve(rhs,dk);double dw=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();
        std::cout<<"rank="<<rank<<" same_s="<<ds<<" wood_s="<<dw<<" speedup="<<ds/dw<<" rel_err="<<e<<" cached="<<wood.cached_update_count()<<"\n";
    }
}
