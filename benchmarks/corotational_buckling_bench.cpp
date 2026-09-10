#include "quake/frame3d.hpp"
#include "quake/stability.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace quake;

static CompiledFrame3D column(int elements,double P){
    Frame3DBuilder b;const double L=4.0;
    for(int i=0;i<=elements;++i)b.add_node(i+1,0,0,L*i/elements);
    b.fix(1);
    for(int i=1;i<=elements;++i){const int n=i+1;b.fix_dof(n,Dof3D::UY);b.fix_dof(n,Dof3D::RX);b.fix_dof(n,Dof3D::RZ);}
    for(int i=0;i<elements;++i)b.add_elastic_frame(i+1,i+1,i+2,30000,12000,2.0,0.3,0.5,0.5,1,0,0,P);
    b.set_response(elements+1,Dof3D::UX);b.set_story_nodes({elements+1},Dof3D::UX);b.set_corotational(true);
    return b.compile();
}
static bool stable(int elements,double P){return assess_positive_definiteness(column(elements,P).K_initial(),1e-9,0).positive_definite;}
int main(){
    const double E=30000.0,I=0.5,L=4.0,exact=M_PI*M_PI*E*I/(4.0*L*L);
    std::cout<<std::fixed<<std::setprecision(6)<<"Euler_Pcr="<<exact<<"\n";
    for(int n:{1,2,4,8,16}){
        double lo=0.0,hi=1.5*exact;
        for(int k=0;k<28;++k){const double mid=0.5*(lo+hi);if(stable(n,mid))lo=mid;else hi=mid;}
        const double p=0.5*(lo+hi);
        std::cout<<"elements="<<n<<" predicted_Pcr="<<p<<" ratio="<<p/exact<<" error_pct="<<100.0*(p/exact-1.0)<<"\n";
    }
}
