#include "quake/soil2d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace quake {
namespace {
void positive_finite(double x,const char* name){
    if(!std::isfinite(x)||x<=0.0)throw std::invalid_argument(std::string(name)+" must be finite and positive");
}
SoilBilinearBackbone backbone(double capacity,double displacement_50){
    positive_finite(capacity,"nodal soil capacity");
    positive_finite(displacement_50,"soil displacement_50");
    // A transparent bilinear approximation: the elastic branch reaches half
    // the nodal ultimate resistance at the supplied 50-percent displacement.
    return {capacity,0.5*capacity/displacement_50};
}
}

AsymmetricElasticPerfectlyPlasticSpring::AsymmetricElasticPerfectlyPlasticSpring(
    double k,double positive,double negative)
    :stiffness_(k),positive_capacity_(positive),negative_capacity_(negative){
    positive_finite(k,"asymmetric spring stiffness");
    if(!std::isfinite(positive)||!std::isfinite(negative)||positive<0.0||negative<0.0||
       (positive==0.0&&negative==0.0))
        throw std::invalid_argument("asymmetric spring capacities must be finite, nonnegative, and not both zero");
}

AsymmetricElasticPerfectlyPlasticTrial AsymmetricElasticPerfectlyPlasticSpring::trial(
    double q,double committed_plastic) const {
    if(!std::isfinite(q)||!std::isfinite(committed_plastic))
        throw std::invalid_argument("nonfinite asymmetric spring state");
    const double elastic=stiffness_*(q-committed_plastic);
    if(elastic>positive_capacity_)
        return {positive_capacity_,0.0,q-positive_capacity_/stiffness_};
    if(elastic<-negative_capacity_)
        return {-negative_capacity_,0.0,q+negative_capacity_/stiffness_};
    return {elastic,stiffness_,committed_plastic};
}

SoilBilinearBackbone py_bilinear_backbone(double p,double L,double y50){
    positive_finite(p,"ultimate resistance per length");positive_finite(L,"tributary length");
    return backbone(p*L,y50);
}
SoilBilinearBackbone tz_bilinear_backbone(double t,double perimeter,double L,double z50){
    positive_finite(t,"ultimate interface stress");positive_finite(perimeter,"pile perimeter");
    positive_finite(L,"tributary length");return backbone(t*perimeter*L,z50);
}
SoilBilinearBackbone qz_bilinear_backbone(double q,double area,double z50){
    positive_finite(q,"ultimate bearing stress");positive_finite(area,"toe area");
    return backbone(q*area,z50);
}

std::vector<double> nodal_tributary_lengths(const std::vector<double>& depths){
    if(depths.size()<2)throw std::invalid_argument("tributary depths require at least two nodes");
    for(std::size_t i=0;i<depths.size();++i){
        if(!std::isfinite(depths[i])||depths[i]<0.0||(i&&depths[i]<=depths[i-1]))
            throw std::invalid_argument("tributary depths must be finite, nonnegative, and strictly increasing");
    }
    std::vector<double> out(depths.size());
    out.front()=0.5*(depths[1]-depths[0]);out.back()=0.5*(depths.back()-depths[depths.size()-2]);
    for(std::size_t i=1;i+1<depths.size();++i)out[i]=0.5*(depths[i+1]-depths[i-1]);
    return out;
}

} // namespace quake
