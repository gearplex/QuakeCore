#pragma once

#include <vector>

namespace quake {

struct AsymmetricElasticPerfectlyPlasticTrial {
    double force{};
    double tangent{};
    double plastic_deformation{};
};

// One-dimensional return-mapping law with independently bounded positive and
// negative resistance. A zero capacity on either side creates a no-resistance
// opening branch and accumulates the corresponding permanent gap in the single
// committed plastic-deformation state variable.
class AsymmetricElasticPerfectlyPlasticSpring {
public:
    AsymmetricElasticPerfectlyPlasticSpring(double stiffness,
                                            double positive_capacity,
                                            double negative_capacity);
    double initial_stiffness() const { return stiffness_; }
    double positive_capacity() const { return positive_capacity_; }
    double negative_capacity() const { return negative_capacity_; }
    AsymmetricElasticPerfectlyPlasticTrial trial(double deformation,
                                                  double committed_plastic) const;
private:
    double stiffness_{};
    double positive_capacity_{};
    double negative_capacity_{};
};

struct SoilBilinearBackbone {
    double nodal_capacity{};
    double initial_stiffness{};
};

// The adapters intentionally perform only transparent force conversion. They
// do not infer geotechnical parameters and are not PySimple1/TzSimple1/QzSimple1.
SoilBilinearBackbone py_bilinear_backbone(double ultimate_resistance_per_length,
                                          double tributary_length,
                                          double displacement_50);
SoilBilinearBackbone tz_bilinear_backbone(double ultimate_interface_stress,
                                          double pile_perimeter,
                                          double tributary_length,
                                          double displacement_50);
SoilBilinearBackbone qz_bilinear_backbone(double ultimate_bearing_stress,
                                          double toe_area,
                                          double displacement_50);

// Tributary lengths for nodes on a line. End nodes receive half an adjacent
// segment and interior nodes receive half of each adjacent segment.
std::vector<double> nodal_tributary_lengths(const std::vector<double>& depths);

} // namespace quake
