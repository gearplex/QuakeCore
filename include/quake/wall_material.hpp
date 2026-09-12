#pragma once
#include "quake/bilinear.hpp"
#include <array>
#include <memory>
#include <vector>

namespace quake {
// Stress/strain laws used by wall fibers. Parameters use consistent stress units.
// All histories are explicit values: trial evaluation never mutates the material.
class WallUniaxial {
public:
    enum class Kind { Elastic, SteelBilinear, Concrete01, MinMax, Parallel };
    static WallUniaxial elastic(double E);
    static WallUniaxial steel(double E, double fy, double b);
    static WallUniaxial concrete01(double fc, double epsc, double fcu, double epsu);
    static WallUniaxial minmax(WallUniaxial material, double min_strain, double max_strain);
    static WallUniaxial parallel(std::vector<WallUniaxial> materials);
    struct Result { double stress{}, tangent{}; };
    Kind kind() const { return kind_; }
    double initial_tangent() const;
    int state_size() const;
    void initialize(double* s) const;
    Result trial(double strain, const double* committed, double* trial_state) const;
private:
    Kind kind_{Kind::Elastic};
    double E_{1}, fy_{1}, b_{0}, fc_{-1}, ec_{-.002}, fu_{0}, eu_{-.006};
    double min_strain_{}, max_strain_{};
    std::shared_ptr<const std::vector<WallUniaxial>> children_;
};

// General plane-stress panel: elastic isotropic background plus oriented
// uniaxial layers. A fixed-angle RC panel is represented by two orthogonal
// concrete struts and horizontal/vertical steel layers. This is NOT FSAM.
struct WallPanelLayer {
    double angle_rad{};
    double weight{};
    WallUniaxial material{WallUniaxial::elastic(1)};
};
struct WallPanelResult {
    std::array<double,3> stress{}; // sigma_x, sigma_y, tau_xy
    std::array<double,9> tangent{}; // derivatives w.r.t. ex,ey,engineering gamma
    std::vector<double> state;
};
class WallPanel {
public:
    WallPanel(double background_E, double poisson, std::vector<WallPanelLayer> layers = {});
    int state_size() const { return state_size_; }
    std::vector<double> initial_state() const;
    std::array<double,9> initial_tangent() const;
    WallPanelResult trial(const std::array<double,3>& strain, const double* committed) const;
private:
    std::array<double,9> background_{};
    std::vector<WallPanelLayer> layers_;
    std::vector<std::array<double,3>> directions_;
    std::vector<int> offsets_;
    int state_size_{};
};
} // namespace quake
