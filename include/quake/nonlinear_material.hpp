#pragma once

#include "quake/bilinear.hpp"
#include "quake/imk_peak_oriented.hpp"
#include "quake/asce41_hinge.hpp"
#include "quake/soil2d.hpp"
#include <variant>

namespace quake {

enum class NonlinearMaterialKind { Bilinear, IMKPeakOriented, ASCE41Hinge, AsymmetricElasticPerfectlyPlastic };

struct MaterialEvalDiagnostics {
    bool fast_path{};
    bool tangent_active{};
    bool transition{};
    bool reversal{};
    bool deterioration{};
    bool failed{};
    bool at_or_beyond_io{};
    bool at_or_beyond_ls{};
    bool at_or_beyond_cp{};
    bool beyond_cp{};
    bool lateral_resistance_lost{};
};

struct MaterialTrialResult {
    double force{};
    double tangent{};
    MaterialEvalDiagnostics diagnostics{};
};

// Read-only response of one compiled nonlinear component at a committed global
// state. The response is evaluated without mutating the supplied constitutive
// history and is intended for validation/comparison instrumentation.
struct NonlinearComponentSnapshot {
    int component_id{};
    int component_index{};
    double deformation{};
    double force{};
    double tangent{};
    MaterialEvalDiagnostics diagnostics{};
};

class NonlinearMaterial {
public:
    explicit NonlinearMaterial(BilinearSpring law);
    explicit NonlinearMaterial(IMKPeakOrientedMaterial law);
    explicit NonlinearMaterial(ASCE41HingeMaterial law);
    explicit NonlinearMaterial(AsymmetricElasticPerfectlyPlasticSpring law);

    NonlinearMaterialKind kind() const;
    int state_size() const;
    double initial_stiffness() const;
    const ASCE41HingeParams* asce41_params() const;
    void initialize_state(double* state) const;
    MaterialTrialResult trial(double deformation, const double* committed,
                              double* trial_state) const;

private:
    std::variant<BilinearSpring, IMKPeakOrientedMaterial, ASCE41HingeMaterial,
                 AsymmetricElasticPerfectlyPlasticSpring> law_;
};

} // namespace quake
