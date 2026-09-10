#pragma once

#include <cstddef>

namespace quake {

// Independent peak-oriented IMK-family material used for QuakeCore research.
// Parameter semantics follow the common Ibarra-Medina-Krawinkler/OpenSees
// convention: Up = pre-capping plastic deformation, Upc = post-capping
// deformation, Uu = ultimate deformation, FcapFy/FresFy are force ratios,
// and Lambda/c control cyclic deterioration energy/cyclic exponents.
//
// This implementation is intentionally kept separate from OpenSees source and
// is validated from envelope/cyclic invariants until external parity tests are
// available.
struct IMKPeakOrientedParams {
    double Ke{};
    double posUp{}, posUpc{}, posUu{}, posFy{}, posFcapFy{}, posFresFy{};
    double negUp{}, negUpc{}, negUu{}, negFy{}, negFcapFy{}, negFresFy{};
    double lambdaS{}, lambdaC{}, lambdaA{}, lambdaK{};
    double cS{1.0}, cC{1.0}, cA{1.0}, cK{1.0};
    double Dpos{1.0}, Dneg{1.0};
};

enum IMKEvent : unsigned {
    IMK_EVENT_NONE          = 0u,
    IMK_EVENT_FAST_ELASTIC  = 1u << 0,
    IMK_EVENT_YIELD         = 1u << 1,
    IMK_EVENT_CAP           = 1u << 2,
    IMK_EVENT_REVERSAL      = 1u << 3,
    IMK_EVENT_DETERIORATION = 1u << 4,
    IMK_EVENT_FAILURE       = 1u << 5,
    IMK_EVENT_RELOAD        = 1u << 6
};

struct IMKTrialResult {
    double force{};
    double tangent{};
    unsigned events{IMK_EVENT_NONE};
};

class IMKPeakOrientedMaterial {
public:
    static constexpr int kStateSize = 16;

    explicit IMKPeakOrientedMaterial(IMKPeakOrientedParams params);

    double initial_stiffness() const { return p_.Ke; }
    int state_size() const { return kStateSize; }
    const IMKPeakOrientedParams& params() const { return p_; }

    void initialize_state(double* state) const;
    IMKTrialResult trial(double deformation, const double* committed,
                         double* trial_state) const;

private:
    struct BackbonePoint { double force{}, tangent{}; unsigned events{}; bool failed{}; };
    BackbonePoint backbone(double u, bool positive,
                           double strength_scale, double cap_scale) const;
    double deterioration(double excursion_energy, double cumulative_energy,
                         double lambda, double exponent) const;

    IMKPeakOrientedParams p_;
};

} // namespace quake
