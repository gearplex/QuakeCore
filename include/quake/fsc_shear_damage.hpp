#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace quake {

// State-only diagnostic for flexure-shear-critical RC columns.
//
// This object intentionally does NOT modify structural resistance.  It is used
// to replay accepted response histories before promoting a shear-failure model
// into the Newton system.  The initiation surface follows the Ghannoum-Moehle
// rotation formulation (plastic-rotation form), with an optional nominal-shear
// trigger.  Post-initiation cycle degradation follows the calibrated
// PinchingLimitStateMaterial cyclic-strength coefficient, but is reported as a
// normalized damage/retention metric rather than claimed as a constitutive
// shear spring.
struct FSCShearDamageParams {
    double b_in{6.0};
    double d_in{4.8};
    double h_in{6.0};
    double clear_length_in{39.0};
    double tie_spacing_in{4.0};
    double longitudinal_steel_area_in2{0.88};
    double confined_concrete_area_in2{18.0};
    double fc_ksi{3.57};
    double fy_ksi{64.0};
    double nominal_shear_limit_kip{-1.0}; // <=0 disables force trigger
    double residual_strength_ratio{0.20}; // OpenSees default when Vr=-1
    double sign_crossing_fraction{0.05};  // suppress zero-force chatter
};

struct FSCShearDamageState {
    bool initiated{false};
    bool residual_reached{false};
    std::size_t initiation_step{static_cast<std::size_t>(-1)};
    std::size_t residual_step{static_cast<std::size_t>(-1)};
    double initiation_time_s{0.0};
    double initiation_rotation_rad{0.0};
    double initiation_rotation_limit_rad{0.0};
    double initiation_axial_kip{0.0};
    double initiation_shear_kip{0.0};
    int initiation_cause{0}; // 1 rotation, 2 shear, 3 simultaneous
    std::size_t opposite_sign_crossings{0};
    double cyclic_strength_coefficient{0.0};
    double cumulative_strength_decrement_kip{0.0};
    double retained_strength_ratio{1.0};
    double minimum_retained_strength_ratio{1.0};
};

class FSCShearDamageDiagnostic {
public:
    explicit FSCShearDamageDiagnostic(FSCShearDamageParams params = {})
        : p_(params) {
        state_.cyclic_strength_coefficient = cyclic_strength_coefficient(p_);
    }

    static double plastic_rotation_limit_rad(const FSCShearDamageParams& p,
                                             double axial_compression_kip,
                                             double shear_kip) {
        const double P = std::max(0.0, axial_compression_kip);
        const double V = std::abs(shear_kip);
        const double Ag = p.b_in * p.h_in;
        const double fc_psi = p.fc_ksi * 1000.0;
        const double v_psi = V / (p.b_in * p.d_in) * 1000.0;
        const double s_over_d = p.tie_spacing_in / p.d_in;
        const double p_ratio = P / (Ag * p.fc_ksi);
        const double vsqrt = v_psi / std::sqrt(fc_psi);
        return std::max(0.0, 0.032 - 0.014*s_over_d - 0.017*p_ratio - 0.0016*vsqrt);
    }

    static double cyclic_strength_coefficient(const FSCShearDamageParams& p) {
        const double Ag = p.b_in * p.h_in;
        const double a = 0.5 * p.clear_length_in;
        const double raw = 0.037133
            + 0.251204 * (p.fy_ksi * p.longitudinal_steel_area_in2 / (p.fc_ksi * Ag))
            - 0.354989 * (p.confined_concrete_area_in2 / Ag)
            + 0.056569 * (a / p.d_in);
        return std::max(0.0, raw);
    }

    const FSCShearDamageState& update(std::size_t step, double time_s,
                                      double bottom_rotation_rad,
                                      double top_rotation_rad,
                                      double axial_compression_kip,
                                      double signed_shear_kip) {
        const double local_rot = std::max(std::abs(bottom_rotation_rad), std::abs(top_rotation_rad));
        const double rot_lim = plastic_rotation_limit_rad(p_, axial_compression_kip, signed_shear_kip);
        const bool rot_hit = rot_lim > 0.0 && local_rot >= rot_lim;
        const bool shear_hit = p_.nominal_shear_limit_kip > 0.0 &&
                               std::abs(signed_shear_kip) >= p_.nominal_shear_limit_kip;

        if (!state_.initiated && (rot_hit || shear_hit)) {
            state_.initiated = true;
            state_.initiation_step = step;
            state_.initiation_time_s = time_s;
            state_.initiation_rotation_rad = local_rot;
            state_.initiation_rotation_limit_rad = rot_lim;
            state_.initiation_axial_kip = std::max(0.0, axial_compression_kip);
            state_.initiation_shear_kip = signed_shear_kip;
            state_.initiation_cause = rot_hit && shear_hit ? 3 : (rot_hit ? 1 : 2);
            failure_strength_kip_ = std::max(std::abs(signed_shear_kip),
                                             p_.nominal_shear_limit_kip > 0.0 ? p_.nominal_shear_limit_kip : 0.0);
            failure_strength_kip_ = std::max(failure_strength_kip_, 1e-12);
            residual_strength_kip_ = std::clamp(p_.residual_strength_ratio,0.0,1.0)*failure_strength_kip_;
            const int s = force_sign(signed_shear_kip, failure_strength_kip_);
            if (s != 0) {
                lobe_sign_ = s;
                lobe_peak_kip_ = std::abs(signed_shear_kip);
            }
        }

        if (!state_.initiated) return state_;

        const int s = force_sign(signed_shear_kip, failure_strength_kip_);
        if (s != 0) {
            if (lobe_sign_ == 0) {
                lobe_sign_ = s;
                lobe_peak_kip_ = std::abs(signed_shear_kip);
            } else if (s == lobe_sign_) {
                lobe_peak_kip_ = std::max(lobe_peak_kip_, std::abs(signed_shear_kip));
            } else {
                // A completed excursion into the opposite force sign.  The
                // OpenSees calibrated material reduces the global envelope by
                // |stress at unload|*dmgStrengthCyclic on this transition.
                state_.cumulative_strength_decrement_kip +=
                    lobe_peak_kip_ * state_.cyclic_strength_coefficient;
                ++state_.opposite_sign_crossings;
                lobe_sign_ = s;
                lobe_peak_kip_ = std::abs(signed_shear_kip);

                const double raw_remaining = failure_strength_kip_ - state_.cumulative_strength_decrement_kip;
                state_.retained_strength_ratio =
                    std::max(residual_strength_kip_, raw_remaining) / failure_strength_kip_;
                state_.minimum_retained_strength_ratio =
                    std::min(state_.minimum_retained_strength_ratio, state_.retained_strength_ratio);
                if (!state_.residual_reached && raw_remaining <= residual_strength_kip_) {
                    state_.residual_reached = true;
                    state_.residual_step = step;
                }
            }
        }
        return state_;
    }

    const FSCShearDamageState& state() const noexcept { return state_; }

private:
    int force_sign(double v, double ref) const {
        const double eps = std::max(1e-9, p_.sign_crossing_fraction * ref);
        if (v > eps) return 1;
        if (v < -eps) return -1;
        return 0;
    }

    FSCShearDamageParams p_{};
    FSCShearDamageState state_{};
    double failure_strength_kip_{0.0};
    double residual_strength_kip_{0.0};
    int lobe_sign_{0};
    double lobe_peak_kip_{0.0};
};

} // namespace quake
