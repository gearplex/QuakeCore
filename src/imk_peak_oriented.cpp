#include "quake/imk_peak_oriented.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace quake {
namespace {
constexpr double eps = 1e-14;
inline double sgn(double x) { return x >= 0.0 ? 1.0 : -1.0; }
inline double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }
}

IMKPeakOrientedMaterial::IMKPeakOrientedMaterial(IMKPeakOrientedParams params)
    : p_(params) {
    if (p_.Ke <= 0.0 || p_.posUp <= 0.0 || p_.posUpc <= 0.0 || p_.posUu <= 0.0 ||
        p_.negUp <= 0.0 || p_.negUpc <= 0.0 || p_.negUu <= 0.0 ||
        p_.posFy <= 0.0 || p_.negFy <= 0.0 || p_.posFcapFy <= 0.0 ||
        p_.negFcapFy <= 0.0 || p_.posFresFy < 0.0 || p_.negFresFy < 0.0 ||
        p_.posFresFy > p_.posFcapFy || p_.negFresFy > p_.negFcapFy ||
        p_.cS <= 0.0 || p_.cC <= 0.0 || p_.cA <= 0.0 || p_.cK <= 0.0 ||
        p_.Dpos < 0.0 || p_.Dneg < 0.0) {
        throw std::invalid_argument("invalid IMK peak-oriented parameters");
    }
}

void IMKPeakOrientedMaterial::initialize_state(double* s) const {
    std::fill(s, s + kStateSize, 0.0);
    // 0 u, 1 f, 2 cumulative energy, 3 excursion energy, 4 unload K,
    // 5/6 positive peak u/f, 7/8 negative peak u/f,
    // 9/10 +/- strength scale, 11/12 +/- cap scale,
    // 13 last direction, 14 envelope flag, 15 failed flag.
    s[4] = p_.Ke;
    s[9] = s[10] = s[11] = s[12] = 1.0;
    s[14] = 1.0;
}

double IMKPeakOrientedMaterial::deterioration(double excursion_energy,
                                               double cumulative_energy,
                                               double lambda,
                                               double exponent) const {
    if (lambda <= 0.0 || excursion_energy <= 0.0) return 0.0;
    // The common IMK deterioration reference energy is Lambda*Fy. Use the
    // positive initial yield force as the reference normalization, consistent
    // with the standard symmetric parameterization.
    const double ref = lambda * p_.posFy;
    const double remaining = ref - cumulative_energy;
    if (remaining <= eps) return 1.0;
    return clamp01(std::pow(std::max(0.0, excursion_energy / remaining), exponent));
}

IMKPeakOrientedMaterial::BackbonePoint
IMKPeakOrientedMaterial::backbone(double u, bool positive,
                                  double strength_scale, double cap_scale) const {
    BackbonePoint r;
    const double sign = positive ? 1.0 : -1.0;
    const double x = std::abs(u);
    const double fy0 = positive ? p_.posFy : p_.negFy;
    const double up = positive ? p_.posUp : p_.negUp;
    const double upc = positive ? p_.posUpc : p_.negUpc;
    const double uu = positive ? p_.posUu : p_.negUu;
    const double fcap_ratio = positive ? p_.posFcapFy : p_.negFcapFy;
    const double fres_ratio = positive ? p_.posFresFy : p_.negFresFy;

    const double fy = std::max(1e-12 * fy0, fy0 * strength_scale);
    const double uy = fy / p_.Ke;
    const double fcap = std::max(fy * 1e-9, fy0 * fcap_ratio * cap_scale);
    const double fres = fy0 * fres_ratio;
    const double kp = (fcap - fy) / up;
    const double ucap = uy + up;
    const double kpc = (fres - fcap) / upc;
    const double ures = ucap + upc;

    if (x >= uu) {
        r.force = 0.0; r.tangent = 0.0; r.events |= IMK_EVENT_FAILURE; r.failed = true;
    } else if (x <= uy) {
        r.force = sign * p_.Ke * x; r.tangent = p_.Ke;
    } else if (x <= ucap) {
        r.force = sign * (fy + kp * (x - uy)); r.tangent = kp; r.events |= IMK_EVENT_YIELD;
    } else if (x <= ures) {
        const double mag = std::max(fres, fcap + kpc * (x - ucap));
        r.force = sign * mag;
        r.tangent = (mag <= fres + 1e-12 * std::max(1.0, fy0)) ? 0.0 : kpc;
        r.events |= IMK_EVENT_YIELD | IMK_EVENT_CAP;
    } else {
        r.force = sign * fres; r.tangent = 0.0; r.events |= IMK_EVENT_YIELD | IMK_EVENT_CAP;
    }
    return r;
}

IMKTrialResult IMKPeakOrientedMaterial::trial(double u, const double* c, double* s) const {
    std::copy(c, c + kStateSize, s);
    IMKTrialResult out;
    if (c[15] > 0.5) {
        s[0] = u; s[1] = 0.0; out.events = IMK_EVENT_FAILURE; return out;
    }

    const double uprev = c[0], fprev = c[1];
    const double du = u - uprev;
    if (std::abs(du) <= eps) {
        out.force = fprev;
        // Preserve a consistent local tangent at zero increment. Virgin
        // elastic components remain on the active-front fast path even when
        // their generalized deformation is exactly unchanged.
        const bool virgin = c[2] <= eps && c[15] < 0.5 &&
            std::abs(fprev-p_.Ke*uprev) <= 1e-10*std::max(1.0,p_.posFy);
        out.tangent = virgin ? p_.Ke : ((std::abs(c[13]) < 0.5) ? p_.Ke : std::max(0.0, c[4]));
        if(virgin) out.events |= IMK_EVENT_FAST_ELASTIC;
        return out;
    }

    const double dir = sgn(du);
    const double last_dir = c[13];
    const bool reversal = std::abs(last_dir) > 0.5 && dir * last_dir < 0.0;

    // Active-front fast path: a virgin component that remains inside the
    // initial elastic domain needs no peak/reversal/deterioration machinery.
    const double fy_dir = u >= 0.0 ? p_.posFy : p_.negFy;
    const double uy_dir = fy_dir / p_.Ke;
    const bool virgin_elastic = c[2] <= eps && c[15] < 0.5 &&
        std::abs(fprev - p_.Ke*uprev) <= 1e-10*std::max(1.0,fy_dir) &&
        std::abs(u) <= uy_dir + 1e-14;
    if (virgin_elastic) {
        out.force = p_.Ke*u; out.tangent = p_.Ke; out.events = IMK_EVENT_FAST_ELASTIC;
        s[0]=u; s[1]=out.force; s[4]=p_.Ke; s[13]=dir; s[14]=1.0;
        if(u>s[5] && out.force>=0.0){s[5]=u;s[6]=out.force;}
        if(u<s[7] && out.force<=0.0){s[7]=u;s[8]=out.force;}
        return out;
    }
    double unload_k = std::clamp(c[4], 1e-9 * p_.Ke, p_.Ke);
    double pos_strength = c[9], neg_strength = c[10];
    double pos_cap = c[11], neg_cap = c[12];
    double pos_peak_u = c[5], pos_peak_f = c[6];
    double neg_peak_u = c[7], neg_peak_f = c[8];
    double excursion_energy = c[3];
    unsigned events = IMK_EVENT_NONE;

    if (reversal) {
        events |= IMK_EVENT_REVERSAL;
        const double betaK = deterioration(excursion_energy, c[2], p_.lambdaK, p_.cK);
        if (betaK > 0.0) {
            unload_k = std::max(1e-9 * p_.Ke, unload_k * (1.0 - betaK));
            events |= IMK_EVENT_DETERIORATION;
        }
        const double betaS = deterioration(excursion_energy, c[2], p_.lambdaS, p_.cS);
        const double betaC = deterioration(excursion_energy, c[2], p_.lambdaC, p_.cC);
        const double betaA = deterioration(excursion_energy, c[2], p_.lambdaA, p_.cA);
        const double D = dir > 0.0 ? p_.Dpos : p_.Dneg;
        if (dir > 0.0) {
            pos_strength = std::max(1e-6, pos_strength * (1.0 - betaS * D));
            pos_cap = std::max(1e-6, pos_cap * (1.0 - betaC * D));
            if (pos_peak_u > 0.0) pos_peak_u *= (1.0 + betaA * D);
        } else {
            neg_strength = std::max(1e-6, neg_strength * (1.0 - betaS * D));
            neg_cap = std::max(1e-6, neg_cap * (1.0 - betaC * D));
            if (neg_peak_u < 0.0) neg_peak_u *= (1.0 + betaA * D);
        }
        if (betaS > 0.0 || betaC > 0.0 || betaA > 0.0) events |= IMK_EVENT_DETERIORATION;
        excursion_energy = 0.0;
    }

    const auto env = backbone(u, u >= 0.0, u >= 0.0 ? pos_strength : neg_strength,
                              u >= 0.0 ? pos_cap : neg_cap);

    // Peak-oriented cyclic rule. Continue on the envelope while extending the
    // current peak. Otherwise unload from the committed point and reload toward
    // the previously reached peak in the opposite direction.
    bool envelope_loading = false;
    if (dir > 0.0) envelope_loading = (u >= pos_peak_u - 1e-14 && fprev >= -1e-12);
    else envelope_loading = (u <= neg_peak_u + 1e-14 && fprev <= 1e-12);

    if (std::abs(last_dir) < 0.5 || (!reversal && c[14] > 0.5 && envelope_loading)) {
        out.force = env.force; out.tangent = env.tangent; events |= env.events;
        s[14] = 1.0;
    } else {
        const double u0 = uprev - fprev / unload_k;
        const bool before_zero = (dir > 0.0) ? (u <= u0) : (u >= u0);
        if ((dir > 0.0 && fprev < 0.0 && before_zero) ||
            (dir < 0.0 && fprev > 0.0 && before_zero)) {
            out.force = fprev + unload_k * du; out.tangent = unload_k;
            events |= IMK_EVENT_RELOAD; s[14] = 0.0;
        } else {
            double target_u = dir > 0.0 ? pos_peak_u : neg_peak_u;
            double target_f = dir > 0.0 ? pos_peak_f : neg_peak_f;
            if ((dir > 0.0 && target_u <= u0 + eps) || (dir < 0.0 && target_u >= u0 - eps) || std::abs(target_f) <= eps) {
                // No prior peak in this direction: use the current envelope as
                // the peak-oriented target and grow onto the backbone.
                target_u = u;
                target_f = env.force;
            }
            const bool beyond_peak = dir > 0.0 ? (u >= target_u) : (u <= target_u);
            if (beyond_peak) {
                out.force = env.force; out.tangent = env.tangent; events |= env.events;
                s[14] = 1.0;
            } else {
                const double denom = target_u - u0;
                double kreload = std::abs(denom) > eps ? target_f / denom : unload_k;
                // A peak-oriented reload should not be stiffer than the initial elastic branch.
                kreload = std::clamp(kreload, 1e-9 * p_.Ke, p_.Ke);
                out.force = kreload * (u - u0); out.tangent = kreload;
                events |= IMK_EVENT_RELOAD; s[14] = 0.0;
            }
        }
    }

    // Update peaks only when the trial point extends the attained response in
    // that direction. This state is committed only after Newton convergence.
    if (u > pos_peak_u && out.force >= 0.0) { pos_peak_u = u; pos_peak_f = out.force; }
    if (u < neg_peak_u && out.force <= 0.0) { neg_peak_u = u; neg_peak_f = out.force; }

    // Approximate dissipated (rather than recoverable elastic) energy. This
    // keeps virgin elastic cycling from consuming deterioration capacity.
    const double work = std::abs(0.5 * (fprev + out.force) * du);
    const double elastic_change = std::abs((out.force*out.force - fprev*fprev)/(2.0*p_.Ke));
    const double dE = std::max(0.0, work - elastic_change);
    s[0] = u; s[1] = out.force; s[2] = c[2] + dE; s[3] = excursion_energy + dE;
    s[4] = unload_k; s[5] = pos_peak_u; s[6] = pos_peak_f; s[7] = neg_peak_u; s[8] = neg_peak_f;
    s[9] = pos_strength; s[10] = neg_strength; s[11] = pos_cap; s[12] = neg_cap;
    s[13] = dir;
    if (env.failed && s[14] > 0.5) { s[15] = 1.0; out.force = 0.0; out.tangent = 0.0; events |= IMK_EVENT_FAILURE; }

    out.events = events;
    return out;
}

} // namespace quake
