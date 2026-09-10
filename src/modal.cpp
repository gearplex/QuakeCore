#include "quake/modal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

extern "C" {
void dggev_(const char* jobvl, const char* jobvr, const int* n,
            double* a, const int* lda, double* b, const int* ldb,
            double* alphar, double* alphai, double* beta,
            double* vl, const int* ldvl, double* vr, const int* ldvr,
            double* work, const int* lwork, int* info);
}

namespace quake {
namespace {

std::vector<double> dense_from_csc(const SparseMatrixCSC& a) {
    if (a.rows() != a.cols()) throw std::invalid_argument("modal stiffness must be square");
    const int n = a.rows();
    std::vector<double> out(static_cast<std::size_t>(n*n), 0.0); // column-major
    for (int col = 0; col < n; ++col) {
        for (int p = a.col_ptr()[static_cast<std::size_t>(col)];
             p < a.col_ptr()[static_cast<std::size_t>(col+1)]; ++p) {
            const int row = a.row_ind()[static_cast<std::size_t>(p)];
            out[static_cast<std::size_t>(col*n + row)] += a.values()[static_cast<std::size_t>(p)];
        }
    }
    return out;
}

std::vector<double> dense_mass(const NonlinearDynamicModel& model) {
    const int n = model.dof();
    std::vector<double> out(static_cast<std::size_t>(n*n), 0.0); // column-major
    std::vector<double> e(static_cast<std::size_t>(n), 0.0);
    for (int col = 0; col < n; ++col) {
        std::fill(e.begin(), e.end(), 0.0);
        e[static_cast<std::size_t>(col)] = 1.0;
        const auto me = model.mass_multiply(e);
        if (static_cast<int>(me.size()) != n) throw std::runtime_error("mass multiply returned wrong size");
        for (int row = 0; row < n; ++row)
            out[static_cast<std::size_t>(col*n + row)] = me[static_cast<std::size_t>(row)];
    }
    return out;
}

} // namespace

std::vector<ModeShape> modal_analysis_from_stiffness(const NonlinearDynamicModel& model,
                                                     const SparseMatrixCSC& stiffness,
                                                     int requested_modes,
                                                     double beta_tolerance,
                                                     double imaginary_tolerance) {
    const int n = model.dof();
    if (n <= 0 || requested_modes <= 0) return {};
    if (beta_tolerance <= 0.0 || imaginary_tolerance <= 0.0)
        throw std::invalid_argument("modal tolerances must be positive");

    if(stiffness.rows()!=n||stiffness.cols()!=n) throw std::invalid_argument("modal tangent stiffness size");
    auto a = dense_from_csc(stiffness);
    auto b = dense_mass(model);
    std::vector<double> ar(static_cast<std::size_t>(n));
    std::vector<double> ai(static_cast<std::size_t>(n));
    std::vector<double> be(static_cast<std::size_t>(n));
    std::vector<double> vr(static_cast<std::size_t>(n*n));
    double vl_dummy = 0.0;
    const int ldvl = 1, ldvr = n, lda = n, ldb = n;
    const char jobvl = 'N', jobvr = 'V';
    int info = 0, lwork = -1;
    double work_query = 0.0;
    dggev_(&jobvl, &jobvr, &n, a.data(), &lda, b.data(), &ldb,
           ar.data(), ai.data(), be.data(), &vl_dummy, &ldvl,
           vr.data(), &ldvr, &work_query, &lwork, &info);
    if (info != 0) throw std::runtime_error("LAPACK dggev workspace query failed");
    lwork = std::max(8*n, static_cast<int>(std::ceil(work_query)));
    std::vector<double> work(static_cast<std::size_t>(lwork));

    // DGGEV overwrites A/B, and the workspace query is permitted to inspect them;
    // rebuild for the actual solve to keep behavior independent of LAPACK version.
    a = dense_from_csc(stiffness);
    b = dense_mass(model);
    dggev_(&jobvl, &jobvr, &n, a.data(), &lda, b.data(), &ldb,
           ar.data(), ai.data(), be.data(), &vl_dummy, &ldvl,
           vr.data(), &ldvr, work.data(), &lwork, &info);
    if (info != 0) throw std::runtime_error("LAPACK dggev failed, info=" + std::to_string(info));

    struct Candidate { double lambda; int index; };
    std::vector<Candidate> candidates;
    candidates.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double scale = std::max({1.0, std::abs(ar[static_cast<std::size_t>(i)]),
                                      std::abs(be[static_cast<std::size_t>(i)])});
        if (std::abs(be[static_cast<std::size_t>(i)]) <= beta_tolerance * scale) continue; // infinite mode
        if (std::abs(ai[static_cast<std::size_t>(i)]) > imaginary_tolerance * scale) continue;
        const double lambda = ar[static_cast<std::size_t>(i)] / be[static_cast<std::size_t>(i)];
        if (!std::isfinite(lambda) || lambda <= 0.0) continue;
        candidates.push_back({lambda, i});
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& x, const Candidate& y){return x.lambda < y.lambda;});
    if (static_cast<int>(candidates.size()) > requested_modes)
        candidates.resize(static_cast<std::size_t>(requested_modes));

    constexpr double two_pi = 6.283185307179586476925286766559;
    std::vector<ModeShape> modes;
    modes.reserve(candidates.size());
    for (const auto& c : candidates) {
        ModeShape mode;
        mode.eigenvalue = c.lambda;
        mode.omega = std::sqrt(c.lambda);
        mode.frequency_hz = mode.omega / two_pi;
        mode.period = two_pi / mode.omega;
        mode.shape.resize(static_cast<std::size_t>(n));
        double max_abs = 0.0;
        for (int row = 0; row < n; ++row) {
            const double v = vr[static_cast<std::size_t>(c.index*n + row)];
            mode.shape[static_cast<std::size_t>(row)] = v;
            max_abs = std::max(max_abs, std::abs(v));
        }
        if (max_abs > 0.0) for (double& v : mode.shape) v /= max_abs;
        modes.push_back(std::move(mode));
    }
    return modes;
}

std::vector<ModeShape> modal_analysis(const NonlinearDynamicModel& model,
                                      int requested_modes,
                                      double beta_tolerance,
                                      double imaginary_tolerance) {
    return modal_analysis_from_stiffness(model,model.K_initial(),requested_modes,
                                         beta_tolerance,imaginary_tolerance);
}

} // namespace quake
