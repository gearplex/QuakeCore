#include "quake/shear_building.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quake {

static SparseMatrixCSC story_stiffness_matrix(int n, double k,
                                               const std::vector<int>& stories = {},
                                               bool include_all_if_empty = true) {
    std::vector<bool> include(static_cast<std::size_t>(n), stories.empty() && include_all_if_empty);
    for (int s : stories) {
        if (s < 0 || s >= n) throw std::out_of_range("story index");
        include[static_cast<std::size_t>(s)] = true;
    }
    std::vector<Triplet> t;
    for (int s = 0; s < n; ++s) {
        if (!include[static_cast<std::size_t>(s)]) continue;
        const int top = s;
        t.push_back({top, top, k});
        if (s > 0) {
            const int bot = s - 1;
            t.push_back({bot, bot, k});
            t.push_back({top, bot, -k});
            t.push_back({bot, top, -k});
        }
    }
    return SparseMatrixCSC::from_triplets(n, n, t, 1e-18);
}

ShearBuilding::ShearBuilding(int stories, double mass_per_story, double linear_story_stiffness,
                             std::vector<int> nonlinear_story_indices,
                             double nonlinear_initial_stiffness, double yield_force,
                             double post_yield_ratio, double rayleigh_alpha_m,
                             double rayleigh_beta_k)
    : n_(stories), mass_(static_cast<std::size_t>(stories), mass_per_story),
      k_linear_(story_stiffness_matrix(stories, linear_story_stiffness)),
      nonlinear_stories_(std::move(nonlinear_story_indices)),
      alpha_m_(rayleigh_alpha_m), beta_k_(rayleigh_beta_k) {
    if (stories <= 0 || mass_per_story <= 0.0 || linear_story_stiffness <= 0.0)
        throw std::invalid_argument("Invalid building parameters");
    std::sort(nonlinear_stories_.begin(), nonlinear_stories_.end());
    nonlinear_stories_.erase(std::unique(nonlinear_stories_.begin(), nonlinear_stories_.end()),
                             nonlinear_stories_.end());
    const int m = static_cast<int>(nonlinear_stories_.size());
    std::vector<std::vector<std::pair<int,double>>> basis_columns(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) {
        const int s = nonlinear_stories_[static_cast<std::size_t>(j)];
        if (s < 0 || s >= n_) throw std::out_of_range("nonlinear story index");
        basis_columns[static_cast<std::size_t>(j)].push_back({s,1.0});
        if (s > 0) basis_columns[static_cast<std::size_t>(j)].push_back({s-1,-1.0});
        springs_.emplace_back(nonlinear_initial_stiffness, yield_force, post_yield_ratio);
        initial_tangents_.push_back(nonlinear_initial_stiffness);
    }
    basis_=SparseUpdateBasis::from_columns(n_,basis_columns);
    auto knl = story_stiffness_matrix(n_, nonlinear_initial_stiffness, nonlinear_stories_, false);
    k_initial_ = add(k_linear_, knl);
}

std::vector<double> ShearBuilding::mass_multiply(const std::vector<double>& a) const {
    if (static_cast<int>(a.size()) != n_) throw std::invalid_argument("mass vector size");
    std::vector<double> out(static_cast<std::size_t>(n_));
    for (int i=0;i<n_;++i) out[static_cast<std::size_t>(i)] = mass_[static_cast<std::size_t>(i)] * a[static_cast<std::size_t>(i)];
    return out;
}

std::vector<double> ShearBuilding::damping_multiply(const std::vector<double>& v) const {
    if (static_cast<int>(v.size()) != n_) throw std::invalid_argument("damping vector size");
    auto kv = k_initial_.multiply(v);
    std::vector<double> out(static_cast<std::size_t>(n_));
    for (int i = 0; i < n_; ++i)
        out[static_cast<std::size_t>(i)] = alpha_m_ * mass_[static_cast<std::size_t>(i)] * v[static_cast<std::size_t>(i)] +
                                           beta_k_ * kv[static_cast<std::size_t>(i)];
    return out;
}

SparseMatrixCSC ShearBuilding::effective_initial_matrix(double a0, double a1) const {
    auto k = add(k_initial_, k_initial_, 1.0, a1 * beta_k_); // (1+a1 betaK) Kinit
    std::vector<double> d(static_cast<std::size_t>(n_));
    for (int i = 0; i < n_; ++i)
        d[static_cast<std::size_t>(i)] = (a0 + a1 * alpha_m_) * mass_[static_cast<std::size_t>(i)];
    return add_diagonal(k, d);
}

SparseMatrixCSC ShearBuilding::effective_tangent_matrix(const std::vector<double>& tangents,
                                                        double a0, double a1) const {
    if (static_cast<int>(tangents.size()) != nonlinear_count())
        throw std::invalid_argument("tangent vector size");
    std::vector<Triplet> t;
    // Build Klinear + current nonlinear tangents.
    for (int c = 0; c < k_linear_.cols(); ++c)
        for (int p = k_linear_.col_ptr()[static_cast<std::size_t>(c)];
             p < k_linear_.col_ptr()[static_cast<std::size_t>(c+1)]; ++p)
            t.push_back({k_linear_.row_ind()[static_cast<std::size_t>(p)], c,
                         k_linear_.values()[static_cast<std::size_t>(p)]});
    for (int j = 0; j < nonlinear_count(); ++j) {
        const double k = tangents[static_cast<std::size_t>(j)];
        const int s = nonlinear_stories_[static_cast<std::size_t>(j)];
        t.push_back({s,s,k});
        if (s > 0) {
            t.push_back({s-1,s-1,k}); t.push_back({s,s-1,-k}); t.push_back({s-1,s,-k});
        }
    }
    auto kt = SparseMatrixCSC::from_triplets(n_, n_, t, 1e-18);
    // C = alphaM M + betaK Kinitial, so Keff = Kt + a0 M + a1 C.
    auto out = add(kt, k_initial_, 1.0, a1 * beta_k_);
    std::vector<double> d(static_cast<std::size_t>(n_));
    for (int i = 0; i < n_; ++i)
        d[static_cast<std::size_t>(i)] = (a0 + a1*alpha_m_) * mass_[static_cast<std::size_t>(i)];
    return add_diagonal(out, d);
}

void ShearBuilding::evaluate_nonlinear_deformations(const std::vector<double>& q,
                                                      const std::vector<double>& committed_state,
                                                      std::vector<double>& component_forces,
                                                      std::vector<double>& tangents,
                                                      std::vector<double>& trial_state) const {
    if (static_cast<int>(q.size()) != nonlinear_count() || static_cast<int>(committed_state.size()) != nonlinear_state_size())
        throw std::invalid_argument("nonlinear deformation/state dimension mismatch");
    component_forces.resize(static_cast<std::size_t>(nonlinear_count()));
    tangents.resize(static_cast<std::size_t>(nonlinear_count()));
    trial_state.resize(static_cast<std::size_t>(nonlinear_state_size()));
    for (int j=0;j<nonlinear_count();++j) {
        const BilinearState c{committed_state[static_cast<std::size_t>(2*j)], committed_state[static_cast<std::size_t>(2*j+1)]};
        const auto tr=springs_[static_cast<std::size_t>(j)].trial(q[static_cast<std::size_t>(j)],c);
        component_forces[static_cast<std::size_t>(j)]=tr.force;
        tangents[static_cast<std::size_t>(j)]=tr.tangent;
        trial_state[static_cast<std::size_t>(2*j)]=tr.state.plastic;
        trial_state[static_cast<std::size_t>(2*j+1)]=tr.state.backstress;
    }
}

void ShearBuilding::internal_force_and_tangent(const std::vector<double>& u,
                                               const std::vector<double>& committed_state,
                                               std::vector<double>& force,
                                               std::vector<double>& tangents,
                                               std::vector<double>& trial_state) const {
    if (static_cast<int>(u.size()) != n_ || static_cast<int>(committed_state.size()) != nonlinear_state_size())
        throw std::invalid_argument("state dimension mismatch");
    force = k_linear_.multiply(u);
    tangents.resize(static_cast<std::size_t>(nonlinear_count()));
    trial_state.resize(static_cast<std::size_t>(nonlinear_state_size()));
    for (int j = 0; j < nonlinear_count(); ++j) {
        const double q=basis_.column_dot(j,u);
        const BilinearState c{committed_state[static_cast<std::size_t>(2*j)], committed_state[static_cast<std::size_t>(2*j+1)]};
        const auto tr = springs_[static_cast<std::size_t>(j)].trial(q, c);
        tangents[static_cast<std::size_t>(j)] = tr.tangent;
        trial_state[static_cast<std::size_t>(2*j)] = tr.state.plastic;
        trial_state[static_cast<std::size_t>(2*j+1)] = tr.state.backstress;
        basis_.axpy_column(j,tr.force,force);
    }
}

std::vector<double> ShearBuilding::base_excitation(double ground_accel) const {
    std::vector<double> p(static_cast<std::size_t>(n_));
    for (int i = 0; i < n_; ++i) p[static_cast<std::size_t>(i)] = -mass_[static_cast<std::size_t>(i)] * ground_accel;
    return p;
}

double ShearBuilding::response_value(const std::vector<double>& u) const {
    if (static_cast<int>(u.size()) != n_) throw std::invalid_argument("response vector size");
    return u.back();
}

double ShearBuilding::max_drift_measure(const std::vector<double>& u) const {
    if (static_cast<int>(u.size()) != n_) throw std::invalid_argument("response vector size");
    double out = 0.0;
    double lower = 0.0;
    for (double ui : u) {
        out = std::max(out, std::abs(ui - lower));
        lower = ui;
    }
    return out;
}

} // namespace quake
