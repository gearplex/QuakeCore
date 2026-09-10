#pragma once
#include "quake/dynamic_model.hpp"
#include "quake/bilinear.hpp"
#include <vector>

namespace quake {

class ShearBuilding final : public NonlinearDynamicModel {
public:
    ShearBuilding(int stories, double mass_per_story, double linear_story_stiffness,
                  std::vector<int> nonlinear_story_indices,
                  double nonlinear_initial_stiffness, double yield_force,
                  double post_yield_ratio,
                  double rayleigh_alpha_m = 0.0,
                  double rayleigh_beta_k = 0.002);

    int dof() const override { return n_; }
    int nonlinear_count() const override { return static_cast<int>(springs_.size()); }
    const SparseMatrixCSC& K_linear() const { return k_linear_; }
    const SparseMatrixCSC& K_initial() const override { return k_initial_; }
    const std::vector<double>& mass() const override { return mass_; }
    std::vector<double> mass_multiply(const std::vector<double>& a) const override;
    const SparseUpdateBasis& nonlinear_basis() const override { return basis_; }
    const std::vector<BilinearSpring>& springs() const { return springs_; }
    int nonlinear_state_size() const override { return 2 * nonlinear_count(); }
    const std::vector<double>& initial_nonlinear_tangents() const override { return initial_tangents_; }
    const std::vector<int>& nonlinear_stories() const { return nonlinear_stories_; }

    std::vector<double> damping_multiply(const std::vector<double>& v) const override;
    SparseMatrixCSC effective_initial_matrix(double a0, double a1) const override;
    SparseMatrixCSC effective_tangent_matrix(const std::vector<double>& tangents,
                                             double a0, double a1) const override;
    void evaluate_nonlinear_deformations(const std::vector<double>& deformations,
                                          const std::vector<double>& committed_state,
                                          std::vector<double>& component_forces,
                                          std::vector<double>& tangents,
                                          std::vector<double>& trial_state) const override;
    void internal_force_and_tangent(const std::vector<double>& u,
                                    const std::vector<double>& committed_state,
                                    std::vector<double>& force,
                                    std::vector<double>& tangents,
                                    std::vector<double>& trial_state) const override;
    std::vector<double> base_excitation(double ground_accel) const override;
    double response_value(const std::vector<double>& u) const override;
    double max_drift_measure(const std::vector<double>& u) const override;

private:
    int n_{};
    std::vector<double> mass_;
    SparseMatrixCSC k_linear_;
    SparseMatrixCSC k_initial_;
    std::vector<int> nonlinear_stories_;
    std::vector<BilinearSpring> springs_;
    std::vector<double> initial_tangents_;
    SparseUpdateBasis basis_;
    double alpha_m_{};
    double beta_k_{};
};

} // namespace quake
