#pragma once

#include "quake/dynamic_model.hpp"

#include <vector>

namespace quake {

struct ReductionBuildInfo {
    int full_dof{};
    int retained_dof{};
    int fixed_interface_modes{};
    int reduced_dof{};
};

// Generic linear-coordinate reduction wrapper. Nonlinear constitutive state is
// preserved exactly; only the kinematic space is reduced by u_full = T q.
class ReducedDynamicModel final : public NonlinearDynamicModel {
public:
    ReducedDynamicModel(const NonlinearDynamicModel& parent,
                        std::vector<double> transform_column_major,
                        int reduced_dof);

    int dof() const override { return nr_; }
    int nonlinear_count() const override { return parent_->nonlinear_count(); }
    const SparseMatrixCSC& K_initial() const override { return k_initial_; }
    const std::vector<double>& mass() const override { return mass_diag_; }
    std::vector<double> mass_multiply(const std::vector<double>& a) const override;
    const SparseUpdateBasis& nonlinear_basis() const override { return basis_; }
    int nonlinear_state_size() const override { return parent_->nonlinear_state_size(); }
    std::vector<double> initial_nonlinear_state() const override { return parent_->initial_nonlinear_state(); }
    const std::vector<double>& initial_nonlinear_tangents() const override {
        return parent_->initial_nonlinear_tangents();
    }

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
    void internal_force_and_tangent_diagnostics(const std::vector<double>& u,
                                    const std::vector<double>& committed_state,
                                    std::vector<double>& force,
                                    std::vector<double>& tangents,
                                    std::vector<double>& trial_state,
                                    NonlinearEvalDiagnostics* diagnostics) const override;
    std::vector<double> base_excitation(double ground_accel) const override;
    double response_value(const std::vector<double>& u) const override;
    double max_drift_measure(const std::vector<double>& u) const override;

    std::vector<double> expand(const std::vector<double>& reduced) const;
    const std::vector<double>& transform() const { return transform_; }
    int full_dof() const { return nf_; }

private:
    std::vector<double> project(const std::vector<double>& full) const;
    SparseMatrixCSC project_matrix(const SparseMatrixCSC& full) const;

    const NonlinearDynamicModel* parent_{};
    int nf_{};
    int nr_{};
    std::vector<double> transform_; // column-major, nf x nr
    SparseMatrixCSC k_initial_;
    SparseMatrixCSC k_linear_;
    SparseMatrixCSC mass_matrix_;
    SparseMatrixCSC damping_matrix_;
    std::vector<double> mass_diag_;
    SparseUpdateBasis basis_;
};

// Returns every global DOF touched by at least one nonlinear generalized
// deformation direction. Useful as the interface set for substructuring.
std::vector<int> nonlinear_support_dofs(const NonlinearDynamicModel& model);

// Craig-Bampton/Guyan prototype. Retained DOFs are represented exactly. The
// interior is represented by static constraint modes plus the requested number
// of lowest finite fixed-interface modes. This is a validation/research utility;
// a sparse Lanczos implementation is planned for production-scale reduction.
ReducedDynamicModel craig_bampton_reduce(const NonlinearDynamicModel& model,
                                         std::vector<int> retained_dofs,
                                         int fixed_interface_modes = 0,
                                         ReductionBuildInfo* info = nullptr);

} // namespace quake
