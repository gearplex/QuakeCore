#pragma once

#include "quake/dynamic_model.hpp"
#include "quake/modal.hpp"

#include <vector>

namespace quake {

// Fixed-elastic-mode viscous damping wrapper.  Modal damping forces are
// evaluated from the modes of the wrapped model's initial elastic stiffness:
//   f_modal = sum_i 2*zeta_i*omega_i*(phi_i^T M v) M phi_i
// with mass-normalized phi_i.  The modal force is deliberately NOT inserted
// into the effective Newton tangent, matching the force-side implementation
// used by PERFORM/OpenSees and avoiding a dense damping matrix.  Any Rayleigh
// damping already present in the wrapped model remains active and does enter
// its effective tangent normally.
class FixedModalDampingModel final : public NonlinearDynamicModel {
public:
    FixedModalDampingModel(const NonlinearDynamicModel& base,
                           double modal_damping_ratio,
                           int requested_modes = 1000,
                           bool exact_low_rank_newton_tangent = false);
    FixedModalDampingModel(const NonlinearDynamicModel& base,
                           std::vector<double> modal_damping_ratios,
                           int requested_modes = 1000,
                           bool exact_low_rank_newton_tangent = false);

    int dof() const override { return base_->dof(); }
    int nonlinear_count() const override { return base_->nonlinear_count(); }
    const SparseMatrixCSC& K_initial() const override { return base_->K_initial(); }
    const std::vector<double>& mass() const override { return base_->mass(); }
    std::vector<double> mass_multiply(const std::vector<double>& a) const override { return base_->mass_multiply(a); }
    const SparseUpdateBasis& nonlinear_basis() const override { return base_->nonlinear_basis(); }
    int nonlinear_state_size() const override { return base_->nonlinear_state_size(); }
    const std::vector<double>& initial_nonlinear_tangents() const override { return base_->initial_nonlinear_tangents(); }
    std::vector<double> initial_nonlinear_state() const override { return base_->initial_nonlinear_state(); }

    std::vector<double> damping_multiply(const std::vector<double>& v) const override;
    bool has_additional_effective_low_rank_update() const override { return exact_low_rank_newton_tangent_ && !modes_.empty(); }
    const SparseUpdateBasis& additional_effective_low_rank_basis() const override { return modal_tangent_basis_; }
    std::vector<double> additional_effective_low_rank_coefficients(double a0,double a1) const override;
    SparseMatrixCSC effective_initial_matrix(double a0,double a1) const override { return base_->effective_initial_matrix(a0,a1); }
    SparseMatrixCSC effective_tangent_matrix(const std::vector<double>& t,double a0,double a1) const override { return base_->effective_tangent_matrix(t,a0,a1); }
    bool has_state_dependent_global_tangent() const override { return base_->has_state_dependent_global_tangent(); }
    SparseMatrixCSC effective_state_tangent_matrix(const std::vector<double>& u,const std::vector<double>& t,double a0,double a1) const override {
        return base_->effective_state_tangent_matrix(u,t,a0,a1);
    }
    SparseMatrixCSC effective_state_tangent_matrix_with_state(const std::vector<double>& u,const std::vector<double>& t,
                                                              const std::vector<double>& s,double a0,double a1) const override {
        return base_->effective_state_tangent_matrix_with_state(u,t,s,a0,a1);
    }
    bool has_generalized_state_update() const override { return base_->has_generalized_state_update(); }
    const SparseUpdateBasis& generalized_state_update_basis() const override { return base_->generalized_state_update_basis(); }
    std::vector<double> generalized_state_update_coefficients(const std::vector<double>& u,const std::vector<double>& s) const override {
        return base_->generalized_state_update_coefficients(u,s);
    }
    void evaluate_nonlinear_deformations(const std::vector<double>& q,const std::vector<double>& s,
                                         std::vector<double>& f,std::vector<double>& t,std::vector<double>& ts) const override {
        base_->evaluate_nonlinear_deformations(q,s,f,t,ts);
    }
    void internal_force_and_tangent(const std::vector<double>& u,const std::vector<double>& s,
                                    std::vector<double>& f,std::vector<double>& t,std::vector<double>& ts) const override {
        base_->internal_force_and_tangent(u,s,f,t,ts);
    }
    void internal_force_and_tangent_diagnostics(const std::vector<double>& u,const std::vector<double>& s,
                                    std::vector<double>& f,std::vector<double>& t,std::vector<double>& ts,
                                    NonlinearEvalDiagnostics* d) const override {
        base_->internal_force_and_tangent_diagnostics(u,s,f,t,ts,d);
    }
    std::vector<double> base_excitation(double ag) const override { return base_->base_excitation(ag); }
    double response_value(const std::vector<double>& u) const override { return base_->response_value(u); }
    double max_drift_measure(const std::vector<double>& u) const override { return base_->max_drift_measure(u); }
    double max_drift_ratio(const std::vector<double>& u) const override { return base_->max_drift_ratio(u); }

    int modal_count() const { return static_cast<int>(modes_.size()); }
    const std::vector<ModeShape>& modes() const { return raw_modes_; }
    const std::vector<double>& damping_ratios() const { return zeta_; }
    bool exact_low_rank_newton_tangent() const { return exact_low_rank_newton_tangent_; }

private:
    struct DampedMode {
        double omega{};
        std::vector<double> phi;   // mass-normalized
        std::vector<double> mphi;  // M*phi
    };
    const NonlinearDynamicModel* base_{};
    std::vector<ModeShape> raw_modes_;
    std::vector<DampedMode> modes_;
    std::vector<double> zeta_;
    SparseUpdateBasis modal_tangent_basis_; // columns are mass-normalized M*phi
    bool exact_low_rank_newton_tangent_{};
};

} // namespace quake
