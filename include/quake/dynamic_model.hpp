#pragma once
#include "quake/sparse.hpp"
#include "quake/update_basis.hpp"
#include <vector>
#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace quake {

// A recoverable local constitutive integration failure. The transient driver
// may backtrack/subdivide these failures; invalid model definitions still throw
// their original exceptions and are never hidden as convergence failures.
class ConstitutiveIntegrationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct NonlinearEvalDiagnostics {
    std::size_t component_evaluations{0};
    std::size_t active_tangent_evaluations{0};
    std::size_t fast_path_evaluations{0};
    std::size_t full_state_evaluations{0};
    std::size_t transition_events{0};
    std::size_t reversal_events{0};
    std::size_t deterioration_events{0};
    std::size_t failure_events{0};
    std::size_t io_or_beyond_evaluations{0};
    std::size_t ls_or_beyond_evaluations{0};
    std::size_t cp_or_beyond_evaluations{0};
    std::size_t beyond_cp_evaluations{0};
    std::size_t lateral_loss_evaluations{0};
};

// Runtime interface consumed by the nonlinear transient solver. Model compilation
// happens outside the time loop; all topology is frozen behind this interface.
class NonlinearDynamicModel {
public:
    virtual ~NonlinearDynamicModel() = default;

    virtual int dof() const = 0;
    virtual int nonlinear_count() const = 0;
    virtual const SparseMatrixCSC& K_initial() const = 0;
    virtual const std::vector<double>& mass() const = 0;
    // Generalized mass action. For unconstrained/lumped models this is diagonal;
    // MPC condensation may create off-diagonal generalized mass terms.
    virtual std::vector<double> mass_multiply(const std::vector<double>& a) const = 0;
    virtual const SparseUpdateBasis& nonlinear_basis() const = 0;
    // Material-neutral nonlinear interface. Solver state is a flat model-owned
    // vector so the transient algorithm does not depend on any constitutive law.
    virtual int nonlinear_state_size() const = 0;
    virtual const std::vector<double>& initial_nonlinear_tangents() const = 0;
    virtual std::vector<double> initial_nonlinear_state() const {
        return std::vector<double>(static_cast<std::size_t>(nonlinear_state_size()), 0.0);
    }

    virtual std::vector<double> damping_multiply(const std::vector<double>& v) const = 0;
    // Optional constant/known low-rank contribution to the effective Newton
    // tangent that is intentionally omitted from the sparse effective matrix.
    // This is useful for fixed-mode modal damping, whose exact velocity
    // Jacobian is dense but very low rank:
    //   K_eff,extra = U C(a0,a1) U^T.
    // The nonlinear solver may apply this contribution exactly through a
    // Woodbury correction without densifying the global sparse matrix.
    virtual bool has_additional_effective_low_rank_update() const { return false; }
    virtual const SparseUpdateBasis& additional_effective_low_rank_basis() const { return nonlinear_basis(); }
    virtual std::vector<double> additional_effective_low_rank_coefficients(double, double) const { return {}; }
    virtual SparseMatrixCSC effective_initial_matrix(double a0, double a1) const = 0;
    virtual SparseMatrixCSC effective_tangent_matrix(const std::vector<double>& tangents,
                                                      double a0, double a1) const = 0;
    // Optional state-dependent global tangent (e.g., updated P-Delta). Such a
    // model may require the direct same-pattern path rather than pure Woodbury.
    virtual bool has_state_dependent_global_tangent() const { return false; }
    virtual SparseMatrixCSC effective_state_tangent_matrix(const std::vector<double>& u,
                                                            const std::vector<double>& tangents,
                                                            double a0,double a1) const {
        (void)u; return effective_tangent_matrix(tangents,a0,a1);
    }
    // State-aware overload for coupled constitutive updates whose consistent
    // tangent depends on the committed hysteretic branch in addition to the
    // current displacement.  Legacy models inherit the state-independent path.
    virtual SparseMatrixCSC effective_state_tangent_matrix_with_state(
                                                            const std::vector<double>& u,
                                                            const std::vector<double>& tangents,
                                                            const std::vector<double>& committed_state,
                                                            double a0,double a1) const {
        (void)committed_state;
        return effective_state_tangent_matrix(u,tangents,a0,a1);
    }
    // Velocity-dependent devices require their consistent damping Jacobian in
    // the Newton matrix.  Under Newmark average acceleration, dv/du=a1, so a
    // device tangent dFd/dv contributes a1*dFd/dv.  The default preserves the
    // displacement-only behavior of existing models.
    virtual bool has_velocity_dependent_global_tangent() const { return false; }
    virtual SparseMatrixCSC effective_state_tangent_matrix_with_state_and_velocity(
                                                            const std::vector<double>& u,
                                                            const std::vector<double>& tangents,
                                                            const std::vector<double>& committed_state,
                                                            const std::vector<double>& velocity,
                                                            double a0,double a1) const {
        (void)velocity;
        return effective_state_tangent_matrix_with_state(u,tangents,committed_state,a0,a1);
    }
    // Optional exact generalized low-rank representation of the state-
    // dependent tangent relative to effective_initial_matrix:
    //   K_eff(u) = A_eff + U C(u) U^T.
    // This extends the diagonal hinge Woodbury update to coupled localized
    // blocks such as changing geometric stiffness.
    virtual bool has_generalized_state_update() const { return false; }
    virtual const SparseUpdateBasis& generalized_state_update_basis() const { return nonlinear_basis(); }
    virtual std::vector<double> generalized_state_update_coefficients(
        const std::vector<double>&, const std::vector<double>&) const { return {}; }
    // Constitutive bank interface. Deformations are generalized component
    // coordinates q_j = b_j^T u. Keeping this separate from global assembly
    // allows reduced models and future GPU material banks to avoid expanding
    // back to the parent/global coordinate space.
    virtual void evaluate_nonlinear_deformations(const std::vector<double>& deformations,
                                                  const std::vector<double>& committed_state,
                                                  std::vector<double>& component_forces,
                                                  std::vector<double>& tangents,
                                                  std::vector<double>& trial_state) const = 0;
    virtual void internal_force_and_tangent(const std::vector<double>& u,
                                             const std::vector<double>& committed_state,
                                             std::vector<double>& force,
                                             std::vector<double>& tangents,
                                             std::vector<double>& trial_state) const = 0;
    virtual void internal_force_and_tangent_diagnostics(const std::vector<double>& u,
                                             const std::vector<double>& committed_state,
                                             std::vector<double>& force,
                                             std::vector<double>& tangents,
                                             std::vector<double>& trial_state,
                                             NonlinearEvalDiagnostics* diagnostics) const {
        internal_force_and_tangent(u, committed_state, force, tangents, trial_state);
        if (!diagnostics) return;
        diagnostics->component_evaluations += static_cast<std::size_t>(nonlinear_count());
        const auto& init = initial_nonlinear_tangents();
        for (int j=0;j<nonlinear_count();++j) {
            const bool active = std::abs(tangents[static_cast<std::size_t>(j)] - init[static_cast<std::size_t>(j)]) > 1e-14;
            diagnostics->active_tangent_evaluations += static_cast<std::size_t>(active);
            diagnostics->fast_path_evaluations += static_cast<std::size_t>(!active);
            diagnostics->full_state_evaluations += static_cast<std::size_t>(active);
        }
    }
    virtual std::vector<double> base_excitation(double ground_accel) const = 0;

    // Lightweight response reducers used by the benchmark/transient driver.
    virtual double response_value(const std::vector<double>& u) const = 0;
    virtual double max_drift_measure(const std::vector<double>& u) const = 0;
    // Story drift ratio when geometry is available. NaN means unsupported.
    virtual double max_drift_ratio(const std::vector<double>&) const { return std::numeric_limits<double>::quiet_NaN(); }
};

} // namespace quake
