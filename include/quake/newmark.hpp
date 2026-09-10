#pragma once
#include "quake/dynamic_model.hpp"
#include "quake/low_rank_solver.hpp"
#include "quake/superlu_solver.hpp"
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>
#include <string>

namespace quake {

enum class LinearStrategy { FullFactorization, SamePatternRefactorization, ModifiedNewton, Woodbury, Adaptive };

struct AnalysisStats {
    std::size_t steps{0};
    std::size_t newton_iterations{0};
    std::size_t global_factorizations{0};
    std::size_t linear_solves{0};
    std::size_t failed_steps{0};
    std::size_t reduced_update_solves{0};
    std::size_t active_rank_sum{0};
    std::size_t max_active_rank{0};
    std::size_t cached_low_rank_columns{0};
    std::size_t residual_evaluations{0};
    std::size_t constitutive_integration_failures{0};
    std::size_t nonfinite_evaluations{0};
    double last_residual_norm{0.0};
    double last_residual_tolerance{0.0};
    std::size_t line_search_backtracks{0};
    std::size_t subdivided_steps{0};
    std::size_t internal_substeps{0};
    std::size_t direct_fallbacks{0};
    std::size_t linear_solve_failures{0};
    std::size_t line_search_failures{0};
    std::size_t max_iteration_failures{0};
    std::size_t adaptive_woodbury_solves{0};
    std::size_t adaptive_direct_solves{0};
    std::size_t generalized_update_solves{0};
    std::size_t max_generalized_update_rank{0};
    // Diagnostics for event-driven/material-frontier research. These count
    // constitutive evaluations requested by the transient driver and the
    // subset whose tangent differs from the initial elastic tangent.
    std::size_t nonlinear_component_evaluations{0};
    std::size_t nonlinear_active_component_evaluations{0};
    std::size_t nonlinear_fast_path_evaluations{0};
    std::size_t nonlinear_full_state_evaluations{0};
    std::size_t nonlinear_transition_events{0};
    std::size_t nonlinear_reversal_events{0};
    std::size_t nonlinear_deterioration_events{0};
    std::size_t nonlinear_failure_events{0};
    std::size_t nonlinear_io_or_beyond{0};
    std::size_t nonlinear_ls_or_beyond{0};
    std::size_t nonlinear_cp_or_beyond{0};
    std::size_t nonlinear_beyond_cp{0};
    std::size_t nonlinear_lateral_loss{0};
    double minimum_dt{0.0};
    double elapsed_seconds{0.0};
    double max_roof_abs{0.0};
    double max_story_drift_abs{0.0};
    double minimum_tangent_ratio{1.0};
    double minimum_tangent_eigenvalue{0.0};
    double minimum_reduced_pivot_ratio{1.0};
    double initial_spd_minimum_pivot_ratio{0.0};
    double initial_symmetry_relative_error{0.0};
    // Robust-driver energy ledger. Internal work includes recoverable plus
    // hysteretic work; balance = input - internal - damping - kinetic.
    double input_energy{0.0};
    double internal_work{0.0};
    double damping_energy{0.0};
    double kinetic_energy{0.0};
    double energy_balance_error{0.0};
    double energy_balance_relative_error{0.0};
};

enum class AnalysisTermination { Completed, PhysicalCollapse, NumericalFailure, InitialInstability };
enum class CollapseMechanism {
    None,
    InitialInstability,
    TangentInstability,
    DriftLimit,
    LateralResistanceLoss,
    GravityResistanceLoss,
    BeyondCPConfigured,
    UnrecoverableLateralLoss,
    NumericalFailure
};

struct AnalysisResult {
    AnalysisStats stats;
    std::vector<double> roof_history;
    std::vector<double> final_displacement;
    std::vector<double> final_state;
    AnalysisTermination termination{AnalysisTermination::Completed};
    CollapseMechanism collapse_mechanism{CollapseMechanism::None};
    std::size_t termination_step{0};
    double termination_time{0.0};
    double max_story_drift_ratio{0.0};
    std::size_t first_io_step{static_cast<std::size_t>(-1)};
    std::size_t first_ls_step{static_cast<std::size_t>(-1)};
    std::size_t first_cp_step{static_cast<std::size_t>(-1)};
    std::size_t first_lateral_loss_step{static_cast<std::size_t>(-1)}; // ASCE 41 E
    std::size_t first_component_failure_step{static_cast<std::size_t>(-1)}; // ASCE 41 F / configured failure
    std::size_t max_components_io_or_beyond{0};
    std::size_t max_components_ls_or_beyond{0};
    std::size_t max_components_cp_or_beyond{0};
    std::size_t max_components_beyond_cp{0};
    std::size_t max_failed_components{0};
    std::size_t max_components_lateral_loss{0};
    std::string termination_reason;
    std::string last_integration_error;
};

AnalysisResult run_newmark(const NonlinearDynamicModel& model,
                           const std::vector<double>& ground_accel,
                           double dt, LinearStrategy strategy,
                           double tolerance = 1e-8,
                           int max_iterations = 20);

// Reuses the model/dt-dependent A factorization and A^-1 B influence matrix
// across many independent records. The model must outlive this object.
class PreparedWoodburyNewmark {
public:
    PreparedWoodburyNewmark(const NonlinearDynamicModel& model, double dt);

    AnalysisResult run(const std::vector<double>& ground_accel,
                       double tolerance = 1e-8,
                       int max_iterations = 20) const;

    double setup_seconds() const { return setup_seconds_; }
    std::size_t setup_factorizations() const { return 1; }
    double dt() const { return dt_; }

private:
    const NonlinearDynamicModel* model_{};
    double dt_{};
    double setup_seconds_{};
    std::unique_ptr<LowRankSolverBase> solver_;
};


// Prepared hybrid solver using a model-calibrated active-rank limit. Both
// exact paths are available: Woodbury below the limit and same-pattern numeric
// refactorization above it. Intended for suite reuse after calibration.
class PreparedAdaptiveNewmark {
public:
    PreparedAdaptiveNewmark(const NonlinearDynamicModel& model, double dt, int woodbury_rank_limit);
    AnalysisResult run(const std::vector<double>& ground_accel,
                       double tolerance=1e-8, int max_iterations=20) const;
    int woodbury_rank_limit() const { return rank_limit_; }
    double setup_seconds() const { return setup_seconds_; }
private:
    const NonlinearDynamicModel* model_{};
    double dt_{}; int rank_limit_{}; double setup_seconds_{};
    std::unique_ptr<LowRankSolverBase> wood_;
    std::unique_ptr<SuperLUSamePatternSolver> direct_;
};

struct CollapseCriteria {
    // Certify that the starting structural tangent is positive definite before
    // integrating. This is independent of the optional near-zero tangent
    // tracking used later to detect mechanism formation.
    bool check_initial_stability{true};
    // 0 disables a criterion. CP is tracked but is not collapse by default.
    double max_story_drift_ratio{0.0};
    std::size_t max_failed_components{0};
    std::size_t max_beyond_cp_components{0};
    std::size_t max_lateral_loss_components{0};
    // Optional near-singularity check. Every N accepted steps, estimate the
    // tangent eigenvalue nearest zero and compare to the initial tangent.
    int tangent_check_interval{0};
    double min_tangent_ratio{0.0};
    // If repeated time-step subdivision cannot find equilibrium and the failed
    // attempts enter ASCE41 E-point lateral-resistance loss, classify the
    // event as physical loss-of-equilibrium rather than a generic numerical failure.
    bool classify_unrecoverable_lateral_loss_as_collapse{false};
};

struct RobustNewmarkOptions {
    double tolerance{1e-8};
    // false selects a fixed absolute infinity-norm residual tolerance, useful
    // for exact convergence-criterion parity with an external solver.
    bool relative_force_tolerance{true};
    // OpenSees' displacement-form Newmark starts Newton at the previous u.
    // The default kinematic predictor is retained for existing callers.
    bool kinematic_initial_guess{true};
    int max_iterations{25};
    bool line_search{true};
    int max_backtracks{7};
    double backtrack_ratio{0.5};
    // Optional nonmonotone merit allowance for nonsmooth degradation events.
    // 1.0 requires strict decrease; values >1 permit the best finite trial
    // to grow modestly before subsequent Newton iterations recover.
    double nonmonotone_line_search_factor{1.0};
    int max_subdivisions{4};
    // For Woodbury, switch to a direct tangent solve if instantaneous reduced
    // rank exceeds this fraction of global DOF. <=0 disables proactive fallback.
    double direct_rank_fraction{0.35};
    CollapseCriteria collapse{};
    bool return_numerical_failure{false};
    // Optional callback at each accepted original output step (after any
    // internal subdivision). Intended for validation/EDP recorders without
    // forcing the production solver to retain full state histories.
    std::function<void(std::size_t step, double time, double ground_accel,
                       const std::vector<double>& displacement,
                       const std::vector<double>& velocity,
                       const std::vector<double>& acceleration)> accepted_step_observer{};
    // State-aware companion used by validation/component instrumentation that
    // needs the committed constitutive history at an accepted output step.
    // Kept separate so existing observers and the hot path remain unchanged.
    std::function<void(std::size_t step, double time, double ground_accel,
                       const std::vector<double>& displacement,
                       const std::vector<double>& velocity,
                       const std::vector<double>& acceleration,
                       const std::vector<double>& committed_state)> accepted_state_observer{};
    // Called for every integration substep that belongs to an accepted original
    // output step. Substep snapshots are buffered until the full original step
    // succeeds, so rejected Newton trials and rolled-back subdivision branches
    // are never exposed. This is the correct hook for exact committed-state
    // demand envelopes used by code-parameter regeneration.
    std::function<void(std::size_t output_step, std::size_t substep_index,
                       int subdivision_depth, double time, double ground_accel,
                       const std::vector<double>& displacement,
                       const std::vector<double>& velocity,
                       const std::vector<double>& acceleration,
                       const std::vector<double>& committed_state)> accepted_substep_state_observer{};
};

AnalysisResult run_newmark_robust(const NonlinearDynamicModel& model,
                                  const std::vector<double>& ground_accel,
                                  double dt, LinearStrategy strategy,
                                  const RobustNewmarkOptions& options = {});

// Prepared robust driver for repeated analyses of one compiled topology. The
// model may replace nonlinear material parameters between runs. The prepared
// linear contexts are reused only while every subdivision-depth effective
// initial matrix remains exactly unchanged; otherwise they are invalidated and
// rebuilt automatically before the next run. This makes comparison/calibration
// loops fast without allowing stale factorizations after a Ke/material-field
// change.
class PreparedRobustNewmark {
public:
    PreparedRobustNewmark(const NonlinearDynamicModel& model, double dt,
                          LinearStrategy strategy, int prepared_subdivision_depth = 4);
    ~PreparedRobustNewmark();
    PreparedRobustNewmark(PreparedRobustNewmark&&) noexcept;
    PreparedRobustNewmark& operator=(PreparedRobustNewmark&&) noexcept;
    PreparedRobustNewmark(const PreparedRobustNewmark&) = delete;
    PreparedRobustNewmark& operator=(const PreparedRobustNewmark&) = delete;

    AnalysisResult run(const std::vector<double>& ground_accel,
                       const RobustNewmarkOptions& options = {});

    double dt() const { return dt_; }
    LinearStrategy strategy() const { return strategy_; }
    int prepared_subdivision_depth() const { return prepared_subdivision_depth_; }
    std::size_t preparation_count() const;
    std::size_t setup_factorizations() const;
    bool last_run_reused_preparation() const;

private:
    struct Impl;
    const NonlinearDynamicModel* model_{};
    double dt_{};
    LinearStrategy strategy_{};
    int prepared_subdivision_depth_{};
    std::unique_ptr<Impl> impl_;
};

std::vector<double> synthetic_ground_motion(int steps, double dt,
                                            double amplitude = 6.0);

} // namespace quake
