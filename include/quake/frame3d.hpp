#pragma once
#include "quake/dynamic_model.hpp"
#include "quake/bilinear.hpp"
#include "quake/nonlinear_material.hpp"
#include "quake/corotational3d.hpp"
#include "quake/fsc_shear_spring.hpp"
#include "quake/pm_interaction_hinge.hpp"

#include <array>
#include <unordered_map>
#include <optional>
#include <utility>
#include <vector>

namespace quake {

enum class Dof3D : int { UX=0, UY=1, UZ=2, RX=3, RY=4, RZ=5 };

struct Node3D {
    int id{};
    double x{}, y{}, z{};
    std::array<double,6> mass{};
};

struct ElasticFrame3D {
    int id{};
    int node_i{}, node_j{};
    double E{}, G{}, A{}, J{}, Iy{}, Iz{};
    std::array<double,3> reference{}; // defines local +y after projection normal to local x
    double axial_compression{};       // positive compression, constant preload
};

struct MpcTerm3D {
    int node{};
    Dof3D dof{};
    double coefficient{};
};

struct GeneralizedSpring3D {
    int id{};
    // Generalized deformation q = sum(term.coefficient * u_term).
    std::vector<MpcTerm3D> terms;
    NonlinearMaterial material;
};

// Local member-force recovery at a committed analysis state. Compression is
// positive. Shear_y/z are local section shears at each end; max vector shear
// is therefore hypot(shear_y, shear_z).
struct ElasticFrame3DResponse {
    int element_id{};
    double axial_compression{};
    double shear_y_i{}, shear_y_j{};
    double shear_z_i{}, shear_z_j{};
    double torsion_i{}, torsion_j{};
    double moment_y_i{}, moment_y_j{};
    double moment_z_i{}, moment_z_j{};
};

// Intermediate Phase-8 P->M coupling: the existing ASCE41 rotational hinge
// remains the hysteretic law, while its yield/capacity magnitude is evaluated
// continuously from the current elastic member axial force.  Axial response is
// not duplicated here; axial_stiffness and axial_deformation_terms only recover
// the force already carried by the associated elastic frame member.
struct AxialCoupledASCE41Params {
    ASCE41HingeParams hinge;
    double axial_preload{};       // positive compression
    double axial_stiffness{};     // EA/L, force per generalized axial deformation
    std::vector<double> axial_force_points;
    std::vector<double> moment_capacity_points;
};

class CompiledFrame3D final : public NonlinearDynamicModel {
public:
    int dof() const override { return n_; }
    int nonlinear_count() const override { return static_cast<int>(materials_.size()); }
    const SparseMatrixCSC& K_linear() const { return k_linear_; }
    const SparseMatrixCSC& K_initial() const override { return k_initial_; }
    const std::vector<double>& mass() const override { return mass_; }
    std::vector<double> mass_multiply(const std::vector<double>& a) const override;
    const SparseMatrixCSC& mass_matrix() const { return mass_matrix_; }
    const SparseUpdateBasis& nonlinear_basis() const override { return basis_; }
    const std::vector<NonlinearMaterial>& materials() const { return materials_; }
    int nonlinear_state_size() const override { return nonlinear_state_size_; }
    std::vector<double> initial_nonlinear_state() const override;
    const std::vector<double>& initial_nonlinear_tangents() const override { return initial_tangents_; }

    std::vector<double> damping_multiply(const std::vector<double>& v) const override;
    SparseMatrixCSC effective_initial_matrix(double a0, double a1) const override;
    SparseMatrixCSC effective_tangent_matrix(const std::vector<double>& tangents,
                                             double a0, double a1) const override;
    bool has_state_dependent_global_tangent() const override { return updated_pdelta_ || corotational_ || !coupled_pm_hinges_.empty() || !pm_return_hinges_.empty() || !fsc_shear_springs_.empty(); }
    bool has_generalized_state_update() const override { return updated_pdelta_ && !corotational_ && coupled_pm_hinges_.empty() && pm_return_hinges_.empty() && generalized_state_update_basis_.cols()>0; }
    const SparseUpdateBasis& generalized_state_update_basis() const override { return generalized_state_update_basis_; }
    std::vector<double> generalized_state_update_coefficients(const std::vector<double>& u,
                                                               const std::vector<double>& tangents) const override;
    SparseMatrixCSC effective_state_tangent_matrix(const std::vector<double>& u,
                                                    const std::vector<double>& tangents,
                                                    double a0,double a1) const override;
    SparseMatrixCSC effective_state_tangent_matrix_with_state(const std::vector<double>& u,
                                                    const std::vector<double>& tangents,
                                                    const std::vector<double>& committed_state,
                                                    double a0,double a1) const override;
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
    double max_drift_ratio(const std::vector<double>& u) const override;
    // Evaluate the compiled story response coordinates for any reduced vector.
    // Useful for validation recorders (displacement, velocity, acceleration).
    std::vector<double> story_response_values(const std::vector<double>& x) const;
    const std::vector<double>& story_elevations() const { return story_z_; }

    int reduced_dof(int node_id, Dof3D dof) const;
    int nonlinear_component_index(int spring_id) const;
    int nonlinear_component_id(int index) const;
    NonlinearComponentSnapshot nonlinear_component_snapshot(int spring_id,
                                                             const std::vector<double>& u,
                                                             const std::vector<double>& committed_state) const;
    PMInteractionTrialResult pm_interaction_snapshot(int spring_id,
                                                             const std::vector<double>& u,
                                                             const std::vector<double>& committed_state) const;
    // Material-field replacement preserves the compiled structural topology and
    // nonlinear state layout. It is intended for synchronous ASCE/ACI parameter
    // regeneration without rebuilding nodes, elements, MPCs, sparse patterns, or
    // nonlinear update directions.
    void replace_nonlinear_material(int spring_id, NonlinearMaterial material);
    void replace_nonlinear_material_field(const std::vector<std::pair<int,NonlinearMaterial>>& replacements);
    struct FSCShearStateSnapshot {
        bool valid{};
        bool initiated{};
        bool residual_reached{};
        double deformation{};
        double force{};
        double failure_deformation{};
        double failure_strength{};
        double cumulative_strength_decrement{};
        double retained_strength_ratio{1.0};
    };
    FSCShearStateSnapshot fsc_shear_state(int spring_index, const std::vector<double>& state) const;
    int node_count() const { return static_cast<int>(node_ids_.size()); }
    int elastic_element_count() const { return elastic_element_count_; }
    int elastic_element_index(int element_id) const;
    ElasticFrame3DResponse elastic_element_response(int element_id,
                                                    const std::vector<double>& u) const;
    const std::vector<double>& reduced_dof_elevations() const { return reduced_dof_z_; }

private:
    friend class Frame3DBuilder;
    int n_{};
    int elastic_element_count_{};
    std::vector<int> node_ids_;
    std::unordered_map<int,int> node_slot_;
    struct ElementResponseData {
        ElasticFrame3D properties;
        std::array<double,3> node_i_xyz{};
        std::array<double,3> node_j_xyz{};
        std::array<std::vector<std::pair<int,double>>,12> terms;
    };
    std::vector<ElementResponseData> element_response_data_;
    std::unordered_map<int,int> element_index_by_id_;
    std::vector<int> full_to_reduced_; // direct one-to-one mapping only; -1 for fixed/MPC rows
    std::vector<double> reduced_dof_z_; // representative elevation; NaN if a generalized DOF spans elevations
    std::vector<std::vector<std::pair<int,double>>> full_to_terms_; // u_full = T q
    std::vector<double> mass_; // generalized diagonal (diagnostic/backward-compatible)
    SparseMatrixCSC mass_matrix_;
    std::vector<double> base_mass_;
    SparseMatrixCSC k_linear_;
    SparseMatrixCSC k_initial_;
    // Physical stiffness reference used for Rayleigh beta*K damping. Numerical
    // penalty coordinates (notably the rigid-plastic P-M axial hinge penalty)
    // are intentionally excluded so they cannot create artificial damping.
    SparseMatrixCSC k_damping_reference_;
    std::vector<NonlinearMaterial> materials_;
    std::vector<int> component_ids_;
    std::unordered_map<int,int> component_index_by_id_;
    std::vector<int> material_state_offsets_;
    int nonlinear_state_size_{};
    std::vector<double> initial_tangents_;
    SparseUpdateBasis basis_;
    SparseMatrixCSC system_pattern_; // union of K and M sparsity, values unused
    std::vector<int> k_to_system_;
    std::vector<int> kd_to_system_;
    std::vector<int> m_to_system_;
    std::vector<std::vector<std::pair<int,double>>> spring_scatter_;
    // Same b*b^T coefficients indexed into k_initial_ rather than the union
    // system pattern. This lets a new material Ke update Kinitial in-place.
    std::vector<std::vector<std::pair<int,double>>> spring_initial_scatter_;
    std::vector<std::vector<std::pair<int,double>>> spring_damping_scatter_;
    struct CoupledPMCrossEntry { int system_pos{-1}; double coefficient{}; };
    struct CoupledPMHinge {
        int spring_index{-1};
        double axial_preload{};
        double axial_stiffness{};
        ASCE41HingeParams base_params;
        std::vector<std::pair<int,double>> axial_terms;
        std::vector<double> axial_force_points;
        std::vector<double> moment_capacity_points;
        std::vector<double> capacity_slopes;
        std::vector<CoupledPMCrossEntry> cross_scatter;
    };
    std::vector<CoupledPMHinge> coupled_pm_hinges_;
    std::vector<int> coupled_pm_by_spring_;
    struct CoupledPMEval { MaterialTrialResult material; double axial_force{}, capacity{}, dcapacity_dP{}, dmoment_dP{}; };
    CoupledPMEval evaluate_coupled_pm(int spring_index,double rotation,
                                      const std::vector<double>& u,
                                      const double* committed,double* trial_state) const;
    struct PMReturnScatterEntry { int system_pos{-1}; int local_row{},local_col{}; double coefficient{}; };
    struct PMReturnHinge {
        int spring_index{-1};
        int state_offset{-1};
        std::optional<PMInteractionHinge2D> law;
        std::vector<std::pair<int,double>> axial_terms;
        std::vector<PMReturnScatterEntry> block_scatter;
    };
    std::vector<PMReturnHinge> pm_return_hinges_;
    std::vector<int> pm_return_by_spring_;
    PMInteractionTrialResult evaluate_pm_return(int spring_index,double rotation,
                                                const std::vector<double>& u,
                                                const double* committed,double* trial_state) const;
    struct FSCShearSpring {
        int spring_index{-1};
        int state_offset{-1};
        double axial_preload{}, axial_stiffness{};
        FSCShearSpringLaw law;
        std::vector<std::pair<int,double>> bottom_rotation_terms;
        std::vector<std::pair<int,double>> top_rotation_terms;
        std::vector<std::pair<int,double>> axial_terms;
    };
    std::vector<FSCShearSpring> fsc_shear_springs_;
    std::vector<int> fsc_shear_by_spring_;
    FSCShearSpringTrial evaluate_fsc_shear(int spring_index,double shear_deformation,
                                           const std::vector<double>& u,
                                           const double* committed,double* trial_state) const;
    struct GeometricEntry { int row{}, col{}, system_pos{-1}; double value{}; };
    struct GeometricCouplingEntry { int row{}, col{}, system_pos{-1}, transpose_pos{-1}; double axial_coefficient{}; };
    struct GeometricUpdate {
        double preload{}, axial_k{};
        std::vector<std::pair<int,double>> axial_terms;
        std::vector<GeometricEntry> entries;          // unit geometric matrix G
        std::vector<GeometricCouplingEntry> coupling; // pattern for -(EA/L)(G u) a^T
    };
    std::vector<GeometricUpdate> geometric_updates_;
    struct GeometricLowRankBlock {
        int start_col{};
        int eigen_count{};
        int axial_col{-1};
        double axial_k{};
        std::vector<double> eigenvalues;
    };
    SparseUpdateBasis generalized_state_update_basis_;
    std::vector<GeometricLowRankBlock> geometric_low_rank_blocks_;
    struct CorotScatterEntry { int local_row{}, local_col{}, system_pos{-1}; double coefficient{}; };
    struct CorotElement {
        CorotationalFrame3DProperties properties;
        std::array<std::vector<std::pair<int,double>>,12> terms;
        std::vector<CorotScatterEntry> scatter;
    };
    std::vector<CorotElement> corotational_elements_;
    bool updated_pdelta_{false};
    bool corotational_{false};
    std::vector<std::pair<int,double>> response_terms_;
    std::vector<std::vector<std::pair<int,double>>> story_terms_;
    std::vector<double> story_z_;
    double alpha_m_{};
    double beta_k_{};
};

class Frame3DBuilder {
public:
    void add_node(int id, double x, double y, double z,
                  double mass_x=0.0, double mass_y=0.0, double mass_z=0.0,
                  double mass_rx=0.0, double mass_ry=0.0, double mass_rz=0.0);
    void add_elastic_frame(int id, int node_i, int node_j,
                           double E, double G, double A, double J,
                           double Iy, double Iz,
                           double ref_x, double ref_y, double ref_z,
                           double axial_compression=0.0);
    void add_bilinear_spring(int id, int node_i, Dof3D dof_i,
                             int node_j, Dof3D dof_j,
                             double k0, double yield_force,
                             double post_yield_ratio);
    void add_linear_bilinear_spring(int id, std::vector<MpcTerm3D> deformation_terms,
                                    double k0, double yield_force,
                                    double post_yield_ratio);
    void add_linear_imk_peak_oriented_spring(int id, std::vector<MpcTerm3D> deformation_terms,
                                             IMKPeakOrientedParams params);
    void add_linear_asce41_hinge(int id, std::vector<MpcTerm3D> deformation_terms,
                                 ASCE41HingeParams params);
    void add_linear_fsc_shear_spring(
                                 int id,
                                 std::vector<MpcTerm3D> shear_deformation_terms,
                                 std::vector<MpcTerm3D> bottom_rotation_terms,
                                 std::vector<MpcTerm3D> top_rotation_terms,
                                 std::vector<MpcTerm3D> axial_deformation_terms,
                                 double axial_preload, double axial_stiffness,
                                 FSCShearSpringParams params);
    void add_linear_axial_coupled_asce41_hinge(
                                 int id,
                                 std::vector<MpcTerm3D> rotation_terms,
                                 std::vector<MpcTerm3D> axial_deformation_terms,
                                 AxialCoupledASCE41Params params);
    // Phase 9F associative two-force P-M hinge. axial_deformation_terms
    // describe the zero-length hinge axial deformation, so the component
    // participates directly in axial equilibrium and plastic flow.
    void add_linear_pm_interaction_hinge(
                                 int id,
                                 std::vector<MpcTerm3D> rotation_terms,
                                 std::vector<MpcTerm3D> axial_deformation_terms,
                                 PMInteractionHingeParams params);
    void add_imk_peak_oriented_spring(int id, int node_i, Dof3D dof_i,
                                      int node_j, Dof3D dof_j,
                                      IMKPeakOrientedParams params);
    void add_asce41_hinge(int id, int node_i, Dof3D dof_i,
                          int node_j, Dof3D dof_j, ASCE41HingeParams params);
    // Relative rotation projected on an arbitrary global unit vector.
    void add_rotational_vector_spring(int id, int node_i, int node_j,
                                      double axis_x, double axis_y, double axis_z,
                                      double k0, double yield_moment,
                                      double post_yield_ratio);

    void fix(int node_id, bool ux=true, bool uy=true, bool uz=true,
             bool rx=true, bool ry=true, bool rz=true);
    void fix_dof(int node_id, Dof3D dof);
    void equal_dof(int master_node, int slave_node, Dof3D dof);
    // General linear MPC: u_slave = sum(coeff_i * u_master_i).
    void linear_constraint(int slave_node, Dof3D slave_dof,
                           std::vector<MpcTerm3D> masters);
    // Horizontal rigid diaphragm in global XY with master RZ coupling:
    // ux_s = ux_m - dy*rz_m, uy_s = uy_m + dx*rz_m, rz_s = rz_m.
    void rigid_diaphragm_z(int master_node, const std::vector<int>& slave_nodes);

    void set_rayleigh(double alpha_m, double beta_k);
    void set_ground_direction(Dof3D translational_dof);
    void set_response(int node_id, Dof3D dof);
    void set_story_nodes(std::vector<int> story_node_ids, Dof3D dof);
    // Research/reference path: update small-displacement geometric stiffness
    // from current member axial force. This is not yet a full corotational
    // formulation and intentionally forces a same-pattern direct solve.
    void set_updated_pdelta(bool enabled=true) { updated_pdelta_=enabled; }
    // Finite-rotation corotational reference path. This uses an objective
    // energy-based 3D element with a numerically differentiated consistent
    // tangent. It is intentionally correctness-first and currently uses the
    // direct state-tangent solve path.
    void set_corotational(bool enabled=true) { corotational_=enabled; }

    CompiledFrame3D compile() const;

private:
    struct FixConstraint { int node{}; Dof3D dof{}; };
    struct EqualConstraint { int master{}; int slave{}; Dof3D dof{}; };
    struct LinearConstraint { int slave{}; Dof3D dof{}; std::vector<MpcTerm3D> masters; };
    struct CoupledPMRaw { int spring_index{-1}; std::vector<MpcTerm3D> axial_terms; AxialCoupledASCE41Params params; };
    struct PMReturnRaw { int spring_index{-1}; std::vector<MpcTerm3D> axial_terms; PMInteractionHingeParams params; };
    struct FSCShearRaw {
        int spring_index{-1};
        std::vector<MpcTerm3D> bottom_rotation_terms,top_rotation_terms,axial_terms;
        double axial_preload{},axial_stiffness{};
        FSCShearSpringParams params;
    };
    std::vector<Node3D> nodes_;
    std::vector<ElasticFrame3D> elements_;
    std::vector<GeneralizedSpring3D> springs_;
    std::vector<CoupledPMRaw> coupled_pm_;
    std::vector<PMReturnRaw> pm_return_;
    std::vector<FSCShearRaw> fsc_shear_;
    std::vector<FixConstraint> fixed_;
    std::vector<EqualConstraint> equal_;
    std::vector<LinearConstraint> linear_;
    int response_node_id_{-1};
    Dof3D response_dof_{Dof3D::UX};
    std::vector<int> story_node_ids_;
    Dof3D story_dof_{Dof3D::UX};
    Dof3D ground_dof_{Dof3D::UX};
    double alpha_m_{0.0};
    double beta_k_{0.002};
    bool updated_pdelta_{false};
    bool corotational_{false};
};

std::array<double,144> frame3d_global_stiffness(
    double xi,double yi,double zi,
    double xj,double yj,double zj,
    double E,double G,double A,double J,double Iy,double Iz,
    const std::array<double,3>& reference,
    double axial_compression=0.0);

} // namespace quake
