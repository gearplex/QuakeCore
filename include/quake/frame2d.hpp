#pragma once
#include "quake/dynamic_model.hpp"
#include "quake/bilinear.hpp"
#include "quake/nonlinear_material.hpp"
#include "quake/steel2d.hpp"
#include "quake/wall2d.hpp"

#include <array>
#include <unordered_map>
#include <utility>
#include <vector>

namespace quake {

enum class Dof2D : int { UX = 0, UY = 1, RZ = 2 };

struct Node2D {
    int id{};
    double x{};
    double y{};
    double mass_x{};
    double mass_y{};
    double mass_r{};
};

struct ElasticFrame2D {
    int id{};
    int node_i{};
    int node_j{};
    double E{};
    double A{};
    double I{};
    // Positive means compression. This is a constant preload geometric stiffness
    // used for P-Delta-capable Phase 1 benchmarking; dynamic axial-force updates
    // are intentionally deferred to a later milestone.
    double axial_compression{};
    bool pdelta_transformation{}; // OpenSees-style sway P/L contribution
};

enum class NonlinearComponent2DKind { RotationalHinge, PanelZone, BRB, SoilSpring };
enum class SteelMember2DRole { Beam, Column };

struct ScalarComponent2D {
    int id{};
    int node_i{};
    int node_j{};
    NonlinearComponent2DKind kind{NonlinearComponent2DKind::RotationalHinge};
    NonlinearMaterial material{BilinearSpring(1.0,1.0,0.0)};
    bool explicit_direction{};
    double direction_x{};
    double direction_y{};
    double direction_rz{};
};

// Committed-state member force recovery used by native RC-column demand
// recording. Compression is positive. End shears/moments follow the element
// local sign convention.
struct ElasticFrame2DResponse {
    int element_id{};
    double axial_compression{};
    double shear_i{};
    double shear_j{};
    double moment_i{};
    double moment_j{};
};

class CompiledFrame2D final : public NonlinearDynamicModel {
public:
    int dof() const override { return n_; }
    int nonlinear_count() const override { return static_cast<int>(materials_.size()); }
    const SparseMatrixCSC& K_linear() const { return k_linear_; }
    const SparseMatrixCSC& K_initial() const override { return k_initial_; }
    const std::vector<double>& mass() const override { return mass_; }
    std::vector<double> mass_multiply(const std::vector<double>& a) const override;
    const SparseUpdateBasis& nonlinear_basis() const override { return basis_; }
    const std::vector<NonlinearMaterial>& materials() const { return materials_; }
    int nonlinear_state_size() const override { return nonlinear_state_size_; }
    std::vector<double> initial_nonlinear_state() const override;
    const std::vector<double>& initial_nonlinear_tangents() const override { return initial_tangents_; }

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
    double max_drift_ratio(const std::vector<double>& u) const override;
    std::vector<double> story_response_values(const std::vector<double>& x) const;
    const std::vector<double>& story_elevations() const { return story_y_; }

    int wall_count() const { return static_cast<int>(walls_.size()); }
    std::vector<int> wall_ids() const;
    Wall2DResponse wall_response(int id,const std::vector<double>& u,const std::vector<double>& committed) const;
    int steel_member_count() const { return static_cast<int>(steel_members_.size()); }
    std::vector<int> steel_member_ids() const;
    SteelMember2DRole steel_member_role(int id) const;
    SteelMember2DResponse steel_member_response(int id,const std::vector<double>& u,const std::vector<double>& committed) const;
    int viscous_damper_count() const { return static_cast<int>(dampers_.size()); }
    std::vector<int> viscous_damper_ids() const;
    ViscousDamper2DTrial viscous_damper_response(int id,const std::vector<double>& v) const;
    NonlinearComponent2DKind nonlinear_component_kind(int id) const;
    bool has_state_dependent_global_tangent() const override { return !walls_.empty()||!steel_members_.empty(); }
    bool has_velocity_dependent_global_tangent() const override { return !dampers_.empty(); }
    SparseMatrixCSC effective_state_tangent_matrix(const std::vector<double>& u,const std::vector<double>& tangents,double a0,double a1) const override;
    SparseMatrixCSC effective_state_tangent_matrix_with_state(const std::vector<double>& u,const std::vector<double>& tangents,const std::vector<double>& committed,double a0,double a1) const override;
    SparseMatrixCSC effective_state_tangent_matrix_with_state_and_velocity(const std::vector<double>& u,const std::vector<double>& tangents,const std::vector<double>& committed,const std::vector<double>& velocity,double a0,double a1) const override;

    int reduced_dof(int node_id, Dof2D dof) const;
    int node_count() const { return static_cast<int>(node_ids_.size()); }
    int elastic_element_count() const { return elastic_element_count_; }
    const std::vector<int>& node_ids() const { return node_ids_; }
    int elastic_element_index(int element_id) const;
    int nonlinear_component_index(int spring_id) const;
    int nonlinear_component_id(int index) const;
    NonlinearComponentSnapshot nonlinear_component_snapshot(int spring_id,
                                                             const std::vector<double>& u,
                                                             const std::vector<double>& committed_state) const;
    // Replace only the nonlinear material field while retaining compiled DOFs,
    // sparse topology, update basis, element maps, and MPC reduction. Replacement
    // materials must preserve the existing state size so rollback/state layouts
    // remain immutable across ASCE parameter iterations.
    void replace_nonlinear_material(int spring_id, NonlinearMaterial material);
    void replace_nonlinear_material_field(const std::vector<std::pair<int,NonlinearMaterial>>& replacements);
    ElasticFrame2DResponse elastic_element_response(int element_id,
                                                    const std::vector<double>& u) const;

private:
    friend class Frame2DBuilder;
    struct CompiledWall {
        int id{}, state_offset{};
        Wall2D element;
        std::array<int,6> dofs{};
        std::array<int,36> scatter{};
        std::array<double,36> k0{};
    };
    std::vector<CompiledWall> walls_;
    struct CompiledSteelMember {
        int id{}, state_offset{};
        SteelMember2DRole role{SteelMember2DRole::Beam};
        SteelMember2D element;
        std::array<int,6> dofs{};
        std::array<int,36> scatter{};
        std::array<double,36> k0{};
    };
    std::vector<CompiledSteelMember> steel_members_;
    struct CompiledDamper {
        int id{};
        ViscousDamper2D element;
        std::array<int,6> dofs{};
        std::array<double,6> b{};
        std::array<int,36> scatter{};
        double c0{};
    };
    std::vector<CompiledDamper> dampers_;
    int n_{};
    int elastic_element_count_{};
    std::vector<int> node_ids_;
    std::unordered_map<int, int> node_slot_;
    struct ElementResponseData {
        ElasticFrame2D properties;
        double xi{}, yi{}, xj{}, yj{};
        std::array<int,6> reduced_dofs{};
    };
    std::vector<ElementResponseData> element_response_data_;
    std::unordered_map<int,int> element_index_by_id_;
    std::unordered_map<int,int> spring_index_by_id_;
    std::vector<int> spring_ids_;
    std::vector<NonlinearComponent2DKind> component_kinds_;
    std::vector<int> full_to_reduced_; // 3*node_count
    std::vector<double> mass_;
    std::vector<double> base_mass_x_;
    SparseMatrixCSC k_linear_;
    SparseMatrixCSC k_initial_;
    std::vector<NonlinearMaterial> materials_;
    std::vector<int> material_state_offsets_;
    int nonlinear_state_size_{};
    std::vector<double> initial_tangents_;
    SparseUpdateBasis basis_;
    std::vector<std::vector<std::pair<int,double>>> spring_scatter_; // CSC value index, b_i*b_j
    std::vector<int> diagonal_positions_;
    int response_dof_{-1};
    std::vector<int> story_dofs_;
    std::vector<double> story_y_;
    double alpha_m_{};
    double beta_k_{};
};

class Frame2DBuilder {
public:
    void add_node(int id, double x, double y,
                  double mass_x = 0.0, double mass_y = 0.0, double mass_r = 0.0);
    void add_elastic_frame(int id, int node_i, int node_j,
                           double E, double A, double I,
                           double axial_compression = 0.0,
                           bool pdelta_transformation = false);
    void add_rotational_spring(int id, int node_i, int node_j,
                               double k0, double yield_moment,
                               double post_yield_ratio);
    void add_asce41_hinge(int id, int node_i, int node_j, ASCE41HingeParams params);
    void add_imk_peak_oriented_hinge(int id, int node_i, int node_j, IMKPeakOrientedParams params);

    void add_panel_zone(int id,int node_i,int node_j,NonlinearMaterial material);
    void add_brb(int id,int node_i,int node_j,double k0,double yield_force,double post_yield_ratio);
    void add_soil_spring(int id,int node_i,int node_j,double direction_x,double direction_y,
                         NonlinearMaterial material);
    void add_rotational_soil_spring(int id,int node_i,int node_j,NonlinearMaterial material);
    void add_steel_member(int id,int node_i,int node_j,SteelMember2DRole role,SteelMember2DProperties properties);
    void add_viscous_damper(int id,int node_i,int node_j,ViscousDamper2DProperties properties);
    void add_directional_viscous_damper(int id,int node_i,int node_j,double direction_x,
                                        double direction_y,ViscousDamper2DProperties properties);
    void add_rotational_soil_dashpot(int id,int node_i,int node_j,ViscousDamper2DProperties properties);

    void add_mvlem(int id,int node_i,int node_j,MVLEMProperties properties);
    void add_sfi_mvlem(int id,int node_i,int node_j,SFIMVLEMProperties properties);

    void fix(int node_id, bool ux = true, bool uy = true, bool rz = true);
    void equal_dof(int master_node, int slave_node, Dof2D dof);
    void rigid_floor_x(int master_node, const std::vector<int>& slave_nodes);

    void set_rayleigh(double alpha_m, double beta_k);
    void set_response_node(int node_id);
    void set_story_nodes(std::vector<int> story_node_ids);

    CompiledFrame2D compile() const;

private:
    struct EqualConstraint { int master{}; int slave{}; Dof2D dof{}; };
    struct WallDefinition { int id{}, i{}, j{}; std::variant<MVLEMProperties,SFIMVLEMProperties> properties; };
    struct SteelMemberDefinition { int id{},i{},j{};SteelMember2DRole role{};SteelMember2DProperties properties; };
    struct DamperDefinition { int id{},i{},j{};ViscousDamper2DProperties properties;bool explicit_direction{};double direction_x{},direction_y{},direction_rz{}; };
    std::vector<WallDefinition> walls_;
    std::vector<SteelMemberDefinition> steel_members_;
    std::vector<DamperDefinition> dampers_;
    std::vector<Node2D> nodes_;
    std::vector<ElasticFrame2D> elements_;
    std::vector<ScalarComponent2D> springs_;
    std::vector<std::array<bool,3>> fixed_;
    std::vector<int> fixed_node_ids_;
    std::vector<EqualConstraint> equal_;
    int response_node_id_{-1};
    std::vector<int> story_node_ids_;
    double alpha_m_{0.0};
    double beta_k_{0.002};
};

// 6x6 global tangent for a prismatic Euler-Bernoulli 2D frame member.
// axial_compression > 0 subtracts the standard small-displacement beam-column
// geometric stiffness associated with a constant compressive preload.
std::array<double,36> frame2d_global_stiffness(double xi, double yi,
                                               double xj, double yj,
                                               double E, double A, double I,
                                               double axial_compression = 0.0,
                                               bool pdelta_transformation = false);

} // namespace quake
