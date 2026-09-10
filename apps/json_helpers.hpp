#pragma once
#include <nlohmann/json.hpp>
#include "quake/newmark.hpp"
#include <fstream>
using Json=nlohmann::json;
inline Json read_json(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("cannot open "+path);Json j;f>>j;return j;}
inline void write_json(const std::string& path,const Json& j){std::ofstream f(path);if(!f)throw std::runtime_error("cannot write "+path);f<<j.dump(2)<<'\n';if(!f)throw std::runtime_error("failed to write "+path);}
inline std::string termination_name(quake::AnalysisTermination t){switch(t){case quake::AnalysisTermination::Completed:return "completed";case quake::AnalysisTermination::PhysicalCollapse:return "configured_collapse_criterion";case quake::AnalysisTermination::InitialInstability:return "initial_instability";default:return "numerical_failure";}}
inline Json analysis_json(const quake::AnalysisResult& a){Json j={{"termination",termination_name(a.termination)},{"reason",a.termination_reason},{"last_integration_error",a.last_integration_error},{"collapse_mechanism_code",static_cast<int>(a.collapse_mechanism)},{"termination_time_s",a.termination_time},{"peak_story_drift_ratio",a.max_story_drift_ratio},{"final_displacement",a.final_displacement}};
#define STAT(x) j["stats"][#x]=a.stats.x
STAT(last_residual_norm);STAT(last_residual_tolerance);STAT(steps);STAT(newton_iterations);STAT(global_factorizations);STAT(linear_solves);STAT(residual_evaluations);STAT(constitutive_integration_failures);STAT(nonfinite_evaluations);STAT(line_search_backtracks);STAT(subdivided_steps);STAT(internal_substeps);STAT(direct_fallbacks);STAT(linear_solve_failures);STAT(line_search_failures);STAT(max_iteration_failures);STAT(max_active_rank);STAT(minimum_dt);STAT(elapsed_seconds);STAT(max_roof_abs);STAT(max_story_drift_abs);STAT(input_energy);STAT(internal_work);STAT(damping_energy);STAT(kinetic_energy);STAT(energy_balance_relative_error);
#undef STAT
return j;}
inline quake::LinearStrategy parse_strategy(const std::string& s){if(s=="woodbury")return quake::LinearStrategy::Woodbury;if(s=="same_pattern")return quake::LinearStrategy::SamePatternRefactorization;if(s=="full")return quake::LinearStrategy::FullFactorization;throw std::invalid_argument("strategy must be full, same_pattern, or woodbury");}
