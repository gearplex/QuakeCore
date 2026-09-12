cmake_minimum_required(VERSION 3.20)
project(quake_core VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(QUAKE_ENABLE_SANITIZERS "Enable ASan/UBSan in debug builds" OFF)

find_path(SUPERLU_INCLUDE_DIR slu_ddefs.h PATH_SUFFIXES superlu)
find_library(SUPERLU_LIBRARY NAMES superlu)
find_package(LAPACK REQUIRED)
find_package(Threads REQUIRED)

if(SUPERLU_INCLUDE_DIR AND SUPERLU_LIBRARY)
  set(QUAKE_DIRECT_SOLVER_SOURCE src/superlu_solver.cpp)
  set(QUAKE_DIRECT_SOLVER_LIB ${SUPERLU_LIBRARY})
  message(STATUS "QuakeCore direct solver: SuperLU")
else()
  set(QUAKE_DIRECT_SOLVER_SOURCE src/dense_lapack_solver.cpp)
  set(QUAKE_DIRECT_SOLVER_LIB "")
  message(WARNING "SuperLU not found: using LAPACK dense fallback (verification only)")
endif()

add_library(quake_core
    src/bilinear.cpp
    src/imk_peak_oriented.cpp
    src/asce41_hinge.cpp
    src/pm_interaction_hinge.cpp
    src/rc_column_asce41.cpp
    src/rc_building_asce41.cpp
    src/rc_frame_asce41.cpp
    src/model_comparison.cpp
    src/nonlinear_material.cpp
    src/sparse.cpp
    src/update_basis.cpp
    ${QUAKE_DIRECT_SOLVER_SOURCE}
    src/low_rank_solver.cpp
  src/generalized_woodbury.cpp
    src/shear_building.cpp
    src/frame2d.cpp
    src/wall_material.cpp
    src/wall2d.cpp
    src/steel2d.cpp
    src/soil2d.cpp
    src/frame3d.cpp
  src/corotational3d.cpp
    src/newmark.cpp
    src/newmark_robust.cpp
    src/modal.cpp
    src/modal_damping.cpp
    src/reduction.cpp
    src/suite.cpp
  src/ida.cpp
    src/calibration.cpp
    src/substructure.cpp
    src/stability.cpp
    src/static_analysis.cpp
)
target_include_directories(quake_core PUBLIC include)
if(SUPERLU_INCLUDE_DIR)
  target_include_directories(quake_core PRIVATE ${SUPERLU_INCLUDE_DIR})
endif()
target_link_libraries(quake_core PUBLIC ${QUAKE_DIRECT_SOLVER_LIB} ${LAPACK_LIBRARIES} Threads::Threads)
target_compile_options(quake_core PRIVATE -Wall -Wextra -Wpedantic)

if(QUAKE_ENABLE_SANITIZERS)
  target_compile_options(quake_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(quake_core PUBLIC -fsanitize=address,undefined)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

add_executable(quake_tests tests/test_main.cpp)
target_link_libraries(quake_tests PRIVATE quake_core)

add_executable(low_rank_bench benchmarks/low_rank_bench.cpp)
target_link_libraries(low_rank_bench PRIVATE quake_core)

add_executable(nrha_bench benchmarks/nrha_bench.cpp)
target_link_libraries(nrha_bench PRIVATE quake_core)

add_executable(frame_nrha_bench benchmarks/frame_nrha_bench.cpp)
target_link_libraries(frame_nrha_bench PRIVATE quake_core)

add_executable(frame_suite_bench benchmarks/frame_suite_bench.cpp)
target_link_libraries(frame_suite_bench PRIVATE quake_core)

add_executable(frame3d_nrha_bench benchmarks/frame3d_nrha_bench.cpp)
target_link_libraries(frame3d_nrha_bench PRIVATE quake_core)

add_executable(frame3d_diaphragm_bench benchmarks/frame3d_diaphragm_bench.cpp)
target_link_libraries(frame3d_diaphragm_bench PRIVATE quake_core)

enable_testing()
add_test(NAME quake_tests COMMAND quake_tests)
add_executable(quake_reliability_tests tests/test_reliability.cpp)
target_link_libraries(quake_reliability_tests PRIVATE quake_core)
add_test(NAME quake_reliability_tests COMMAND quake_reliability_tests)

add_executable(reduction_bench benchmarks/reduction_bench.cpp)
target_link_libraries(reduction_bench PRIVATE quake_core)

add_executable(rank_crossover_bench benchmarks/rank_crossover_bench.cpp)
target_link_libraries(rank_crossover_bench PRIVATE quake_core)

add_executable(frame3d_suite_bench benchmarks/frame3d_suite_bench.cpp)
target_link_libraries(frame3d_suite_bench PRIVATE quake_core)

add_executable(calibration_bench benchmarks/calibration_bench.cpp)
target_link_libraries(calibration_bench PRIVATE quake_core)

add_executable(substructure_bench benchmarks/substructure_bench.cpp)
target_link_libraries(substructure_bench PRIVATE quake_core)

add_executable(imk_active_front_bench benchmarks/imk_active_front_bench.cpp)
target_link_libraries(imk_active_front_bench PRIVATE quake_core)

add_executable(collapse_bench benchmarks/collapse_bench.cpp)
target_link_libraries(collapse_bench PRIVATE quake_core)

add_executable(stability_bench benchmarks/stability_bench.cpp)
target_link_libraries(stability_bench PRIVATE quake_core)

add_executable(corotational_bench benchmarks/corotational_bench.cpp)
target_link_libraries(corotational_bench PRIVATE quake_core)
add_executable(geometric_update_bench benchmarks/geometric_update_bench.cpp)
target_link_libraries(geometric_update_bench PRIVATE quake_core)
add_executable(ida_bench benchmarks/ida_bench.cpp)
target_link_libraries(ida_bench PRIVATE quake_core)

add_executable(corotational_buckling_bench benchmarks/corotational_buckling_bench.cpp)
target_link_libraries(corotational_buckling_bench PRIVATE quake_core)

add_executable(ucb_validation_run validation/uc_berkeley_3story/run_validation.cpp)
target_link_libraries(ucb_validation_run PRIVATE quake_core)
add_executable(ucb_validation_refined validation/uc_berkeley_3story/run_validation_refined.cpp)
target_link_libraries(ucb_validation_refined PRIVATE quake_core)

add_executable(ucb_validation_phase8 validation/uc_berkeley_3story/run_validation_phase8.cpp)
target_link_libraries(ucb_validation_phase8 PRIVATE quake_core)

add_executable(ucb_phase9f_pm_experiment validation/uc_berkeley_3story/run_phase9f_pm_experiment.cpp)
target_link_libraries(ucb_phase9f_pm_experiment PRIVATE quake_core)
add_executable(ucb_phase9f2_mroz_experiment validation/uc_berkeley_3story/run_phase9f2_mroz_experiment.cpp)
target_link_libraries(ucb_phase9f2_mroz_experiment PRIVATE quake_core)

add_executable(ucb_phase9g_damping_experiment validation/uc_berkeley_3story/run_phase9g_damping_experiment.cpp)
target_link_libraries(ucb_phase9g_damping_experiment PRIVATE quake_core)

add_executable(ucb_phase9g_hybrid_experiment validation/uc_berkeley_3story/run_phase9g_hybrid_experiment.cpp)
target_link_libraries(ucb_phase9g_hybrid_experiment PRIVATE quake_core)

add_executable(ucb_phase9g_modal3_experiment validation/uc_berkeley_3story/run_phase9g_modal3_experiment.cpp)
target_link_libraries(ucb_phase9g_modal3_experiment PRIVATE quake_core)

add_executable(ucb_phase9h_statepath_audit validation/uc_berkeley_3story/run_phase9h_statepath_audit.cpp)
target_link_libraries(ucb_phase9h_statepath_audit PRIVATE quake_core)

add_executable(quake_pm_safety_tests tests/test_pm_safety.cpp)
target_link_libraries(quake_pm_safety_tests PRIVATE quake_core)
add_test(NAME quake_pm_safety_tests COMMAND quake_pm_safety_tests)
add_executable(pm_audit validation/component_audit/pm_audit.cpp)
target_link_libraries(pm_audit PRIVATE quake_core)

add_executable(quake_run apps/quake_run.cpp)
target_include_directories(quake_run PRIVATE third_party)
target_link_libraries(quake_run PRIVATE quake_core)
add_executable(ucb_phase9i validation/uc_berkeley_3story/run_phase9i_validation.cpp)
target_include_directories(ucb_phase9i PRIVATE third_party apps)
target_link_libraries(ucb_phase9i PRIVATE quake_core)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME quake_job_contracts COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_job_cli.py $<TARGET_FILE:quake_run> ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_executable(quake_wall_tests tests/test_walls.cpp)
target_link_libraries(quake_wall_tests PRIVATE quake_core)
add_test(NAME quake_wall_tests COMMAND quake_wall_tests)
add_executable(quake_wall_material_wrapper_tests tests/test_wall_material_wrappers.cpp)
target_link_libraries(quake_wall_material_wrapper_tests PRIVATE quake_core)
add_test(NAME quake_wall_material_wrapper_tests COMMAND quake_wall_material_wrapper_tests)
add_executable(quake_wall_wrapper_lifecycle_tests tests/test_wall_wrapper_lifecycle.cpp)
target_link_libraries(quake_wall_wrapper_lifecycle_tests PRIVATE quake_core)
add_test(NAME quake_wall_wrapper_lifecycle_tests COMMAND quake_wall_wrapper_lifecycle_tests)
set_tests_properties(quake_wall_wrapper_lifecycle_tests PROPERTIES LABELS "wall;material;gate4;lifecycle")
add_executable(quake_wall2d_wrapper_tests tests/test_wall2d_wrappers.cpp)
target_link_libraries(quake_wall2d_wrapper_tests PRIVATE quake_core)
add_test(NAME quake_wall2d_wrapper_tests COMMAND quake_wall2d_wrapper_tests)
add_executable(quake_steel_component_tests tests/test_steel_components.cpp)
target_link_libraries(quake_steel_component_tests PRIVATE quake_core)
add_test(NAME quake_steel_component_tests COMMAND quake_steel_component_tests)
add_executable(quake_soil_component_tests tests/test_soil_components.cpp)
target_link_libraries(quake_soil_component_tests PRIVATE quake_core)
add_test(NAME quake_soil_component_tests COMMAND quake_soil_component_tests)
add_executable(wall_probe validation/walls/wall_probe.cpp)
target_include_directories(wall_probe PRIVATE third_party apps)
target_link_libraries(wall_probe PRIVATE quake_core)
add_executable(steel_probe validation/steel/steel_probe.cpp)
target_include_directories(steel_probe PRIVATE third_party apps)
target_link_libraries(steel_probe PRIVATE quake_core)
if(Python3_Interpreter_FOUND)
  add_test(NAME quake_wall_job_contracts COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wall_jobs.py $<TARGET_FILE:quake_run> ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME quake_steel_job_contracts COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_steel_jobs.py $<TARGET_FILE:quake_run> ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME quake_soil_job_contracts COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_soil_jobs.py $<TARGET_FILE:quake_run> ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME quake_yori_material_oracle_integrity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/validation/yori_wall_10001/verify_material_oracle.py
    ${CMAKE_CURRENT_SOURCE_DIR}/validation/yori_wall_10001/material_oracle_story1.json)
endif()
