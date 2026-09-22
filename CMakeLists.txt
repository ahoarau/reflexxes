cmake_minimum_required(VERSION 3.22)
project(reflexxes VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(REFLEXXES_BUILD_EXAMPLES "Build examples" ON)
option(REFLEXXES_BUILD_TESTS "Build tests" ON)
option(REFLEXXES_BUILD_BENCHMARKS "Build Google Benchmark comparison against Ruckig" OFF)

find_package(absl REQUIRED CONFIG)

set(REFLEXXES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/intrinsic/icon/reflexxes)

add_library(reflexxes
  ${REFLEXXES_DIR}/calc_min_execution_time_utils.cc
  ${REFLEXXES_DIR}/input_output_utility_functions.cc
  ${REFLEXXES_DIR}/inputs.cc
  ${REFLEXXES_DIR}/motion_polynomials.cc
  ${REFLEXXES_DIR}/outputs.cc
  ${REFLEXXES_DIR}/polynomial_solvers.cc
  ${REFLEXXES_DIR}/position_calc_exec_time_funcs.cc
  ${REFLEXXES_DIR}/position_calc_inoperative_zone_end.cc
  ${REFLEXXES_DIR}/position_calc_inoperative_zone_start.cc
  ${REFLEXXES_DIR}/position_calc_min_execution_time.cc
  ${REFLEXXES_DIR}/position_calc_min_execution_time_per_dof.cc
  ${REFLEXXES_DIR}/position_calc_traj_funcs.cc
  ${REFLEXXES_DIR}/position_inputs.cc
  ${REFLEXXES_DIR}/position_synchronize_dofs.cc
  ${REFLEXXES_DIR}/reflexxes_api.cc
  ${REFLEXXES_DIR}/synchronize_dofs_base.cc
  ${REFLEXXES_DIR}/velocity_calc_min_execution_time.cc
  ${REFLEXXES_DIR}/velocity_inputs.cc
  ${REFLEXXES_DIR}/velocity_synchronize_dofs.cc
)
add_library(reflexxes::reflexxes ALIAS reflexxes)

target_include_directories(reflexxes PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(reflexxes PUBLIC
  absl::check
  absl::log
  absl::inlined_vector
  absl::algorithm_container
  absl::bind_front
  absl::function_ref
  absl::strings
  absl::str_format
  absl::time
  absl::span
)
target_compile_features(reflexxes PUBLIC cxx_std_20)
set_target_properties(reflexxes PROPERTIES
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
)

if(REFLEXXES_BUILD_EXAMPLES)
  add_executable(reflexxes_position_example examples/position_example.cc)
  target_link_libraries(reflexxes_position_example PRIVATE reflexxes::reflexxes)
  set_target_properties(reflexxes_position_example PROPERTIES CXX_EXTENSIONS OFF)
endif()

if(REFLEXXES_BUILD_TESTS)
  enable_testing()
  if(REFLEXXES_BUILD_EXAMPLES)
    add_test(NAME position_example COMMAND reflexxes_position_example)
  endif()
endif()

if(REFLEXXES_BUILD_BENCHMARKS)
  find_package(benchmark REQUIRED CONFIG)
  find_package(ruckig REQUIRED CONFIG)

  add_executable(reflexxes_ruckig_benchmark benchmarks/ruckig_comparison.cc)
  target_link_libraries(reflexxes_ruckig_benchmark PRIVATE
    reflexxes::reflexxes ruckig::ruckig benchmark::benchmark)
  set_target_properties(reflexxes_ruckig_benchmark PROPERTIES CXX_EXTENSIONS OFF)
endif()

# Install
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS reflexxes EXPORT reflexxesTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY src/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
install(EXPORT reflexxesTargets NAMESPACE reflexxes::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/reflexxes)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/reflexxesConfig.cmake.in [=[
@PACKAGE_INIT@
include(CMakeFindDependencyMacro)
find_dependency(absl CONFIG)
include("${CMAKE_CURRENT_LIST_DIR}/reflexxesTargets.cmake")
]=])
configure_package_config_file(
  ${CMAKE_CURRENT_BINARY_DIR}/reflexxesConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/reflexxesConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/reflexxes)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/reflexxesConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/reflexxesConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/reflexxesConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/reflexxes)
