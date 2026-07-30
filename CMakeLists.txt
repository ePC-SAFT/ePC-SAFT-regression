cmake_minimum_required(VERSION 3.24)
project(epcsaft_regression LANGUAGES CXX)

find_package(Python 3.13 REQUIRED COMPONENTS Interpreter Development.Module)
find_path(Unwind_INCLUDE_DIR NAMES unwind.h PATH_SUFFIXES libunwind)
find_package(glog CONFIG REQUIRED)
find_package(Ceres 2.2 REQUIRED)

if(NOT DEFINED ENV{EPCSAFT_INCLUDE_DIR} OR "$ENV{EPCSAFT_INCLUDE_DIR}" STREQUAL "")
  message(FATAL_ERROR "EPCSAFT_INCLUDE_DIR must name the installed EOS public include directory")
endif()

set(EPCSAFT_INCLUDE_DIR "$ENV{EPCSAFT_INCLUDE_DIR}")
if(NOT EXISTS "${EPCSAFT_INCLUDE_DIR}/epcsaft/native_sdk_v1.h")
  message(FATAL_ERROR "installed epcsaft/native_sdk_v1.h was not found under EPCSAFT_INCLUDE_DIR")
endif()

Python_add_library(_native MODULE WITH_SOABI
  src/epcsaft_regression/native/born_diameter_fit.cpp
  src/epcsaft_regression/native/ceres_core.cpp
  src/epcsaft_regression/native/evaluator_fit.cpp
  src/epcsaft_regression/native/figiel_kij_fit.cpp
  src/epcsaft_regression/native/figiel_water_factor_fit.cpp
  src/epcsaft_regression/native/general_fit.cpp
  src/epcsaft_regression/native/pure_saturation_contract.cpp
  src/epcsaft_regression/native/pure_saturation_fit.cpp
  src/epcsaft_regression/native/module.cpp
)
target_compile_features(_native PRIVATE cxx_std_17)
target_compile_options(_native PRIVATE -Wall -Wextra -Wpedantic -Werror)
target_include_directories(_native PRIVATE "${EPCSAFT_INCLUDE_DIR}")
target_link_libraries(_native PRIVATE Ceres::ceres)
set_target_properties(_native PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
)
install(TARGETS _native LIBRARY DESTINATION epcsaft_regression)

include(CTest)
if(BUILD_TESTING)
  add_executable(ceres_core_test EXCLUDE_FROM_ALL
    tests/native/test_ceres_core.cpp
    tests/native/test_evaluator_fit.cpp
    src/epcsaft_regression/native/ceres_core.cpp
    src/epcsaft_regression/native/evaluator_fit.cpp
  )
  target_compile_features(ceres_core_test PRIVATE cxx_std_17)
  target_compile_definitions(ceres_core_test PRIVATE
    EPCSAFT_REGRESSION_EVALUATOR_CORE_ONLY
  )
  target_compile_options(ceres_core_test PRIVATE
    -Wall -Wextra -Wpedantic -Werror
  )
  target_include_directories(ceres_core_test PRIVATE
    src/epcsaft_regression/native
  )
  target_link_libraries(ceres_core_test PRIVATE Ceres::ceres)
  add_test(NAME ceres_core COMMAND ceres_core_test)
endif()
