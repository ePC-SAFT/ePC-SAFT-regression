cmake_minimum_required(VERSION 3.24)
project(epcsaft_regression LANGUAGES CXX)

find_package(Python 3.13 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(Ceres 2.2 REQUIRED)

if(NOT DEFINED ENV{EPCSAFT_INCLUDE_DIR} OR "$ENV{EPCSAFT_INCLUDE_DIR}" STREQUAL "")
  message(FATAL_ERROR "EPCSAFT_INCLUDE_DIR must name the installed provider public include directory")
endif()

set(EPCSAFT_INCLUDE_DIR "$ENV{EPCSAFT_INCLUDE_DIR}")
if(NOT EXISTS "${EPCSAFT_INCLUDE_DIR}/epcsaft/native_sdk_v1.h")
  message(FATAL_ERROR "installed epcsaft/native_sdk_v1.h was not found under EPCSAFT_INCLUDE_DIR")
endif()

Python_add_library(_native MODULE WITH_SOABI
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

option(EPCSAFT_BUILD_BINARY_KIJ_PREFLIGHT "Build the non-installed binary-kij preflight" OFF)
if(EPCSAFT_BUILD_BINARY_KIJ_PREFLIGHT)
  Python_add_library(_binary_kij_preflight MODULE WITH_SOABI tools/binary_kij_preflight.cpp)
  target_compile_features(_binary_kij_preflight PRIVATE cxx_std_17)
  target_compile_options(_binary_kij_preflight PRIVATE -Wall -Wextra -Wpedantic -Werror)
  target_include_directories(_binary_kij_preflight PRIVATE "${EPCSAFT_INCLUDE_DIR}")
  target_link_libraries(_binary_kij_preflight PRIVATE Ceres::ceres)
  set_target_properties(_binary_kij_preflight PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/preflight"
  )
endif()
