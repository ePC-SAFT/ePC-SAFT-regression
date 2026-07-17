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
  src/epcsaft_regression/native/methane_contract.cpp
  src/epcsaft_regression/native/methane_fit.cpp
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
