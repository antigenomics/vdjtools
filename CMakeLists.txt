cmake_minimum_required(VERSION 3.20)

# Single source of truth for the version: parse it from pyproject.toml's [project] table so the
# C++ core, the pybind11 module and the C++ test can never drift from the Python package version.
# (A past release bumped src/core.cpp but not the test's hand-copied literal, turning CI red on a
# no-op — VDJTOOLS_VERSION below removes the literal entirely.)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" _vdjtools_pyproject)
string(REGEX MATCH "\nversion = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" _ "\n${_vdjtools_pyproject}")
if(NOT CMAKE_MATCH_1)
  message(FATAL_ERROR "could not parse [project].version from pyproject.toml")
endif()
project(vdjtools VERSION ${CMAKE_MATCH_1} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()  # MSVC Release already uses /O2; -O3 is not a valid MSVC flag

option(VDJTOOLS_TESTS  "Build C++ tests" OFF)
option(VDJTOOLS_PYTHON "Build the pybind11 module" OFF)

# Native core: Pgen DP, generation sampler, EM E-step inner loop land here (Phase 1).
# hamming is a minimal hot-path primitive; inext.cpp is the iNEXT size-based
# diversity kernel (curve + bootstrap + parallel batch).
add_library(vdjtools_core STATIC
  src/core.cpp
  src/inext.cpp
  src/pgen.cpp
)
target_include_directories(vdjtools_core PUBLIC include PRIVATE src)
# Inject the single-sourced version; PUBLIC so the pybind11 module and the C++ test (both link
# vdjtools_core) see the same VDJTOOLS_VERSION macro.
target_compile_definitions(vdjtools_core PUBLIC VDJTOOLS_VERSION="${PROJECT_VERSION}")

find_package(Threads REQUIRED)
target_link_libraries(vdjtools_core PUBLIC Threads::Threads)

if(VDJTOOLS_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core src/_bindings.cpp)
  target_link_libraries(_core PRIVATE vdjtools_core)
  install(TARGETS _core DESTINATION vdjtools)
endif()

if(VDJTOOLS_TESTS)
  enable_testing()
  add_executable(vdjtools_tests tests/cpp/test_core.cpp)
  target_link_libraries(vdjtools_tests PRIVATE vdjtools_core)
  add_test(NAME vdjtools_tests COMMAND vdjtools_tests)
endif()
