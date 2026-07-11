cmake_minimum_required(VERSION 3.20)
project(vdjtools LANGUAGES CXX)

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
# Currently a minimal hot-path primitive (hamming) so the toolchain + seam are testable.
add_library(vdjtools_core STATIC
  src/core.cpp
  src/pgen.cpp
)
target_include_directories(vdjtools_core PUBLIC include PRIVATE src)

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
