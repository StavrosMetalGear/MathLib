cmake_minimum_required(VERSION 3.20)
project(MathLib VERSION 0.1.0 LANGUAGES CXX)

option(MATHLIB_BUILD_TESTS "Build MathLib tests" ON)
option(MATHLIB_BUILD_EXAMPLES "Build MathLib examples" ON)

add_library(MathLib INTERFACE)
add_library(MathLib::MathLib ALIAS MathLib)

target_compile_features(MathLib INTERFACE cxx_std_20)
target_include_directories(MathLib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# -----------------------
# Examples
# -----------------------
if(MATHLIB_BUILD_EXAMPLES)
  add_executable(mathlib_demo examples/demo.cpp)
  target_link_libraries(mathlib_demo PRIVATE MathLib::MathLib)
endif()

# -----------------------
# Tests (GoogleTest via FetchContent)
# -----------------------
if(MATHLIB_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE) # helps on Windows
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_executable(mathlib_tests tests/test_linalg.cpp
  tests/test_core.cpp
  tests/test_solve.cpp
  tests/test_calculus.cpp
  )
  target_link_libraries(mathlib_tests PRIVATE MathLib::MathLib GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(mathlib_tests)
endif()
