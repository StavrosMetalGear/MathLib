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
  tests/test_next_math.cpp
  )
  target_link_libraries(mathlib_tests PRIVATE MathLib::MathLib GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(mathlib_tests)
endif()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Your library target
# (Assuming you already have: add_library(MathLib INTERFACE) and ALIAS MathLib::MathLib)

# Install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install/export the target
install(TARGETS MathLib
        EXPORT MathLibTargets)

install(EXPORT MathLibTargets
        NAMESPACE MathLib::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MathLib)

# Create a Config file for find_package(MathLib CONFIG REQUIRED)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MathLibConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MathLibConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/MathLibConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MathLib
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/MathLibConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/MathLibConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MathLib
)
