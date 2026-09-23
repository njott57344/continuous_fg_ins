cmake_minimum_required(VERSION 3.22)
project(ContinuousBSplineIns VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Find Dependencies
# ------------------------------------------------------------------------------

find_package(Eigen3 3.4 REQUIRED)
find_package(Ceres REQUIRED)

# ------------------------------------------------------------------------------
# Build Custom Shared/Static Library
# ------------------------------------------------------------------------------

add_library(spline_tools
    src/spline_tools.cpp
    src/so3_tools.cpp
)

target_include_directories(spline_tools PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(spline_tools PUBLIC
    Eigen3::Eigen
    Ceres::ceres
)

target_compile_options(spline_tools PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -O3
)

# ------------------------------------------------------------------------------
# Build Executable
# ------------------------------------------------------------------------------

add_executable(solve_basis_spline
    src/main.cpp
)

target_link_libraries(solve_basis_spline PRIVATE
    spline_tools
)