#include <ceres/ceres.h>
#include <ceres/gradient_checker.h>
#include <ceres/manifold.h>
#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <sophus/so3.hpp>
#include <vector>

#include "spline_tools/spline_factors.hpp"

namespace CubicBasisSpline {
namespace testing {}  // namespace testing
}  // namespace CubicBasisSpline

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}