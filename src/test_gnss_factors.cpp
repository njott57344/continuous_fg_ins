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

/*
Unit Tests:

1) We want to make sure that we can recover the actual pos,vel values given
truth ctrl points,time (i.e. with zero noise,a known model our residuals are zero)
This breaks down into 2 unit tests, one for position and one for velocity


*/

namespace CubicBasisSpline {
namespace testing {}  // namespace testing
}  // namespace CubicBasisSpline

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}