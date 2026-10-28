#ifndef CUBIC_SPLINE_FACTORS_HPP
#define CUBIC_SPLINE_FACTORS_HPP

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "gtsam/geometry/Point3.h"
#include "gtsam/geometry/Pose3.h"
#include "gtsam/geometry/Rot3.h"
#include "iostream"

namespace CubicBasisSpline {

struct SplineEvaluation {
  gtsam::Point3 pos;
  gtsam::Rot3 orientation;
  gtsam::Vector3 lin_vel;
  gtsam::Vector3 lin_acc;
  gtsam::Vector3 ang_vel;
};

SplineEvaluation evaluateCubicBSpline(gtsam::Pose3& T_im1, gtsam::Pose3& T_i, gtsam::Pose3& T_ip1,
                                      gtsam::Pose3& T_ip2, double u);

}  // namespace CubicBasisSpline
#endif