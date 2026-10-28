#ifndef CUBIC_SPLINE_FACTORS_HPP
#define CUBIC_SPLINE_FACTORS_HPP

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "gtsam/geometry/Point3.h"
#include "gtsam/geometry/Pose3.h"
#include "gtsam/geometry/Rot3.h"
#include "iostream"

namespace CubicBasisSpline {

inline const Eigen::Matrix4d C =
    (Eigen::Matrix4d() << 6, 0, 0, 0, 5, 3, -3, 1, 1, 3, 3, -2, 0, 0, 0, 1).finished();

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

Eigen::Matrix<double, 4, 3> calculateWeights(double u);

#endif