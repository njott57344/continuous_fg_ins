#ifndef CUBIC_SPLINE_FACTORS_HPP
#define CUBIC_SPLINE_FACTORS_HPP

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "gtsam/geometry/Point3.h"
#include "gtsam/geometry/Pose3.h"
#include "gtsam/geometry/Rot3.h"
#include "gtsam/geometry/SO3.h"
#include "iostream"

namespace CubicBasisSpline {

inline const Eigen::Matrix4d C =
    (Eigen::Matrix4d() << 6, 0, 0, 0, 5, 3, -3, 1, 1, 3, 3, -2, 0, 0, 0, 1).finished();

struct SplineEvaluation {
  gtsam::Point3 pos = gtsam::Point3::Zero();
  gtsam::Rot3 orientation = gtsam::Rot3::Identity();
  gtsam::Vector3 lin_vel = gtsam::Vector3::Zero();
  gtsam::Vector3 lin_acc = gtsam::Vector3::Zero();
  gtsam::Vector3 ang_vel = gtsam::Vector3::Zero();
};

struct SplineKnotSet {
  gtsam::Pose3 T_im1;
  gtsam::Pose3 T_i;
  gtsam::Pose3 T_ip1;
  gtsam::Pose3 T_ip2;
  double u;
  double dT;
};

bool gnssFactor(SplineKnotSet& ctrl_pts, gtsam::Point3& pos, gtsam::Vector3& vel,
                gtsam::Vector3& l_a, gtsam::Point3& residual_pos, gtsam::Vector3& residual_vel);

bool evaluateCubicBSpline(SplineKnotSet& ctrl_pts, SplineEvaluation& spline_eval);

}  // namespace CubicBasisSpline

Eigen::Matrix<double, 4, 3> calculateWeights(double u);

#endif