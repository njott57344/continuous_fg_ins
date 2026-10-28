#include "spline_tools/cubic_spline_factors.hpp"

namespace CubicBasisSpline {

SplineEvaluation evaluateCubicBSpline(gtsam::Pose3& T_im1, gtsam::Pose3& T_i, gtsam::Pose3& T_ip1,
                                      gtsam::Pose3& T_ip2, double u) {
  Eigen::Matrix<double, 4, 3> U = calculateWeights(u);

  Eigen::Vector3d p_im1 = T_im1.translation();
  Eigen::Vector3d p_i = T_i.translation();
  Eigen::Vector3d p_ip1 = T_ip1.translation();
  Eigen::Vector3d p_ip2 = T_ip2.translation();

  Eigen::Matrix3d R_im1 = T_im1.rotation().matrix();
  Eigen::Matrix3d R_i = T_i.rotation().matrix();
  Eigen::Matrix3d R_ip1 = T_ip1.rotation().matrix();
  Eigen::Matrix3d R_ip2 = T_ip2.rotation().matrix();

  Eigen::Vector4d B = (1 / 6) * C * U.block<4, 1>(0, 0);
  Eigen::Vector4d dB = (1 / 6) * C * U.block<4, 1>(0, 1);
  Eigen::Vector4d ddB = (1 / 6) * C * U.block<4, 1>(0, 2);
}

Eigen::Matrix<double, 4, 3> calculateWeights(double u) {
  Eigen::Matrix<double, 4, 3> U;
  Eigen::Vector4d b, db, ddb;

  b << 1, u, u * u, u * u * u;
  db << 0, 1, 2 * u, 3 * u * u;
  ddb << 0, 0, 2, 6 * u;

  U.block<4, 1>(0, 0) = b;
  U.block<4, 1>(0, 1) = db;
  U.block<4, 1>(0, 2) = ddb;

  return U;
}
}  // namespace CubicBasisSpline