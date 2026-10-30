#include "spline_tools/cubic_spline_factors.hpp"

namespace CubicBasisSpline {

bool evaluateCubicBSpline(SplineKnotSet& ctrl_pts, SplineEvaluation& spline_eval) {
  try {
    // parsing out the current set of ctrl points
    Eigen::Vector3d p_im1 = ctrl_pts.T_im1.translation();
    Eigen::Vector3d p_i = ctrl_pts.T_i.translation();
    Eigen::Vector3d p_ip1 = ctrl_pts.T_ip1.translation();
    Eigen::Vector3d p_ip2 = ctrl_pts.T_ip2.translation();

    Eigen::Matrix3d R_im1 = ctrl_pts.T_im1.rotation().matrix();
    Eigen::Matrix3d R_i = ctrl_pts.T_i.rotation().matrix();
    Eigen::Matrix3d R_ip1 = ctrl_pts.T_ip1.rotation().matrix();
    Eigen::Matrix3d R_ip2 = ctrl_pts.T_ip2.rotation().matrix();

    gtsam::Rot3 R_1(R_im1);
    gtsam::Rot3 R_rel_1(R_im1.inverse() * R_i);
    gtsam::Rot3 R_rel_2(R_i.inverse() * R_ip1);
    gtsam::Rot3 R_rel_3(R_ip1.inverse() * R_ip2);

    // Making Weight Vectors
    Eigen::Matrix<double, 4, 3> U = calculateWeights(ctrl_pts.u);

    Eigen::Vector4d B = (1 / 6) * C * U.block<4, 1>(0, 0);
    Eigen::Vector4d dB = (1 / 6) * C * U.block<4, 1>(0, 1);
    Eigen::Vector4d ddB = (1 / 6) * C * U.block<4, 1>(0, 2);

    // ===== R3 Spline Evaluation =====

    // position
    gtsam::Point3 dp1 = B[1] * (p_i - p_im1);
    gtsam::Point3 dp2 = B[2] * (p_ip1 - p_i);
    gtsam::Point3 dp3 = B[3] * (p_ip2 - p_ip1);

    spline_eval.pos = p_im1 + dp1 + dp2 + dp3;

    // Velocity
    gtsam::Vector3 dv1 = dB[1] * (p_i - p_im1);
    gtsam::Vector3 dv2 = dB[2] * (p_ip1 - p_i);
    gtsam::Vector3 dv3 = dB[3] * (p_ip2 - p_ip1);

    spline_eval.lin_vel = (1 / ctrl_pts.dT) * (dv1 + dv2 + dv3);

    // acceleration
    gtsam::Vector3 da1 = ddB[1] * (p_i - p_im1);
    gtsam::Vector3 da2 = ddB[2] * (p_ip1 - p_i);
    gtsam::Vector3 da3 = ddB[3] * (p_ip2 - p_ip1);

    spline_eval.lin_acc = (1 / std::pow(ctrl_pts.dT, 2)) * (da1 + da2 + da3);

    // ===== SO(3) Spline Evaluation =====

    // Orientation
    gtsam::Vector3 Omega1 = gtsam::Rot3::Logmap(R_rel_1);
    gtsam::Vector3 Omega2 = gtsam::Rot3::Logmap(R_rel_2);
    gtsam::Vector3 Omega3 = gtsam::Rot3::Logmap(R_rel_3);

    gtsam::Rot3 A1 = gtsam::Rot3::Expmap(B[1] * Omega1);
    gtsam::Rot3 A2 = gtsam::Rot3::Expmap(B[2] * Omega2);
    gtsam::Rot3 A3 = gtsam::Rot3::Expmap(B[3] * Omega3);

    spline_eval.orientation = R_1 * A1 * A2 * A3;

    // Angular Velocity
    gtsam::Vector3 term1 = dB[1] * (A2 * A3).transpose() * Omega1;
    gtsam::Vector3 term2 = dB[2] * A3.transpose() * Omega2;
    gtsam::Vector3 term3 = dB[3] * Omega3;

    spline_eval.ang_vel = (1 / ctrl_pts.dT) * (term1 + term2 + term3);

    return true;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    return false;
  }
}

bool gnssFactor(SplineKnotSet& ctrl_pts, gtsam::Point3& pos, gtsam::Vector3& vel,
                gtsam::Vector3& l_a, gtsam::Point3& residual_pos, gtsam::Vector3& residual_vel) {
  try {
    SplineEvaluation X;
    evaluateCubicBSpline(ctrl_pts, X);

    // position residual
    gtsam::Point3 pos_hat = X.pos + X.orientation * l_a;
    residual_pos = pos_hat - pos;

    // velocity residual
    gtsam::Vector3 vel_hat = X.lin_vel + X.orientation * (X.ang_vel.cross(l_a));
    residual_vel = vel_hat - vel;

    return true;

  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    return false;
  }
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