#include "spline_tools/cubic_spline_factors.hpp"

#include "spline_tools/nav_utils.hpp"

namespace CubicBasisSpline {

namespace {

// Must match imuFactor():  f_hat = a + kGravSign*g(p) + kCorSign*2*Omega*v
//                          r_w   = w + kEarthSign*R^T w_e + b_g
// Groves-consistent values shown. If your imuFactor still has the original signs
// use kGravSign=+1, kCorSign=-1, kEarthSign=-1.
constexpr double kGravSign = -1.0, kCorSign = +1.0, kEarthSign = +1.0;
}  // namespace

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

    Eigen::Vector4d B = (1.0 / 6.0) * C * U.block<4, 1>(0, 0);
    Eigen::Vector4d dB = (1.0 / 6.0) * C * U.block<4, 1>(0, 1);
    Eigen::Vector4d ddB = (1.0 / 6.0) * C * U.block<4, 1>(0, 2);

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

bool evaluateBiasSpline(FullBiasSet& ctrl_pts, BiasEvaluation& spline_eval) {
  try {
    gtsam::Vector3 w_im1 = ctrl_pts.knots_b_g.P_im1;
    gtsam::Vector3 w_i = ctrl_pts.knots_b_g.P_i;
    gtsam::Vector3 w_ip1 = ctrl_pts.knots_b_g.P_ip1;
    gtsam::Vector3 w_ip2 = ctrl_pts.knots_b_g.P_ip2;

    gtsam::Vector3 f_im1 = ctrl_pts.knots_b_a.P_im1;
    gtsam::Vector3 f_i = ctrl_pts.knots_b_a.P_i;
    gtsam::Vector3 f_ip1 = ctrl_pts.knots_b_a.P_ip1;
    gtsam::Vector3 f_ip2 = ctrl_pts.knots_b_a.P_ip2;

    Eigen::Matrix<double, 4, 3> U = calculateWeights(ctrl_pts.u);

    Eigen::Vector4d B = (1.0 / 6.0) * C * U.block<4, 1>(0, 0);

    // gyro bias
    gtsam::Point3 dw1 = B[1] * (w_i - w_im1);
    gtsam::Point3 dw2 = B[2] * (w_ip1 - w_i);
    gtsam::Point3 dw3 = B[3] * (w_ip2 - w_ip1);

    spline_eval.b_g = w_im1 + dw1 + dw2 + dw3;

    gtsam::Point3 df1 = B[1] * (f_i - f_im1);
    gtsam::Point3 df2 = B[2] * (f_ip1 - f_i);
    gtsam::Point3 df3 = B[3] * (f_ip2 - f_ip1);

    spline_eval.b_a = f_im1 + df1 + df2 + df3;

    return true;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    return false;
  }
}

bool gnssFactor(SplineKnotSet& ctrl_pts, const gtsam::Point3& pos, const gtsam::Vector3& vel,
                const gtsam::Vector3& l_a, gtsam::Point3& residual_pos,
                gtsam::Vector3& residual_vel) {
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

bool imuFactor(SplineKnotSet& ctrl_pts, FullBiasSet& bias_pts, const gtsam::Vector3& w_ib_b,
               const gtsam::Vector3& f_ib_b, gtsam::Vector3& residual_w,
               gtsam::Vector3& residual_f) {
  try {
    SplineEvaluation X;
    BiasEvaluation B;

    if (!evaluateCubicBSpline(ctrl_pts, X)) return false;
    if (!evaluateBiasSpline(bias_pts, B)) return false;

    const Eigen::Vector3d ecef_gravity = CubicBasisSpline::ecefSomiglianaGravity(X.pos);
    const Eigen::Matrix3d C_ie_e = gtsam::skewSymmetric(CubicBasisSpline::w_e_);

    // specific force in ECEF, then rotate to body
    const gtsam::Vector3 f_hat =
        X.lin_acc + kGravSign * ecef_gravity + kCorSign * 2.0 * C_ie_e * X.lin_vel;
    residual_f = X.orientation.transpose() * f_hat + B.b_a - f_ib_b;

    // angular rate
    residual_w = X.ang_vel + kEarthSign * (X.orientation.transpose() * CubicBasisSpline::w_e_) +
                 B.b_g - w_ib_b;

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

namespace {

struct SplineJacobians {
  Eigen::Vector4d c, d, e;    // position / velocity / acceleration weights (dT included)
  Eigen::Matrix3d Dphi[4];    // d(orientation, right perturbation) / d(knot rotation tangent)
  Eigen::Matrix3d Domega[4];  // d(angular velocity) / d(knot rotation tangent)
};

void splineJacobians(const SplineKnotSet& K, SplineJacobians& J) {
  const Eigen::Matrix<double, 4, 3> U = calculateWeights(K.u);
  const Eigen::Vector4d B = (1.0 / 6.0) * C * U.col(0);
  const Eigen::Vector4d dB = (1.0 / 6.0) * C * U.col(1);
  const Eigen::Vector4d ddB = (1.0 / 6.0) * C * U.col(2);
  const double iT = 1.0 / K.dT;

  auto weights = [](const Eigen::Vector4d& b) {
    Eigen::Vector4d w;
    w << b[0] - b[1], b[1] - b[2], b[2] - b[3], b[3];
    return w;
  };
  J.c = weights(B);
  J.d = weights(dB) * iT;
  J.e = weights(ddB) * iT * iT;

  const Eigen::Matrix3d R[4] = {K.T_im1.rotation().matrix(), K.T_i.rotation().matrix(),
                                K.T_ip1.rotation().matrix(), K.T_ip2.rotation().matrix()};

  Eigen::Matrix3d dR[4], A[4], Jinv[4], Gam[4];
  Eigen::Vector3d Om[4];
  for (int j = 1; j <= 3; j++) {
    dR[j] = R[j - 1].transpose() * R[j];
    Om[j] = gtsam::Rot3::Logmap(gtsam::Rot3(dR[j]));
    A[j] = gtsam::Rot3::Expmap(B[j] * Om[j]).matrix();
    Jinv[j] = gtsam::Rot3::LogmapDerivative(Om[j]);               // Jr^{-1}
    Gam[j] = B[j] * gtsam::Rot3::ExpmapDerivative(B[j] * Om[j]);  // B * Jr
  }

  const Eigen::Matrix3d M0 = (A[1] * A[2] * A[3]).transpose();
  const Eigen::Matrix3d M1 = (A[2] * A[3]).transpose();
  const Eigen::Matrix3d M2 = A[3].transpose();

  // orientation
  const Eigen::Matrix3d K1 = M1 * Gam[1] * Jinv[1];
  const Eigen::Matrix3d K2 = M2 * Gam[2] * Jinv[2];
  const Eigen::Matrix3d K3 = Gam[3] * Jinv[3];
  J.Dphi[0] = M0 - K1 * dR[1].transpose();
  J.Dphi[1] = K1 - K2 * dR[2].transpose();
  J.Dphi[2] = K2 - K3 * dR[3].transpose();
  J.Dphi[3] = K3;

  // angular velocity
  const double s1 = dB[1] * iT, s2 = dB[2] * iT, s3 = dB[3] * iT;
  const Eigen::Matrix3d V1 = gtsam::skewSymmetric(Eigen::Vector3d(M1 * Om[1]));
  const Eigen::Matrix3d V2 = gtsam::skewSymmetric(Eigen::Vector3d(M2 * Om[2]));
  const Eigen::Matrix3d H1 = (s1 * M1) * Jinv[1];
  const Eigen::Matrix3d H2 = (s2 * M2 + s1 * V1 * M2 * Gam[2]) * Jinv[2];
  const Eigen::Matrix3d H3 =
      (s3 * Eigen::Matrix3d::Identity() + (s1 * V1 + s2 * V2) * Gam[3]) * Jinv[3];
  J.Domega[0] = -H1 * dR[1].transpose();
  J.Domega[1] = H1 - H2 * dR[2].transpose();
  J.Domega[2] = H2 - H3 * dR[3].transpose();
  J.Domega[3] = H3;
}

Eigen::Matrix3d gravityJacobian(const Eigen::Vector3d& p) {
  Eigen::Matrix3d G;
  const double h = 1.0;  // meters
  for (int i = 0; i < 3; i++) {
    const Eigen::Vector3d dp = Eigen::Vector3d::Unit(i) * h;
    G.col(i) = (ecefSomiglianaGravity(p + dp) - ecefSomiglianaGravity(p - dp)) / (2.0 * h);
  }
  return G;
}

std::array<Eigen::Matrix3d, 4> knotRotations(const SplineKnotSet& K) {
  return {K.T_im1.rotation().matrix(), K.T_i.rotation().matrix(), K.T_ip1.rotation().matrix(),
          K.T_ip2.rotation().matrix()};
}

}  // namespace

void gnssJacobians(const SplineKnotSet& knots, const gtsam::Vector3& l_a,
                   std::array<gtsam::Matrix6, 4>& H) {
  SplineKnotSet k = knots;
  SplineEvaluation X;
  evaluateCubicBSpline(k, X);
  SplineJacobians J;
  splineJacobians(k, J);

  const Eigen::Matrix3d Rs = X.orientation.matrix();
  const Eigen::Matrix3d Lx = gtsam::skewSymmetric(l_a);
  const Eigen::Matrix3d WLx = gtsam::skewSymmetric(Eigen::Vector3d(X.ang_vel.cross(l_a)));
  const auto Rk = knotRotations(k);

  for (int i = 0; i < 4; i++) {
    H[i].block<3, 3>(0, 0) = -Rs * Lx * J.Dphi[i];
    H[i].block<3, 3>(0, 3) = J.c[i] * Rk[i];
    H[i].block<3, 3>(3, 0) = -Rs * WLx * J.Dphi[i] - Rs * Lx * J.Domega[i];
    H[i].block<3, 3>(3, 3) = J.d[i] * Rk[i];
  }
}

void imuPoseJacobians(const SplineKnotSet& knots, std::array<gtsam::Matrix6, 4>& H) {
  SplineKnotSet k = knots;
  SplineEvaluation X;
  evaluateCubicBSpline(k, X);
  SplineJacobians J;
  splineJacobians(k, J);

  const Eigen::Matrix3d Rs = X.orientation.matrix();
  const Eigen::Vector3d we = CubicBasisSpline::w_e_;
  const Eigen::Matrix3d Om = gtsam::skewSymmetric(we);
  const Eigen::Matrix3d G = gravityJacobian(X.pos);
  const Eigen::Vector3d g = ecefSomiglianaGravity(X.pos);

  const Eigen::Vector3d f_hat = X.lin_acc + kGravSign * g + kCorSign * 2.0 * Om * X.lin_vel;
  const Eigen::Matrix3d Fx = gtsam::skewSymmetric(Eigen::Vector3d(Rs.transpose() * f_hat));
  const Eigen::Matrix3d Wx = gtsam::skewSymmetric(Eigen::Vector3d(Rs.transpose() * we));
  const auto Rk = knotRotations(k);

  for (int i = 0; i < 4; i++) {
    const Eigen::Matrix3d dpos = J.e[i] * Eigen::Matrix3d::Identity() + kGravSign * J.c[i] * G +
                                 2.0 * kCorSign * J.d[i] * Om;
    H[i].block<3, 3>(0, 0) = Fx * J.Dphi[i];
    H[i].block<3, 3>(0, 3) = Rs.transpose() * dpos * Rk[i];
    H[i].block<3, 3>(3, 0) = J.Domega[i] + kEarthSign * Wx * J.Dphi[i];
    H[i].block<3, 3>(3, 3).setZero();
  }
}

}  // namespace CubicBasisSpline