#include "spline_tools/spline_tools.hpp"

#include <iostream>

#include "spline_tools/so3_tools.hpp"

namespace CubicBasisSplines {

double calculateNormalizedTime(const double &t_i, const double &t_eval, const double &dt) {
  return (t_eval - t_i) / dt;
}

Weights calculateWeights(const double &u) {
  Eigen::Vector4d U, dU, ddU;

  Weights calc_weights;

  U << 1, u, std::pow(u, 2), std::pow(u, 3);
  dU << 0, 1.0, 2 * u, 3 * std::pow(u, 2);
  ddU << 0.0, 0.0, 2.0, 6 * u;

  calc_weights.Btilde = (1.0 / 6.0) * C * U;
  calc_weights.dBtilde = (1.0 / 6.0) * C * dU;
  calc_weights.ddBtilde = (1.0 / 6.0) * C * ddU;

  return calc_weights;
}

bool evaluateTranslationSpline(const std::vector<CtrlPt> &ctrl_pts, double &t_eval,
                               SplineEval &eval_spline) {
  int num_pts = ctrl_pts.size();

  if (num_pts != 4) {
    return false;
  } else {
    CtrlPt P_im1, P, P_ip1, P_ip2;
    Weights current_weights;

    P_im1 = ctrl_pts[0];
    P = ctrl_pts[1];
    P_ip1 = ctrl_pts[2];
    P_ip2 = ctrl_pts[3];

    double t_i = P.t;
    double t_ip1 = P_ip1.t;
    double dt = t_ip1 - t_i;

    if (dt <= 0.0) {
      return false;
    }

    double u = calculateNormalizedTime(t_i, t_eval, dt);
    current_weights = calculateWeights(u);

    Eigen::Vector3d p_im1, p, p_ip1, p_ip2;
    Eigen::Vector3d dp_i, dp_ip1, dp_ip2;

    p_im1 = P_im1.p;
    p = P.p;
    p_ip1 = P_ip1.p;
    p_ip2 = P_ip2.p;

    dp_i = p - p_im1;
    dp_ip1 = p_ip1 - p;
    dp_ip2 = p_ip2 - p_ip1;

    eval_spline.pos = p_im1 + current_weights.Btilde[1] * dp_i +
                      current_weights.Btilde[2] * dp_ip1 + current_weights.Btilde[3] * dp_ip2;

    eval_spline.vel =
        (1.0 / dt) * (current_weights.dBtilde[1] * dp_i + current_weights.dBtilde[2] * dp_ip1 +
                      current_weights.dBtilde[3] * dp_ip2);

    eval_spline.acc = (1.0 / std::pow(dt, 2)) *
                      (current_weights.ddBtilde[1] * dp_i + current_weights.ddBtilde[2] * dp_ip1 +
                       current_weights.ddBtilde[3] * dp_ip2);

    return true;
  }
}

bool evaluateRotationSpline(const std::vector<CtrlPt> &ctrl_pts, double &t_eval,
                            SplineEval &eval_spline) {
  int num_pts = ctrl_pts.size();

  if (num_pts != 4) {
    return false;
  } else {
    CtrlPt P_im1, P, P_ip1, P_ip2;
    Weights current_weights;

    P_im1 = ctrl_pts[0];
    P = ctrl_pts[1];
    P_ip1 = ctrl_pts[2];
    P_ip2 = ctrl_pts[3];

    double t_i = P.t;
    double t_ip1 = P_ip1.t;
    double dt = t_ip1 - t_i;

    if (dt <= 0.0) {
      return false;
    }

    double u = calculateNormalizedTime(t_i, t_eval, dt);
    current_weights = calculateWeights(u);

    Eigen::Matrix3d R_im1, R, R_ip1, R_ip2;
    Eigen::Matrix3d delR_i, delR_ip1, delR_ip2;
    Eigen::Vector3d Omega_i, Omega_ip1, Omega_ip2;
    Eigen::Matrix3d A1, A2, A3;

    R_im1 = P_im1.R;
    R = P.R;
    R_ip1 = P_ip1.R;
    R_ip2 = P_ip2.R;

    delR_i = (R_im1.transpose()) * R;
    delR_ip1 = (R.transpose()) * R_ip1;
    delR_ip2 = (R_ip1.transpose()) * R_ip2;

    Omega_i = logm(delR_i);
    Omega_ip1 = logm(delR_ip1);
    Omega_ip2 = logm(delR_ip2);

    A1 = expm(current_weights.Btilde[1] * Omega_i);
    A2 = expm(current_weights.Btilde[2] * Omega_ip1);
    A3 = expm(current_weights.Btilde[3] * Omega_ip2);

    eval_spline.C_bw = R_im1 * A1 * A2 * A3;
    eval_spline.w_bw = (current_weights.dBtilde[1] * ((A2 * A3).transpose()) * Omega_i +
                        current_weights.dBtilde[2] * A3.transpose() * Omega_ip1 +
                        current_weights.dBtilde[3] * Omega_ip2) *
                       (1.0 / dt);

    return true;
  }
}

bool evaluateSplitSpline(const std::vector<CtrlPt> &ctrl_pts, double &t_eval,
                         SplineEval &eval_spline) {
  try {
    evaluateTranslationSpline(ctrl_pts, t_eval, eval_spline);
    evaluateRotationSpline(ctrl_pts, t_eval, eval_spline);

    return true;
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return false;
  }
}
}  // namespace CubicBasisSplines