#ifndef SPLINE_EVAL_HPP
#define SPLINE_EVAL_HPP

#include <sophus/so3.hpp>

#include "spline_tools/spline_basis.hpp"

namespace CubicBasisSpline {

template <typename T>
struct SO3SplineEval {
  Sophus::SO3<T> C_bw;
  Eigen::Matrix<T, 3, 1> w_bw;
};

template <typename T>
struct R3SplineEval {
  Eigen::Matrix<T, 3, 1> p;  // position
  Eigen::Matrix<T, 3, 1> v;  // velocity
  Eigen::Matrix<T, 3, 1> a;  // acceleration
};

/**
 @brief Evaluates continuous SO3 Spline for Orientation and Angular Velocity
 @param R_im1 Orientation Control Point at i minus 1
 @param R_i Orientation Control Point at i
 @param R_ip1 Orientation Control Point at i plus 1
 @param R_ip2 Orientation Control Point at i plus 2
 @param u Normalized Time [0,1]
 @param dt Time spacing between control points
*/
template <typename T>
SO3SplineEval<T> evaluateSO3(const Sophus::SO3<T>& R_im1, const Sophus::SO3<T>& R_i,
                             const Sophus::SO3<T>& R_ip1, const Sophus::SO3<T>& R_ip2, double u,
                             double dt) {
  SO3SplineEval<T> spline_eval;

  SplineWeights w = calculateSplineWeights(u);

  Sophus::SO3<T> dR1 = R_im1.inverse() * R_i;
  Sophus::SO3<T> dR2 = R_i.inverse() * R_ip1;
  Sophus::SO3<T> dR3 = R_ip1.inverse() * R_ip2;

  Eigen::Matrix<T, 3, 1> Omega1 = dR1.log();
  Eigen::Matrix<T, 3, 1> Omega2 = dR2.log();
  Eigen::Matrix<T, 3, 1> Omega3 = dR3.log();

  // calculating our orientation
  Sophus::SO3<T> A1 = Sophus::SO3<T>::exp(T(w.Btilde[1]) * Omega1);
  Sophus::SO3<T> A2 = Sophus::SO3<T>::exp(T(w.Btilde[2]) * Omega2);
  Sophus::SO3<T> A3 = Sophus::SO3<T>::exp(T(w.Btilde[3]) * Omega3);

  spline_eval.C_bw = R_im1 * A1 * A2 * A3;

  // calculating our angular velocity

  T inv_dt = T(1.0 / dt);

  Eigen::Matrix<T, 3, 3> R_A2_A3_T = (A2 * A3).matrix().transpose();
  Eigen::Matrix<T, 3, 3> R_A3_T = A3.matrix().transpose();

  spline_eval.w_bw = inv_dt * (T(w.dBtilde[1]) * (R_A2_A3_T * Omega1) +
                               T(w.dBtilde[2]) * (R_A3_T * Omega2) + T(w.dBtilde[3]) * Omega3);

  return spline_eval;
}

/**
 @brief Evaluates continuous R3 Spline for Position, Velocity, and Acceleration
 @param p_im1 Translation Control Point at i minus 1
 @param p_i   Translation Control Point at i
 @param p_ip1 Translation Control Point at i plus 1
 @param p_ip2 Translation Control Point at i plus 2
 @param u     Normalized Time [0,1]
 @param dt    Time spacing between control points
*/
template <typename T>
R3SplineEval<T> evaluateR3(const Eigen::Matrix<T, 3, 1>& p_im1, const Eigen::Matrix<T, 3, 1>& p_i,
                           const Eigen::Matrix<T, 3, 1>& p_ip1, const Eigen::Matrix<T, 3, 1>& p_ip2,
                           double u, double dt) {
  R3SplineEval<T> spline_eval;

  SplineWeights w = calculateSplineWeights(u);

  Eigen::Matrix<T, 3, 1> dp1 = p_i - p_im1;
  Eigen::Matrix<T, 3, 1> dp2 = p_ip1 - p_i;
  Eigen::Matrix<T, 3, 1> dp3 = p_ip2 - p_ip1;

  T inv_dt = T(1.0 / dt);
  T inv_dt2 = inv_dt * inv_dt;

  spline_eval.p = p_im1 + T(w.Btilde[1]) * dp1 + T(w.Btilde[2]) * dp2 + T(w.Btilde[3]) * dp3;
  spline_eval.v = inv_dt * (T(w.dBtilde[1]) * dp1 + T(w.dBtilde[2]) * dp2 + T(w.dBtilde[3]) * dp3);
  spline_eval.a =
      inv_dt2 * (T(w.ddBtilde[1]) * dp1 + T(w.ddBtilde[2]) * dp2 + T(w.ddBtilde[3]) * dp3);

  return spline_eval;
}
}  // namespace CubicBasisSpline

#endif