#ifndef SPLINE_FACTORS_HPP
#define SPLINE_FACTORS_HPP

#include <ceres/ceres.h>

#include <Eigen/Core>
#include <sophus/so3.hpp>

#include "spline_tools/spline_eval.hpp"

namespace CubicBasisSpline {

/*
p = r + C * l
*/

struct GnssPositionFactor {
  GnssPositionFactor(const Eigen::Vector3d& z_gnss, const Eigen::Vector3d& l_b,
                     const Eigen::Vector3d& sigmas, double u, double dt)
      : z_gnss_(z_gnss), l_b_(l_b), sigmas_(sigmas), u_(u), dt_(dt) {}

  template <typename T>
  bool operator()(const T* const p_im1, const T* const p_i, const T* const p_ip1,
                  const T* const p_ip2, const T* const R_im1, const T* const R_i,
                  const T* const R_ip1, const T* const R_ip2, T* residuals) const {
    Eigen::Matrix<T, 3, 1> p0 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_im1);
    Eigen::Matrix<T, 3, 1> p1 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_i);
    Eigen::Matrix<T, 3, 1> p2 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_ip1);
    Eigen::Matrix<T, 3, 1> p3 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_ip2);

    Sophus::SO3<T> R0((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_im1)));
    Sophus::SO3<T> R1((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_i)));
    Sophus::SO3<T> R2((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_ip1)));
    Sophus::SO3<T> R3((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_ip2)));

    R3SplineEval<T> r3 = evaluateR3(p0, p1, p2, p3, u_, dt_);
    SO3SplineEval<T> so3 = evaluateSO3(R0, R1, R2, R3, u_, dt_);

    Eigen::Matrix<T, 3, 1> p_antenna = r3.p + (so3.C_bw * l_b_.template cast<T>());

    Eigen::Map<Eigen::Matrix<T, 3, 1>> residual(residuals);
    Eigen::Matrix<T, 3, 1> inv_sigmas = sigmas_.cwiseInverse().template cast<T>();

    residual = inv_sigmas.cwiseProduct(p_antenna - z_gnss_.template cast<T>());

    return true;
  }

 private:
  Eigen::Vector3d z_gnss_;
  Eigen::Vector3d l_b_;
  Eigen::Vector3d sigmas_;
  double u_, dt_;
};

/*

p_dot = v + C(w.cross(l))

*/
struct GnssVelocityFactor {
  GnssVelocityFactor(const Eigen::Vector3d& z_gnss, const Eigen::Vector3d& l_b,
                     const Eigen::Vector3d& sigmas, double u, double dt)
      : z_gnss_(z_gnss), l_b_(l_b), sigmas_(sigmas), u_(u), dt_(dt) {}

  template <typename T>
  bool operator()(const T* const p_im1, const T* const p_i, const T* const p_ip1,
                  const T* const p_ip2, const T* const R_im1, const T* const R_i,
                  const T* const R_ip1, const T* const R_ip2, T* residuals) const {
    Eigen::Matrix<T, 3, 1> p0 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_im1);
    Eigen::Matrix<T, 3, 1> p1 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_i);
    Eigen::Matrix<T, 3, 1> p2 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_ip1);
    Eigen::Matrix<T, 3, 1> p3 = Eigen::Map<const Eigen::Matrix<T, 3, 1>>(p_ip2);

    Sophus::SO3<T> R0((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_im1)));
    Sophus::SO3<T> R1((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_i)));
    Sophus::SO3<T> R2((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_ip1)));
    Sophus::SO3<T> R3((Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>>(R_ip2)));

    R3SplineEval<T> r3 = evaluateR3(p0, p1, p2, p3, u_, dt_);
    SO3SplineEval<T> so3 = evaluateSO3(R0, R1, R2, R3, u_, dt_);

    Eigen::Matrix<T, 3, 1> v_antenna = r3.v + so3.C_bw * (so3.w_bw.cross(l_b_.template cast<T>()));

    Eigen::Map<Eigen::Matrix<T, 3, 1>> residual(residuals);
    Eigen::Matrix<T, 3, 1> inv_sigmas = sigmas_.cwiseInverse().template cast<T>();

    residual = inv_sigmas.cwiseProduct(v_antenna - z_gnss_.template cast<T>());

    return true;
  }

 private:
  Eigen::Vector3d z_gnss_;
  Eigen::Vector3d l_b_;
  Eigen::Vector3d sigmas_;
  double u_, dt_;
};

}  // namespace CubicBasisSpline
#endif