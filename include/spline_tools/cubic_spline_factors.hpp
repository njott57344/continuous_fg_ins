#ifndef CUBIC_SPLINE_FACTORS_HPP
#define CUBIC_SPLINE_FACTORS_HPP

#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "gtsam/base/numericalDerivative.h"
#include "gtsam/geometry/Point3.h"
#include "gtsam/geometry/Pose3.h"
#include "gtsam/geometry/Rot3.h"
#include "gtsam/geometry/SO3.h"
#include "gtsam/nonlinear/NonlinearFactor.h"
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

Eigen::Matrix<double, 4, 3> calculateWeights(double u);

struct SplineKnotSet {
  gtsam::Pose3 T_im1;
  gtsam::Pose3 T_i;
  gtsam::Pose3 T_ip1;
  gtsam::Pose3 T_ip2;
  double u;
  double dT;
};

struct BiasEvaluation {
  gtsam::Vector3 b_g = gtsam::Vector3::Zero();
  gtsam::Vector3 b_a = gtsam::Vector3::Zero();
};

struct BiasKnotSet {
  gtsam::Vector3 P_im1;
  gtsam::Vector3 P_i;
  gtsam::Vector3 P_ip1;
  gtsam::Vector3 P_ip2;
};

struct FullBiasSet {
  BiasKnotSet knots_b_g;
  BiasKnotSet knots_b_a;
  double u;
};

bool gnssFactor(SplineKnotSet& ctrl_pts, const gtsam::Point3& pos, const gtsam::Vector3& vel,
                const gtsam::Vector3& l_a, gtsam::Point3& residual_pos,
                gtsam::Vector3& residual_vel);

bool imuFactor(SplineKnotSet& ctrl_pts, FullBiasSet& bias_pts, const gtsam::Vector3& w_ib_b,
               const gtsam::Vector3& f_ib_b, gtsam::Vector3& residual_w,
               gtsam::Vector3& residual_f);

bool evaluateCubicBSpline(SplineKnotSet& ctrl_pts, SplineEvaluation& spline_eval);
bool evaluateBiasSpline(FullBiasSet& ctrl_pts, BiasEvaluation& spline_eval);

class CubicGnssFactor
    : public gtsam::NoiseModelFactorN<gtsam::Pose3, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3> {
 public:
  using Base = gtsam::NoiseModelFactorN<gtsam::Pose3, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3>;

  CubicGnssFactor(gtsam::Key key_im1, gtsam::Key key_i, gtsam::Key key_ip1, gtsam::Key key_ip2,
                  const gtsam::Point3& z_pos, const gtsam::Vector3& z_vel, const gtsam::Point3& l_a,
                  double u, double dT, const gtsam::SharedNoiseModel& noise_model)
      : Base(noise_model, key_im1, key_i, key_ip1, key_ip2),
        z_pos_(z_pos),
        z_vel_(z_vel),
        l_a_(l_a),
        u_(u),
        dT_(dT) {}

  gtsam::Vector evaluateError(const gtsam::Pose3& T_im1, const gtsam::Pose3& T_i,
                              const gtsam::Pose3& T_ip1, const gtsam::Pose3& T_ip2,
                              gtsam::OptionalMatrixType H1 = nullptr,
                              gtsam::OptionalMatrixType H2 = nullptr,
                              gtsam::OptionalMatrixType H3 = nullptr,
                              gtsam::OptionalMatrixType H4 = nullptr) const override {
    auto residualFn = [this](const gtsam::Pose3& a, const gtsam::Pose3& b, const gtsam::Pose3& c,
                             const gtsam::Pose3& d) -> gtsam::Vector6 {
      SplineKnotSet knots{a, b, c, d, u_, dT_};
      gtsam::Point3 rp;
      gtsam::Vector3 rv;
      gnssFactor(knots, z_pos_, z_vel_, l_a_, rp, rv);
      gtsam::Vector6 res;
      res.head<3>() = rp;
      res.tail<3>() = rv;
      return res;
    };

    if (H1) {
      *H1 = gtsam::numericalDerivative41<gtsam::Vector6, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3,
                                         gtsam::Pose3>(residualFn, T_im1, T_i, T_ip1, T_ip2, 1e-5);
    }
    if (H2) {
      *H2 = gtsam::numericalDerivative42<gtsam::Vector6, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3,
                                         gtsam::Pose3>(residualFn, T_im1, T_i, T_ip1, T_ip2, 1e-5);
    }
    if (H3) {
      *H3 = gtsam::numericalDerivative43<gtsam::Vector6, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3,
                                         gtsam::Pose3>(residualFn, T_im1, T_i, T_ip1, T_ip2, 1e-5);
    }
    if (H4) {
      *H4 = gtsam::numericalDerivative44<gtsam::Vector6, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3,
                                         gtsam::Pose3>(residualFn, T_im1, T_i, T_ip1, T_ip2, 1e-5);
    }

    return residualFn(T_im1, T_i, T_ip1, T_ip2);
  }

 private:
  gtsam::Point3 z_pos_;
  gtsam::Vector3 z_vel_;
  gtsam::Point3 l_a_;
  double u_;
  double dT_;
};

class CubicImuFactor : public gtsam::NoiseModelFactorN<gtsam::Pose3, gtsam::Pose3, gtsam::Pose3,
                                                       gtsam::Pose3, gtsam::Vector6, gtsam::Vector6,
                                                       gtsam::Vector6, gtsam::Vector6> {
 public:
  using Base =
      gtsam::NoiseModelFactorN<gtsam::Pose3, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3,
                               gtsam::Vector6, gtsam::Vector6, gtsam::Vector6, gtsam::Vector6>;

  CubicImuFactor(gtsam::Key p_im1, gtsam::Key p_i, gtsam::Key p_ip1, gtsam::Key p_ip2,
                 gtsam::Key b_im1, gtsam::Key b_i, gtsam::Key b_ip1, gtsam::Key b_ip2,
                 const gtsam::Vector3& z_w, const gtsam::Vector3& z_f, double u_p, double dT_p,
                 double u_b, const gtsam::SharedNoiseModel& noise_model)
      : Base(noise_model, p_im1, p_i, p_ip1, p_ip2, b_im1, b_i, b_ip1, b_ip2),
        z_w_(z_w),
        z_f_(z_f),
        u_p_(u_p),
        dT_p_(dT_p),
        u_b_(u_b) {}

  gtsam::Vector evaluateError(
      const gtsam::Pose3& T0, const gtsam::Pose3& T1, const gtsam::Pose3& T2,
      const gtsam::Pose3& T3, const gtsam::Vector6& b0, const gtsam::Vector6& b1,
      const gtsam::Vector6& b2, const gtsam::Vector6& b3, gtsam::OptionalMatrixType H1 = nullptr,
      gtsam::OptionalMatrixType H2 = nullptr, gtsam::OptionalMatrixType H3 = nullptr,
      gtsam::OptionalMatrixType H4 = nullptr, gtsam::OptionalMatrixType H5 = nullptr,
      gtsam::OptionalMatrixType H6 = nullptr, gtsam::OptionalMatrixType H7 = nullptr,
      gtsam::OptionalMatrixType H8 = nullptr) const override {
    // bias knots -> FullBiasSet (b = [b_a; b_g])
    FullBiasSet bias_set;
    bias_set.knots_b_a = {b0.head<3>(), b1.head<3>(), b2.head<3>(), b3.head<3>()};
    bias_set.knots_b_g = {b0.tail<3>(), b1.tail<3>(), b2.tail<3>(), b3.tail<3>()};
    bias_set.u = u_b_;  // normalized time on the BIAS grid

    auto residualFn = [&](const gtsam::Pose3& a, const gtsam::Pose3& b, const gtsam::Pose3& c,
                          const gtsam::Pose3& d) -> gtsam::Vector6 {
      SplineKnotSet knots{a, b, c, d, u_p_, dT_p_};  // u on the POSE grid
      gtsam::Vector3 rw = gtsam::Vector3::Zero(), rf = gtsam::Vector3::Zero();
      FullBiasSet bs = bias_set;  // imuFactor takes non-const ref
      if (!imuFactor(knots, bs, z_w_, z_f_, rw, rf)) throw std::runtime_error("imuFactor failed");
      gtsam::Vector6 res;
      res << rf, rw;
      return res;
    };

    using P = gtsam::Pose3;
    using V = gtsam::Vector6;
    if (H1) *H1 = gtsam::numericalDerivative41<V, P, P, P, P>(residualFn, T0, T1, T2, T3, 1e-5);
    if (H2) *H2 = gtsam::numericalDerivative42<V, P, P, P, P>(residualFn, T0, T1, T2, T3, 1e-5);
    if (H3) *H3 = gtsam::numericalDerivative43<V, P, P, P, P>(residualFn, T0, T1, T2, T3, 1e-5);
    if (H4) *H4 = gtsam::numericalDerivative44<V, P, P, P, P>(residualFn, T0, T1, T2, T3, 1e-5);

    // Bias is linear in the knots: b(u) = (1-B1)P0 + (B1-B2)P1 + (B2-B3)P2 + B3 P3
    if (H5 || H6 || H7 || H8) {
      Eigen::Vector4d B = (1.0 / 6.0) * C * calculateWeights(u_b_).col(0);
      const gtsam::Matrix6 I6 = gtsam::Matrix6::Identity();
      if (H5) *H5 = (1.0 - B[1]) * I6;
      if (H6) *H6 = (B[1] - B[2]) * I6;
      if (H7) *H7 = (B[2] - B[3]) * I6;
      if (H8) *H8 = B[3] * I6;
    }

    return residualFn(T0, T1, T2, T3);
  }

 private:
  gtsam::Vector3 z_w_, z_f_;
  double u_p_, dT_p_, u_b_;
};

void testGnssTrajectory();

}  // namespace CubicBasisSpline

#endif