#ifndef CUBIC_SPLINE_FACTORS_HPP
#define CUBIC_SPLINE_FACTORS_HPP

#include <gtsam/inference/Symbol.h>
#include <gtsam/linear/GaussianFactorGraph.h>
#include <gtsam/linear/HessianFactor.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>

#include "array"
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

void gnssJacobians(const SplineKnotSet& knots, const gtsam::Vector3& l_a,
                   std::array<gtsam::Matrix6, 4>& H);
void imuPoseJacobians(const SplineKnotSet& knots, std::array<gtsam::Matrix6, 4>& H);

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
    SplineKnotSet knots{T_im1, T_i, T_ip1, T_ip2, u_, dT_};
    gtsam::Point3 rp;
    gtsam::Vector3 rv;
    if (!gnssFactor(knots, z_pos_, z_vel_, l_a_, rp, rv))
      throw std::runtime_error("gnssFactor failed");
    gtsam::Vector6 res;
    res << rp, rv;

    if (H1 || H2 || H3 || H4) {
      std::array<gtsam::Matrix6, 4> H;
      gnssJacobians(knots, l_a_, H);
      if (H1) *H1 = H[0];
      if (H2) *H2 = H[1];
      if (H3) *H3 = H[2];
      if (H4) *H4 = H[3];
    }
    return res;
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

    SplineKnotSet knots{T0, T1, T2, T3, u_p_, dT_p_};
    gtsam::Vector3 rw = gtsam::Vector3::Zero(), rf = gtsam::Vector3::Zero();
    FullBiasSet bs = bias_set;
    if (!imuFactor(knots, bs, z_w_, z_f_, rw, rf)) throw std::runtime_error("imuFactor failed");
    gtsam::Vector6 res;
    res << rf, rw;

    if (H1 || H2 || H3 || H4) {
      std::array<gtsam::Matrix6, 4> H;
      imuPoseJacobians(knots, H);
      if (H1) *H1 = H[0];
      if (H2) *H2 = H[1];
      if (H3) *H3 = H[2];
      if (H4) *H4 = H[3];
    }

    // Bias is linear in the knots: b(u) = (1-B1)P0 + (B1-B2)P1 + (B2-B3)P2 + B3 P3
    if (H5 || H6 || H7 || H8) {
      Eigen::Vector4d B = (1.0 / 6.0) * C * calculateWeights(u_b_).col(0);
      const gtsam::Matrix6 I6 = gtsam::Matrix6::Identity();
      if (H5) *H5 = (1.0 - B[1]) * I6;
      if (H6) *H6 = (B[1] - B[2]) * I6;
      if (H7) *H7 = (B[2] - B[3]) * I6;
      if (H8) *H8 = B[3] * I6;
    }

    return res;
  }

 private:
  gtsam::Vector3 z_w_, z_f_;
  double u_p_, dT_p_, u_b_;
};

class ImuWindowFactor : public gtsam::NonlinearFactor {
 public:
  using Sample = CubicBasisSpline::CubicImuFactor;
  const std::vector<std::shared_ptr<Sample>>& samples() const { return samples_; }
  explicit ImuWindowFactor(std::vector<std::shared_ptr<Sample>> samples)
      : gtsam::NonlinearFactor(samples.front()->keys()), samples_(std::move(samples)) {}

  double error(const gtsam::Values& v) const override {
    double e = 0;
    for (auto& s : samples_) e += s->error(v);
    return e;
  }
  size_t dim() const override { return 6 * samples_.size(); }

  std::shared_ptr<gtsam::GaussianFactor> linearize(const gtsam::Values& v) const override {
    gtsam::GaussianFactorGraph gfg;
    for (auto& s : samples_) gfg.push_back(s->linearize(v));
    return std::make_shared<gtsam::HessianFactor>(gfg);  // sums J^T J, J^T r, r^T r
  }

  gtsam::NonlinearFactor::shared_ptr clone() const override {
    return std::make_shared<ImuWindowFactor>(*this);
  }

 private:
  std::vector<std::shared_ptr<Sample>> samples_;
};

}  // namespace CubicBasisSpline

#endif