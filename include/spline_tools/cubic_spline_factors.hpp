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

bool gnssFactor(SplineKnotSet& ctrl_pts, const gtsam::Point3& pos, const gtsam::Vector3& vel,
                const gtsam::Vector3& l_a, gtsam::Point3& residual_pos,
                gtsam::Vector3& residual_vel);

bool evaluateCubicBSpline(SplineKnotSet& ctrl_pts, SplineEvaluation& spline_eval);

class CubicGnssFactor
    : public gtsam::NoiseModelFactor4<gtsam::Pose3, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3> {
 public:
  typedef gtsam::NoiseModelFactor4<gtsam::Pose3, gtsam::Pose3, gtsam::Pose3, gtsam::Pose3> Base;

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
                              boost::optional<gtsam::Matrix&> H1 = boost::none,
                              boost::optional<gtsam::Matrix&> H2 = boost::none,
                              boost::optional<gtsam::Matrix&> H3 = boost::none,
                              boost::optional<gtsam::Matrix&> H4 = boost::none) const override {
    SplineKnotSet ctrl_pts{T_im1, T_i, T_ip1, T_ip2, u_, dT_};
    gtsam::Point3 res_pos;
    gtsam::Vector3 res_vel;

    gnssFactor(ctrl_pts, z_pos_, z_vel_, l_a_, res_pos, res_vel);
    gtsam::Vector6 residual;
    residual.head<3>() = res_pos;
    residual.tail<3>() = res_vel;

    if (H1) {
      *H1 = gtsam::numericalDerivative11<gtsam::Vector6, gtsam::Pose3>(
          [this, &T_i, &T_ip1, &T_ip2](const gtsam::Pose3& T_im1_val) {
            SplineKnotSet knots{T_im1_val, T_i, T_ip1, T_ip2, u_, dT_};
            gtsam::Point3 rp;
            gtsam::Vector3 rv;
            gnssFactor(knots, z_pos_, z_vel_, l_a_, rp, rv);
            gtsam::Vector6 res;
            res.head<3>() = rp;
            res.tail<3>() = rv;
            return res;
          },
          T_im1, 1e-5);
    }

    if (H2) {
      *H2 = gtsam::numericalDerivative11<gtsam::Vector6, gtsam::Pose3>(
          [this, &T_im1, &T_ip1, &T_ip2](const gtsam::Pose3& T_i_val) {
            SplineKnotSet knots{T_im1, T_i_val, T_ip1, T_ip2, u_, dT_};
            gtsam::Point3 rp;
            gtsam::Vector3 rv;
            gnssFactor(knots, z_pos_, z_vel_, l_a_, rp, rv);
            gtsam::Vector6 res;
            res.head<3>() = rp;
            res.tail<3>() = rv;
            return res;
          },
          T_i, 1e-5);
    }

    if (H3) {
      *H3 = gtsam::numericalDerivative11<gtsam::Vector6, gtsam::Pose3>(
          [this, &T_im1, &T_i, &T_ip2](const gtsam::Pose3& T_ip1_val) {
            SplineKnotSet knots{T_im1, T_i, T_ip1_val, T_ip2, u_, dT_};
            gtsam::Point3 rp;
            gtsam::Vector3 rv;
            gnssFactor(knots, z_pos_, z_vel_, l_a_, rp, rv);
            gtsam::Vector6 res;
            res.head<3>() = rp;
            res.tail<3>() = rv;
            return res;
          },
          T_ip1, 1e-5);
    }

    if (H4) {
      *H4 = gtsam::numericalDerivative11<gtsam::Vector6, gtsam::Pose3>(
          [this, &T_im1, &T_i, &T_ip1](const gtsam::Pose3& T_ip2_val) {
            SplineKnotSet knots{T_im1, T_i, T_ip1, T_ip2_val, u_, dT_};
            gtsam::Point3 rp;
            gtsam::Vector3 rv;
            gnssFactor(knots, z_pos_, z_vel_, l_a_, rp, rv);
            gtsam::Vector6 res;
            res.head<3>() = rp;
            res.tail<3>() = rv;
            return res;
          },
          T_ip2, 1e-5);
    }

    return residual;
  }

 private:
  gtsam::Point3 z_pos_;
  gtsam::Vector3 z_vel_;
  gtsam::Point3 l_a_;
  double u_;
  double dT_;
};

void optimizeGnssTrajectory();

}  // namespace CubicBasisSpline

#endif