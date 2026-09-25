#include <ceres/ceres.h>
#include <ceres/gradient_checker.h>
#include <ceres/manifold.h>
#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <sophus/so3.hpp>
#include <vector>

#include "spline_tools/spline_factors.hpp"

namespace CubicBasisSpline {
namespace testing {

/**
 * @brief Helper to convert double or ceres::Jet scalars to pure double
 */
template <typename T>
inline double scalarToDouble(const T& val) {
  return static_cast<double>(val);
}

template <typename T, int N>
inline double scalarToDouble(const ceres::Jet<T, N>& val) {
  return val.a;
}

/**
 * @brief Safely constructs a Sophus::SO3 object from a raw 3x3 matrix block,
 * ensuring strict orthogonality even under finite-difference perturbations.
 */
template <typename T>
Sophus::SO3<T> matrixToSO3Safely(const T* const R_ptr) {
  Eigen::Map<const Eigen::Matrix<T, 3, 3, Eigen::RowMajor>> R_raw(R_ptr);

  // Extract pure double values for SVD orthogonalization
  Eigen::Matrix3d R_double;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      R_double(r, c) = scalarToDouble(R_raw(r, c));
    }
  }

  Eigen::JacobiSVD<Eigen::Matrix3d> svd(R_double, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d R_ortho = svd.matrixU() * svd.matrixV().transpose();

  if (R_ortho.determinant() < 0.0) {
    Eigen::Matrix3d U = svd.matrixU();
    U.col(2) *= -1.0;
    R_ortho = U * svd.matrixV().transpose();
  }

  // Cast back to scalar type T
  Eigen::Matrix<T, 3, 3> R_clean;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      R_clean(r, c) = T(R_ortho(r, c));
    }
  }

  return Sophus::SO3<T>(R_clean);
}

/**
 * @brief Ceres Manifold for a 3x3 Row-Major Rotation Matrix parameter block
 * (ambient size 9, tangent size 3)
 */
class SO3MatrixManifold : public ceres::Manifold {
 public:
  int AmbientSize() const override { return 9; }
  int TangentSize() const override { return 3; }

  // Plus: R_new = R * exp(tau)
  bool Plus(const double* x, const double* delta, double* x_plus_delta) const override {
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> R(x);
    Eigen::Map<const Eigen::Vector3d> tau(delta);

    Sophus::SO3d dR = Sophus::SO3d::exp(tau);

    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> R_new(x_plus_delta);
    R_new = R * dR.matrix();
    return true;
  }

  // PlusJacobian: d(R * exp(tau)) / d(tau) at tau = 0
  bool PlusJacobian(const double* x, double* jacobian) const override {
    Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> J(jacobian);
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> R(x);

    J.setZero();
    for (int i = 0; i < 3; ++i) {
      Eigen::Vector3d e_i = Eigen::Vector3d::Unit(i);
      Eigen::Matrix3d dR = R * Sophus::SO3d::hat(e_i);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          J(r * 3 + c, i) = dR(r, c);
        }
      }
    }
    return true;
  }

  // Minus: tau = log(R_x^T * R_y)
  bool Minus(const double* y, const double* x, double* y_minus_x) const override {
    Sophus::SO3d SO3_y = matrixToSO3Safely(y);
    Sophus::SO3d SO3_x = matrixToSO3Safely(x);

    Eigen::Map<Eigen::Vector3d> tau(y_minus_x);
    tau = (SO3_x.inverse() * SO3_y).log();
    return true;
  }

  bool MinusJacobian(const double* x, double* jacobian) const override {
    Eigen::Map<Eigen::Matrix<double, 3, 9, Eigen::RowMajor>> J(jacobian);
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> R(x);

    J.setZero();
    for (int i = 0; i < 3; ++i) {
      Eigen::Vector3d e_i = Eigen::Vector3d::Unit(i);
      Eigen::Matrix3d dR = R * Sophus::SO3d::hat(e_i);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          J(i, r * 3 + c) = dR(r, c);
        }
      }
    }
    return true;
  }
};

class GnssFactorsTest : public ::testing::Test {
 protected:
  void SetUp() override {
    dt_ = 0.2;  // 5 Hz knot interval
    u_ = 0.4;   // Sample point at 40% between knot i and i+1

    // GNSS antenna offset in body frame
    l_b_ << 0.2, -0.1, 0.5;

    // Measurement noise standard deviations
    sigmas_pos_ << 0.05, 0.05, 0.10;
    sigmas_vel_ << 0.02, 0.02, 0.05;

    // 1. Define 4 R3 control points
    p_gt_[0] << 0.0, 0.0, 0.0;
    p_gt_[1] << 1.0, 0.2, 0.05;
    p_gt_[2] << 2.1, 0.5, 0.12;
    p_gt_[3] << 3.3, 0.9, 0.20;

    // 2. Define 4 SO(3) control points
    SO3_gt_[0] = Sophus::SO3d::exp(Eigen::Vector3d(0.0, 0.0, 0.0));
    SO3_gt_[1] = Sophus::SO3d::exp(Eigen::Vector3d(0.01, 0.02, 0.05));
    SO3_gt_[2] = Sophus::SO3d::exp(Eigen::Vector3d(0.02, 0.04, 0.11));
    SO3_gt_[3] = Sophus::SO3d::exp(Eigen::Vector3d(0.03, 0.07, 0.18));
  }

  double dt_;
  double u_;
  Eigen::Vector3d l_b_;
  Eigen::Vector3d sigmas_pos_;
  Eigen::Vector3d sigmas_vel_;

  Eigen::Vector3d p_gt_[4];
  Sophus::SO3d SO3_gt_[4];
};

// -----------------------------------------------------------------------------
// TEST 1: Direct Residual Evaluation at Ground Truth
// -----------------------------------------------------------------------------
TEST_F(GnssFactorsTest, ZeroResidualAtGroundTruth) {
  const std::vector<double> u_samples = {0.0, 0.25, 0.5, 0.75, 0.99};

  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R_mats[4];
  for (int i = 0; i < 4; ++i) {
    R_mats[i] = SO3_gt_[i].matrix();
  }

  for (double u : u_samples) {
    R3SplineEval<double> r3_eval = evaluateR3(p_gt_[0], p_gt_[1], p_gt_[2], p_gt_[3], u, dt_);
    SO3SplineEval<double> so3_eval =
        evaluateSO3(SO3_gt_[0], SO3_gt_[1], SO3_gt_[2], SO3_gt_[3], u, dt_);

    Eigen::Vector3d z_pos_gt = r3_eval.p + so3_eval.C_bw * l_b_;
    Eigen::Vector3d z_vel_gt = r3_eval.v + so3_eval.C_bw * (so3_eval.w_bw.cross(l_b_));

    GnssPositionFactor pos_factor(z_pos_gt, l_b_, sigmas_pos_, u, dt_);
    GnssVelocityFactor vel_factor(z_vel_gt, l_b_, sigmas_vel_, u, dt_);

    Eigen::Vector3d pos_res, vel_res;
    EXPECT_TRUE(pos_factor(p_gt_[0].data(), p_gt_[1].data(), p_gt_[2].data(), p_gt_[3].data(),
                           R_mats[0].data(), R_mats[1].data(), R_mats[2].data(), R_mats[3].data(),
                           pos_res.data()));

    EXPECT_TRUE(vel_factor(p_gt_[0].data(), p_gt_[1].data(), p_gt_[2].data(), p_gt_[3].data(),
                           R_mats[0].data(), R_mats[1].data(), R_mats[2].data(), R_mats[3].data(),
                           vel_res.data()));

    EXPECT_LT(pos_res.norm(), 1e-9);
    EXPECT_LT(vel_res.norm(), 1e-9);
  }
}

// -----------------------------------------------------------------------------
// TEST 2: Finite-Difference Gradient Checking for AutoDiff
// -----------------------------------------------------------------------------
TEST_F(GnssFactorsTest, AutoDiffGradientCheck) {
  R3SplineEval<double> r3_eval = evaluateR3(p_gt_[0], p_gt_[1], p_gt_[2], p_gt_[3], u_, dt_);
  SO3SplineEval<double> so3_eval =
      evaluateSO3(SO3_gt_[0], SO3_gt_[1], SO3_gt_[2], SO3_gt_[3], u_, dt_);

  Eigen::Vector3d z_pos_gt = r3_eval.p + so3_eval.C_bw * l_b_;
  Eigen::Vector3d z_vel_gt = r3_eval.v + so3_eval.C_bw * (so3_eval.w_bw.cross(l_b_));

  ceres::CostFunction* pos_cost =
      new ceres::AutoDiffCostFunction<GnssPositionFactor, 3, 3, 3, 3, 3, 9, 9, 9, 9>(
          new GnssPositionFactor(z_pos_gt, l_b_, sigmas_pos_, u_, dt_));

  ceres::CostFunction* vel_cost =
      new ceres::AutoDiffCostFunction<GnssVelocityFactor, 3, 3, 3, 3, 3, 9, 9, 9, 9>(
          new GnssVelocityFactor(z_vel_gt, l_b_, sigmas_vel_, u_, dt_));

  SO3MatrixManifold so3_manifold;

  // Parameter block map: 4 R3 blocks (nullptr) followed by 4 SO3 blocks (&so3_manifold)
  std::vector<const ceres::Manifold*> manifolds = {nullptr,       nullptr,       nullptr,
                                                   nullptr,       &so3_manifold, &so3_manifold,
                                                   &so3_manifold, &so3_manifold};

  ceres::NumericDiffOptions num_diff_options;
  num_diff_options.relative_step_size = 1e-6;

  ceres::GradientChecker pos_checker(pos_cost, &manifolds, num_diff_options);
  ceres::GradientChecker vel_checker(vel_cost, &manifolds, num_diff_options);

  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R_mats[4];
  for (int i = 0; i < 4; ++i) {
    R_mats[i] = SO3_gt_[i].matrix();
  }

  const double* params[8] = {p_gt_[0].data(),  p_gt_[1].data(),  p_gt_[2].data(),
                             p_gt_[3].data(),  R_mats[0].data(), R_mats[1].data(),
                             R_mats[2].data(), R_mats[3].data()};

  ceres::GradientChecker::ProbeResults pos_results, vel_results;
  const double relative_precision_tolerance = 1e-4;

  EXPECT_TRUE(pos_checker.Probe(params, relative_precision_tolerance, &pos_results))
      << "Position Factor AutoDiff failed gradient check:\n"
      << pos_results.error_log;

  EXPECT_TRUE(vel_checker.Probe(params, relative_precision_tolerance, &vel_results))
      << "Velocity Factor AutoDiff failed gradient check:\n"
      << vel_results.error_log;
}

// -----------------------------------------------------------------------------
// TEST 3: Observable Translation Control Point Optimization
// -----------------------------------------------------------------------------
TEST_F(GnssFactorsTest, OptimizationTranslationRecovery) {
  const std::vector<double> u_samples = {0.0, 0.25, 0.5, 0.75};
  const size_t num_meas = u_samples.size();

  std::vector<Eigen::Vector3d> z_pos_gt(num_meas);
  std::vector<Eigen::Vector3d> z_vel_gt(num_meas);

  for (size_t k = 0; k < num_meas; ++k) {
    R3SplineEval<double> r3_eval =
        evaluateR3(p_gt_[0], p_gt_[1], p_gt_[2], p_gt_[3], u_samples[k], dt_);
    SO3SplineEval<double> so3_eval =
        evaluateSO3(SO3_gt_[0], SO3_gt_[1], SO3_gt_[2], SO3_gt_[3], u_samples[k], dt_);

    z_pos_gt[k] = r3_eval.p + so3_eval.C_bw * l_b_;
    z_vel_gt[k] = r3_eval.v + so3_eval.C_bw * (so3_eval.w_bw.cross(l_b_));
  }

  // Corrupt interior R3 control point initial estimates
  Eigen::Vector3d p_opt[4] = {p_gt_[0], p_gt_[1] + Eigen::Vector3d(0.5, -0.3, 0.2),
                              p_gt_[2] + Eigen::Vector3d(-0.4, 0.2, -0.1), p_gt_[3]};

  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R_mats[4];
  for (int i = 0; i < 4; ++i) {
    R_mats[i] = SO3_gt_[i].matrix();
  }

  ceres::Problem problem;

  for (size_t k = 0; k < num_meas; ++k) {
    ceres::CostFunction* pos_cost =
        new ceres::AutoDiffCostFunction<GnssPositionFactor, 3, 3, 3, 3, 3, 9, 9, 9, 9>(
            new GnssPositionFactor(z_pos_gt[k], l_b_, sigmas_pos_, u_samples[k], dt_));

    ceres::CostFunction* vel_cost =
        new ceres::AutoDiffCostFunction<GnssVelocityFactor, 3, 3, 3, 3, 3, 9, 9, 9, 9>(
            new GnssVelocityFactor(z_vel_gt[k], l_b_, sigmas_vel_, u_samples[k], dt_));

    problem.AddResidualBlock(pos_cost, nullptr, p_opt[0].data(), p_opt[1].data(), p_opt[2].data(),
                             p_opt[3].data(), R_mats[0].data(), R_mats[1].data(), R_mats[2].data(),
                             R_mats[3].data());

    problem.AddResidualBlock(vel_cost, nullptr, p_opt[0].data(), p_opt[1].data(), p_opt[2].data(),
                             p_opt[3].data(), R_mats[0].data(), R_mats[1].data(), R_mats[2].data(),
                             R_mats[3].data());
  }

  // Pin boundary points and rotation states constant
  problem.SetParameterBlockConstant(p_opt[0].data());
  problem.SetParameterBlockConstant(p_opt[3].data());
  for (int i = 0; i < 4; ++i) {
    problem.SetParameterBlockConstant(R_mats[i].data());
  }

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  EXPECT_TRUE(summary.IsSolutionUsable());
  EXPECT_LT(summary.final_cost, 1e-10);
  EXPECT_LT((p_opt[1] - p_gt_[1]).norm(), 1e-3);
  EXPECT_LT((p_opt[2] - p_gt_[2]).norm(), 1e-3);
}

}  // namespace testing
}  // namespace CubicBasisSpline

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}