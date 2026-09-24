#include "spline_tools/so3_tools.hpp"

namespace CubicBasisSplines {

Eigen::Matrix3d skew(const Eigen::Vector3d& vec) {
  Eigen::Matrix3d skew_out = Eigen::Matrix3d::Zero();

  skew_out(1, 0) = vec[2];
  skew_out(0, 1) = -1 * vec[2];

  skew_out(2, 0) = -1 * vec[1];
  skew_out(0, 2) = vec[1];

  skew_out(2, 1) = vec[0];
  skew_out(1, 2) = -1 * vec[0];

  return skew_out;
}

Eigen::Vector3d vee(const Eigen::Matrix3d& mat) {
  Eigen::Vector3d vector_out;

  vector_out[0] = mat(2, 1);
  vector_out[1] = mat(0, 2);
  vector_out[2] = mat(1, 0);

  return vector_out;
}

Eigen::Matrix3d expm(const Eigen::Vector3d& vec) {
  const double phi_sq = vec.squaredNorm();
  const double phi = std::sqrt(phi_sq);
  const Eigen::Matrix3d skew_vec = skew(vec);

  double A, B;
  if (phi < 1e-4) {
    // Taylor expansion for small angles
    A = 1.0 - phi_sq / 6.0;
    B = 0.5 - phi_sq / 24.0;
  } else {
    A = std::sin(phi) / phi;
    B = (1.0 - std::cos(phi)) / phi_sq;
  }

  return Eigen::Matrix3d::Identity() + A * skew_vec + B * (skew_vec * skew_vec);
}

Eigen::Vector3d logm(const Eigen::Matrix3d& mat) {
  // Clamp trace input to prevent acos domain errors
  const double cos_theta = std::clamp((mat.trace() - 1.0) * 0.5, -1.0, 1.0);
  const double theta = std::acos(cos_theta);
  const double theta_sq = theta * theta;

  const Eigen::Matrix3d skew_sym = 0.5 * (mat - mat.transpose());

  double scale;
  if (theta < 1e-4) {
    // Taylor expansion: theta / sin(theta) ~ 1 + theta^2 / 6
    scale = 1.0 + theta_sq / 6.0;
  } else {
    scale = theta / std::sin(theta);
  }

  return vee(scale * skew_sym);
}
}  // namespace CubicBasisSplines