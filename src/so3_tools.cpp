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
  Eigen::Matrix3d output;
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();

  double eps = 1e-4;

  double phi = vec.norm();

  Eigen::Matrix3d skew_vec;

  skew_vec = skew(vec);

  if (phi < eps) {
    output = I3 + skew_vec;
  } else {
    output = I3 + (std::sin(phi) / phi) * skew_vec +
             ((1 - std::cos(phi)) / (phi * phi)) * (skew_vec * skew_vec);
  }

  return output;
}

Eigen::Matrix3d logm(const Eigen::Matrix3d& mat) {}
}  // namespace CubicBasisSplines