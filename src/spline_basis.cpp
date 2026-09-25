#include "spline_tools/spline_basis.hpp"

namespace CubicBasisSpline {

SplineWeights calculateSplineWeights(double u) {
  SplineWeights w;

  static const Eigen::Matrix4d C =
      (1.0 / 6.0) * (Eigen::Matrix4d() << 6.0, 0.0, 0.0, 0.0, 5.0, 3.0, -3.0, 1.0, 1.0, 3.0, 3.0,
                     -2.0, 0.0, 0.0, 0.0, 1.0)
                        .finished();

  double u2 = u * u;
  double u3 = u2 * u;

  Eigen::Vector4d U(1.0, u, u2, u3);
  Eigen::Vector4d dU(0.0, 1.0, 2.0 * u, 3.0 * u2);
  Eigen::Vector4d ddU(0.0, 0.0, 2.0, 6.0 * u);

  Eigen::Vector4d Btilde = C * U;
  Eigen::Vector4d dBtilde = C * dU;
  Eigen::Vector4d ddBtilde = C * ddU;

  for (int i = 0; i < 4; i++) {
    w.Btilde[i] = Btilde[i];
    w.dBtilde[i] = dBtilde[i];
    w.ddBtilde[i] = ddBtilde[i];
  }

  return w;
}

}  // namespace CubicBasisSpline