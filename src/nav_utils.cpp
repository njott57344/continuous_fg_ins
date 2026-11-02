#include "spline_tools/nav_utils.hpp"

namespace CubicBasisSpline {

Eigen::Vector3d lla2Ecef(const Eigen::Vector3d& lla) {
  // see groves 2.112
  double lat, lon, alt;
  double x, y, z;
  Eigen::Vector3d ecef;

  lat = lla(0) * M_PI / 180.0;
  lon = lla(1) * M_PI / 180.0;
  alt = lla(2);

  double Re = transverseRadiusOfCurvature(lat);
  double e2 = pow(0.0818191908425, 2);

  x = (Re + alt) * std::cos(lat) * std::cos(lon);
  y = (Re + alt) * std::cos(lat) * std::sin(lon);
  z = ((1 - e2) * Re + alt) * sin(lat);

  ecef << x, y, z;
  return ecef;
}

Eigen::Matrix3d rotateLLaToEcef(const Eigen::Vector3d& lla, const Eigen::Matrix3d& Rlla) {
  Eigen::Matrix3d Recef;                            // output covariance matrix
  Eigen::Matrix3d C = Eigen::Matrix3d::Identity();  // lla to ecef rotation matrix

  double lat, lon;

  lat = lla(0) * M_PI / 180.0;
  lon = lla(1) * M_PI / 180.0;

  C(0, 0) = -std::sin(lat) * std::cos(lon);
  C(0, 1) = -std::sin(lon);
  C(0, 2) = -std::cos(lat) * std::cos(lon);

  C(1, 0) = -std::sin(lat) * std::sin(lon);
  C(1, 1) = std::cos(lon);
  C(1, 2) = -std::cos(lat) * std::sin(lon);

  C(2, 0) = std::cos(lat);
  C(2, 1) = 0;
  C(2, 2) = -std::sin(lat);

  Recef = C * Rlla * C.transpose();

  return Recef;
}

double transverseRadiusOfCurvature(const double& lat) {
  double denom = 1 - pow(0.0818191908425, 2) * pow(sin(lat), 2);
  return 6378137.0 / sqrt(denom);
}
}  // namespace CubicBasisSpline
