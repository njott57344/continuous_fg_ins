#ifndef NAV_UTILS_HPP
#define NAV_UTILS_HPP

#include "cmath"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "iostream"

namespace CubicBasisSpline {

Eigen::Vector3d lla2Ecef(const Eigen::Vector3d& lla);
Eigen::Matrix3d rotateLLaToEcef(const Eigen::Vector3d& lla, const Eigen::Matrix3d& Rlla);
double transverseRadiusOfCurvature(const double& lat);

}  // namespace CubicBasisSpline
#endif