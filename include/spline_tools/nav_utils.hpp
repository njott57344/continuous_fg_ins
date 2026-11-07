#ifndef NAV_UTILS_HPP
#define NAV_UTILS_HPP

#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "cmath"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "iostream"

namespace CubicBasisSpline {

struct SensorMeasurement {
  double timestamp_ns;
  double timestamp_s;
  double timestamp;

  double publish_time_ns;
  uint32_t header_stamp_sec;
  uint32_t header_stamp_nanosec;
  std::string header_frame_id;
  std::string child_frame_id;

  // Pose
  double pose_x, pose_y, pose_z;
  double orient_x, orient_y, orient_z, orient_w;
  std::array<double, 36> pose_covariance;

  // Twist
  double linear_x, linear_y, linear_z;
  double angular_x, angular_y, angular_z;
  std::array<double, 36> twist_covariance;
};

Eigen::Vector3d lla2Ecef(const Eigen::Vector3d& lla);
Eigen::Vector3d lla2EcefVel(const Eigen::Vector3d& v_ned, const Eigen::Vector3d& lla);

Eigen::Matrix3d rotateLLaToEcef(const Eigen::Vector3d& lla, const Eigen::Matrix3d& Rlla);

double transverseRadiusOfCurvature(const double& lat);

std::vector<SensorMeasurement> readOdometryCsv(const std::string& csv_filename);
std::vector<std::string> parseCSVLine(const std::string& line);
std::array<double, 36> parseCovarianceString(std::string cov_str);

}  // namespace CubicBasisSpline
#endif