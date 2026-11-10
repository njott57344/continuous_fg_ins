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

struct OdometryMeasurement {
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

struct ImuMeasurement {
  double timestamp_ns;
  double timestamp_s;
  double timestamp;

  double publish_time_ns;
  uint32_t header_stamp_sec;
  uint32_t header_stamp_nanosec;
  std::string header_frame_id;
  std::string child_frame_id;

  double w_x, w_y, w_z;
  double f_x, f_y, f_z;
};

// GPS Constants
const double gps_pi = M_PI;  // PI
const double GM = 3.986005 * pow(10, 14);
const double F = -4.442807633 * pow(10, -10);
const int half_week = 302400;                   // [s]
const double c = 299792458.0;                   // [m/s] speed of light
const double J2 = 1.082627 * pow(10, -3);       // J2 gravity constant (from Groves p. 72);
const double mu_g = 3.986004418 * pow(10, 14);  // Earth gravity constant (from Goves p. 71)
const double w_e = 7.292115 * pow(10, -5);      // rotation rate of earth [rad/s]
inline const Eigen::Vector3d w_e_ = (Eigen::Vector3d() << 0, 0, w_e).finished();

// Parameters of WGS84 Ellipsoid
const double Ro = 6378137.0;       // WGS84 equatorial radius [m]
const double e = 0.0818191908425;  // ecentricity of earth
const double e2 = std::pow(e, 2);
const double flattening = 1 / 298.257223563;
const double Rp = 6356752.31425;  // WGS84 polar radius [m]
const double A = 6378137.0;

Eigen::Vector3d lla2Ecef(const Eigen::Vector3d& lla);
Eigen::Vector3d lla2EcefVel(const Eigen::Vector3d& v_ned, const Eigen::Vector3d& lla);

Eigen::Matrix3d rotateLLaToEcef(const Eigen::Vector3d& lla, const Eigen::Matrix3d& Rlla);

Eigen::Vector3d ecefSomiglianaGravity(const Eigen::Vector3d& ecef_pos);

double transverseRadiusOfCurvature(const double& lat);

std::vector<OdometryMeasurement> readOdometryCsv(const std::string& csv_filename);
std::vector<ImuMeasurement> readImuCsv(const std::string& csv_filename);

std::vector<std::string> parseCSVLine(const std::string& line);
std::array<double, 36> parseCovarianceString(std::string cov_str);

}  // namespace CubicBasisSpline
#endif