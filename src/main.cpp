#include <iostream>

#include "spline_tools/cubic_spline_factors.hpp"
#include "spline_tools/nav_utils.hpp"

int main() {
  // gemini script testing a random toy example
  // CubicBasisSpline::testGnssTrajectory();

  /*
I have merged pos,vel at putnam, we will test first on just gnss at putnam,
then bring in the IMU
*/
  std::string gnss_a_merged =
      "/home/njo0004/atr_data/putnam_run2_gnss_merged_0/navigation_merged_rcvr_a.csv";

  std::string gnss_b_merged =
      "/home/njo0004/atr_data/putnam_run2_gnss_merged_0/navigation_merged_rcvr_b.csv";

  std::vector<CubicBasisSpline::SensorMeasurement> rcvr_a =
      CubicBasisSpline::readOdometryCsv(gnss_a_merged);

  std::vector<CubicBasisSpline::SensorMeasurement> rcvr_b =
      CubicBasisSpline::readOdometryCsv(gnss_b_merged);

  // lever-arms from localization (nav point to rcvr's)
  Eigen::Vector3d l_a;
  Eigen::Vector3d l_b;

  l_a << 1.54, -0.0751, -0.473;
  l_b << 0.419, 0.422, -0.454;

  // setting up factor graph
  gtsam::NonlinearFactorGraph graph;
  gtsam::Values initial_values;

  double dT = 1 / 5;  // 5 Hz knot spacing

  int num_a_measurements = rcvr_a.size();

  for (int i = 0; i < num_a_measurements; i++) {
    CubicBasisSpline::SensorMeasurement current_a = rcvr_a[i];

    Eigen::Vector3d lla_pos{current_a.pose_x, current_a.pose_y, current_a.pose_z};
    Eigen::Vector3d ecef_pos = CubicBasisSpline::lla2Ecef(lla_pos);

    Eigen::Vector3d ned_vel{current_a.linear_x, current_a.linear_y, current_a.linear_z};
  }

  return 0;
}