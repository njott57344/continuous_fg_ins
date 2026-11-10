#include <algorithm>
#include <iostream>

#include "spline_tools/cubic_spline_factors.hpp"
#include "spline_tools/nav_utils.hpp"

int findCtrlPtIdx(const std::vector<double>& vec, double target) {
  auto it = std::upper_bound(vec.begin(), vec.end(), target);
  return static_cast<int>(std::distance(vec.begin(), it)) - 1;
}

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

  std::vector<CubicBasisSpline::OdometryMeasurement> rcvr_a =
      CubicBasisSpline::readOdometryCsv(gnss_a_merged);

  std::vector<CubicBasisSpline::OdometryMeasurement> rcvr_b =
      CubicBasisSpline::readOdometryCsv(gnss_b_merged);

  std::string imu_a_file =
      "/home/njo0004/atr_data/putnam_run2_gnss_merged_0/sensors_novatel_a_imu_data_raw.csv";
  std::string imu_b_file =
      "/home/njo0004/data/putnam_run2_gnss_merged_0/sensors_novatel_b_imu_data_raw.csv";

  std::vector<CubicBasisSpline::ImuMeasurement> imu_a = CubicBasisSpline::readImuCsv(imu_a_file);
  // std::vector<CubicBasisSpline::ImuMeasurement> imu_b = CubicBasisSpline::readImuCsv(imu_b_file);

  int num_a_measurements = rcvr_a.size();
  int num_b_measurements = rcvr_b.size();

  int num_a_imu = imu_a.size();
  // int num_b_imu = imu_b.size();

  // lever-arms from localization (nav point to rcvr's)
  Eigen::Vector3d l_a;
  Eigen::Vector3d l_b;

  l_a << 1.54, -0.0751, -0.473;
  l_b << 0.419, 0.422, -0.454;

  // setting up factor graph
  gtsam::NonlinearFactorGraph graph;
  gtsam::Values initial_values;

  // pose b spline
  double dT = 1.0 / 5.0;  // 5 Hz knot spacing

  /*
    We need 2 control points before the first timestamp
    We need 2 control points after the last timestamp
  */

  double start_time_a = rcvr_a[0].timestamp;                      // start time of data
  double stop_time_a = rcvr_a[num_a_measurements - 1].timestamp;  // stop time of data

  int num_ctrl_pts = static_cast<int>(std::floor((stop_time_a - start_time_a) / dT)) + 4;

  std::vector<double> ctrl_pt_times;
  int time_idx = -1;  // this is the index we are at w.r.t. the start time

  // before the knot loop: precompute ECEF positions and a timestamp vector
  std::vector<double> meas_times;
  std::vector<Eigen::Vector3d> meas_ecef;
  for (auto& m : rcvr_a) {
    meas_times.push_back(m.timestamp);
    meas_ecef.push_back(CubicBasisSpline::lla2Ecef({m.pose_x, m.pose_y, m.pose_z}));
  }

  auto nedToEcef = [](double lat, double lon) {  // radians
    Eigen::Matrix3d R;
    R << -sin(lat) * cos(lon), -sin(lon), -cos(lat) * cos(lon), -sin(lat) * sin(lon), cos(lon),
        -cos(lat) * sin(lon), cos(lat), 0.0, -sin(lat);
    return R;
  };

  for (int i = 0; i < num_ctrl_pts; i++) {
    gtsam::Key key = gtsam::Symbol('T', i);
    // in the knot loop, replace the origin pose:
    double t_knot = start_time_a + time_idx * dT;
    auto it = std::lower_bound(meas_times.begin(), meas_times.end(), t_knot);
    size_t k = std::min<size_t>(std::distance(meas_times.begin(), it), meas_times.size() - 1);
    // optional: step back one if the previous sample is closer

    const auto& m = rcvr_a[k];
    double lat = m.pose_x * M_PI / 180.0,
           lon = m.pose_y * M_PI / 180.0;  // drop the conversion if already radians
    gtsam::Pose3 initial_pose(gtsam::Rot3(nedToEcef(lat, lon)), gtsam::Point3(meas_ecef[k]));
    initial_values.insert(key, initial_pose);

    ctrl_pt_times.push_back(start_time_a + time_idx * dT);
    time_idx++;
  }

  time_idx = -1;

  // bias b splines
  double dT_bias = 5.0;  // fairly sparse (biases don't change quickly)
  int num_imu_pts = static_cast<int>(std::floor((stop_time_a - start_time_a) / dT_bias)) + 4;
  std::vector<double> bias_times;

  for (int i = 0; i < num_imu_pts; i++) {
    gtsam::Key key_w = gtsam::Symbol('w', i);
    gtsam::Key key_f = gtsam::Symbol('f', i);

    double t_knot = start_time_a + time_idx * dT_bias;

    // imu biases to zero
    gtsam::Vector3 Z3 = gtsam::Vector3::Zero();

    // initial_values.insert(key_w, Z3);
    // initial_values.insert(key_f, Z3);

    bias_times.push_back(start_time_a + time_idx * dT_bias);
    time_idx++;
  }

  // adding gnss A
  for (int i = 0; i < num_a_measurements; i++) {
    CubicBasisSpline::OdometryMeasurement current_a = rcvr_a[i];

    // position
    Eigen::Vector3d lla_pos{current_a.pose_x, current_a.pose_y, current_a.pose_z};
    Eigen::Vector3d ecef_pos = CubicBasisSpline::lla2Ecef(lla_pos);

    // position covariance
    Eigen::Matrix3d lla_cov = Eigen::Matrix3d::Zero();
    lla_cov(0, 0) = current_a.pose_covariance[0];
    lla_cov(1, 1) = current_a.pose_covariance[3];
    lla_cov(2, 2) = current_a.pose_covariance[6];

    Eigen::Matrix3d ecef_cov_p = CubicBasisSpline::rotateLLaToEcef(lla_pos, lla_cov);

    // velocity
    Eigen::Vector3d ned_vel{current_a.linear_x, current_a.linear_y, current_a.linear_z};
    Eigen::Vector3d ecef_vel = CubicBasisSpline::lla2EcefVel(ned_vel, lla_pos);

    // velocity covariance
    Eigen::Matrix3d ned_cov = Eigen::Matrix3d::Zero();
    ned_cov(0, 0) = current_a.twist_covariance[0];
    ned_cov(1, 1) = current_a.twist_covariance[3];
    ned_cov(2, 2) = current_a.twist_covariance[6];

    Eigen::Matrix3d ecef_cov_v = CubicBasisSpline::rotateLLaToEcef(lla_pos, ned_cov);

    auto gnss_noise = gtsam::noiseModel::Diagonal::Sigmas(
        (gtsam::Vector(6) << std::sqrt(ecef_cov_p(0, 0)), std::sqrt(ecef_cov_p(1, 1)),
         std::sqrt(ecef_cov_p(2, 2)), 0.4, 0.4, 0.4)
            .finished());

    // adding measurement to graph
    /*
      1) We need to identify the set of control points we are working with (based on time of
      measurement) 2) We need to take the data and put it in custom cubic b spline gnss factor
      3) Solve?
    */

    double t_meas = current_a.timestamp;

    // indices for keys for graph
    int idx_t = findCtrlPtIdx(ctrl_pt_times, t_meas);
    int idx_tm1 = idx_t - 1;
    int idx_tp1 = idx_t + 1;
    int idx_tp2 = idx_t + 2;

    gtsam::Key key_im1 = gtsam::Symbol('T', idx_tm1);
    gtsam::Key key_i = gtsam::Symbol('T', idx_t);
    gtsam::Key key_ip1 = gtsam::Symbol('T', idx_tp1);
    gtsam::Key key_ip2 = gtsam::Symbol('T', idx_tp2);

    double time_t = ctrl_pt_times[idx_t];
    double u = (t_meas - time_t) / dT;

    graph.add(CubicBasisSpline::CubicGnssFactor(key_im1, key_i, key_ip1, key_ip2, ecef_pos,
                                                ecef_vel, l_a, u, dT, gnss_noise));
  }

  // adding gnss b
  for (int i = 0; i < num_b_measurements; i++) {
    CubicBasisSpline::OdometryMeasurement current_b = rcvr_b[i];

    // position
    Eigen::Vector3d lla_pos{current_b.pose_x, current_b.pose_y, current_b.pose_z};
    Eigen::Vector3d ecef_pos = CubicBasisSpline::lla2Ecef(lla_pos);

    // position covariance
    Eigen::Matrix3d lla_cov = Eigen::Matrix3d::Zero();
    lla_cov(0, 0) = current_b.pose_covariance[0];
    lla_cov(1, 1) = current_b.pose_covariance[3];
    lla_cov(2, 2) = current_b.pose_covariance[6];

    Eigen::Matrix3d ecef_cov_p = CubicBasisSpline::rotateLLaToEcef(lla_pos, lla_cov);

    // velocity
    Eigen::Vector3d ned_vel{current_b.linear_x, current_b.linear_y, current_b.linear_z};
    Eigen::Vector3d ecef_vel = CubicBasisSpline::lla2EcefVel(ned_vel, lla_pos);

    // velocity covariance
    Eigen::Matrix3d ned_cov = Eigen::Matrix3d::Zero();
    ned_cov(0, 0) = current_b.twist_covariance[0];
    ned_cov(1, 1) = current_b.twist_covariance[3];
    ned_cov(2, 2) = current_b.twist_covariance[6];

    Eigen::Matrix3d ecef_cov_v = CubicBasisSpline::rotateLLaToEcef(lla_pos, ned_cov);

    auto gnss_noise = gtsam::noiseModel::Diagonal::Sigmas(
        (gtsam::Vector(6) << std::sqrt(ecef_cov_p(0, 0)), std::sqrt(ecef_cov_p(1, 1)),
         std::sqrt(ecef_cov_p(2, 2)), 0.4, 0.4, 0.4)
            .finished());

    // adding measurement to graph
    /*
      1) We need to identify the set of control points we are working with (based on time of
      measurement) 2) We need to take the data and put it in custom cubic b spline gnss factor
      3) Solve?
    */

    double t_meas = current_b.timestamp;

    // indices for keys for graph
    int idx_t = findCtrlPtIdx(ctrl_pt_times, t_meas);
    int idx_tm1 = idx_t - 1;
    int idx_tp1 = idx_t + 1;
    int idx_tp2 = idx_t + 2;

    gtsam::Key key_im1 = gtsam::Symbol('T', idx_tm1);
    gtsam::Key key_i = gtsam::Symbol('T', idx_t);
    gtsam::Key key_ip1 = gtsam::Symbol('T', idx_tp1);
    gtsam::Key key_ip2 = gtsam::Symbol('T', idx_tp2);

    double time_t = ctrl_pt_times[idx_t];
    double u = (t_meas - time_t) / dT;

    graph.add(CubicBasisSpline::CubicGnssFactor(key_im1, key_i, key_ip1, key_ip2, ecef_pos,
                                                ecef_vel, l_b, u, dT, gnss_noise));
  }

  // adding imu a
  for (int i = 0; i < num_a_imu; i++) {
    CubicBasisSpline::ImuMeasurement curr_imu = imu_a[i];

    // putting in vehicle body frame

    gtsam::Vector3 w_ib_b{curr_imu.w_y, curr_imu.w_x, -1 * curr_imu.w_z};
    gtsam::Vector3 f_ib_b{curr_imu.f_y, curr_imu.f_x, -1 * curr_imu.f_z};

    double curr_time = curr_imu.timestamp;
  }

  gtsam::LevenbergMarquardtParams params;
  params.setVerbosityLM("SUMMARY");  // prints each lambda trial and its error
  params.setVerbosity("ERROR");
  params.setMaxIterations(50);

  gtsam::LevenbergMarquardtOptimizer optimizer(graph, initial_values, params);
  gtsam::Values result = optimizer.optimize();

  double pos_cost = 0, vel_cost = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    auto nf = std::dynamic_pointer_cast<gtsam::NoiseModelFactor>(graph[i]);
    gtsam::Vector w = nf->whitenedError(result);
    pos_cost += 0.5 * w.head<3>().squaredNorm();
    vel_cost += 0.5 * w.tail<3>().squaredNorm();
  }
  std::cout << "pos: " << pos_cost << "  vel: " << vel_cost << std::endl;

  {
    std::ofstream g("ctrl_pts.g2o");
    g << std::setprecision(17);
    for (int i = 0; i < num_ctrl_pts; i++) {
      gtsam::Pose3 T = result.at<gtsam::Pose3>(gtsam::Symbol('T', i));
      Eigen::Quaterniond q(T.rotation().matrix());  // body -> ECEF
      g << "VERTEX_SE3:QUAT " << i << " " << T.x() << " " << T.y() << " " << T.z() << " " << q.x()
        << " " << q.y() << " " << q.z() << " " << q.w() << "\n";  // g2o order: x y z w
    }

    std::ofstream meta("ctrl_pts_meta.txt");
    meta << std::setprecision(17) << "dT " << dT << "\n"
         << "l_a " << l_a.transpose() << "\n"
         << "start_time " << start_time_a << "\n";
  }

  return 0;
}