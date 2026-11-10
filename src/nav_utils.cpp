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

  x = (Re + alt) * std::cos(lat) * std::cos(lon);
  y = (Re + alt) * std::cos(lat) * std::sin(lon);
  z = ((1 - CubicBasisSpline::e2) * Re + alt) * sin(lat);

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

Eigen::Vector3d lla2EcefVel(const Eigen::Vector3d& v_ned, const Eigen::Vector3d& lla) {
  Eigen::Matrix3d C = Eigen::Matrix3d::Identity();

  double lat, lon;

  lat = lla[0] * M_PI / 180.0;
  lon = lla[1] * M_PI / 180.0;

  C(0, 0) = -std::sin(lat) * std::cos(lon);
  C(0, 1) = -std::sin(lon);
  C(0, 2) = -std::cos(lat) * std::cos(lon);

  C(1, 0) = -std::sin(lat) * std::sin(lon);
  C(1, 1) = std::cos(lon);
  C(1, 2) = -std::cos(lat) * std::sin(lon);

  C(2, 0) = std::cos(lat);
  C(2, 1) = 0;
  C(2, 2) = -std::sin(lat);

  return C * v_ned;
}

double transverseRadiusOfCurvature(const double& lat) {
  double denom = 1 - CubicBasisSpline::e2 * pow(sin(lat), 2);
  return CubicBasisSpline::Ro / sqrt(denom);
}

std::vector<std::string> parseCSVLine(const std::string& line) {
  std::vector<std::string> fields;
  std::string field;
  bool in_quotes = false;

  for (size_t i = 0; i < line.length(); ++i) {
    char c = line[i];
    if (c == '"') {
      in_quotes = !in_quotes;  // Toggle quote flag
    } else if (c == ',' && !in_quotes) {
      fields.push_back(field);
      field.clear();
    } else {
      field += c;
    }
  }
  fields.push_back(field);
  return fields;
}

Eigen::Vector3d ecefSomiglianaGravity(const Eigen::Vector3d& ecef_pos) {
  Eigen::Vector3d inner;
  Eigen::Vector3d ecef_gravitation;
  Eigen::Vector3d ecef_gravity;
  Eigen::Matrix3d C;

  double term_x, term_y, term_z, common_term;
  double abs_pos = ecef_pos.norm();
  double pos_3, R02, pos_2;
  double outer_lump, inner_lump;

  common_term = 5 * pow((ecef_pos[2] / abs_pos), 2);
  term_x = (1 - common_term) * ecef_pos[0];
  term_y = (1 - common_term) * ecef_pos[1];
  term_z = (3 - common_term) * ecef_pos[2];

  inner << term_x, term_y, term_z;

  pos_3 = pow(abs_pos, 3);
  pos_2 = pow(abs_pos, 2);
  R02 = pow(Ro, 2);

  outer_lump = mu_g / pos_3;            // lumped together constants
  inner_lump = 1.5 * J2 * R02 / pos_2;  // lumped together constants

  ecef_gravitation = -outer_lump * (ecef_pos + inner_lump * inner);

  C.setZero();
  C(0, 0) = 1;
  C(1, 1) = 1;

  if (ecef_pos.norm() == 0) {
    ecef_gravity.setZero();
  } else {
    ecef_gravity = ecef_gravitation + w_e * w_e * C * ecef_pos;
  }
  return ecef_gravity;
}

// Helper to parse 36 elements from a covariance string (strips brackets and splits by whitespace)
std::array<double, 36> parseCovarianceString(std::string cov_str) {
  std::array<double, 36> cov{};
  // Remove brackets and quotes if present
  cov_str.erase(std::remove(cov_str.begin(), cov_str.end(), '['), cov_str.end());
  cov_str.erase(std::remove(cov_str.begin(), cov_str.end(), ']'), cov_str.end());
  cov_str.erase(std::remove(cov_str.begin(), cov_str.end(), '"'), cov_str.end());

  std::stringstream ss(cov_str);
  for (int i = 0; i < 36; ++i) {
    std::string val;
    if (ss >> val) {
      try {
        cov[i] = std::stod(val);
      } catch (...) {
        cov[i] = 0.0;
      }
    }
  }
  return cov;
}

std::vector<OdometryMeasurement> readOdometryCsv(const std::string& csv_filename) {
  std::vector<OdometryMeasurement> measurements;
  std::ifstream file(csv_filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << csv_filename << std::endl;
    return measurements;
  }

  std::string line;
  // Skip header row
  if (!std::getline(file, line)) {
    return measurements;
  }

  std::string accumulated_record = "";
  int open_brackets = 0;

  while (std::getline(file, line)) {
    // Count brackets in the current line to track if we are inside a multiline covariance array
    for (char c : line) {
      if (c == '[') open_brackets++;
      if (c == ']') open_brackets--;
    }

    if (!accumulated_record.empty()) {
      accumulated_record += "\n" + line;
    } else {
      accumulated_record = line;
    }

    // Only process the record once all brackets are closed and we aren't splitting mid-row
    if (open_brackets <= 0) {
      open_brackets = 0;  // safety clamp
      std::vector<std::string> fields = parseCSVLine(accumulated_record);
      accumulated_record.clear();  // Reset for the next record

      if (fields.size() < 21) {
        continue;  // Skip malformed or empty padding lines
      }

      OdometryMeasurement m{};
      try {
        m.timestamp_ns = std::stod(fields[0]);
        m.timestamp_s = std::stod(fields[1]);
        m.publish_time_ns = std::stod(fields[2]);
        m.header_stamp_sec = static_cast<uint32_t>(std::stoul(fields[3]));
        m.header_stamp_nanosec = static_cast<uint32_t>(std::stoul(fields[4]));

        m.timestamp = static_cast<double>(m.header_stamp_sec) +
                      1e-9 * static_cast<double>(m.header_stamp_nanosec);

        m.header_frame_id = fields[5];
        m.child_frame_id = fields[6];

        m.pose_x = std::stod(fields[7]);
        m.pose_y = std::stod(fields[8]);
        m.pose_z = std::stod(fields[9]);

        m.orient_x = std::stod(fields[10]);
        m.orient_y = std::stod(fields[11]);
        m.orient_z = std::stod(fields[12]);
        m.orient_w = std::stod(fields[13]);

        m.pose_covariance = parseCovarianceString(fields[14]);

        m.linear_x = std::stod(fields[15]);
        m.linear_y = std::stod(fields[16]);
        m.linear_z = std::stod(fields[17]);

        m.angular_x = std::stod(fields[18]);
        m.angular_y = std::stod(fields[19]);
        m.angular_z = std::stod(fields[20]);

        measurements.push_back(m);
      } catch (const std::exception& e) {
        // Silently catch or handle parsing exceptions
      }
    }
  }

  file.close();
  return measurements;
}

std::vector<ImuMeasurement> readImuCsv(const std::string& csv_filename) {
  std::vector<ImuMeasurement> measurements;
  std::ifstream file(csv_filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << csv_filename << std::endl;
    return measurements;
  }

  std::string line;
  // Skip header row
  if (!std::getline(file, line)) {
    return measurements;
  }

  std::string accumulated_record = "";
  int open_brackets = 0;

  while (std::getline(file, line)) {
    // Count brackets in the current line to track if we are inside a multiline covariance array
    for (char c : line) {
      if (c == '[') open_brackets++;
      if (c == ']') open_brackets--;
    }

    if (!accumulated_record.empty()) {
      accumulated_record += "\n" + line;
    } else {
      accumulated_record = line;
    }

    // Only process the record once all brackets are closed and we aren't splitting mid-row
    if (open_brackets <= 0) {
      open_brackets = 0;  // safety clamp
      std::vector<std::string> fields = parseCSVLine(accumulated_record);
      accumulated_record.clear();  // Reset for the next record

      if (fields.size() < 21) {
        continue;  // Skip malformed or empty padding lines
      }

      ImuMeasurement m{};
      try {
        m.timestamp_ns = std::stod(fields[0]);
        m.timestamp_s = std::stod(fields[1]);
        m.publish_time_ns = std::stod(fields[2]);
        m.header_stamp_sec = static_cast<uint32_t>(std::stoul(fields[3]));
        m.header_stamp_nanosec = static_cast<uint32_t>(std::stoul(fields[4]));

        m.timestamp = static_cast<double>(m.header_stamp_sec) +
                      1e-9 * static_cast<double>(m.header_stamp_nanosec);

        m.header_frame_id = fields[5];
        m.child_frame_id = fields[6];

        m.w_x = std::stod(fields[11]);
        m.w_y = std::stod(fields[12]);
        m.w_z = std::stod(fields[13]);

        m.f_x = std::stod(fields[15]);
        m.f_y = std::stod(fields[16]);
        m.f_z = std::stod(fields[17]);

        measurements.push_back(m);
      } catch (const std::exception& e) {
        // Silently catch or handle parsing exceptions
      }
    }
  }

  file.close();
  return measurements;
}

}  // namespace CubicBasisSpline
