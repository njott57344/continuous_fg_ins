#ifndef SO3_TOOLS_HPP
#define SO3_TOOLS_HPP

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "cmath"

namespace CubicBasisSplines {

Eigen::Matrix3d skew(const Eigen::Vector3d& vec);

Eigen::Vector3d vee(const Eigen::Matrix3d& mat);

Eigen::Matrix3d expm(const Eigen::Vector3d& vec);

Eigen::Matrix3d logm(const Eigen::Matrix3d& mat);

}  // namespace CubicBasisSplines

#endif