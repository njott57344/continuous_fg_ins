#ifndef CUBIC_BASIS_SPLINE_LIB_HPP
#define CUBIC_BASIS_SPLINE_LIB_HPP

#include <ceres/ceres.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

namespace CubicBasisSplines {

const Eigen::Matrix4d C{
    {6.0, 0.0, 0.0, 0.0}, {5.0, 3.0, -3.0, 1.0}, {1.0, 3.0, 3.0, -2.0}, {0.0, 0.0, 0.0, 1.0}};

// struct for a single control point
// p->R3 position ctrl pnt
// R->SO(3) orientation ctrl pnt
// t->time for that control point
struct CtrlPt {
  Eigen::Vector3d p;
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  double t;
};

struct Weights {
  Eigen::Vector4d Btilde;
  Eigen::Vector4d dBtilde;
  Eigen::Vector4d ddBtilde;
};

struct SplineEval {
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
  Eigen::Vector3d acc;

  Eigen::Matrix3d C_bw;
  Eigen::Vector3d w_bw;
};

double calculateNormalizedTime(const double& t_i, const double& t_eval, const double& dt);

Weights calculateWeights(const double& u);

bool evaluateTranslationSpline(const std::vector<CtrlPt>& ctrl_pts, double& t_eval,
                               SplineEval& eval_spline);

bool evaluateRotationSpline(const std::vector<CtrlPt>& ctrl_pts, double& t_eval,
                            SplineEval& eval_spline);

bool evaluateSplitSpline(const std::vector<CtrlPt>& ctrl_pts, double& t_eval,
                         SplineEval& eval_spline);

}  // namespace CubicBasisSplines

#endif