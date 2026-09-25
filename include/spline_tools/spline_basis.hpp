#ifndef SPLINE_BASIS_HPP
#define SPLINE_BASIS_HPP

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "iostream"

namespace CubicBasisSpline {

struct SplineWeights {
  double Btilde[4];
  double dBtilde[4];
  double ddBtilde[4];
};

SplineWeights calculateSplineWeights(double u);

}  // namespace CubicBasisSpline

#endif