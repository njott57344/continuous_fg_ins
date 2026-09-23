#ifndef CUBIC_BASIS_SPLINE_LIB_HPP
#define CUBIC_BASIS_SPLINE_LIB_HPP

#include <Eigen/Core>
#include <ceres/ceres.h>

namespace CubicBasisSplines {

// Simple Ceres cost functor using Eigen types
struct CostFunctor {
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        residual[0] = x[0] - T(10.0);
        return true;
    }
};

void run_optimization();

} 

#endif