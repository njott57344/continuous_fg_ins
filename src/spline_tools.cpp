#include "spline_tools/spline_tools.hpp"
#include <iostream>

namespace CubicBasisSplines {

void run_optimization() {
    double initial_x = 5.0;
    double x = initial_x;

    ceres::Problem problem;
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor),
        nullptr,
        &x
    );

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.minimizer_progress_to_stdout = false;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << "  Initial x: " << initial_x << "\n"
              << "  Final x  : " << x << "\n";
}

}