#include "spline_tools/spline_tools.hpp"
#include <iostream>

int main() {
    std::cout << "Starting Application..." << std::endl;
    CubicBasisSplines::run_optimization();
    return 0;
}