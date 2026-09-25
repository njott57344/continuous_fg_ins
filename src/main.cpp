#include <iostream>

#include "spline_tools/spline_eval.hpp"

/*
TODO:

1) Need to make a parser for ATR INS Data
    what is most principled way to setup a csv parser?
2) Need to calculate the number of control points
    even spaced at a tuned dt
3) Need to parse ecef positions,velocities from novatel data
4) Need to build the ceres problem
5) Need to solve (hopefully it just goes)

Unknowns:
1) How to setup automatic differentiation
(so I don't have to manually plug in jacobians)

2) I think the evaluation is correct, but I don't know for sure

*/

int main() { return 0; }