# continuous_fg_ins

Continuous Time INS Factor graph based on Cubic Basis Splines. This is my attempt at using GTSAM develop banches new CUDA capability to estimate parameters for a Cubic Basis Spline around a closed circuit using an IMU and multipl Novatek Pwr-Pak 7D receivers as aiding sources. Much of this work will be based off of a recent ION Navigation Journal: [Spline-Based Factor Graph Estimation with High-Grade Inertial Sensors](https://navi.ion.org/content/73/1/navi.742). This paper gives a formulation for how to parameterize a vehicle's trajectory as a set of Cubic Basis Splines (Cubic B Splines) parameterized by a set of **Control Points**. I will define what these control points are in following sections. This formulation has several advantages over discrete-time state estimation techniques. Primarily, by describing the problem in continuous time, we are allowed to apply constraints at an arbitrary point in time. This helps alleviate some of the need for systems like PTP timing synchronization to enable multi-sensor fusion into a single navigator. This approach is based off of factor graphs and in particular a Cubic B Spline introduces some ***non-sparsity*** that would make this more difficult to solve in a real-time manner than something like an EKF or fixed lag smoother. As such I am implementing this in c++ and not initially requiring super-real time solve times. My hope is that eventually we can make this solution super-real time, making that happen comes later. Now I introduce what the cubic B Splines are, and how we can use them to introduce IMU constraints.

# Cubic Basis Splines

## Cumulative Splines

# Dependencies

As solving the factor graph is outside the scope of this work, we will use some external tools to actually solve the optimization problem. Initially, I was going to solve with GTSAM, as their ISAM2 tool is good and it is very useful for a variety of slam problems. However, I might try out ceres as it is a more general optimization tool and this is all pretty custom anyways. As such the following is the list of external dependencies for using ceres to do the CCBS Optimization piece:

To use Ceres autodiff functionality, I guess it is easier if the SO(3) lie group operations are done w/ Sophus as there are some restrictions on the "setup" to get Ceres autodiff to work properly. Regardless it might just be easier to use Sophus as opposed to my own so(3) tools. After all, I don't really feel the need to do the entire factor graph by hand.

SuiteSparse: git clone https://github.com/DrTimothyAldenDavis/SuiteSparse.git on branch v7.5.0 Build options: cmake ..     -DCMAKE_BUILD_TYPE=Releasee     -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.6     -DSUITESPARSE_USE_CUDA=ON     -DBUILD_SHARED_LIBS=ON

Ceres: git clone https://github.com/ceres-solver/ceres-solver.git on branch 2.2.0. Build options: cmake ..   -DCMAKE_BUILD_TYPE=Release   -DCMAKE_INSTALL_PREFIX=/usr/local   -DCMAKE_PREFIX_PATH=/usr/local   -DSUITESPARSE=ON   -DBUILD_EXAMPLES=OFF   -DBUILD_TESTING=OFF   -DBUILD_SHARED_LIBS=ON

Sophus: git clone https://github.com/strasdat/Sophus.git on main
Build options: cmake .. -DCMAKE_BUILD_TYPE=Release

all of these should install to /usr/local with >> sudo make install 
all should be sudo make install (ed)