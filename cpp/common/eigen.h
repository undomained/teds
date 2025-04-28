// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include <Eigen/Dense>

// This array is most often used through-out TEDS. Row-major access is
// optimal in many cases and, most importantly, this makes it easy to
// read and write NetCDF variables because NetCDF is row-major (by
// default at least). For linear algebra operations, where needed,
// we'll make use of the normal Eigen::Matrix class.
using ArrayXXd =
  Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <int N>
using ArrayXNd = Eigen::Array<double, Eigen::Dynamic, N, Eigen::RowMajor>;

using ArrayXb = Eigen::Array<bool, Eigen::Dynamic, 1>;
