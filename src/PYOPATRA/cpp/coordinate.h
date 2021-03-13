//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_COORDINATE_H
#define PYOPATRA_COORDINATE_H

#include <Eigen/Dense>

using Eigen::Array3d;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Vector4d;
using Eigen::ArrayXXd;

// Latitude = 0 index
// Longitude = 1 index
// Depth = 2 index
using Coordinate3D = Vector3d;

// Latitude = 0 index
// Longitude = 1 index
using Coordinate2D = Vector2d;

#endif //PYOPATRA_COORDINATE_H
