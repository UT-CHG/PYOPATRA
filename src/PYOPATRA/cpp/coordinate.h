//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_COORDINATE_H
#define PYOPATRA_COORDINATE_H

#include <Eigen/Dense>
#include <cmath>

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

// Convert meters to degrees latitude and longitude. Not accurate near the poles.
inline double meters_to_latitude(double meters) { return meters / 111111.0; }
inline double meters_to_longitude(double meters, double latitude) { return meters / (111111.0 * cos(latitude * M_PI / 180.0)); }

#endif //PYOPATRA_COORDINATE_H
