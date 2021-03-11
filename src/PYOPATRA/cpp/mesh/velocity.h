//
// Created by Georgia Stuart on 3/11/21.
//

#ifndef PYOPATRA_VELOCITY_H
#define PYOPATRA_VELOCITY_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "../coordinate.h"

class Velocity {
private:
    Eigen::Tensor<double, 3> velocity;

public:
    Velocity() : velocity(1, 3, 1) {}
    Velocity(size_t num_depth_steps, size_t num_time_steps) : velocity(num_depth_steps, 3, num_time_steps) {}


    const Eigen::Tensor<double, 3>& get_velocity() const { return velocity; }

    Coordinate3D get_velocity(size_t depth_index, size_t time_index) const {
        return {velocity(depth_index, 0, time_index), velocity(depth_index, 1, time_index), velocity(depth_index, 2, time_index)};
    }

    void set_velocity(size_t depth_index, size_t time_index, double latitude, double longitude, double vertical) {
        velocity(depth_index, 0, time_index) = latitude;
        velocity(depth_index, 1, time_index) = longitude;
        velocity(depth_index, 2, time_index) = vertical;
    }
};

#endif //PYOPATRA_VELOCITY_H
