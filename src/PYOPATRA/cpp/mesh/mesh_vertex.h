//
// Created by georgia on 3/10/21.
//

#ifndef PYOPATRA_MESH_VERTEX_H
#define PYOPATRA_MESH_VERTEX_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "../coordinate.h"
#include "velocity.h"

class MeshVertex {
private:
    Coordinate3D location;
    Velocity velocity;
    ArrayXXd density, temperature, water_viscosity, viscosity;

public:
    MeshVertex();
    MeshVertex(double latitude, double longitude);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_depth_layers);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_depth_layers, int num_time_steps);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, Eigen::Ref<const Eigen::ArrayXXd>& density, Eigen::Ref<const Eigen::ArrayXXd>& temperature);

    double get_latitude() const { return location[0]; }
    double get_longitude() const { return location[1]; }
    double get_bathymetric_depth() const { return location[2]; }
    const Eigen::ArrayXXd& get_density() const { return density; }
    const Eigen::ArrayXXd& get_temperature() const { return temperature; }
    const Eigen::ArrayXXd& get_water_viscosity() const { return water_viscosity; }
    const Eigen::ArrayXXd& get_viscosity() const { return viscosity; }
    const Eigen::Tensor<double, 3>& get_velocity() const { return velocity.get_velocity(); }
    Coordinate3D get_velocity(size_t depth_index, size_t time_index) const { return velocity.get_velocity(depth_index, time_index); }

    void set_temperature(Eigen::Ref<const Eigen::ArrayXXd>& temp) { temperature = temp; }
    void set_density(Eigen::Ref<const Eigen::ArrayXXd>& dense) { density = dense; }

    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};

#endif //PYOPATRA_MESH_VERTEX_H
