//
// Created by georgia on 3/10/21.
//

#ifndef PYOPATRA_MESH_VERTEX_H
#define PYOPATRA_MESH_VERTEX_H

#include <Eigen/Dense>

class MeshVertex {
private:
    const double latitude, longitude, bathymetric_depth;
    Eigen::ArrayXXd density, temperature, water_viscosity, viscosity;

public:
    MeshVertex();
    MeshVertex(double latitude, double longitude);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_depth_layers);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_depth_layers, int num_time_steps);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, Eigen::Ref<const Eigen::ArrayXXd>& density, Eigen::Ref<const Eigen::ArrayXXd>& temperature);

    double get_latitude() const { return latitude; }
    double get_longitude() const { return longitude; }
    double get_bathymetric_depth() const { return bathymetric_depth; }
    const Eigen::ArrayXXd& get_density() const { return density; }
    const Eigen::ArrayXXd& get_temperature() const { return temperature; }
    const Eigen::ArrayXXd& get_water_viscosity() const { return water_viscosity; }
    const Eigen::ArrayXXd& get_viscosity() const { return viscosity; }

    void set_temperature(Eigen::Ref<const Eigen::ArrayXXd>& temp) { temperature = temp; }
    void set_density(Eigen::Ref<const Eigen::ArrayXXd>& dense) { density = dense; }

    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};

#endif //PYOPATRA_MESH_VERTEX_H
