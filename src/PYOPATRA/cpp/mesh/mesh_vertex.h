//
// Created by georgia on 3/10/21.
//

#ifndef PYOPATRA_MESH_VERTEX_H
#define PYOPATRA_MESH_VERTEX_H

#include <Eigen/Dense>

#include "../coordinate.h"
#include "velocity.h"

class MeshVertex {
private:
    Coordinate3D location;
    Vector3d velocity;
    double density, temperature, water_viscosity, viscosity;

public:
    MeshVertex();
    MeshVertex(double latitude, double longitude);
    MeshVertex(double latitude, double longitude, double bathymetric_depth);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature);

    Coordinate3D get_location() const { return location; }
    double get_latitude() const { return location[0]; }
    double get_longitude() const { return location[1]; }
    double get_depth() const { return location[2]; }
    double get_density() const { return density; }
    double get_temperature() const { return temperature; }
    double get_water_viscosity() const { return water_viscosity; }
    double get_viscosity() const { return viscosity; }
    Vector3d get_velocity() const { return velocity; }

    void set_temperature(double new_temperature);
    void set_density(double new_density);

    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};

bool operator==(const MeshVertex& lhs, const MeshVertex& rhs);

#endif //PYOPATRA_MESH_VERTEX_H
