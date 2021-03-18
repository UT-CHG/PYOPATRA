//
// Created by georgia on 3/10/21.
//

#ifndef PYOPATRA_MESH_VERTEX_H
#define PYOPATRA_MESH_VERTEX_H

#include <Eigen/Dense>
#include <initializer_list>

#include "../coordinate.h"

template <int dimension>
class MeshVertexBase {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

protected:
    Vector location;
    Vector velocity;

public:
    MeshVertexBase();
    MeshVertexBase(double latitude, double longitude, std::initializer_list<double> velocity);

    Vector get_location() const { return location; }
    double get_latitude() const { return location[0]; }
    double get_longitude() const { return location[1]; }
    Vector get_velocity() const { return velocity; }
};

template<int dimension>
bool operator==(const MeshVertexBase<dimension>& lhs, const MeshVertexBase<dimension>& rhs);

template <int dimension>
class MeshVertex : MeshVertexBase<dimension> {};

template <>
class MeshVertex<3> : public MeshVertexBase<3> {
private:
    double density, temperature, water_viscosity, viscosity;

public:
    MeshVertex();
    MeshVertex(double latitude, double longitude, double bathymetric_depth);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature);

    double get_density() const { return density; }
    double get_temperature() const { return temperature; }
    double get_water_viscosity() const { return water_viscosity; }
    double get_viscosity() const { return viscosity; }
    double get_depth() const { return location[2]; }

    void set_temperature(double new_temperature);
    void set_density(double new_density);

    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};


typedef MeshVertex<3> MeshVertex3D;
typedef MeshVertex<2> MeshVertex2D;

#include "mesh_vertex.inl"

#endif //PYOPATRA_MESH_VERTEX_H
