//
// Created by georgia on 3/10/21.
//

#ifndef PYOPATRA_MESH_VERTEX_H
#define PYOPATRA_MESH_VERTEX_H

#include <Eigen/Dense>
#include <initializer_list>
#include <array>

#include "../coordinate.h"

template <int dimension>
class MeshVertexBase {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

protected:
    Vector location;
    int velocity;
    int diffusion_coefficient;

public:
    explicit MeshVertexBase(int index, int num_timesteps);
    MeshVertexBase(double latitude, double longitude, int index, int num_timesteps);
//    MeshVertexBase(double latitude, double longitude, Vector velocity);
//    MeshVertexBase(double latitude, double longitude, Vector velocity, Vector diffusion_coefficient);

    Vector get_location() const { return location; }
    double get_latitude() const { return location[0]; }
    double get_longitude() const { return location[1]; }
    int get_velocity() const { return velocity; }
    int get_diffusion_coefficient() const { return diffusion_coefficient; }

    void set_diffusion_coefficient(int new_diffusion_coefficient) { diffusion_coefficient = new_diffusion_coefficient; }
    void set_velocity(int new_velocity) { velocity = new_velocity; }
    void set_location(Vector new_location) { location = new_location; }
};

template<int dimension>
bool operator==(const MeshVertexBase<dimension>& lhs, const MeshVertexBase<dimension>& rhs);

template <int dimension>
class MeshVertex : MeshVertexBase<dimension> {};

template <>
class MeshVertex<2> : public MeshVertexBase<2> {
public:
    MeshVertex(int index, int num_timesteps) : MeshVertexBase<2>(index, num_timesteps) {}
    MeshVertex(double latitude, double longitude, int index, int num_timesteps) : MeshVertexBase<2>(latitude, longitude, index, num_timesteps) {}
//    MeshVertex(double latitude, double longitude, Vector velocity, Vector diffusion_coefficient)
//        : MeshVertexBase<2>(latitude, longitude, velocity, diffusion_coefficient)
//    {}
};

template <>
class MeshVertex<3> : public MeshVertexBase<3> {
private:
    std::vector<double> density, temperature, water_viscosity, viscosity;

public:
    explicit MeshVertex(int index, int num_timesteps);
    MeshVertex(double latitude, double longitude, double bathymetric_depth, int index, int num_timesteps);
//    MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature, int num_timesteps);
//    MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature, Vector velocity);
//    MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature, Vector velocity, Vector diffusion_coefficient);

    const std::vector<double>& get_density() const { return density; }
    const std::vector<double>& get_temperature() const { return temperature; }
    const std::vector<double>& get_water_viscosity() const { return water_viscosity; }
    const std::vector<double>& get_viscosity() const { return viscosity; }
    double get_depth() const { return location[2]; }

    void set_temperature_and_density(double new_temperature, double new_density, int time_index) {
        temperature[time_index] = new_temperature;
        density[time_index] = new_density;

        water_viscosity[time_index] = calculate_pure_water_viscosity(new_temperature);
        viscosity[time_index] = calculate_fluid_viscosity(new_temperature, new_density);
    }

    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};


using MeshVertexBase2D = MeshVertexBase<2>;
using MeshVertexBase3D = MeshVertexBase<3>;

using MeshVertex3D = MeshVertex<3> ;
using MeshVertex2D = MeshVertex<2> ;

#include "mesh_vertex.inl"

#endif //PYOPATRA_MESH_VERTEX_H
