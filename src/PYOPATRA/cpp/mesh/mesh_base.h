//
// Created by Georgia Stuart on 3/10/21.
//

#ifndef PYOPATRA_MESH_BASE_H
#define PYOPATRA_MESH_BASE_H

#include <ctime>

#include "mesh_element.h"
#include "mesh_vertex.h"
#include "mesh_water_column.h"
#include "../particle.h"

template <int num_vertices_per_element, int dimension>
class MeshBase {
protected:
    time_t current_time;
    int current_time_step, total_time_steps, time_step_size;
    std::vector<time_t> measured_times;
    std::vector<std::vector<WaterColumn<num_vertices_per_element, dimension>>> water_columns;
    std::vector<MeshVertex<dimension>> vertices;

public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

    MeshBase()
        : current_time(0)
        , current_time_step(0)
        , total_time_steps(0)
        , time_step_size(0)
        , measured_times()
        , water_columns()
        , vertices()
    {}

    MeshBase(int total_time_steps, int time_step_size, int num_water_columns, int total_vertices, std::vector<time_t>&& measured_times)
        : current_time(0)
        , current_time_step(0)
        , total_time_steps(total_time_steps)
        , time_step_size(time_step_size)
        , measured_times(measured_times)
        , water_columns(num_water_columns, std::vector<WaterColumn<num_vertices_per_element, dimension>>(measured_times.size(), WaterColumn<num_vertices_per_element, dimension>()))
        , vertices(total_vertices, MeshVertex<dimension>())
    {}

    virtual ~MeshBase() = default;
    virtual WaterColumn<num_vertices_per_element, dimension>* find_particle_location(ParticleBase<dimension> &particle) = 0;
    time_t get_current_time() { return current_time; }

    void setup_vertices(std::vector<double>& latitude, std::vector<double>& longitude, std::vector<Vector>& velocity, std::vector<Vector>& diffusion_coefficient);
};

using TriangularMesh2D = MeshBase<3, 2>;
using TriangularMesh3D = MeshBase<3, 3>;

#endif //PYOPATRA_MESH_BASE_H
