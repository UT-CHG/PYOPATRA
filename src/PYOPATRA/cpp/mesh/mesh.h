//
// Created by Georgia Stuart on 3/10/21.
//

#ifndef PYOPATRA_MESH_H
#define PYOPATRA_MESH_H

#include <ctime>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mesh_element.h"
#include "mesh_vertex.h"
#include "mesh_water_column.h"
#include "../particle.h"

namespace py = pybind11;

template <int num_vertices_per_element, int dimension>
class MeshFrame {
protected:
    time_t current_time;
    int current_time_step, total_time_steps, time_step_size;
    std::vector<time_t> measured_times;
    std::vector<WaterColumn<num_vertices_per_element, dimension>> water_columns;
    std::vector<MeshVertexBase<dimension>> vertices;

public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

    Mesh()
        : current_time(0)
        , current_time_step(0)
        , total_time_steps(0)
        , time_step_size(0)
        , measured_times()
        , water_columns()
        , vertices()
    {}

    Mesh(int total_time_steps, int time_step_size, int num_water_columns, int num_vertices, std::vector<time_t>&& measured_times)
        : current_time(0)
        , current_time_step(0)
        , total_time_steps(total_time_steps)
        , time_step_size(time_step_size)
        , measured_times(measured_times)
        , water_columns(num_water_columns, WaterColumn<num_vertices_per_element, dimension>())
        , vertices(num_vertices, MeshVertexBase<dimension>(measured_times.size()))
    {}

    WaterColumn<num_vertices_per_element, dimension>* find_particle_location(ParticleBase<dimension> &particle);
    time_t get_current_time() { return current_time; }

    void set_vertices(const py::list& new_vertices) { vertices = new_vertices; }
};

using TriangularMesh2D = Mesh<3, 2>;
using TriangularMesh3D = Mesh<3, 3>;

#endif //PYOPATRA_MESH_H
