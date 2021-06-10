//
// Created by Georgia Stuart on 3/10/21.
//

#ifndef PYOPATRA_MESH_H
#define PYOPATRA_MESH_H

#include <ctime>

#include "mesh_element.h"
#include "mesh_vertex.h"
#include "mesh_water_column.h"
#include "../particle_list.h"

template <int num_vertices_per_element, int dimension>
class Mesh {
public:
    using WaterCol = WaterColumn<num_vertices_per_element, dimension>;
    using Vertex = MeshVertex<dimension>;

protected:
    time_t current_time;
    int current_time_step, total_time_steps, time_step_size;
    std::vector<time_t> measured_times;
    std::vector<WaterCol> water_columns;
    std::vector<Vertex> vertices;
    ParticleList<dimension> particles;

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

    Mesh(int num_water_columns, int num_vertices, std::vector<time_t>&& measured_times)
        : current_time(0)
        , current_time_step(0)
        , total_time_steps(0)
        , time_step_size(0)
        , measured_times(measured_times)
        , water_columns(num_water_columns, WaterCol())
        , vertices(num_vertices, Vertex(measured_times.size()))
    {
        std::cout << water_columns.size() << std::endl;
    }

    time_t get_current_time() { return current_time; }
    std::vector<Vertex>& get_vertices() { return vertices; }
    std::vector<WaterCol>& get_water_columns() { return water_columns; }
    void set_vertex_location(int vertex_index, Vector new_location) { vertices[vertex_index].set_location(new_location); }
    void set_vertex_velocity(int vertex_index, int time_index, Vector new_velocity) { vertices[vertex_index].set_velocity(new_velocity, time_index); }
    void set_vertex_diffusion(int vertex_index, int time_index, Vector new_diffusion) { vertices[vertex_index].set_diffusion_coefficient(new_diffusion, time_index); }
    void set_water_column_adjacency(int water_column_index, int adjacent_index, int position) { water_columns[water_column_index].set_adjacent_column(&water_columns[adjacent_index], position); }
    void set_element_vertex(int water_column_index, int element_depth_index, int position, int vertex_index) {
        water_columns[water_column_index].set_element_vertex(element_depth_index, position, &vertices[vertex_index]);
    }

    size_t get_water_columns_size() { return water_columns.size(); }
    Vertex* get_vertex_pointer(int vertex_index) { return &vertices[vertex_index]; }
    WaterCol* get_water_column_pointer(int water_column_index) { return &water_columns[water_column_index]; }
    void add_particle() { particles.push(); }
};

using TriangularMesh2D = Mesh<3, 2>;
using TriangularMesh3D = Mesh<3, 3>;

#endif //PYOPATRA_MESH_H
