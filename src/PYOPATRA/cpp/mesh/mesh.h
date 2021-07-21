//
// Created by Georgia Stuart on 3/10/21.
//

#ifndef PYOPATRA_MESH_H
#define PYOPATRA_MESH_H

#include <ctime>
#include "Eigen/Geometry"

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
        for (size_t i = 0; i < water_columns.size(); i++) {
            water_columns[i].set_index(i);
        }
    }

    time_t get_current_time() { return current_time; }
    std::vector<Vertex>& get_vertices() { return vertices; }
    Eigen::MatrixXd get_vertex_locations() {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(vertices.size(), dimension);

        for (size_t index = 0; index < vertices.size(); index++) {
            temp.row(index) = vertices[index].get_location();
        }
        return temp;
    }
    Eigen::MatrixXd get_velocities(int time_index) {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(vertices.size(), dimension);
        for (size_t index = 0; index < vertices.size(); index++) {
            temp.row(index) = vertices[index].get_velocity()[time_index];
        }
        return temp;
    }
    std::vector<WaterCol>& get_water_columns() { return water_columns; }
    void set_vertex_location(int vertex_index, Vector new_location) { vertices[vertex_index].set_location(new_location); }
    void set_vertex_velocity(int vertex_index, int time_index, Vector new_velocity) { vertices[vertex_index].set_velocity(new_velocity, time_index); }
    void set_vertex_diffusion(int vertex_index, int time_index, Vector new_diffusion) { vertices[vertex_index].set_diffusion_coefficient(new_diffusion, time_index); }
    void set_water_column_adjacency(int water_column_index, int adjacent_index, int position) {
        water_columns[water_column_index].set_adjacent_column(&water_columns[adjacent_index], position);
    }
    void set_element_vertex(int water_column_index, int element_depth_index, int position, int vertex_index) {
        water_columns[water_column_index].set_element_vertex(element_depth_index, position, &vertices[vertex_index]);
    }
    bool check_water_column_adjacency(int origin_index, int destination_index, int side) {
        return water_columns[origin_index].get_adjacencies()[side] == &water_columns[destination_index];
    }
    bool check_mesh_element_vertex(int water_column_index, int element_index, int vertex_index, int position) {
        return water_columns[water_column_index].get_mesh_elements()[element_index].get_vertices()[position] == &vertices[vertex_index];
    }
    size_t get_water_columns_size() { return water_columns.size(); }
    Vertex* get_vertex_pointer(int vertex_index) { return &vertices[vertex_index]; }
    const WaterCol* get_water_column_pointer(int water_column_index) const { return &water_columns[water_column_index]; }
    const std::array<WaterCol*, num_vertices_per_element>& get_water_column_adjacencies(int water_column_index) const { return  water_columns[water_column_index].get_adjacencies(); }
    void add_particle(Vector& location) {
        particles.create_particle(location);
        update_particle_mesh_location(*particles.get_tail());
    }
    Eigen::MatrixXd get_all_particle_locations() const { return particles.get_all_particle_locations(); };
    Eigen::VectorXi get_all_particle_column_indices() const { return particles.get_all_particle_column_indices(); }

    void update_particle_locations() {
        auto current = particles.get_head();

        while (current) {
            update_particle_mesh_location(*current);
            current = current->get_next();
        }
    }

    void update_particle_mesh_location(typename ParticleList<dimension>::ParticleN& particle) {
        auto last_known_water_column = water_columns[particle.get_last_known_water_column_index()];

        particle.set_water_column_index(locate_new_water_column(last_known_water_column, particle.get_location()));
    }

    // From https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    size_t locate_new_water_column(const WaterCol& starting_water_col, const Vector& location) {
        if constexpr (dimension == 2 && num_vertices_per_element == 3) {
            Vector A = starting_water_col.get_mesh_elements()[0].get_vertices()[0]->get_location();
            Vector B = starting_water_col.get_mesh_elements()[0].get_vertices()[1]->get_location();
            Vector C = starting_water_col.get_mesh_elements()[0].get_vertices()[2]->get_location();

            int side_1 = sgn((B[0] - A[0]) * (location[1] - A[1]) - (B[1] - A[1]) * (location[0] - A[0]));
            std::cout << side_1 << std::endl;
            int side_2 = sgn((C[0] - B[0]) * (location[1] - B[1]) - (C[1] - B[1]) * (location[0] - B[0]));
            std::cout << side_2 << std::endl;
            int side_3 = sgn((A[0] - C[0]) * (location[1] - C[1]) - (A[1] - C[1]) * (location[0] - C[0]));
            std::cout << side_3 << std::endl;

            if (side_1 <= 0 && side_2 <= 0 && side_3 <=0) {
                return starting_water_col.get_index();
            } else if (side_1 > 0) {
                return locate_new_water_column(*starting_water_col.get_adjacencies()[0], location);
            } else if (side_2 > 0) {
                return locate_new_water_column(*starting_water_col.get_adjacencies()[1], location);
            } else {
                return locate_new_water_column(*starting_water_col.get_adjacencies()[2], location);
            }
        }

        return 0;
    }

};

using TriangularMesh2D = Mesh<3, 2>;
using TriangularMesh3D = Mesh<3, 3>;

#endif //PYOPATRA_MESH_H
