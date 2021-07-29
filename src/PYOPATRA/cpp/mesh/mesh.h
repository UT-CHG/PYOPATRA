//
// Created by Georgia Stuart on 3/10/21.
//

#ifndef PYOPATRA_MESH_H
#define PYOPATRA_MESH_H

#include <string>
#include <ctime>
#include "mpi.h"
#include "mesh_element.h"
#include "mesh_vertex.h"
#include "mesh_water_column.h"
#include "../particle_list.h"
#include "../inversion_tools/objective_functions.h"
#include "../util.h"

template <int num_vertices_per_element, int dimension>
class Mesh {
public:
    using WaterCol = WaterColumn<num_vertices_per_element, dimension>;
    using Vertex = MeshVertex<dimension>;
    using MeshElement = MeshElementT<num_vertices_per_element, dimension>;

protected:
    WaterCol *water_columns;
    Vertex *vertices;
    MeshElement *elements;
    MPI_Win columns_win, vertices_win, elements_win;
    PointerWrapper<Mesh<num_vertices_per_element, dimension>> ptr_wrapper;
    int rank, num_water_columns, num_vertices, num_mesh_elements, num_depths;
    MPI_Comm node_comm;

public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

    Mesh()
        : water_columns(nullptr)
        , vertices(nullptr)
        , elements(nullptr)
        , num_water_columns(0)
        , num_vertices(0)
        , num_mesh_elements(0)
        , num_depths(0)
    {
        ptr_wrapper.set_pointer(this);
    }

    Mesh(int num_water_columns, int num_vertices, std::vector<double>&& measured_times)
        : num_water_columns(num_water_columns)
        , num_vertices(num_vertices)
        , num_mesh_elements(num_water_columns)
        , num_depths(1)
    {
        int world_size;

        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
        MPI_Comm_rank(node_comm, &rank);
        MPI_Comm_size(node_comm, &world_size);


        MPI_Aint size;

        if (rank == 0) {
            size = num_vertices * sizeof(Vertex);
            MPI_Win_allocate_shared(size, sizeof(Vertex), MPI_INFO_NULL, node_comm, &vertices, &vertices_win);

            for (int i = 0; i < num_vertices; i++) {
                new (vertices + i) Vertex(measured_times.size());
            }

            size = num_water_columns * sizeof(MeshElement);
            MPI_Win_allocate_shared(size, sizeof(MeshElement), MPI_INFO_NULL, node_comm, &elements, &elements_win);

            for (int i = 0; i < num_vertices; i++) {
                new (elements + i) MeshElement();
            }

            size = num_water_columns * sizeof(WaterCol);
            MPI_Win_allocate_shared(size, sizeof(WaterCol), MPI_INFO_NULL, node_comm, &water_columns, &columns_win);

            for (int i = 0; i < num_water_columns; i++) {
                new (water_columns + i) WaterCol();
                water_columns[i].set_index(i);
                water_columns[i].set_element_head(&elements[i * num_depths]);
            }


        } else {
            MPI_Aint ssize;
            int disp_unit;

            size = 0;

            MPI_Win_allocate_shared(size, sizeof(WaterCol), MPI_INFO_NULL, node_comm, &water_columns, &columns_win);
            MPI_Win_allocate_shared(size, sizeof(Vertex), MPI_INFO_NULL, node_comm, &vertices, &vertices_win);
            MPI_Win_allocate_shared(size, sizeof(MeshElement), MPI_INFO_NULL, node_comm, &elements, &elements_win);

            MPI_Win_shared_query(columns_win, 0, &ssize, &disp_unit, &water_columns);
            MPI_Win_shared_query(vertices_win, 0, &ssize, &disp_unit, &vertices);
            MPI_Win_shared_query(elements_win, 0, &ssize, &disp_unit, &elements);
        }

        MPI_Barrier(node_comm);

        ptr_wrapper.set_pointer(this);
    }

    ~Mesh() {
        MPI_Barrier(node_comm);

        if (rank == 0) {
            for (int i = 0; i < num_water_columns; i++) {
                water_columns[i].~WaterCol();
            }

            for (int i = 0; i < num_vertices; i++) {
                vertices[i].~Vertex();
            }

            for (int i = 0; i < num_mesh_elements; i++) {
                elements[i].~MeshElement();
            }
        }

        MPI_Barrier(node_comm);
        MPI_Win_free(&columns_win);
        MPI_Win_free(&vertices_win);
        MPI_Win_free(&elements_win);
    }

    PointerWrapper<Mesh<num_vertices_per_element, dimension>> get_pointer_wrapper() { return ptr_wrapper; }
    Vertex* get_vertices() { return vertices; }
    int get_num_vertices() { return num_vertices; }
    int get_num_columns() { return num_water_columns; }
    int get_num_elements() { return num_mesh_elements; }
    Eigen::MatrixXd get_vertex_locations() {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(num_vertices, dimension);

        for (int index = 0; index < num_vertices; index++) {
            temp.row(index) = vertices[index].get_location();
        }
        return temp;
    }
    Eigen::MatrixXd get_velocities(int time_index) {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(num_vertices, dimension);
        for (int index = 0; index < num_vertices; index++) {
            temp.row(index) = vertices[index].get_velocity()[time_index];
        }
        return temp;
    }
    WaterCol* get_water_columns() { return water_columns; }
    void set_vertex_location(int vertex_index, Vector new_location) { vertices[vertex_index].set_location(new_location); }
    void set_vertex_velocity(int vertex_index, int time_index, Vector new_velocity) { vertices[vertex_index].set_velocity(new_velocity, time_index); }
    void set_vertex_diffusion(int vertex_index, int time_index, Vector new_diffusion) { vertices[vertex_index].set_diffusion_coefficient(new_diffusion, time_index); }
    void set_water_column_adjacency(int water_column_index, int adjacent_index, int position) {
        water_columns[water_column_index].set_adjacent_column(&water_columns[adjacent_index], position);
    }
    void set_element_vertex(int water_column_index, int element_depth_index, int position, int vertex_index) {
        elements[water_column_index * num_depths + element_depth_index].set_vertex(&vertices[vertex_index], position);
    }
    bool check_water_column_adjacency(int origin_index, int destination_index, int side) {
        return water_columns[origin_index].get_adjacencies()[side] == &water_columns[destination_index];
    }
    bool check_mesh_element_vertex(int water_column_index, int element_index, int vertex_index, int position) {
        return elements[water_column_index * num_depths + element_index].get_vertices()[position] == &vertices[vertex_index];
    }
    int get_water_columns_size() { return num_water_columns; }
    Vertex* get_vertex_pointer(int vertex_index) { return &vertices[vertex_index]; }
    const WaterCol* get_water_column_pointer(int water_column_index) const { return &water_columns[water_column_index]; }
    const std::array<WaterCol*, num_vertices_per_element>& get_water_column_adjacencies(int water_column_index) const { return  water_columns[water_column_index].get_adjacencies(); }

    // From https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    int locate_new_water_column(const WaterCol* starting_water_col, const Vector& location) {
        if constexpr (dimension == 2 && num_vertices_per_element == 3) {
            if (starting_water_col) {
                Vector A = starting_water_col->get_mesh_elements()[0].get_vertices()[0]->get_location();
                Vector B = starting_water_col->get_mesh_elements()[0].get_vertices()[1]->get_location();
                Vector C = starting_water_col->get_mesh_elements()[0].get_vertices()[2]->get_location();

                int side_1 = sgn((B[0] - A[0]) * (location[1] - A[1]) - (B[1] - A[1]) * (location[0] - A[0]));
                int side_2 = sgn((C[0] - B[0]) * (location[1] - B[1]) - (C[1] - B[1]) * (location[0] - B[0]));
                int side_3 = sgn((A[0] - C[0]) * (location[1] - C[1]) - (A[1] - C[1]) * (location[0] - C[0]));

                if (side_1 <= 0 && side_2 <= 0 && side_3 <= 0) {
                    return starting_water_col->get_index();
                } else if (side_1 > 0) {
                    return locate_new_water_column(starting_water_col->get_adjacencies()[0], location);
                } else if (side_2 > 0) {
                    return locate_new_water_column(starting_water_col->get_adjacencies()[1], location);
                } else {
                    return locate_new_water_column(starting_water_col->get_adjacencies()[2], location);
                }
            } else {
                return -1;
            }
        }

        return -1;
    }
};

using TriangularMesh2D = Mesh<3, 2>;
using TriangularMesh3D = Mesh<3, 3>;

using TriangularMesh2DPtrWrapper = PointerWrapper<TriangularMesh2D>;

#endif //PYOPATRA_MESH_H
