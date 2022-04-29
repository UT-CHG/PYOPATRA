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
#include <stdexcept>

template <int num_vertices_per_element, int dimension>
class Mesh {
public:
    using WaterCol = WaterColumn<num_vertices_per_element, dimension>;
    using Vertex = MeshVertex<dimension>;
    using MeshElement = MeshElementT<num_vertices_per_element, dimension>;
    using Vector = Eigen::Matrix<double, dimension, 1>;

protected:
    WaterCol *water_columns; /**< Array of water columns, initialized in MPI shared memory. */
    Vertex *vertices; /**< Array of vertices, initialized in MPI shared memory. */
    MeshElement *elements; /**< Array of mesh elements, initialized in MPI shared memory. */
    Vector *velocities; /**< Array of velocities in time and space, initialized in MPI shared memory. */
    Vector *diffusions; /**< Array of diffusions in time and space, initialized in MPI shared memory. */
    Vector *winds; /**< Array of wind velocities in time and space, initialized in MPI shared memory. */
    MPI_Win columns_win; /**< Water columns MPI window for accessing MPI shared memory. */
    MPI_Win vertices_win; /**< Vertices MPI window for accessing MPI shared memory. */
    MPI_Win elements_win; /**< Mesh element MPI window for accessing MPI shared memory. */
    MPI_Win velocities_win; /**< Velocities MPI window for accessing MPI shared memory. */
    MPI_Win diffusions_win; /**< Diffusions MPI window for accessing MPI shared memory. */
    MPI_Win wind_win; /**< Wind MPI window for accessing MPI shared memory. */
    PointerWrapper<Mesh<num_vertices_per_element, dimension>> ptr_wrapper; /**< Pointer wrapper for passing addresses through the Python layer. */
    size_t  num_water_columns, num_vertices, num_mesh_elements, num_depths, num_time_steps, num_wind_time_steps;
    int rank;
    MPI_Comm node_comm; /**< Node-level communication for accessing and initializing MPI shared memory. */
    double wind_coef; /**< Wind contribution coefficient - mesh-wide. */

public:
    /**
     * The base mesh class for 2D and 3D triangular meshes.
     */
    Mesh()
        : water_columns(nullptr)
        , vertices(nullptr)
        , elements(nullptr)
        , velocities(nullptr)
        , diffusions(nullptr)
        , num_water_columns(0)
        , num_vertices(0)
        , num_mesh_elements(0)
        , num_depths(0)
        , wind_coef(0)
    {
        ptr_wrapper.set_pointer(this);
    }

    Mesh(size_t num_water_columns, size_t num_vertices, std::vector<double>&& measured_times, std::vector<double>&& winds_measured_times, double wind_coef)
        : num_water_columns(num_water_columns)
        , num_vertices(num_vertices)
        , num_mesh_elements(num_water_columns)
        , num_depths(1)
        , num_time_steps(measured_times.size())
        , num_wind_time_steps(winds_measured_times.size())
        , wind_coef(wind_coef)
    {
        int world_size, full_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &full_rank);

        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, full_rank, MPI_INFO_NULL, &node_comm);
        MPI_Comm_rank(node_comm, &rank);
        MPI_Comm_size(node_comm, &world_size);


        MPI_Aint size;
        void *vert_bsptr, *el_bsptr, *col_bsptr, *vel_bsptr, *dif_bsptr, *wind_bsptr;

        if (rank == 0) {
            size = num_vertices * sizeof(Vertex);
            MPI_Win_allocate_shared(size, sizeof(Vertex), MPI_INFO_NULL, node_comm, &vert_bsptr, &vertices_win);

            size = num_mesh_elements * sizeof(MeshElement);
            MPI_Win_allocate_shared(size, sizeof(MeshElement), MPI_INFO_NULL, node_comm, &el_bsptr, &elements_win);

            size = num_water_columns * sizeof(WaterCol);
            MPI_Win_allocate_shared(size, sizeof(WaterCol), MPI_INFO_NULL, node_comm, &col_bsptr, &columns_win);

            size = num_vertices * num_time_steps * sizeof(Vector);
            MPI_Win_allocate_shared(size, sizeof(Vector), MPI_INFO_NULL, node_comm, &vel_bsptr, &velocities_win);
            MPI_Win_allocate_shared(size, sizeof(Vector), MPI_INFO_NULL, node_comm, &dif_bsptr, &diffusions_win);

            size = num_vertices * (num_wind_time_steps == 0 ? 2 : num_wind_time_steps) * sizeof(Vector);
            MPI_Win_allocate_shared(size, sizeof(Vector), MPI_INFO_NULL, node_comm, &wind_bsptr, &wind_win);
        } else {
            size = 0;

            MPI_Win_allocate_shared(size, sizeof(Vertex), MPI_INFO_NULL, node_comm, &vert_bsptr, &vertices_win);
            MPI_Win_allocate_shared(size, sizeof(MeshElement), MPI_INFO_NULL, node_comm, &el_bsptr, &elements_win);
            MPI_Win_allocate_shared(size, sizeof(WaterCol), MPI_INFO_NULL, node_comm, &col_bsptr, &columns_win);
            MPI_Win_allocate_shared(size, sizeof(Vector), MPI_INFO_NULL, node_comm, &vel_bsptr, &velocities_win);
            MPI_Win_allocate_shared(size, sizeof(Vector), MPI_INFO_NULL, node_comm, &dif_bsptr, &diffusions_win);
            MPI_Win_allocate_shared(size, sizeof(Vector), MPI_INFO_NULL, node_comm, &wind_bsptr, &wind_win);
        }

        MPI_Aint ssize;
        int disp_unit;

        MPI_Win_shared_query(vertices_win, MPI_PROC_NULL, &ssize, &disp_unit, &vert_bsptr);
        MPI_Win_shared_query(elements_win, MPI_PROC_NULL, &ssize, &disp_unit, &el_bsptr);
        MPI_Win_shared_query(columns_win, MPI_PROC_NULL, &ssize, &disp_unit, &col_bsptr);
        MPI_Win_shared_query(velocities_win, MPI_PROC_NULL, &ssize, &disp_unit, &vel_bsptr);
        MPI_Win_shared_query(diffusions_win, MPI_PROC_NULL, &ssize, &disp_unit, &dif_bsptr);
        MPI_Win_shared_query(wind_win, MPI_PROC_NULL, &ssize, &disp_unit, &wind_bsptr);

        water_columns = static_cast<WaterCol*>(col_bsptr);
        vertices = static_cast<Vertex*>(vert_bsptr);
        elements = static_cast<MeshElement*>(el_bsptr);
        velocities = static_cast<Vector*>(vel_bsptr);
        diffusions = static_cast<Vector*>(dif_bsptr);
        winds = static_cast<Vector*>(wind_bsptr);

        if (rank == 0) {
            for (size_t i = 0; i < num_vertices; i++) {
                new (vertices + i) Vertex(i, measured_times.size(), winds_measured_times.size());
            }

            for (size_t i = 0; i < num_mesh_elements; i++) {
                new (elements + i) MeshElement();
            }

            for (size_t i = 0; i < num_water_columns; i++) {
                new (water_columns + i) WaterCol();
                water_columns[i].set_index(i);
                water_columns[i].set_element_head(i);
            }

            for (size_t i = 0; i < num_vertices * num_time_steps; i++) {
                new (velocities + i) Vector();
                new (diffusions + i) Vector();
            }

            for (size_t i = 0; i < num_vertices * (num_wind_time_steps == 0 ? 2 : num_wind_time_steps); i++) {
                new (winds + i) Vector();
            }
        }

        MPI_Barrier(node_comm);

        ptr_wrapper.set_pointer(this);
        // Seed the random number generator with the rank to avoid duplicate behavior across processes
        generator.seed(full_rank);
    }

    ~Mesh() {
        MPI_Barrier(node_comm);

        if (rank == 0) {
            for (size_t i = 0; i < num_water_columns; i++) {
                water_columns[i].~WaterCol();
            }

            for (size_t i = 0; i < num_vertices; i++) {
                vertices[i].~Vertex();
            }

            for (size_t i = 0; i < num_mesh_elements; i++) {
                elements[i].~MeshElement();
            }

            for (size_t i = 0; i < num_vertices * num_time_steps; i++) {
                velocities[i].~Vector();
                diffusions[i].~Vector();
            }

            for (size_t i = 0; i < num_vertices * (num_wind_time_steps == 0 ? 2 : num_wind_time_steps); i++) {
                winds[i].~Vector();
            }
        }

        MPI_Barrier(node_comm);
        MPI_Win_free(&columns_win);
        MPI_Win_free(&vertices_win);
        MPI_Win_free(&elements_win);
        MPI_Win_free(&velocities_win);
        MPI_Win_free(&diffusions_win);
        MPI_Win_free(&wind_win);
    }

    PointerWrapper<Mesh<num_vertices_per_element, dimension>> get_pointer_wrapper() { return ptr_wrapper; }
    Vertex* get_vertices() { return vertices; }
    size_t get_num_vertices() { return num_vertices; }
    size_t get_num_columns() { return num_water_columns; }
    size_t get_num_elements() { return num_mesh_elements; }
    void set_wind_coef(double new_wind_coef) { wind_coef = new_wind_coef; }
    double get_wind_coef() { return wind_coef; }
    Eigen::MatrixXd get_vertex_locations() {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(num_vertices, dimension);

        for (size_t index = 0; index < num_vertices; index++) {
            temp.row(index) = vertices[index].get_location();
        }
        return temp;
    }
    Vector* get_velocities_ptr() { return velocities; }
    Vector* get_diffusions_ptr() { return diffusions; }
    Vector* get_winds_ptr() { return winds; }
    Eigen::MatrixXd get_velocities(size_t time_index) {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(num_vertices, dimension);
        for (size_t index = 0; index < num_vertices; index++) {
            temp.row(index) = velocities[vertices[index].get_velocity() + time_index];
        }
        return temp;
    }
    MeshElement* get_elements() { return elements; }
    WaterCol* get_water_columns() { return water_columns; }
    void set_vertex_location(size_t vertex_index, Vector new_location) {
        vertices[vertex_index].set_location(new_location);
    }

    void set_vertex_locations(Eigen::MatrixXd new_locations) {
        size_t num_rows = new_locations.rows();
        size_t num_cols = new_locations.cols();

        if (num_rows != num_vertices || num_cols != dimension) {
            throw std::invalid_argument("Locations must have rows = " + std::to_string(num_vertices) + 
                " and cols = " + std::to_string(dimension));
        }

        for (size_t i = 0; i < num_rows; i++) {
            vertices[i].set_location(new_locations.row(i));
        }
    }

    void set_vertex_velocity(size_t vertex_index, size_t time_index, Vector new_velocity) {
        velocities[vertex_index * num_time_steps + time_index] = new_velocity;
    }

    void set_velocities(Eigen::MatrixXd new_velocities) {
        set_vertices_data(velocities, new_velocities, num_vertices*num_time_steps);
    }

    void set_vertex_diffusion(size_t vertex_index, size_t time_index, Vector new_diffusion) {
        diffusions[vertex_index * num_time_steps + time_index] = new_diffusion;
    }

    void set_diffusions(Eigen::MatrixXd new_diffusions) {
        set_vertices_data(diffusions, new_diffusions, num_vertices*num_time_steps);
    }

    void set_vertex_wind(size_t vertex_index, size_t time_index, Vector new_wind) {
        winds[vertex_index * (num_wind_time_steps == 0 ? 2 : num_wind_time_steps) + time_index] = new_wind;
    }

    void set_winds(Eigen::MatrixXd new_winds) {
        set_vertices_data(winds, new_winds, num_vertices * (num_wind_time_steps == 0 ? 2 : num_wind_time_steps));
    }

    void set_water_column_adjacency(size_t water_column_index, size_t adjacent_index, size_t position) {
        if (rank == 0) {
            water_columns[water_column_index].set_adjacent_column(adjacent_index, position);
        }
    }

    void set_water_column_adjacencies(Eigen::Matrix<long int, Eigen::Dynamic, num_vertices_per_element> new_adjacencies) {
        if (rank != 0) return;

        size_t num_rows = new_adjacencies.rows();
        size_t num_cols = new_adjacencies.cols();

        if (num_rows != num_water_columns) {
            throw std::invalid_argument("Adjacencies must have rows = " + std::to_string(num_water_columns) + 
                " and cols = " + std::to_string(num_vertices_per_element));
        }

        for (size_t i = 0; i < num_rows; i++) {
            for (size_t j = 0; j < num_cols; j++) {
                water_columns[i].set_adjacent_column(new_adjacencies(i,j), j);
            }
        }
    }

    void set_element_vertex(size_t water_column_index, size_t element_depth_index, size_t position, size_t vertex_index) {
        if (rank == 0) {
            elements[water_column_index * num_depths + element_depth_index].set_vertex(vertex_index, position);
        }
    }

    void set_element_vertices(Eigen::Matrix<size_t, Eigen::Dynamic, num_vertices_per_element> new_vertices) {
        if (rank != 0) return;

        size_t num_rows = new_vertices.rows();
        size_t num_cols = new_vertices.cols();

        if (num_rows != num_mesh_elements) {
            throw std::invalid_argument("Vertices must have rows = " + std::to_string(num_water_columns) + 
                " and cols = " + std::to_string(num_vertices_per_element));
        }

        for (size_t i = 0; i < num_rows; i++) {
            for (size_t j = 0; j < num_cols; j++) {
                elements[i].set_vertex(new_vertices(i,j), j);
            }
        }        
    }

    void set_vertices_data(Vector * target, Eigen::MatrixXd source, size_t expected_rows) {
        size_t num_rows = source.rows();
        size_t num_cols = source.cols();

        if (num_rows != expected_rows || num_cols != dimension) {
            throw std::invalid_argument("Input matrix must have rows = " + std::to_string(expected_rows) + 
                "  and cols = " + std::to_string(num_cols));
        }

        for (size_t i = 0; i < num_rows; i++) {
            for (size_t j = 0; j < num_cols; j++) {
                target[i](j) = source(i,j);
            }
        }
    }

    bool check_water_column_adjacency(size_t origin_index, size_t destination_index, int side) {
        return water_columns[origin_index].get_adjacencies()[side] == (long int)destination_index;
    }
    bool check_mesh_element_vertex(size_t water_column_index, size_t element_index, size_t vertex_index, size_t position) {
        return elements[water_column_index * num_depths + element_index].get_vertices()[position] == vertex_index;
    }
    size_t get_water_columns_size() { return num_water_columns; }
    Vertex* get_vertex_pointer(size_t vertex_index) { return &vertices[vertex_index]; }
    const WaterCol* get_water_column_pointer(size_t water_column_index) const { return &water_columns[water_column_index]; }
    const std::array<long int, num_vertices_per_element>& get_water_column_adjacencies(size_t water_column_index) const { return  water_columns[water_column_index].get_adjacencies(); }

    // From https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    long int locate_new_water_column(const WaterCol* starting_water_col, const Vector& location) {
        if constexpr (dimension == 2 && num_vertices_per_element == 3) {
            long int next_col;

            Vector A = vertices[elements[starting_water_col->get_mesh_elements()].get_vertices()[0]].get_location();
            Vector B = vertices[elements[starting_water_col->get_mesh_elements()].get_vertices()[1]].get_location();
            Vector C = vertices[elements[starting_water_col->get_mesh_elements()].get_vertices()[2]].get_location();

            long int side_1 = sgn((B[0] - A[0]) * (location[1] - A[1]) - (B[1] - A[1]) * (location[0] - A[0]));
            long int side_2 = sgn((C[0] - B[0]) * (location[1] - B[1]) - (C[1] - B[1]) * (location[0] - B[0]));
            long int side_3 = sgn((A[0] - C[0]) * (location[1] - C[1]) - (A[1] - C[1]) * (location[0] - C[0]));

            if (side_1 <= 0 && side_2 <= 0 && side_3 <= 0) {
                return starting_water_col->get_index();
            } else if (side_1 > 0) {
                next_col = starting_water_col->get_adjacencies()[0];
            } else if (side_2 > 0) {
                next_col = starting_water_col->get_adjacencies()[1];
            } else {
                next_col = starting_water_col->get_adjacencies()[2];
            }

            if (next_col >= 0) {
                return locate_new_water_column(&water_columns[next_col], location);
            } else {
                return -1;
            }
        }

    }
};

using TriangularMesh2D = Mesh<3, 2>;
using TriangularMesh3D = Mesh<3, 3>;

using TriangularMesh2DPtrWrapper = PointerWrapper<TriangularMesh2D>;

#endif //PYOPATRA_MESH_H
