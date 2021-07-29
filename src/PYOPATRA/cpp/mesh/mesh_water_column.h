//
// Created by Georgia Stuart on 3/14/21.
//

#ifndef PYOPATRA_MESH_WATER_COLUMN_H
#define PYOPATRA_MESH_WATER_COLUMN_H

#include <vector>
#include <tuple>
#include <iostream>
#include <cmath>

#include "mesh_element.h"
#include "../coordinate.h"
#include "../util.h"


template <int vertices, int dimension>
class WaterColumn {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;
private:
    MeshElementT<vertices, dimension> *mesh_elements;
    int num_depths;
    std::array<WaterColumn<vertices, dimension>*, vertices> adjacent_columns;
//    std::tuple<MeshElementT<vertices, dimension>*, MeshElementT<vertices, dimension>*> get_element_depth_bounds(const Vector &location);
    int index;

public:
    WaterColumn()
        : mesh_elements(nullptr)
        , num_depths(1)
        , index(0)
    {}
    explicit WaterColumn(int num_depths)
        : mesh_elements(nullptr)
        , num_depths(num_depths)
        , index(0)
    {}

    Vector interpolate_velocity(const Vector location, size_t time_index, double delta_t,  double time, double lower_time, double upper_time) {
        if constexpr (dimension == 2) {
            double t = ((time - lower_time) / (upper_time - lower_time));
            auto barycentric = mesh_elements[0].calculate_barycentric_coordinate(location);

            Vector lb = mesh_elements[0].sample_velocity(barycentric, time_index);
            Vector ub = mesh_elements[0].sample_velocity(barycentric, time_index + 1);
            Vector interpolated_velocity = (1 - t) * lb + t * ub;

            lb = mesh_elements[0].sample_diffusion_coefficient(barycentric, time_index);
            ub = mesh_elements[0].sample_diffusion_coefficient(barycentric, time_index + 1);
            Vector interpolated_diffusion = (1 - t) * lb + t * ub;

            Vector new_velocity = interpolated_velocity;
            double ra = unif_pi(generator);
            double rn = normal(generator);

            new_velocity(0) += rn * sqrt(4.0 * interpolated_diffusion(0) / (delta_t * 3600.0)) * sin(ra);
            new_velocity(1) += rn * sqrt(4.0 * interpolated_diffusion(1) / (delta_t * 3600.0)) * cos(ra);

            return new_velocity;
        }
    }
    void set_num_depths(int new_num_depths) { num_depths = new_num_depths; }
    void set_element_head(MeshElementT<vertices, dimension>* head) { mesh_elements = head; }
    void set_adjacent_columns(WaterColumn<vertices, dimension>* a, WaterColumn<vertices, dimension>* b, WaterColumn<vertices, dimension>* c) { adjacent_columns = {a, b, c}; }
    void set_adjacent_column(WaterColumn<vertices, dimension>* col, int position) { adjacent_columns[position] = col; }
    const std::array<WaterColumn<vertices, dimension>*, vertices>& get_adjacencies() const { return adjacent_columns; }
    MeshElementT<vertices, dimension>* get_mesh_elements() const { return mesh_elements; }
    int get_index() const { return index; }
    void set_index(int new_index) { index = new_index; }
};

using TriangularWaterColumn2D = WaterColumn<3, 2>;
using TriangularWaterColumn3D = WaterColumn<3, 3>;

#include "mesh_water_column.inl"

#endif //PYOPATRA_MESH_WATER_COLUMN_H
