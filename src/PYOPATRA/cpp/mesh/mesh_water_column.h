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
    size_t mesh_elements;
    int num_depths;
    std::array<long int, vertices> adjacent_columns;
//    std::tuple<MeshElementT<vertices, dimension>*, MeshElementT<vertices, dimension>*> get_element_depth_bounds(const Vector &location);
    size_t index;

public:
    WaterColumn()
        : mesh_elements(0)
        , num_depths(1)
        , index(0)
    {
        for (int i = 0; i < vertices; i++) {
            adjacent_columns[i] = -1;
        }
    }
    explicit WaterColumn(int num_depths)
        : mesh_elements(0)
        , num_depths(num_depths)
        , index(0)
    {
        for (int i = 0; i < vertices; i++) {
            adjacent_columns[i] = -1;
        }
    }

    void interpolate_velocity(const MeshElementT<vertices, dimension>* elements, const MeshVertex<dimension>* vertex_array,
                              Vector* velocities, Vector* diffusions, Vector* winds, const Vector& location, Eigen::Index time_index, Eigen::Index wind_time_index, double delta_t,
                              double time, double lower_time, double upper_time, double lower_wind_time, double upper_wind_time, 
                              double wind_coef, bool constant_diffusion, double constant_diffusion_coef, Vector &out_vec) {
        if constexpr (dimension == 2) {
            double t = ((time - lower_time) / (upper_time - lower_time));
            double wt = 0.0;
            if (upper_wind_time != lower_wind_time)
                wt = ((time - lower_wind_time) / (upper_wind_time - lower_wind_time));
            auto barycentric = elements[mesh_elements].calculate_barycentric_coordinate(vertex_array, location);

            Vector lb = elements[mesh_elements].sample_velocity(vertex_array, velocities, barycentric, time_index);
            Vector ub = elements[mesh_elements].sample_velocity(vertex_array, velocities, barycentric, time_index + 1);
            Vector interpolated_velocity = (1 - t) * lb + t * ub;

            Vector interpolated_diffusion;
            if (constant_diffusion) {
                interpolated_diffusion << constant_diffusion_coef, constant_diffusion_coef;
            }
            else {
                lb = elements[mesh_elements].sample_diffusion_coefficient(vertex_array, diffusions, barycentric, time_index);
                ub = elements[mesh_elements].sample_diffusion_coefficient(vertex_array, diffusions, barycentric, time_index + 1);
                interpolated_diffusion = (1 - t) * lb + t * ub;
            }

            lb = elements[mesh_elements].sample_wind(vertex_array, winds, barycentric, wind_time_index);
            ub = elements[mesh_elements].sample_wind(vertex_array, winds, barycentric, wind_time_index + 1);
            Vector interpolated_wind = (1 - wt) * lb + wt * ub;

            out_vec = interpolated_velocity;
            double ra = unif_pi(generator);
            double rn = normal(generator);

            out_vec(0) += rn * sqrt(4.0 * interpolated_diffusion(0) / (delta_t * 3600.0)) * sin(ra);
            out_vec(1) += rn * sqrt(4.0 * interpolated_diffusion(1) / (delta_t * 3600.0)) * cos(ra);

            double magnitude = sqrt(pow(interpolated_wind(0), 2.0) + pow(interpolated_wind(1), 2.0));
            double theta = 0.0;
            if (magnitude <= 25.0) {
                theta = 40 - 8 * magnitude;
            }

            double cos_theta = cos(theta * M_PI / 180.0);
            double sin_theta = sin(theta * M_PI / 180.0);
            out_vec(0) += wind_coef * (-1.0 * sin_theta * interpolated_wind(1) + cos_theta * interpolated_wind(0));
            out_vec(1) += wind_coef * (cos_theta * interpolated_wind(1) + sin_theta * interpolated_wind(0));
        }
    }
    void set_num_depths(int new_num_depths) { num_depths = new_num_depths; }
    void set_element_head(size_t head) { mesh_elements = head; }
    void set_adjacent_columns(size_t a, size_t b, size_t c) { adjacent_columns = {a, b, c}; }
    void set_adjacent_column(size_t col, size_t position) { adjacent_columns[position] = col; }
    const std::array<long int, vertices>& get_adjacencies() const { return adjacent_columns; }
    size_t get_mesh_elements() const { return mesh_elements; }
    size_t get_index() const { return index; }
    void set_index(size_t new_index) { index = new_index; }
};

using TriangularWaterColumn2D = WaterColumn<3, 2>;
using TriangularWaterColumn3D = WaterColumn<3, 3>;

#include "mesh_water_column.inl"

#endif //PYOPATRA_MESH_WATER_COLUMN_H
