//
// Created by Georgia Stuart on 3/14/21.
//

#ifndef PYOPATRA_MESH_WATER_COLUMN_H
#define PYOPATRA_MESH_WATER_COLUMN_H

#include <vector>
#include <tuple>
#include <iostream>

#include "mesh_element.h"
#include "../coordinate.h"
//#include "../particle.h"

template <int vertices, int dimension>
class WaterColumn {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;
private:
    std::vector<MeshElementT<vertices, dimension>> mesh_elements;
    int num_depths;
    std::array<WaterColumn<vertices, dimension>*, vertices> adjacent_columns;
    std::tuple<MeshElementT<vertices, dimension>*, MeshElementT<vertices, dimension>*> get_element_depth_bounds(const Vector &location);
    size_t index;

public:
    WaterColumn() : mesh_elements(1, MeshElementT<vertices, dimension>()), num_depths(1) {}
    explicit WaterColumn(int num_depths) : mesh_elements(num_depths, MeshElementT<vertices, dimension>()), num_depths(num_depths) {}

    Vector interpolate_velocity(const Vector location, size_t time_index) {
        if constexpr (dimension == 2) {
            auto barycentric = mesh_elements[0].calculate_barycentric_coordinate(location);
            return mesh_elements[0].sample_velocity(barycentric, time_index);
        }
    }
    void set_num_depths(int new_num_depths) { num_depths = new_num_depths; }
    void set_adjacent_columns(WaterColumn<vertices, dimension>* a, WaterColumn<vertices, dimension>* b, WaterColumn<vertices, dimension>* c) { adjacent_columns = {a, b, c}; }
    void set_adjacent_column(WaterColumn<vertices, dimension>* col, int position) { adjacent_columns[position] = col; }
    void set_element_vertex(int depth_index, int position, MeshVertex<dimension>* vert) { mesh_elements[depth_index].set_vertex(vert, position); }
    const std::array<WaterColumn<vertices, dimension>*, vertices>& get_adjacencies() const { return adjacent_columns; }
    const std::vector<MeshElementT<vertices, dimension>>& get_mesh_elements() const { return mesh_elements; }
    size_t get_index() const { return index; }
    void set_index(size_t new_index) { index = new_index; }
};

using TriangularWaterColumn2D = WaterColumn<3, 2>;
using TriangularWaterColumn3D = WaterColumn<3, 3>;

#include "mesh_water_column.inl"

#endif //PYOPATRA_MESH_WATER_COLUMN_H
