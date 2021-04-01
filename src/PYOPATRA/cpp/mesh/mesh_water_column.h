//
// Created by Georgia Stuart on 3/14/21.
//

#ifndef PYOPATRA_MESH_WATER_COLUMN_H
#define PYOPATRA_MESH_WATER_COLUMN_H

#include <vector>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "mesh_element.h"
#include "../coordinate.h"
#include "../particle.h"

namespace py = pybind11;

template <int vertices, int dimension>
class WaterColumn {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;
private:
    std::vector<MeshElementT<vertices, dimension>> mesh_elements;
    int num_depths;
    std::array<WaterColumn<vertices, dimension>*, vertices> adjacent_columns;

    std::tuple<MeshElementT<vertices, dimension>*, MeshElementT<vertices, dimension>*> get_element_depth_bounds(const Vector &location);

public:
    WaterColumn() : mesh_elements(0, MeshElementT<vertices, dimension>()), num_depths(0) {}
    explicit WaterColumn(int num_depths) : mesh_elements(num_depths, MeshElementT<vertices, dimension>()), num_depths(num_depths) {}

    Vector interpolate_velocity(const Particle<dimension>& particle);
    void set_num_depths(int new_num_depths) { num_depths = new_num_depths; }
    void set_adjacent_columns(WaterColumn<vertices, dimension>* a, WaterColumn<vertices, dimension>* b, WaterColumn<vertices, dimension>* c) { adjacent_columns = {a, b, c}; }
    void set_adjacent_columns(WaterColumn<vertices, dimension>* a, WaterColumn<vertices, dimension>* b, WaterColumn<vertices, dimension>* c, WaterColumn<vertices, dimension>* d) { adjacent_columns = {a, b, c, d}; }
};

using TriangularWaterColumn2D = WaterColumn<3, 2>;
using TriangularWaterColumn3D = WaterColumn<3, 3>;

template <int vertices, int dimension>
class WaterColumnOverTime {
private:
    // Must be a python list composed of water columns
    py::list water_columns_over_time;

public:
    explicit WaterColumnOverTime(py::list& column_list) : water_columns_over_time(column_list) {}


};

#include "mesh_water_column.inl"

#endif //PYOPATRA_MESH_WATER_COLUMN_H
