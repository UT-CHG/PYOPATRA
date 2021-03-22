//
// Created by Georgia Stuart on 3/14/21.
//

#ifndef PYOPATRA_MESH_WATER_COLUMN_H
#define PYOPATRA_MESH_WATER_COLUMN_H

#include <vector>
#include <tuple>
#include "mesh_element.h"
#include "../coordinate.h"
#include "../particle.h"

template <int dimension>
class WaterColumn {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;
private:
    std::vector<MeshElementCursor<dimension>*> mesh_elements;
    int num_depths;
    std::tuple<MeshElementCursor<dimension>*, MeshElementCursor<dimension>*> get_element_depth_bounds(const Vector3d &location);

public:
    Vector interpolate_velocity(const Particle<dimension>& particle);
};

using WaterColumn2D = WaterColumn<2>;
using WaterColumn3D = WaterColumn<3>;

#endif //PYOPATRA_MESH_WATER_COLUMN_H
