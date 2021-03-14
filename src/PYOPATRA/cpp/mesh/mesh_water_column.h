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

class WaterColumn {
private:
    std::vector<MeshElementCursor*> mesh_elements;
    int num_depths;
    std::tuple<MeshElementCursor*, MeshElementCursor*> get_element_depth_bounds(const Vector3d &location);

public:
    Vector3d interpolate_velocity(const Particle& particle);
};

#endif //PYOPATRA_MESH_WATER_COLUMN_H
