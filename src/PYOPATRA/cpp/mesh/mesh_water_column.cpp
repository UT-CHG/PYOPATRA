//
// Created by Georgia Stuart on 3/14/21.
//

#include <cmath>
#include "mesh_water_column.h"

std::tuple<MeshElementCursor*, MeshElementCursor*> WaterColumn::get_element_depth_bounds(const Vector3d& location) {
    int above = 0;
    int below = num_depths - 1;

    while (below <= above - 1) {
        int m = (below + above) / 2;
        int direction = mesh_elements[m]->check_halfspace(location);

        if (direction == -1) {
            above = m + 1;
        } else if (direction == 1) {
            below = m + 1;
        } else {
            above = m;
            break;
        }
    }

    return {mesh_elements[above], mesh_elements[above + 1]};
}