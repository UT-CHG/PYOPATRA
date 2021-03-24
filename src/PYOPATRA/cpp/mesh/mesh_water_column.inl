//
// Created by Georgia Stuart on 3/14/21.
//

#include <cmath>
#include <exception>
#include "mesh_water_column.h"

template <int dimension>
std::tuple<MeshElementCursor<dimension>*, MeshElementCursor<dimension>*> WaterColumn<dimension>::get_element_depth_bounds(const Vector& location) {
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


template <int dimension>
Vector WaterColumn<dimension>::interpolate_velocity(const Particle<dimension> &particle) {
    if constexpr (dimension == 2) {
        return mesh_elements[0]->calculate_velocity(particle.get_location());
    } else {
        throw std::logic_error("Velocity interpolation for this dimension is not implemented.\n")
    }
}