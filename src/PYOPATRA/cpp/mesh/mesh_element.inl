//
// Created by Georgia Stuart on 2/4/21.
//

#include <stdexcept>

//template <int num_vertices>
//Vector3d MeshElementT<num_vertices, 3>::calculate_velocity(Vector3d &point) {
//    auto barycentric_coords = calculate_barycentric_coordinate(point);
//    double density = sample_density_at_point(barycentric_coords);
//    double viscosity = sample_viscosity_at_point(barycentric_coords);
//    double water_viscosity = sample_water_viscosity_at_point(barycentric_coords);
//}
//

template <int num_vertices, int dimensions>
typename MeshElementT<num_vertices, dimensions>::VectorTd MeshElementT<num_vertices, dimensions>::calculate_barycentric_coordinate(const Vector3d & point) const {
    // From Real Time Collision Detection by Christer Ericson
    if constexpr (num_vertices == 3) {
        Vector3d pq(0.0, 0.0, 1.0);
        Vector3d pa = vertices[0]->get_location() - point;
        Vector3d pb = vertices[1]->get_location() - point;
        Vector3d pc = vertices[2]->get_location() - point;

        double u = pq.cross(pc).dot(pb);
        double v = pq.cross(pa).dot(pc);
        double w = pq.cross(pb).dot(pa);

        double denom = 1.0 / (u + v + w);

        return {u * denom, v * denom, w * denom};
    } else {
        throw std::logic_error(std::string("Barycentric coordinates are not implemented for polygons with ") + std::to_string(num_vertices) + std::string(" vertices."));
    }
}


template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::calculate_depth_at_point(const Vector3d& point) const {
    auto bc = calculate_barycentric_coordinate(point);
    return (vertices[0]->get_location() * bc[0] + vertices[1]->get_location() * bc[1] + vertices[2]->get_location() * bc[2])[2];
}


template <int num_vertices, int dimensions>
int MeshElementT<num_vertices, dimensions>::check_halfspace(const Vector3d &point) {
    if constexpr (num_vertices == 3) {
        double dot = normals[0].dot(point - vertices[0]->get_location());

        if (dot < 0) {
            return -1;
        } else if (dot > 0) {
            return 1;
        } else {
            return 0;
        }

    } else {
        throw std::logic_error(std::string("Check halfspace is not implemented for polygons with ") + std::to_string(num_vertices) + std::string(" vertices."));
    }
}

template <int num_vertices, int dimensions>
Vector3d MeshElementT<num_vertices, dimensions>::sample_velocity_at_barycentric_coordinate(const VectorTd& barycentric_coordinates) {
    Vector3d velocity = Eigen::Vector3d::Zero();

    for (int i = 0; i < num_vertices; i++) {
        velocity += vertices[i]->get_velocity() * barycentric_coordinates[i];
    }

    return velocity;
}

template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::sample_viscosity_at_barycentric_coordinate(const VectorTd& barycentric_coordinates) {
    double viscosity = 0.0;

    for (int i = 0; i < num_vertices; i++) {
        viscosity += vertices[i]->get_viscosity() * barycentric_coordinates[i];
    }

    return viscosity;
}

template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::sample_density_at_barycentric_coordinate(const VectorTd& barycentric_coordinates) {
    double density = 0.0;

    for (int i = 0; i < num_vertices; i++) {
        density += vertices[i]->get_density() * barycentric_coordinates[i];
    }

    return density;
}

template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::sample_water_viscosity_at_barycentric_coordinate(const VectorTd& barycentric_coordinates) {
    double viscosity = 0.0;

    for (int i = 0; i < num_vertices; i++) {
        viscosity += vertices[i]->get_water_viscosity() * barycentric_coordinates[i];
    }

    return viscosity;
}

//template <int num_vertices, int dimensions>
//Vector3d MeshElementT<num_vertices, dimensions>::calculate_velocity(const Vector3d& point) const {
//    auto bc = calculate_barycentric_coordinate(point);
//}