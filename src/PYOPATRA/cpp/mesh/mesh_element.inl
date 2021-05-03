//
// Created by Georgia Stuart on 2/4/21.
//

#include <stdexcept>


template <int num_vertices, int dimensions>
typename MeshElementT<num_vertices, dimensions>::VectorTd MeshElementT<num_vertices, dimensions>::calculate_barycentric_coordinate(const Vector& point) const {
    // From Real Time Collision Detection by Christer Ericson
    if constexpr ((num_vertices == 3) && (dimensions == 3)) {
        Vector3d pq(0.0, 0.0, 1.0);
        Vector3d pa = vertices[0]->get_location() - point;
        Vector3d pb = vertices[1]->get_location() - point;
        Vector3d pc = vertices[2]->get_location() - point;

        double u = pq.cross(pc).dot(pb);
        double v = pq.cross(pa).dot(pc);
        double w = pq.cross(pb).dot(pa);

        double denom = 1.0 / (u + v + w);

        return {u * denom, v * denom, w * denom};
    } else if constexpr ((num_vertices == 3) && (dimensions == 2)) {
        Vector2d v0 = vertices[1]->get_location() - vertices[0]->get_location();
        Vector2d v1 = vertices[2]->get_location() - vertices[0]->get_location();
        Vector2d v2 = point - vertices[0]->get_location();

        double d00 = v0.dot(v0);
        double d01 = v0.dot(v1);
        double d11 = v1.dot(v1);
        double d20 = v2.dot(v0);
        double d21 = v2.dot(v1);
        double denom = 1.0 / (d00 * d11 - d01 * d01);

        double v = (d11 * d20 - d01 * d21) * denom;
        double w = (d00 * d21 - d01 * d20) * denom;

        return {v, w, 1.0 - v - w};
    } else {
        throw std::logic_error(std::string("Barycentric coordinates are not implemented for polygons with ") + std::to_string(num_vertices) + std::string(" vertices."));
    }
}


template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::calculate_depth_at_point(const Vector& point) const {
    auto bc = calculate_barycentric_coordinate(point);
    return (vertices[0]->get_location() * bc[0] + vertices[1]->get_location() * bc[1] + vertices[2]->get_location() * bc[2])[2];
}


template <typename T>
void ignore(T &&)
{ }

template <int num_vertices, int dimensions>
int MeshElementT<num_vertices, dimensions>::check_halfspace(const typename MeshElementT<num_vertices, dimensions>::Vector& point) const {
    ignore(point);

    if constexpr ((num_vertices == 3) or (num_vertices == 2)) {
        return 0;
    } else {
        throw std::logic_error(std::string("Check halfspace is not implemented for polygons with ") + std::to_string(num_vertices) + std::string(" num_vertices_per_element."));
    }
}

template <int num_vertices>
int MeshElement3DT<num_vertices>::check_halfspace(const Vector3d &point) const {
    double dot = normals[0].dot(point - this->vertices[0]->get_location());

    if (dot < 0) {
        return -1;
    } else if (dot > 0) {
        return 1;
    } else {
        return 0;
    }
}


template <int num_vertices, int dimensions>
typename MeshElementT<num_vertices, dimensions>::Vector MeshElementT<num_vertices, dimensions>::sample_velocity(const VectorTd& barycentric_coordinates, int time_index) const {
    MeshElementT<num_vertices, dimensions>::Vector velocity = MeshElementT<num_vertices, dimensions>::Vector::Zero();

    for (int i = 0; i < num_vertices; i++) {
        velocity += vertices[i]->get_velocity()[time_index] * barycentric_coordinates(i);
    }

    return velocity;
}


template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::sample_viscosity(const VectorTd& barycentric_coordinates, int time_index) const {
    if constexpr (dimensions == 3) {
        double viscosity = 0.0;

        for (int i = 0; i < num_vertices; i++) {
            viscosity += vertices[i]->get_viscosity()[time_index] * barycentric_coordinates(i);
        }

        return viscosity;
    } else {
        return 0.0;
    }
}


template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::sample_density(const VectorTd& barycentric_coordinates, int time_index) const {
    if constexpr (dimensions == 3) {
        double density = 0.0;

        for (int i = 0; i < num_vertices; i++) {
            density += vertices[i]->get_density()[time_index] * barycentric_coordinates(i);
        }

        return density;
    } else {
        return 0.0;
    }
}


template <int num_vertices, int dimensions>
double MeshElementT<num_vertices, dimensions>::sample_water_viscosity(const VectorTd& barycentric_coordinates, int time_index) const {
    if constexpr (dimensions == 3) {
        double viscosity = 0.0;

        for (int i = 0; i < num_vertices; i++) {
            viscosity += vertices[i]->get_water_viscosity()[time_index] * barycentric_coordinates(i);
        }

        return viscosity;
    } else {
        return 0.0;
    }
}

template <int num_vertices, int dimensions>
typename MeshElementT<num_vertices, dimensions>::Vector MeshElementT<num_vertices, dimensions>::sample_diffusion_coefficient(const VectorTd& barycentric_coordinates, int time_index) const {
    Vector diffusion = Vector::Zero();

    for (int i = 0; i < num_vertices; i++) {
        diffusion += vertices[i]->get_diffusion_coefficient()[time_index] * barycentric_coordinates(i);
    }

    return diffusion;
}


