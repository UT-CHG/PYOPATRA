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
typename MeshElementT<num_vertices, dimensions>::VectorTd MeshElementT<num_vertices, dimensions>::calculate_barycentric_coordinate(const Vector3d & point) {
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