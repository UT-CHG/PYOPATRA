//
// Created by Georgia Stuart on 2/4/21.
//

#include "mesh_element.h"

// From Real Time Collision Detection by Christer Ericson
template <>
MeshElementT<3>::VectorTd MeshElementT<3>::calculate_barycentric_coordinate(const MeshElementT<3>::VectorTd& point) {
    Vector3d pq(0.0, 0.0, 1.0);
    Vector3d pa = vertices[0]->get_location() - point;
    Vector3d pb = vertices[1]->get_location() - point;
    Vector3d pc = vertices[2]->get_location() - point;

    double u = pq.cross(pc).dot(pb);
    double v = pq.cross(pa).dot(pc);
    double w = pq.cross(pb).dot(pa);

    double denom = 1.0 / (u + v + w);

    return {u * denom, v * denom, w * denom};
}