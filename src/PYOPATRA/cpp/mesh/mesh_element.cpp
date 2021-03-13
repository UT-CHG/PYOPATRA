//
// Created by Georgia Stuart on 2/4/21.
//

#include "mesh_element.h"

TriangularMeshElement::TriangularMeshElement(MeshVertex *vertex_0, MeshVertex *vertex_1, MeshVertex *vertex_2)
    : vertices({vertex_0, vertex_1, vertex_2})
{}

// From Real Time Collision Detection by Christer Ericson
Vector3d TriangularMeshElement::calculate_barycentric_coordinate(const Coordinate3D &point) {
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