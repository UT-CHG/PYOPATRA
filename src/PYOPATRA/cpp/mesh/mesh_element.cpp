//
// Created by Georgia Stuart on 2/4/21.
//

#include "mesh_element.h"

TriangularMeshElement::TriangularMeshElement(MeshVertex *vertex_0, MeshVertex *vertex_1, MeshVertex *vertex_2)
    : vertices({vertex_0, vertex_1, vertex_2})
    , v0(vertex_1->get_location() - vertex_0->get_location())
    , v1(vertex_2->get_location() - vertex_0->get_location())
    , d00(v0.dot(v0))
    , d01(v0.dot(v1))
    , d11(v1.dot(v1))
    , denom(d00 * d11 - d01 * d01)
{}

// From https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
Vector3d TriangularMeshElement::calculate_barycentric_coordinate(const Coordinate3D &point) {
    Vector3d v2 = point - vertices[0]->get_location();

    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;

    return {v, w, 1.0 - v - w};
}