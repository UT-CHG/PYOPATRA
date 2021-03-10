//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_ELEMENT_H
#define PYOPATRA_MESH_ELEMENT_H

#include <vector>
#include <array>
#include "../coordinate.h"
#include "mesh_vertex.h"

template <size_t T>
class PolygonMeshElement {
private:
    int mesh_index;
    std::array<MeshVertex*, T> vertices;

public:
    PolygonMeshElement() : mesh_index(0), vertices(0) {}
    explicit PolygonMeshElement(int mesh_index) : mesh_index(mesh_index), vertices(0) {}
    explicit PolygonMeshElement(int mesh_index, const std::array<MeshVertex*, T>& vertices) : mesh_index(mesh_index), vertices(vertices) {}

};

typedef PolygonMeshElement<3> TriangularMeshElement;
typedef PolygonMeshElement<4> QuadrilateralMeshElement;

#include "mesh_element.inl"

#endif //PYOPATRA_MESH_ELEMENT_H
