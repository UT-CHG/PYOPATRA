//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_ELEMENT_H
#define PYOPATRA_MESH_ELEMENT_H

#include <vector>
#include <array>
#include "../coordinate.h"
#include "mesh_vertex.h"

class MeshElementBase {
private:
    int mesh_index;

public:
    virtual ~MeshElementBase() = default;

    template<typename vector_type>
    double sample_density_at_point(const vector_type barycentric_coordinates);

    template<typename vector_type>
    double sample_viscosity_at_point(const vector_type barycentric_coordinates);

    template<typename vector_type>
    double sample_water_viscosity_at_point(const vector_type barycentric_coordinates);
};


class TriangularMeshElement: public MeshElementBase {
private:
    std::array<MeshVertex*, 3> vertices;

public:
    TriangularMeshElement(MeshVertex *vertex_0, MeshVertex *vertex_1, MeshVertex *vertex_2);

    Vector3d calculate_barycentric_coordinate(const Coordinate3D& point);
    std::array<MeshVertex*, 3>& get_vertices() { return vertices; }
};

#include "mesh_element.inl"

#endif //PYOPATRA_MESH_ELEMENT_H
