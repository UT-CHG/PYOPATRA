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
    virtual double sample_density_at_point(const Coordinate3D& location);
    virtual double sample_viscosity_at_point(const Coordinate3D& location);
    virtual double sample_water_viscosity_at_point(const Coordinate3D& location);
};


class TriangularMeshElement: public MeshElementBase {
private:
    std::array<MeshVertex*, 3> vertices;
    Vector3d v0, v1;
    double d00, d01, d11, denom;

public:
    TriangularMeshElement(MeshVertex *vertex_0, MeshVertex *vertex_1, MeshVertex *vertex_2);

    Vector3d calculate_barycentric_coordinate(const Coordinate3D& point);
};

#endif //PYOPATRA_MESH_ELEMENT_H
