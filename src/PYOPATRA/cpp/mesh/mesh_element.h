//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_ELEMENT_H
#define PYOPATRA_MESH_ELEMENT_H

#include <vector>
#include <array>
#include <initializer_list>
#include "../coordinate.h"
#include "mesh_vertex.h"

template <int num_vertices, int dimensions>
class MeshElementT {
public:
    using VertexArray = std::array<MeshVertex*, num_vertices>;
    using VectorTd = Eigen::Matrix<double, num_vertices, 1>;
private:
    VertexArray vertices;
    int mesh_index;

public:
    MeshElementT(MeshVertex* a, MeshVertex* b, MeshVertex* c) : vertices{a, b, c} {};
    MeshElementT(MeshVertex* a, MeshVertex* b, MeshVertex* c, MeshVertex* d) : vertices{a, b, c, d} {};

    Vector3d calculate_velocity(Vector3d& point);
    VectorTd calculate_barycentric_coordinate(const Vector3d& point);
    const VertexArray& get_vertices() const { return vertices; }
    double sample_density_at_point(const VectorTd& barycentric_coordinates);
    double sample_viscosity_at_point(const VectorTd& barycentric_coordinates);
    double sample_water_viscosity_at_point(const VectorTd& barycentric_coordinates);
};

class MeshElementCursor {
public:
    virtual ~MeshElementCursor() = default;
    virtual Vector3d calculate_velocity(Vector3d& point) = 0;
};

template <int num_vertices, int num_dimensions>
class MeshElementCursorT : public MeshElementCursor {
private:
    MeshElementT<num_vertices, num_dimensions> *p_impl;

public:
    using MeshElementImpl = MeshElementT<num_vertices, num_dimensions>;
    MeshElementCursorT(MeshElementImpl* mesh_element) : p_impl(mesh_element) {}

    virtual Vector3d calculate_velocity(Vector3d& point) {
        return p_impl->calculate_velocity(point);
    }

    void move(MeshElementT<num_vertices, num_dimensions>* new_p_impl) { p_impl = new_p_impl; }
};

typedef MeshElementT<3, 2> TriangularMeshElement2D;
typedef MeshElementT<3, 3> TriangularMeshElement3D;

#include "mesh_element.inl"

#endif //PYOPATRA_MESH_ELEMENT_H
