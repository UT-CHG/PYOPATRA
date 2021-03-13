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

template <int num_vertices>
class MeshElementT {
private:
    std::vector<MeshVertex*> vertices;

    int mesh_index;

public:
    using VectorTd = Eigen::Matrix<double, num_vertices, 1>;

    MeshElementT(std::initializer_list<MeshVertex*> vertex_list) : vertices{vertex_list} {};

    Vector3d calculate_velocity(Vector3d& point);
    VectorTd calculate_barycentric_coordinate(const VectorTd& point);
    const std::vector<MeshVertex*>& get_vertices() const { return vertices; }
    double sample_density_at_point(const VectorTd& barycentric_coordinates);
    double sample_viscosity_at_point(const VectorTd& barycentric_coordinates);
    double sample_water_viscosity_at_point(const VectorTd& barycentric_coordinates);
};

class MeshElementCursor {
public:
    virtual ~MeshElementCursor() = default;
    virtual Vector3d calculate_velocity(Vector3d& point) = 0;
};

template <int num_vertices>
class MeshElementCursorT : public MeshElementCursor {
private:
    MeshElementT<num_vertices> *p_impl;

public:
    using MeshElementImpl = MeshElementT<num_vertices>;
    MeshElementCursorT(MeshElementImpl* mesh_element) : p_impl(mesh_element) {}

    virtual Vector3d calculate_velocity(Vector3d& point) {
        return p_impl->calculate_velocity(point);
    }

    void move(MeshElementT<num_vertices>* new_p_impl) { p_impl = new_p_impl; }
};

typedef  MeshElementT<3> TriangularMeshElement;

#include "mesh_element.inl"

#endif //PYOPATRA_MESH_ELEMENT_H
