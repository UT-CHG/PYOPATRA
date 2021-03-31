//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_ELEMENT_H
#define PYOPATRA_MESH_ELEMENT_H

#include <vector>
#include <array>
#include <tuple>
#include <initializer_list>
#include "../coordinate.h"
#include "mesh_vertex.h"

template <int num_vertices, int dimensions>
class MeshElementT {
public:
    using VertexArray = std::array<MeshVertex<dimensions>*, num_vertices>;
    using VectorTd = Eigen::Matrix<double, num_vertices, 1>;
    using Vector = Eigen::Matrix<double, dimensions, 1>;

protected:
    VertexArray vertices;
    int mesh_index;

public:
    MeshElementT()
        : vertices()
        , mesh_index(0)
    {}

    MeshElementT(MeshVertex<dimensions>* a, MeshVertex<dimensions>* b, MeshVertex<dimensions>* c, int mesh_index)
        : vertices{a, b, c}
        , mesh_index(mesh_index)
    {}

    MeshElementT(MeshVertex<dimensions>* a, MeshVertex<dimensions>* b, MeshVertex<dimensions>* c, MeshVertex<dimensions>* d, int mesh_index)
        : vertices{a, b, c, d}
        , mesh_index(mesh_index)
    {}

    void set_vertices(MeshVertex<dimensions>* a, MeshVertex<dimensions>* b, MeshVertex<dimensions>* c) { vertices = {a, b, c}; }
//    void set_vertices(MeshVertex<dimensions>* a, MeshVertex<dimensions>* b, MeshVertex<dimensions>* c, MeshVertex<dimensions>* d) { vertices = {a, b, c, d}; }
    void set_mesh_index(int new_mesh_index) { mesh_index = new_mesh_index; }

    VectorTd calculate_barycentric_coordinate(const Vector& point) const;
    double calculate_depth_at_point(const Vector& point) const;
    const VertexArray& get_vertices() const { return vertices; }
    Vector sample_velocity(const VectorTd& barycentric_coordinates) const;
    double sample_density(const VectorTd& barycentric_coordinates) const;
    double sample_viscosity(const VectorTd& barycentric_coordinates) const;
    double sample_water_viscosity(const VectorTd& barycentric_coordinates) const;
    Vector sample_diffusion_coefficient(const VectorTd& barycentric_coordinates) const;
    // 1 if deeper, -1 if shallower, 0 if on the element
    virtual int check_halfspace(const Vector& point) const;
};

template <int num_vertices>
class MeshElement3DT : public MeshElementT<num_vertices, 3> {
public:
    using Vector = Eigen::Matrix<double, 3, 1>;
    using PlaneNormalArray = std::array<Vector3d, num_vertices - 2>;
private:
    PlaneNormalArray normals;
public:
    MeshElement3DT()
        : MeshElementT<num_vertices, 3>()
        , normals()
    {}

    MeshElement3DT(MeshVertex<3>* a, MeshVertex<3>* b, MeshVertex<3>* c, int mesh_index)
        : MeshElementT<num_vertices, 3>(a, b, c, mesh_index)
        , normals{(a->get_location() - b->get_location()).cross(a->get_location() - c->get_location())}
    {}

    MeshElement3DT(MeshVertex<3>* a, MeshVertex<3>* b, MeshVertex<3>* c, MeshVertex<3>* d, int mesh_index)
        : MeshElementT<num_vertices, 3>(a, b, c, d, mesh_index)
        , normals{(a->get_location() - b->get_location()).cross(a->get_location() - c->get_location()), (a->get_location() - c->get_location()).cross(a->get_location() - d->get_location())}
    {}

    int check_halfspace(const Vector& point) const override;
};

template <int dimensions>
class InterpolatedValues {
public:
    using Vector = Eigen::Matrix<double, dimensions, 1>;

    Vector velocity, diffusion_coefficient;
    double density, viscosity, water_viscosity;
};

template <int dimensions>
class MeshElementCursor {
public:
    using Vector = Eigen::Matrix<double, dimensions, 1>;
    virtual ~MeshElementCursor() = default;
    virtual Vector calculate_velocity(const Vector& point) const = 0;
    virtual double get_depth_at_point(const Vector& point) const = 0;
    virtual int check_halfspace(const Vector& point) const = 0;
    virtual void get_interpolated_values(const Vector& point, InterpolatedValues<dimensions>& interpolated_values) const = 0;
};

template <int num_vertices, int num_dimensions>
class MeshElementCursorT : public MeshElementCursor<num_dimensions> {
private:
    MeshElementT<num_vertices, num_dimensions> *p_impl;

public:
    using MeshElementImpl = MeshElementT<num_vertices, num_dimensions>;
    using Vector = Eigen::Matrix<double, num_dimensions, 1>;
    explicit MeshElementCursorT(MeshElementImpl* mesh_element) : p_impl(mesh_element) {}

    Vector calculate_velocity(const Vector3d& point) const override {
        return p_impl->sample_velocity(p_impl->calculate_barycentric_coordinate(point));
    }

    double get_depth_at_point(const Vector3d& point) const override {
        return p_impl->calculate_depth_at_point(point);
    }

    int check_halfspace(const Vector& point) const override {
        if constexpr (num_dimensions == 3) {
            return p_impl->check_halfspace(point);
        } else {
            return 0.0;
        }
    }

    void get_interpolated_values(const Vector3d& point, InterpolatedValues<num_dimensions>& interpolated_values) const override {
        auto bc = p_impl->calculate_barycentric_coordinate(point);
        interpolated_values.velocity = p_impl->sample_velocity(bc);
        interpolated_values.density = p_impl->sample_density(bc);
        interpolated_values.viscosity = p_impl->sample_viscosity(bc);
        interpolated_values.water_viscosity = p_impl->sample_water_viscosity(bc);
    }


    void move(MeshElementT<num_vertices, num_dimensions>* new_p_impl) { p_impl = new_p_impl; }
};

typedef MeshElementT<3, 2> TriangularMeshElement2D;
typedef MeshElement3DT<3> TriangularMeshElement3D;

typedef MeshElementCursorT<3, 2> TriangularMeshCursor2D;
typedef MeshElementCursorT<3, 3> TriangularMeshCursor3D;

#include "mesh_element.inl"

#endif //PYOPATRA_MESH_ELEMENT_H
