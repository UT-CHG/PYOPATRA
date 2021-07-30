//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_ELEMENT_H
#define PYOPATRA_MESH_ELEMENT_H

#include <vector>
#include <array>
#include <tuple>
#include <initializer_list>
#include <iostream>
#include "../coordinate.h"
#include "mesh_vertex.h"

template <int num_vertices, int dimensions>
class MeshElementT {
public:
    using VertexArray = std::array<int, num_vertices>;
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

    MeshElementT(int a, int b, int c, int mesh_index)
        : vertices{a, b, c}
        , mesh_index(mesh_index)
    {}

    MeshElementT(int a, int b, int c, int d, int mesh_index)
        : vertices{a, b, c, d}
        , mesh_index(mesh_index)
    {}

    void set_vertices(int a, int b, int c) { vertices = {a, b, c}; }
//    void set_vertices(MeshVertex<dimensions>* a, MeshVertex<dimensions>* b, MeshVertex<dimensions>* c, MeshVertex<dimensions>* d) { vertices = {a, b, c, d}; }
    void set_vertex(int vert, int position) { vertices[position] = vert; }
    void set_mesh_index(int new_mesh_index) { mesh_index = new_mesh_index; }

    VectorTd calculate_barycentric_coordinate(const MeshVertex<dimensions>* vertex_array, const Vector& point) const;
    double calculate_depth_at_point(const MeshVertex<dimensions>* vertex_array, const Vector& point) const;
    const VertexArray& get_vertices() const { return vertices; }
    Vector sample_velocity(const MeshVertex<dimensions>* vertex_array, const Eigen::Matrix<double, dimensions, 1>* velocities, const VectorTd& barycentric_coordinates, int time_index) const;
//    double sample_density(const MeshVertex<dimensions>* vertex_array, const VectorTd& barycentric_coordinates, int time_index) const;
//    double sample_viscosity(const MeshVertex<dimensions>* vertex_array, const VectorTd& barycentric_coordinates, int time_index) const;
//    double sample_water_viscosity(const MeshVertex<dimensions>* vertex_array, const VectorTd& barycentric_coordinates, int time_index) const;
    Vector sample_diffusion_coefficient(const MeshVertex<dimensions>* vertex_array, const Eigen::Matrix<double, dimensions, 1>* diffusions, const VectorTd& barycentric_coordinates, int time_index) const;
    // 1 if deeper, -1 if shallower, 0 if on the element
//    virtual int check_halfspace(const MeshVertex<3>* vertex_array, const Vector& point) const;
};
//
//template <int num_vertices>
//class MeshElement3DT : public MeshElementT<num_vertices, 3> {
//public:
//    using Vector = Eigen::Matrix<double, 3, 1>;
//    using PlaneNormalArray = std::array<Vector3d, num_vertices - 2>;
//private:
//    PlaneNormalArray normals;
//public:
//    MeshElement3DT()
//        : MeshElementT<num_vertices, 3>()
//        , normals()
//    {}
//
//    MeshElement3DT(int a, int b, int c, int mesh_index, MeshVertex<3>* vertex_array)
//        : MeshElementT<num_vertices, 3>(a, b, c, mesh_index)
//        , normals{(vertex_array[a].get_location() - vertex_array[b].get_location())
//        .cross(vertex_array[a].get_location() - vertex_array[c].get_location())}
//    {}
//
//    MeshElement3DT(int a, int b, int c, int d, int mesh_index, MeshVertex<3>* vertex_array)
//        : MeshElementT<num_vertices, 3>(a, b, c, d, mesh_index)
//        , normals{(vertex_array[a].get_location() - vertex_array[b].get_location())
//        .cross(vertex_array[a].get_location() - vertex_array[c].get_location()),
//        (vertex_array[a].get_location() - vertex_array[c].get_location())
//        .cross(vertex_array[a].get_location() - vertex_array[d].get_location())}
//    {}
//
//    int check_halfspace(const MeshVertex<3>* vertex_array, const Vector& point) const override;
//};
//
//template <int dimensions>
//class InterpolatedValues {
//public:
//    using Vector = Eigen::Matrix<double, dimensions, 1>;
//
//    Vector velocity, diffusion_coefficient;
//    double density, viscosity, water_viscosity;
//};

template <int dimensions>
class MeshElementCursor {
public:
    using Vector = Eigen::Matrix<double, dimensions, 1>;
    virtual ~MeshElementCursor() = default;
    virtual Vector calculate_velocity(const MeshVertex<3>* vertex_array, const Vector& point) const = 0;
//    virtual double get_depth_at_point(const MeshVertex<3>* vertex_array, const Vector& point) const = 0;
//    virtual int check_halfspace(const MeshVertex<3>* vertex_array, const Vector& point) const = 0;
//    virtual void get_interpolated_values(const MeshVertex<3>* vertex_array, const Vector& point, InterpolatedValues<dimensions>& interpolated_values, int time_index) const = 0;
};

template <int num_vertices, int num_dimensions>
class MeshElementCursorT : public MeshElementCursor<num_dimensions> {
private:
    MeshElementT<num_vertices, num_dimensions> *p_impl;

public:
    using MeshElementImpl = MeshElementT<num_vertices, num_dimensions>;
    using Vector = Eigen::Matrix<double, num_dimensions, 1>;
    explicit MeshElementCursorT(MeshElementImpl* mesh_element) : p_impl(mesh_element) {}

    Vector calculate_velocity(const MeshVertex<3>* vertex_array, const Eigen::Matrix<double, num_dimensions, 1>* velocities, const Vector3d& point) const override {
        return p_impl->sample_velocity(vertex_array, velocities, p_impl->calculate_barycentric_coordinate(vertex_array, point), 0);
    }

//    double get_depth_at_point(const MeshVertex<3>* vertex_array, const Vector3d& point) const override {
//        return p_impl->calculate_depth_at_point(vertex_array, point);
//    }
//
//    int check_halfspace(const MeshVertex<3>* vertex_array, const Vector& point) const override {
//        if constexpr (num_dimensions == 3) {
//            return p_impl->check_halfspace(vertex_array, point);
//        } else {
//            return 0.0;
//        }
//    }

//    void get_interpolated_values(const MeshVertex<3>* vertex_array,
//                                 const Vector3d& point,
//                                 InterpolatedValues<num_dimensions>& interpolated_values, int time_index) const override {
//        auto bc = p_impl->calculate_barycentric_coordinate(vertex_array, point);
//        interpolated_values.velocity = p_impl->sample_velocity(vertex_array, bc, time_index);
//        interpolated_values.density = p_impl->sample_density(vertex_array, bc, time_index);
//        interpolated_values.viscosity = p_impl->sample_viscosity(vertex_array, bc, time_index);
//        interpolated_values.water_viscosity = p_impl->sample_water_viscosity(vertex_array, bc, time_index);
//    }
//
//
//    void move(MeshElementT<num_vertices, num_dimensions>* new_p_impl) { p_impl = new_p_impl; }
};

typedef MeshElementT<3, 2> TriangularMeshElement2D;
//typedef MeshElement3DT<3> TriangularMeshElement3D;

typedef MeshElementCursorT<3, 2> TriangularMeshCursor2D;
//typedef MeshElementCursorT<3, 3> TriangularMeshCursor3D;

#include "mesh_element.inl"

#endif //PYOPATRA_MESH_ELEMENT_H
