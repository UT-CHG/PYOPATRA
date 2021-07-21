//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYTHONLPT_BINDINGS_CPP
#define PYTHONLPT_BINDINGS_CPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "mesh/mesh_vertex.h"
#include "mesh/mesh_water_column.h"
#include "mesh/mesh.h"

namespace py = pybind11;

PYBIND11_MODULE(pyopatra_pybind, m) {
    // Two dimensional MeshVertex and Mesh
    py::class_<MeshVertexBase2D>(m, "CppMeshVertexBase2D")
            .def(py::init<int>())
            .def("get_location", &MeshVertexBase2D::get_location)
            .def("get_latitude", &MeshVertexBase2D::get_latitude)
            .def("get_longitude", &MeshVertexBase2D::get_longitude)
            .def("get_velocity", &MeshVertexBase2D::get_velocity)
            .def("get_diffusion_coefficient", &MeshVertexBase2D::get_diffusion_coefficient)
            .def("set_diffusion_coefficient", &MeshVertexBase2D::set_diffusion_coefficient)
            .def("set_velocity", &MeshVertexBase2D::set_velocity)
            .def("set_location", &MeshVertexBase2D::set_location);

    py::class_<MeshVertex2D, MeshVertexBase2D>(m, "CppMeshVertex2D")
            .def(py::init<int>());

    py::class_<TriangularMeshElement2D>(m, "TriangularMeshElement2D")
            .def(py::init<MeshVertex<2>*, MeshVertex<2>*, MeshVertex<2>*, int>())
            .def("calculate_barycentric_coordinate", &TriangularMeshElement2D::calculate_barycentric_coordinate)
            .def("sample_velocity", &TriangularMeshElement2D::sample_velocity)
            .def("sample_diffusion_coefficient", &TriangularMeshElement2D::sample_diffusion_coefficient)
            .def("set_vertices", &TriangularMeshElement2D::set_vertices)
            .def("get_vertices", &TriangularMeshElement2D::get_vertices);

    py::class_<TriangularMesh2D>(m, "CppTriangularMesh2D")
            .def(py::init<int, int, std::vector<time_t>&&>())
            .def("set_vertex_location", &TriangularMesh2D::set_vertex_location)
            .def("set_vertex_velocity", &TriangularMesh2D::set_vertex_velocity)
            .def("set_vertex_diffusion", &TriangularMesh2D::set_vertex_diffusion)
            .def("get_vertex_pointer", &TriangularMesh2D::get_vertex_pointer, py::return_value_policy::reference)
            .def("get_vertex_locations", &TriangularMesh2D::get_vertex_locations)
            .def("get_velocities", &TriangularMesh2D::get_velocities)
            .def("set_water_column_adjacency", &TriangularMesh2D::set_water_column_adjacency)
            .def("get_water_column_adjacencies", &TriangularMesh2D::get_water_column_adjacencies, py::return_value_policy::reference)
            .def("check_water_column_adjacency", &TriangularMesh2D::check_water_column_adjacency)
            .def("check_mesh_element_vertex", &TriangularMesh2D::check_mesh_element_vertex)
            .def("get_water_column_pointer", &TriangularMesh2D::get_water_column_pointer, py::return_value_policy::reference)
            .def("set_element_vertex", &TriangularMesh2D::set_element_vertex)
            .def("get_water_columns_size", &TriangularMesh2D::get_water_columns_size)
            .def("add_particle", &TriangularMesh2D::add_particle)
            .def("get_all_particle_locations", &TriangularMesh2D::get_all_particle_locations);
//            .def()

    py::class_<TriangularWaterColumn2D>(m, "CppTriangleWaterColumn2D")
            .def("set_adjacent_columns", &TriangularWaterColumn2D::set_adjacent_columns);

    py::class_<ParticleList2D>(m, "CppParticleList2D")
            .def("create_particle", &ParticleList2D::create_particle);

}

#endif //PYTHONLPT_BINDINGS_CPP