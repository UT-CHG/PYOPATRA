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
#include "mesh/mesh.h"

namespace py = pybind11;

PYBIND11_MODULE(pyopatra_pybind, m) {
    // Two dimensional MeshVertex and Mesh
    py::class_<MeshVertexBase2D>(m, "CppMeshVertexBase2D")
            .def(py::init<double, double, Eigen::Matrix<double, 2, 1>>())
            .def("get_location", &MeshVertexBase2D::get_location)
            .def("get_latitude", &MeshVertexBase2D::get_latitude)
            .def("get_longitude", &MeshVertexBase2D::get_longitude)
            .def("get_velocity", &MeshVertexBase2D::get_velocity)
            .def("get_diffusion_coefficient", &MeshVertexBase2D::get_diffusion_coefficient)
            .def("set_diffusion_coefficient", &MeshVertexBase2D::set_diffusion_coefficient)
            .def("set_velocity", &MeshVertexBase2D::set_velocity)
            .def("set_location", &MeshVertexBase2D::set_location);

    py::class_<MeshVertex2D, MeshVertexBase2D>(m, "CppMeshVertex2D")
            .def(py::init<double, double, Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 1>>());

    py::class_<TriangularMeshElement2D>(m, "TriangularMeshElement2D")
            .def(py::init<MeshVertex<2>*, MeshVertex<2>*, MeshVertex<2>*, int>())
            .def("calculate_barycentric_coordinate", &TriangularMeshElement2D::calculate_barycentric_coordinate)
            .def("sample_velocity", &TriangularMeshElement2D::sample_velocity)
            .def("sample_diffusion_coefficient", &TriangularMeshElement2D::sample_diffusion_coefficient)
            .def("set_vertices", &TriangularMeshElement2D::set_vertices)
            .def("get_vertices", &TriangularMeshElement2D::get_vertices);

    py::class_<TriangularMesh2D>(m, "CppTriangularMesh2D")
            .def(py::init<int, int, int, std::vector<time_t>&&>());
//            .def()

}

#endif //PYTHONLPT_BINDINGS_CPP