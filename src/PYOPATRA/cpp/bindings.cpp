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

namespace py = pybind11;

PYBIND11_MODULE(pyopatra_pybind, m) {
    py::class_<MeshVertexBase2D>(m, "Mesh")
            .def(py::init<double, double, Eigen::Matrix<double, 2, 1>>())
            .def("get_location", &MeshVertexBase2D::get_location)
            .def("get_latitude", &MeshVertexBase2D::get_latitude)
            .def("get_longitude", &MeshVertexBase2D::get_longitude)
            .def("get_velocity", &MeshVertexBase2D::get_velocity)
            .def("get_diffusion_coefficient", &MeshVertexBase2D::get_diffusion_coefficient)
            .def("set_diffusion_coefficient", &MeshVertexBase2D::set_diffusion_coefficient)
            .def("set_velocity", &MeshVertexBase2D::set_velocity)
            .def("set_location", &MeshVertexBase2D::set_location);

    py::class_<MeshVertex2D, MeshVertexBase2D>(m, "MeshVertex2D")
            .def(py::init<double, double, Eigen::Matrix<double, 2, 1>>());

}

#endif //PYTHONLPT_BINDINGS_CPP