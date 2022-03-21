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
#include "inversion_tools/objective_functions.h"
#include "solver.h"
#include "particle_list.h"

namespace py = pybind11;

PYBIND11_MODULE(pyopatra_pybind, m) {
    // Two dimensional MeshVertex and Mesh
    py::class_<MeshVertexBase2D>(m, "CppMeshVertexBase2D")
            .def(py::init<int, int, int>())
            .def("get_location", &MeshVertexBase2D::get_location)
            .def("get_latitude", &MeshVertexBase2D::get_latitude)
            .def("get_longitude", &MeshVertexBase2D::get_longitude)
            .def("get_velocity", &MeshVertexBase2D::get_velocity)
            .def("get_diffusion_coefficient", &MeshVertexBase2D::get_diffusion_coefficient)
            .def("set_diffusion_coefficient", &MeshVertexBase2D::set_diffusion_coefficient)
            .def("set_velocity", &MeshVertexBase2D::set_velocity)
            .def("set_wind", &MeshVertex2D::set_wind)
            .def("set_location", &MeshVertexBase2D::set_location);

    py::class_<MeshVertex2D, MeshVertexBase2D>(m, "CppMeshVertex2D")
            .def(py::init<int, int, int>());

    py::class_<TriangularMeshElement2D>(m, "TriangularMeshElement2D")
            .def(py::init<>())
            .def("calculate_barycentric_coordinate", &TriangularMeshElement2D::calculate_barycentric_coordinate)
            .def("sample_velocity", &TriangularMeshElement2D::sample_velocity)
            .def("sample_diffusion_coefficient", &TriangularMeshElement2D::sample_diffusion_coefficient)
//            .def("set_vertices", &TriangularMeshElement2D::set_vertices)
            .def("get_vertices", &TriangularMeshElement2D::get_vertices);

//    py::class_<TriangularMesh2D, std::unique_ptr<TriangularMesh2D, py::nodelete>>(m, "CppTriangularMesh2D")
    py::class_<TriangularMesh2D>(m, "CppTriangularMesh2D")
            .def(py::init<int, int, std::vector<double>&&, std::vector<double>&&, double>())
            .def("set_vertex_location", &TriangularMesh2D::set_vertex_location)
            .def("set_vertex_locations", &TriangularMesh2D::set_vertex_locations)
            .def("set_vertex_velocity", &TriangularMesh2D::set_vertex_velocity)
            .def("set_velocities", &TriangularMesh2D::set_velocities)
            .def("set_vertex_diffusion", &TriangularMesh2D::set_vertex_diffusion)
            .def("set_diffusions", &TriangularMesh2D::set_diffusions)
            .def("set_vertex_wind", &TriangularMesh2D::set_vertex_wind)
            .def("set_winds", &TriangularMesh2D::set_winds)
            .def("get_vertex_pointer", &TriangularMesh2D::get_vertex_pointer, py::return_value_policy::reference)
            .def("get_vertex_locations", &TriangularMesh2D::get_vertex_locations)
            .def("get_velocities", &TriangularMesh2D::get_velocities)
            .def("set_water_column_adjacency", &TriangularMesh2D::set_water_column_adjacency)
            .def("set_water_column_adjacencies", &TriangularMesh2D::set_water_column_adjacencies)
            .def("get_water_column_adjacencies", &TriangularMesh2D::get_water_column_adjacencies, py::return_value_policy::reference)
            .def("check_water_column_adjacency", &TriangularMesh2D::check_water_column_adjacency)
            .def("check_mesh_element_vertex", &TriangularMesh2D::check_mesh_element_vertex)
            .def("get_water_column_pointer", &TriangularMesh2D::get_water_column_pointer, py::return_value_policy::reference)
            .def("set_element_vertex", &TriangularMesh2D::set_element_vertex)
            .def("set_element_vertices", &TriangularMesh2D::set_element_vertices)
            .def("get_water_columns_size", &TriangularMesh2D::get_water_columns_size)
            .def("get_pointer_wrapper", &TriangularMesh2D::get_pointer_wrapper, py::return_value_policy::copy)
            .def("set_wind_coef", &TriangularMesh2D::set_wind_coef);

//    py::class_<TriangularMesh2DSolver, std::unique_ptr<TriangularMesh2DSolver, py::nodelete>>(m, "CppTriangularMesh2DSolver")
    py::class_<TriangularMesh2DSolver>(m, "CppTriangularMesh2DSolver")
            .def(py::init<TriangularMesh2DPtrWrapper&, ParticleList2DPtrWrapper&, ObjectiveFunction2DPtrWrapper&, Eigen::VectorXd, Eigen::VectorXd>())
            .def(py::init<TriangularMesh2DPtrWrapper&, ParticleList2DPtrWrapper&, Eigen::VectorXd, Eigen::VectorXd>())
            .def("get_current_time", &TriangularMesh2DSolver::get_current_time)
            .def("time_step", &TriangularMesh2DSolver::time_step)
            .def("reset_solver", &TriangularMesh2DSolver::reset_solver)
            .def("update_particle_location_indices", &TriangularMesh2DSolver::update_particle_location_indices)
            .def("calculate_objective_value", &TriangularMesh2DSolver::calculate_objective_value);

//    py::class_<ParticleList2D, std::unique_ptr<ParticleList2D, py::nodelete>>(m, "CppParticleList2D")
    py::class_<ParticleList2D>(m, "CppParticleList2D")
            .def(py::init<>())
            .def("create_particle", &ParticleList2D::create_particle)
            .def("reset_particles", &ParticleList2D::delete_all_particles)
            .def("get_pointer_wrapper", &ParticleList2D::get_pointer_wrapper, py::return_value_policy::copy)
            .def("get_all_particle_locations", &ParticleList2D::get_all_particle_locations)
            .def("get_all_particle_column_indices", &ParticleList2D::get_all_particle_column_indices)
            .def("get_length", &ParticleList2D::get_length);

//    py::class_<SlicedWassersteinDistance2D, std::unique_ptr<SlicedWassersteinDistance2D, py::nodelete>>(m, "CppSlicedWassersteinDistance2D")
    py::class_<SlicedWassersteinDistance2D>(m, "CppSlicedWassersteinDistance2D")
            .def(py::init<int, int, const Eigen::Vector4d&, int, unsigned int>())
            .def("calculate_value", &SlicedWassersteinDistance2D::calculate_value)
            .def("set_observed_values", &SlicedWassersteinDistance2D::set_observed_values)
            .def("get_pointer_wrapper", &SlicedWassersteinDistance2D::get_pointer_wrapper, py::return_value_policy::copy);

    py::class_<BhattacharyyaDistance2D>(m, "CppBhattacharyyaDistance2D")
            .def(py::init<int, int, const Eigen::Vector4d&>())
            .def("calculate_value", &SlicedWassersteinDistance2D::calculate_value)
            .def("set_observed_values", &SlicedWassersteinDistance2D::set_observed_values)
            .def("get_pointer_wrapper", &SlicedWassersteinDistance2D::get_pointer_wrapper, py::return_value_policy::copy);

    py::class_<TriangularMesh2DPtrWrapper>(m, "CppTriangularMesh2DPtrWrapper");
    py::class_<ParticleList2DPtrWrapper>(m, "CppParticleList2DPtrWrapper");
    py::class_<ObjectiveFunction2DPtrWrapper>(m, "CppObjectiveFunction2DPtrWrapper");
}

#endif //PYTHONLPT_BINDINGS_CPP
