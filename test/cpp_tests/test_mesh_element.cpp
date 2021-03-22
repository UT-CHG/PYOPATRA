//
// Created by Georgia Stuart on 3/13/21.
//

#include <catch2/catch.hpp>
#include <iostream>
#include "PYOPATRA/cpp/mesh/mesh_vertex.h"
#include "PYOPATRA/cpp/mesh/mesh_element.h"


TEST_CASE("Triangle Mesh Element Tests", "[triangle-mesh-element-tests]") {
    MeshVertex3D v1(-89.4717160000, 30.0866470000, 3.0000000000, 998, 273.15, {5, 6, 7});
    MeshVertex3D v2(-89.4966660000, 30.0901150000, 3.9491493700, 1005, 275.15, {8, 9, 10});
    MeshVertex3D v3(-89.4889400000, 30.0744220000, 1.9666023300, 1010, 277.15, {11, 12, 13});

    Coordinate3D p(-89.486445 , 30.0822685, 4.0);
    Coordinate3D p2(-89.486445 , 30.0822685, 1.0);
    Coordinate3D bc(0.27354909, 0.2869026 , 0.43954832);

    TriangularMeshElement3D mesh_element(&v1, &v2, &v3, 0);

    REQUIRE(*(mesh_element.get_vertices()[0]) == v1);
    REQUIRE(*(mesh_element.get_vertices()[1]) == v2);
    REQUIRE(*(mesh_element.get_vertices()[2]) == v3);

    SECTION("Barycentric Coordinate Tests") {
        REQUIRE((mesh_element.calculate_barycentric_coordinate(p) - bc).norm() < 10e-6);
        REQUIRE(mesh_element.calculate_depth_at_point(p) == Approx(2.81808521));

        REQUIRE(mesh_element.sample_density_at_barycentric_coordinate(bc) == Approx(1005.2828979957392));
        REQUIRE(mesh_element.sample_viscosity_at_barycentric_coordinate(bc) == Approx(
                bc[0] * mesh_element.get_vertices()[0]->get_viscosity() +
                bc[1] * mesh_element.get_vertices()[1]->get_viscosity() +
                bc[2] * mesh_element.get_vertices()[2]->get_viscosity()));
        REQUIRE(mesh_element.sample_water_viscosity_at_barycentric_coordinate(bc) == Approx(
                bc[0] * mesh_element.get_vertices()[0]->get_water_viscosity() +
                bc[1] * mesh_element.get_vertices()[1]->get_water_viscosity() +
                bc[2] * mesh_element.get_vertices()[2]->get_water_viscosity()));

        Vector3d bc_velocity = bc[0] * mesh_element.get_vertices()[0]->get_velocity() +
                               bc[1] * mesh_element.get_vertices()[1]->get_velocity() +
                               bc[2] * mesh_element.get_vertices()[2]->get_velocity();

        REQUIRE((mesh_element.sample_velocity_at_barycentric_coordinate(bc) - bc_velocity).norm() < 10e-6);


        REQUIRE(mesh_element.check_halfspace(p) == 1);
        REQUIRE(mesh_element.check_halfspace(p2) == -1);
    }

    SECTION("Cursor Checks") {
        TriangularMeshCursor3D cursor(&mesh_element);
        InterpolatedValues interp;

        REQUIRE(cursor.get_depth_at_point(p) == Approx(2.81808521));
        REQUIRE(cursor.check_halfspace(p) == 1);
        REQUIRE(cursor.check_halfspace(p2) == -1);

        cursor.get_interpolated_values(p, interp);

        REQUIRE(interp.density == Approx(1005.2828979957392));
        REQUIRE(interp.viscosity == Approx(
                bc[0] * mesh_element.get_vertices()[0]->get_viscosity() +
                bc[1] * mesh_element.get_vertices()[1]->get_viscosity() +
                bc[2] * mesh_element.get_vertices()[2]->get_viscosity()));
        REQUIRE(interp.water_viscosity == Approx(
                bc[0] * mesh_element.get_vertices()[0]->get_water_viscosity() +
                bc[1] * mesh_element.get_vertices()[1]->get_water_viscosity() +
                bc[2] * mesh_element.get_vertices()[2]->get_water_viscosity()));

        Vector3d bc_velocity = bc[0] * mesh_element.get_vertices()[0]->get_velocity() +
                               bc[1] * mesh_element.get_vertices()[1]->get_velocity() +
                               bc[2] * mesh_element.get_vertices()[2]->get_velocity();

        REQUIRE((interp.velocity - bc_velocity).norm() < 10e-6);
    }
}