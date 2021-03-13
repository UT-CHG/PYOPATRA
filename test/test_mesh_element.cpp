//
// Created by Georgia Stuart on 3/13/21.
//

#include <catch2/catch.hpp>
#include <iostream>
#include "PYOPATRA/cpp/mesh/mesh_vertex.h"
#include "PYOPATRA/cpp/mesh/mesh_element.h"


TEST_CASE("Triangle Mesh Element Tests", "[triangle-mesh-element-tests]") {
    MeshVertex v1(-89.4717160000, 30.0866470000, 3.0000000000);
    MeshVertex v2(-89.4966660000, 30.0901150000, 3.9491493700);
    MeshVertex v3(-89.4889400000, 30.0744220000, 1.9666023300);

    TriangularMeshElement3D mesh_element(&v1, &v2, &v3);

    REQUIRE(*(mesh_element.get_vertices()[0]) == v1);
    REQUIRE(*(mesh_element.get_vertices()[1]) == v2);
    REQUIRE(*(mesh_element.get_vertices()[2]) == v3);

    SECTION("Barycentric Coordinate Tests") {
        Coordinate3D p(-89.486445 , 30.0822685, 4.0);
        Coordinate3D bc(0.27354909, 0.2869026 , 0.43954832);

        REQUIRE((mesh_element.calculate_barycentric_coordinate(p) - bc).norm() < 10e-6);
    }
}