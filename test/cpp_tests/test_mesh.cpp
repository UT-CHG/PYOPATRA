//
// Created by Georgia Stuart on 5/3/21
//

#include <catch2/catch.hpp>
#include "PYOPATRA/cpp/mesh/mesh.h"

TEST_CASE("Mesh Tests", "[mesh-tests]") {
    TriangularMesh2D mesh(12, 12, {0, 3, 6, 9}, {}, 0.03);

    REQUIRE(mesh.get_num_vertices() == 12);
    REQUIRE(mesh.get_num_columns() == 12);
    REQUIRE(mesh.get_water_column_pointer(1) != nullptr);

    mesh.set_water_column_adjacency(0, 1, 1);
    mesh.set_water_column_adjacency(0, 1, 1);

}