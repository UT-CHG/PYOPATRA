//
// Created by Georgia Stuart on 2/21/21.
//

#include <catch2/catch.hpp>
#include "../src/pythonLPT/cpp/mesh.h"

TEST_CASE("Mesh Tests", "[mesh-tests]") {
    SECTION("Water Viscosity Tests") {
        REQUIRE(MeshNode::calculate_fluid_viscosity(298.15, 998) == Approx(0.000889735100));
        REQUIRE(MeshNode::calculate_fluid_viscosity(298.15, 1200) == Approx(0.001437649467));
        REQUIRE(MeshNode::calculate_fluid_viscosity(373.15, 1000) == Approx(0.000307883622));
    }
}