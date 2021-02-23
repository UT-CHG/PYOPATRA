//
// Created by Georgia Stuart on 2/21/21.
//

#include <catch2/catch.hpp>
#include "../src/PYOPATRA/cpp/mesh.h"

TEST_CASE("Mesh Tests", "[mesh-tests]") {
    SECTION("Water Viscosity Tests") {
        // Values from Huber et al (2009)
        REQUIRE(MeshNode::calculate_fluid_viscosity(298.15, 998) == Approx(0.000889735100));
        REQUIRE(MeshNode::calculate_fluid_viscosity(298.15, 1200) == Approx(0.001437649467));
        REQUIRE(MeshNode::calculate_fluid_viscosity(373.15, 1000) == Approx(0.000307883622));
        REQUIRE(MeshNode::calculate_fluid_viscosity(433.15, 1) == Approx(0.000014538324));
        REQUIRE(MeshNode::calculate_fluid_viscosity(433.15, 1000) == Approx(0.000217685358));
        REQUIRE(MeshNode::calculate_fluid_viscosity(873.15, 1) == Approx(0.000032619287));
        REQUIRE(MeshNode::calculate_fluid_viscosity(873.15, 100) == Approx(0.000035802262));
        REQUIRE(MeshNode::calculate_fluid_viscosity(873.15, 600) == Approx(0.000077430195));
        REQUIRE(MeshNode::calculate_fluid_viscosity(1173.15, 1) == Approx(0.000044217245));
        REQUIRE(MeshNode::calculate_fluid_viscosity(1173.15, 100) == Approx(0.000047640433));
        REQUIRE(MeshNode::calculate_fluid_viscosity(1173.15, 400) == Approx(0.000064154608));
    }

    SECTION("Pure water viscosity tests") {
        REQUIRE(MeshNode::calculate_pure_water_viscosity(10.0 + 273.15) == Approx(1305.90172775e-6));
        REQUIRE(MeshNode::calculate_pure_water_viscosity(298.15) == Approx(889.996773679e-6));

    }

    SECTION("Mesh Node Constructor Test") {
        MeshNode node1(25);

        REQUIRE(node1.num_depth_layers == 25);
        REQUIRE(node1.density.size() == 25);
        REQUIRE(node1.temperature.size() == 25);
        REQUIRE(node1.water_viscosity.size() == 25);
        REQUIRE(node1.viscosity.size() == 25);
        REQUIRE(node1.velocity.x == 0.0);
        REQUIRE(node1.velocity.y == 0.0);
        REQUIRE(node1.velocity.z == 0.0);
        REQUIRE(node1.location.x == 0.0);
        REQUIRE(node1.location.y == 0.0);
        REQUIRE(node1.location.z == 0.0);

        MeshNode node2;

        REQUIRE(node2.num_depth_layers == 0);
        REQUIRE(node2.density.empty());
        REQUIRE(node2.temperature.empty());
        REQUIRE(node2.water_viscosity.empty());
        REQUIRE(node2.viscosity.empty());
        REQUIRE(node1.velocity.x == 0.0);
        REQUIRE(node1.velocity.y == 0.0);
        REQUIRE(node1.velocity.z == 0.0);
        REQUIRE(node1.location.x == 0.0);
        REQUIRE(node1.location.y == 0.0);
        REQUIRE(node1.location.z == 0.0);
    }
}