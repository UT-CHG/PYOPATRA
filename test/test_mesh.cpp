//
// Created by Georgia Stuart on 2/21/21.
//

#include <catch2/catch.hpp>
#include "PYOPATRA/cpp/mesh/mesh_element.h"

TEST_CASE("Mesh Tests", "[mesh-tests]") {
    SECTION("Water Viscosity Tests") {
        // Values from Huber et al (2009)
        REQUIRE(MeshElement::calculate_fluid_viscosity(298.15, 998) == Approx(0.000889735100));
        REQUIRE(MeshElement::calculate_fluid_viscosity(298.15, 1200) == Approx(0.001437649467));
        REQUIRE(MeshElement::calculate_fluid_viscosity(373.15, 1000) == Approx(0.000307883622));
        REQUIRE(MeshElement::calculate_fluid_viscosity(433.15, 1) == Approx(0.000014538324));
        REQUIRE(MeshElement::calculate_fluid_viscosity(433.15, 1000) == Approx(0.000217685358));
        REQUIRE(MeshElement::calculate_fluid_viscosity(873.15, 1) == Approx(0.000032619287));
        REQUIRE(MeshElement::calculate_fluid_viscosity(873.15, 100) == Approx(0.000035802262));
        REQUIRE(MeshElement::calculate_fluid_viscosity(873.15, 600) == Approx(0.000077430195));
        REQUIRE(MeshElement::calculate_fluid_viscosity(1173.15, 1) == Approx(0.000044217245));
        REQUIRE(MeshElement::calculate_fluid_viscosity(1173.15, 100) == Approx(0.000047640433));
        REQUIRE(MeshElement::calculate_fluid_viscosity(1173.15, 400) == Approx(0.000064154608));

        REQUIRE(MeshElement::calculate_fluid_viscosity(10.0 + 273.15, 998.2071) == Approx(0.00130877));
    }

    SECTION("Pure water viscosity tests") {
        REQUIRE(MeshElement::calculate_pure_water_viscosity(10.0 + 273.15) == Approx(1305.90172775e-6));
        REQUIRE(MeshElement::calculate_pure_water_viscosity(298.15) == Approx(889.996773679e-6));

    }

    SECTION("Mesh Node Constructor Test") {
        MeshElement node1(50, 25);

        REQUIRE(node1.mesh_index == 50);
        REQUIRE(node1.num_depth_layers == 25);
        REQUIRE(node1.density.size() == 25);
        REQUIRE(node1.temperature.size() == 25);
        REQUIRE(node1.water_viscosity.size() == 25);
        REQUIRE(node1.viscosity.size() == 25);
        REQUIRE(node1.velocity.x == 0.0);
        REQUIRE(node1.velocity.y == 0.0);
        REQUIRE(node1.velocity.z == 0.0);
        REQUIRE(node1.location.latitude == 0.0);
        REQUIRE(node1.location.longitude == 0.0);

        MeshElement node2;

        REQUIRE(node2.mesh_index == 0);
        REQUIRE(node2.num_depth_layers == 0);
        REQUIRE(node2.density.empty());
        REQUIRE(node2.temperature.empty());
        REQUIRE(node2.water_viscosity.empty());
        REQUIRE(node2.viscosity.empty());
        REQUIRE(node2.velocity.x == 0.0);
        REQUIRE(node2.velocity.y == 0.0);
        REQUIRE(node2.velocity.z == 0.0);
        REQUIRE(node1.location.latitude == 0.0);
        REQUIRE(node1.location.longitude == 0.0);

        MeshElement node3(50);

        REQUIRE(node3.mesh_index == 50);
        REQUIRE(node3.num_depth_layers == 0);
        REQUIRE(node3.density.empty());
        REQUIRE(node3.temperature.empty());
        REQUIRE(node3.water_viscosity.empty());
        REQUIRE(node3.viscosity.empty());
        REQUIRE(node3.velocity.x == 0.0);
        REQUIRE(node3.velocity.y == 0.0);
        REQUIRE(node3.velocity.z == 0.0);
        REQUIRE(node1.location.latitude == 0.0);
        REQUIRE(node1.location.longitude == 0.0);

        MeshElement node4(50, 998, 298.15);
        REQUIRE(node4.mesh_index == 50);
        REQUIRE(node4.num_depth_layers == 1);
        REQUIRE(node4.density[0] == 998);
        REQUIRE(node4.temperature[0] == 298.15);
        REQUIRE(node4.water_viscosity[0] == Approx(889.996773679e-6));
        REQUIRE(node4.viscosity[0] == Approx(0.000889735100));
        REQUIRE(node4.velocity.x == 0.0);
        REQUIRE(node4.velocity.y == 0.0);
        REQUIRE(node4.velocity.z == 0.0);
        REQUIRE(node1.location.latitude == 0.0);
        REQUIRE(node1.location.longitude == 0.0);
    }
}