//
// Created by Georgia Stuart on 2/21/21.
//

#include <catch2/catch.hpp>
#include "PYOPATRA/cpp/mesh/mesh_vertex.h"

TEST_CASE("Mesh Vertex Tests", "[mesh-vertex-tests]") {
    SECTION("Water Viscosity Tests") {
        // Values from Huber et al (2009)
        REQUIRE(MeshVertex::calculate_fluid_viscosity(298.15, 998) == Approx(0.000889735100));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(298.15, 1200) == Approx(0.001437649467));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(373.15, 1000) == Approx(0.000307883622));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(433.15, 1) == Approx(0.000014538324));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(433.15, 1000) == Approx(0.000217685358));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(873.15, 1) == Approx(0.000032619287));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(873.15, 100) == Approx(0.000035802262));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(873.15, 600) == Approx(0.000077430195));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(1173.15, 1) == Approx(0.000044217245));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(1173.15, 100) == Approx(0.000047640433));
        REQUIRE(MeshVertex::calculate_fluid_viscosity(1173.15, 400) == Approx(0.000064154608));

        REQUIRE(MeshVertex::calculate_fluid_viscosity(10.0 + 273.15, 998.2071) == Approx(0.00130877));
    }

    SECTION("Pure water viscosity tests") {
        REQUIRE(MeshVertex::calculate_pure_water_viscosity(10.0 + 273.15) == Approx(1305.90172775e-6));
        REQUIRE(MeshVertex::calculate_pure_water_viscosity(298.15) == Approx(889.996773679e-6));

    }

    SECTION("Mesh Vertex Constructor Test") {
        MeshVertex node1(50, 25);

        REQUIRE(node1.get_density().size() == 1);
        REQUIRE(node1.get_temperature().size() == 1);
        REQUIRE(node1.get_water_viscosity().size() == 1);
        REQUIRE(node1.get_viscosity().size() == 1);
        REQUIRE(node1.get_latitude() == 50);
        REQUIRE(node1.get_longitude() == 25);
        REQUIRE(node1.get_bathymetric_depth() == 0.0);

        MeshVertex node2;

        REQUIRE(node2.get_density().size() == 1);
        REQUIRE(node2.get_temperature().size() == 1);
        REQUIRE(node2.get_water_viscosity().size() == 1);
        REQUIRE(node2.get_viscosity().size() == 1);
        REQUIRE(node2.get_latitude() == 0.0);
        REQUIRE(node2.get_longitude() == 0.0);
        REQUIRE(node2.get_bathymetric_depth() == 0.0);

        MeshVertex node3(50, 25, 33.7, 12);

        REQUIRE(node3.get_density().size() == 12);
        REQUIRE(node3.get_temperature().size() == 12);
        REQUIRE(node3.get_water_viscosity().size() == 12);
        REQUIRE(node3.get_viscosity().size() == 12);
        REQUIRE(node3.get_latitude() == 50.0);
        REQUIRE(node3.get_longitude() == 25.0);
        REQUIRE(node3.get_bathymetric_depth() == 33.7);

//        MeshVertex node4(50, 998, 298.15);
//        REQUIRE(node4.mesh_index == 50);
//        REQUIRE(node4.num_depth_layers == 1);
//        REQUIRE(node4.density[0] == 998);
//        REQUIRE(node4.temperature[0] == 298.15);
//        REQUIRE(node4.water_viscosity[0] == Approx(889.996773679e-6));
//        REQUIRE(node4.viscosity[0] == Approx(0.000889735100));
//        REQUIRE(node4.velocity.x == 0.0);
//        REQUIRE(node4.velocity.y == 0.0);
//        REQUIRE(node4.velocity.z == 0.0);
//        REQUIRE(node1.location.latitude == 0.0);
//        REQUIRE(node1.location.longitude == 0.0);
    }
}