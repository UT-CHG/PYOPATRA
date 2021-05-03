//
// Created by Georgia Stuart on 2/21/21.
//

#include <catch2/catch.hpp>
#include "PYOPATRA/cpp/mesh/mesh_vertex.h"

TEST_CASE("Mesh Vertex Tests", "[mesh-vertex-tests]") {
    SECTION("Water Viscosity Tests") {
        // Values from Huber et al (2009)
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(298.15, 998) == Approx(0.000889735100));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(298.15, 1200) == Approx(0.001437649467));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(373.15, 1000) == Approx(0.000307883622));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(433.15, 1) == Approx(0.000014538324));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(433.15, 1000) == Approx(0.000217685358));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(873.15, 1) == Approx(0.000032619287));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(873.15, 100) == Approx(0.000035802262));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(873.15, 600) == Approx(0.000077430195));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(1173.15, 1) == Approx(0.000044217245));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(1173.15, 100) == Approx(0.000047640433));
        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(1173.15, 400) == Approx(0.000064154608));

        REQUIRE(MeshVertex3D::calculate_fluid_viscosity(10.0 + 273.15, 998.2071) == Approx(0.00130877));
    }

    SECTION("Pure water viscosity tests") {
        REQUIRE(MeshVertex3D::calculate_pure_water_viscosity(10.0 + 273.15) == Approx(1305.90172775e-6));
        REQUIRE(MeshVertex3D::calculate_pure_water_viscosity(298.15) == Approx(889.996773679e-6));

    }

    SECTION("Mesh Vertex Constructor Test") {
        MeshVertex3D node1(50, 25, 0.0, 1);

        REQUIRE(node1.get_density()[0] == 0.0);
        REQUIRE(node1.get_temperature()[0] == 0.0);
        REQUIRE(node1.get_water_viscosity()[0] == 0.0);
        REQUIRE(node1.get_viscosity()[0] == 0.0);
        REQUIRE(node1.get_latitude() == 50);
        REQUIRE(node1.get_longitude() == 25);
        REQUIRE(node1.get_depth() == 0.0);

        node1.set_location({-70, 30, 20.0});
        node1.set_velocity({-0.5, 0.5, 0.25}, 0);
        REQUIRE(node1.get_latitude() == -70);
        REQUIRE(node1.get_longitude() == 30);
        REQUIRE(node1.get_depth() == 20.0);

        MeshVertex3D::Vector v = {-0.5, 0.5, 0.25};

        REQUIRE(node1.get_velocity()[0].isApprox(v));

        MeshVertex3D node2(1);

        REQUIRE(node2.get_density()[0] == 0.0);
        REQUIRE(node2.get_temperature()[0] == 0.0);
        REQUIRE(node2.get_water_viscosity()[0] == 0.0);
        REQUIRE(node2.get_viscosity()[0] == 0.0);
        REQUIRE(node2.get_latitude() == 0.0);
        REQUIRE(node2.get_longitude() == 0.0);
        REQUIRE(node2.get_depth() == 0.0);

        MeshVertex3D node3(-89.4960990000, 30.1932510000, -0.7478574510, 1);
        node3.set_temperature_and_density(283.15, 998.2071, 0);

        REQUIRE(node3.get_density()[0] == 998.2071);
        REQUIRE(node3.get_temperature()[0] == 283.15);
        REQUIRE(node3.get_water_viscosity()[0] == Approx(1305.90172775e-6));
        REQUIRE(node3.get_viscosity()[0] == Approx(0.00130877));
        REQUIRE(node3.get_latitude() == -89.4960990000);
        REQUIRE(node3.get_longitude() == 30.1932510000);
        REQUIRE(node3.get_depth() == -0.7478574510);

        MeshVertex3D node4(-89.4960990000, 30.1932510000, -0.7478574510, 1);

        REQUIRE(node4.get_density()[0] == 0.0);
        REQUIRE(node4.get_temperature()[0] == 0.0);
        REQUIRE(node4.get_water_viscosity()[0] == 0.0);
        REQUIRE(node4.get_viscosity()[0] == 0.0);
        REQUIRE(node4.get_latitude() == -89.4960990000);
        REQUIRE(node4.get_longitude() == 30.1932510000);
        REQUIRE(node4.get_depth() == -0.7478574510);

        node4.set_temperature_and_density(283.15, 998.2071, 0);

        REQUIRE(node4.get_density()[0] == 998.2071);
        REQUIRE(node4.get_temperature()[0] == 283.15);
    }
}