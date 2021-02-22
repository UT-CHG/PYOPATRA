//
// Created by Georgia Stuart on 2/5/21.
//

#include <catch2/catch.hpp>
#include <stdexcept>
#include <iostream>
#include "../src/pythonLPT/cpp/particle.h"
#include "../src/pythonLPT/cpp/mesh.h"

TEST_CASE("Buoyancy Calculated Correctly", "[buoyancy]") {
    SECTION("Reynold's Numbers") {
        REQUIRE(Particle::calculate_reynolds(50.0) == Approx(1.72923));
        REQUIRE(Particle::calculate_reynolds(400.0) == Approx(9.50261));
        REQUIRE(Particle::calculate_reynolds(1000.0) == Approx(19.0142));
//        REQUIRE_THROWS_AS(Particle::calculate_reynolds(1.0e8), std::runtime_error);
    }

    SECTION("Water Viscosity") {
//        REQUIRE(Particle::calculate_fluid_viscosity(10.0) == Approx(0.001308));
    }

    SECTION("Test configuration 1, small particle, oil in water") {
        Particle particle(0.0, 0.0, 0.00005, 858.0, -15.0, 0.023);
        MeshNode node(998.2071, 10.0 + 273.15);
        particle.current_node = &node;

//        REQUIRE(particle.calculate_nd() == Approx(0.227683));
//        REQUIRE(particle.calculate_reynolds(particle.calculate_nd()) == Approx(0.009477679));
//        REQUIRE(particle.terminal_buoyancy_velocity() == Approx(0.000190274));
//        REQUIRE(particle.calculate_morton_number() == Approx(1.142463e-10));
//        REQUIRE(particle.calculate_eotvos_number() == Approx(0.00014935));
//        REQUIRE(particle.calculate_critical_diameter() == Approx())
    }
}