//
// Created by Georgia Stuart on 2/5/21.
//

#include <catch2/catch.hpp>
#include <stdexcept>
#include <iostream>
#include "../src/PYOPATRA/cpp/particle.h"
#include "PYOPATRA/cpp/mesh/mesh_element.h"

TEST_CASE("Particles Initialized Correctly", "[particle-constructors]") {
    SECTION("Default Constructor") {
        Particle p;

        REQUIRE(p.diameter == 0.0);
        REQUIRE(p.density == 0.0);
        REQUIRE(p.interfacial_tension == 0.0);
        REQUIRE(p.depth_index == 0);
        REQUIRE(p.current_mesh_node == nullptr);
        REQUIRE(p.position.x == 0.0);
        REQUIRE(p.position.y == 0.0);
        REQUIRE(p.position.z == 0.0);
    }

    SECTION("Full Constructor") {
    Particle p(10.0, -20.0, 0.0005, 858.0, -15.0, 0.023);
        REQUIRE(p.diameter == 0.0005);
        REQUIRE(p.density == 858.0);
        REQUIRE(p.interfacial_tension == 0.023);
        REQUIRE(p.depth_index == 0);
        REQUIRE(p.current_mesh_node == nullptr);
        REQUIRE(p.position.x == 10.0);
        REQUIRE(p.position.y == -20.0);
        REQUIRE(p.position.z == -15.0);
    }
}

TEST_CASE("Buoyancy Calculated Correctly", "[buoyancy]") {
    SECTION("Reynold's Numbers") {
        REQUIRE(Particle::calculate_reynolds(50.0) == Approx(1.72923));
        REQUIRE(Particle::calculate_reynolds(400.0) == Approx(9.50261));
        REQUIRE(Particle::calculate_reynolds(1000.0) == Approx(19.0142));
        REQUIRE_THROWS_AS(Particle::calculate_reynolds(1.0e8), std::runtime_error);
    }

    SECTION("Test configuration 1, small particle, oil in water") {
        Particle particle(0.0, 0.0, 0.0005, 858.0, -15.0, 0.023);
        MeshElement node(0, 998.2071, 10.0 + 273.15);
        particle.current_mesh_node = &node;


        REQUIRE(particle.calculate_nd() == Approx(133.45677288001536));
        REQUIRE(particle.calculate_reynolds(particle.calculate_nd()) == Approx(4.008095982211944));
        REQUIRE(particle.terminal_buoyancy_velocity() == Approx(0.010510173562365449));
        REQUIRE(particle.calculate_morton_number() == Approx(3.3252247175369584e-10));
        REQUIRE(particle.calculate_eotvos_number(0.0005) == Approx(0.01493510413043478));
        REQUIRE(particle.calculate_eotvos_number(0.015) == Approx(13.441593717391305));
        REQUIRE(particle.calculate_diameter_from_H(59.3, particle.calculate_morton_number()) == Approx(0.005368663147596897));
        REQUIRE(particle.calculate_critical_diameter() == Approx(0.01798749114892274));
    }

    SECTION("Test configuration 2, ellpsoid particle, oil in water") {
        Particle particle(0.0, 0.0, 0.01, 858.0, -15.0, 0.023);
        MeshElement node(0, 998.2071, 10.0 + 273.15);
        particle.current_mesh_node = &node;

        REQUIRE(particle.calculate_morton_number() == Approx(3.3252247175369584e-10));
        REQUIRE(particle.calculate_eotvos_number(0.01) == Approx(5.974041652173913));
        REQUIRE(particle.terminal_buoyancy_velocity() == Approx(0.11846135362594433));
    }

    SECTION("Test configuration 2, spherical cap particle, oil in water") {
        Particle particle(0.0, 0.0, 0.02, 858.0, -15.0, 0.023);
        MeshElement node(0, 998.2071, 10.0 + 273.15);
        particle.current_mesh_node = &node;

        REQUIRE(particle.calculate_morton_number() == Approx(3.3252247175369584e-10));
        REQUIRE(particle.calculate_eotvos_number(0.02) == Approx(23.89616660869565));
        REQUIRE(particle.terminal_buoyancy_velocity() == Approx(0.11797019910948671));
    }
}