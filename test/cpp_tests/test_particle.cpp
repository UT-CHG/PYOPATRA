//
// Created by Georgia Stuart on 2/5/21.
//

#include <catch2/catch.hpp>
#include <stdexcept>
#include <iostream>
#include "PYOPATRA/cpp/particle.h"
#include "PYOPATRA/cpp/mesh/mesh_element.h"

TEST_CASE("Particles Initialized Correctly", "[particle-constructors]") {
    SECTION("Default Constructor") {
        Particle3D p;

        REQUIRE(p.get_diameter() == 0.0);
        REQUIRE(p.get_density() == 0.0);
        REQUIRE(p.get_interfacial_tension() == 0.0);
        REQUIRE(p.get_depth_index() == 0);

        Vector3d location = p.get_location();
        REQUIRE(location(0) == Approx(0.0));
        REQUIRE(location(1) == Approx(0.0));
        REQUIRE(location(2) == Approx(0.0));
    }

    SECTION("Full Constructor") {
        Particle3D p(10.0, -20.0, 0.0005, 858.0, -15.0, 0.023);
        REQUIRE(p.get_diameter() == 0.0005);
        REQUIRE(p.get_density() == 858.0);
        REQUIRE(p.get_interfacial_tension() == 0.023);
        REQUIRE(p.get_depth_index() == 0);

        Vector3d location = p.get_location();
        REQUIRE(location[0] == 10.0);
        REQUIRE(location[1] == -20.0);
        REQUIRE(location[2] == -15.0);
    }
}

TEST_CASE("Buoyancy Calculated Correctly", "[buoyancy]") {

    double fluid_density = 998.2071;
    double fluid_visosity = 1.308767294838539e-03;
    double water_viscosity = 1305.90172775e-6;

    SECTION("Reynold's Numbers") {
        REQUIRE(Particle3D::calculate_reynolds(50.0) == Approx(1.72923));
        REQUIRE(Particle3D::calculate_reynolds(400.0) == Approx(9.50261));
        REQUIRE(Particle3D::calculate_reynolds(1000.0) == Approx(19.0142));
//        REQUIRE_THROWS_AS(Particle::calculate_reynolds(1.0e8), std::runtime_error);
    }

    SECTION("Test configuration 1, small particle, oil in water") {
        Particle3D particle(0.0, 0.0, 0.0005, 858.0, -15.0, 0.023);

        REQUIRE(particle.calculate_nd(fluid_density, fluid_visosity) == Approx(133.45677288001536));
        REQUIRE(particle.calculate_reynolds(particle.calculate_nd(fluid_density, fluid_visosity)) == Approx(4.008095982211944));
        REQUIRE(particle.terminal_buoyancy_velocity(fluid_density, fluid_visosity, water_viscosity) == Approx(0.010510173562365449));
        REQUIRE(particle.calculate_morton_number(fluid_density, fluid_visosity) == Approx(3.3252247175369584e-10));
        REQUIRE(particle.calculate_eotvos_number(0.0005, fluid_density) == Approx(0.01493510413043478));
        REQUIRE(particle.calculate_eotvos_number(0.015, fluid_density) == Approx(13.441593717391305));
        REQUIRE(particle.calculate_diameter_from_H(59.3, particle.calculate_morton_number(fluid_density, fluid_visosity), fluid_density, fluid_visosity, water_viscosity) == Approx(0.005368663147596897));
        REQUIRE(particle.calculate_critical_diameter(fluid_density, fluid_visosity, water_viscosity) == Approx(0.01798749114892274));
    }

    SECTION("Test configuration 2, ellpsoid particle, oil in water") {
        Particle3D particle(0.0, 0.0, 0.01, 858.0, -15.0, 0.023);



        REQUIRE(particle.calculate_morton_number(fluid_density, fluid_visosity) == Approx(3.3252247175369584e-10));
        REQUIRE(particle.calculate_eotvos_number(0.01, fluid_density) == Approx(5.974041652173913));
        REQUIRE(particle.terminal_buoyancy_velocity(fluid_density, fluid_visosity, water_viscosity) == Approx(0.11846135362594433));
    }

    SECTION("Test configuration 2, spherical cap particle, oil in water") {
        Particle3D particle(0.0, 0.0, 0.02, 858.0, -15.0, 0.023);

        REQUIRE(particle.calculate_morton_number(fluid_density, fluid_visosity) == Approx(3.3252247175369584e-10));
        REQUIRE(particle.calculate_eotvos_number(0.02, fluid_density) == Approx(23.89616660869565));
        REQUIRE(particle.terminal_buoyancy_velocity(fluid_density, fluid_visosity, water_viscosity) == Approx(0.11797019910948671));
    }
}