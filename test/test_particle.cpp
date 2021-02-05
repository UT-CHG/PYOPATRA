//
// Created by Georgia Stuart on 2/5/21.
//

#include <catch2/catch.hpp>
#include "../src/pythonLPT/cpp/particle.h"
#include "../src/pythonLPT/cpp/mesh.h"

TEST_CASE("Buoyancy Calculated Correctly", "[buoyancy]") {
    Particle particle(0.0, 0.0, 0.00005, 858.0, -15.0, 0.023);
    MeshNode node(998.2071, 10.0, 0.001002);

    SECTION("Eotvos number check") {

    }
}