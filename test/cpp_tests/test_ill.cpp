//
// Created by Georgia Stuart on 3/8/21.
//

#include <catch2/catch.hpp>

#include "PYOPATRA/cpp/illnode.h"
#include "PYOPATRA/cpp/list.h"
#include "PYOPATRA/cpp/particle.h"

TEST_CASE("Intrusive Linked List", "[ill]") {

    SECTION("Test particle") {
        Particle3D p(10.0, -20.0, 0.0005, 858.0, -15.0, 0.023);
        Particle3D p2(20.0, -30.0, 0.0005, 858.0, -15.0, 0.023);
        Particle3D p3(30.0, -40.0, 0.0005, 858.0, -15.0, 0.023);


        REQUIRE(p.get_node().owner == &p);
        REQUIRE(p.get_node().next == nullptr);
        REQUIRE(p.get_node().prev == nullptr);

        p2.get_node().insert_after(p.get_node());
        p3.get_node().insert_after(p2.get_node());

        REQUIRE(p.get_node().prev == nullptr);
        REQUIRE(p.get_node().next == &p2.get_node());
        REQUIRE(p2.get_node().prev == &p.get_node());
        REQUIRE(p2.get_node().next == &p3.get_node());
        REQUIRE(p3.get_node().prev == &p2.get_node());
        REQUIRE(p3.get_node().next == nullptr);

        ParticleBase<3> *particle = &p;
        int len = 1;

        while (particle->get_node().next != nullptr) {
            particle = particle->get_node().next->owner;
            len++;
        }

        REQUIRE(len == 3);

        p2.get_node().remove();

        REQUIRE(p.get_node().next == &p3.get_node());
        REQUIRE(p3.get_node().prev == &p.get_node());
        REQUIRE(p2.get_node().next == nullptr);
        REQUIRE(p2.get_node().prev == nullptr);

        p2.get_node().insert_after(p.get_node());

        REQUIRE(p.get_node().prev == nullptr);
        REQUIRE(p.get_node().next == &p2.get_node());
        REQUIRE(p2.get_node().prev == &p.get_node());
        REQUIRE(p2.get_node().next == &p3.get_node());
        REQUIRE(p3.get_node().prev == &p2.get_node());
        REQUIRE(p3.get_node().next == nullptr);
    }
}