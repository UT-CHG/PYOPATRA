//
// Created by Georgia Stuart on 3/8/21.
//

#include <catch2/catch.hpp>

#include "../src/PYOPATRA/cpp/illnode.h"
#include "../src/PYOPATRA/cpp/list.h"
#include "../src/PYOPATRA/cpp/particle.h"

TEST_CASE("Intrusive Linked List", "[ill]") {

    SECTION("Test particle") {
        Particle p(10.0, -20.0, 0.0005, 858.0, -15.0, 0.023);
        Particle p2(20.0, -30.0, 0.0005, 858.0, -15.0, 0.023);
        Particle p3(30.0, -40.0, 0.0005, 858.0, -15.0, 0.023);


        REQUIRE(p.node.owner == &p);
        REQUIRE(p.node.next == nullptr);
        REQUIRE(p.node.prev == nullptr);

        p2.node.insert_after(p.node);
        p3.node.insert_after(p2.node);

        REQUIRE(p.node.prev == nullptr);
        REQUIRE(p.node.next == &p2.node);
        REQUIRE(p2.node.prev == &p.node);
        REQUIRE(p2.node.next == &p3.node);
        REQUIRE(p3.node.prev == &p2.node);
        REQUIRE(p3.node.next == nullptr);
    }
}