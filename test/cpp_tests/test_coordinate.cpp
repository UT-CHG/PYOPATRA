//
// Created by Georgia Stuart on 3/16/21.
//

#include <catch2/catch.hpp>
#include "PYOPATRA/cpp/coordinate.h"

TEST_CASE("Meters to degrees latitude and Longitude", "[latlong]") {
    double latitude = 30.266666;

    REQUIRE(meters_to_latitude(100) == Approx(0.0009));
    REQUIRE(meters_to_longitude(100, latitude) == Approx(0.0010420429));
}