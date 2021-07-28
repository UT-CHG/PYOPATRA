//
// Created by Georgia Stuart on 7/26/21.
//

#include <catch2/catch.hpp>
#include "PYOPATRA/cpp/inversion_tools/objective_functions.h"
#include "PYOPATRA/cpp/particle_list.h"

TEST_CASE("Objective Tests", "[objective-tests]") {
    SlicedWassersteinDistance2D obj(5, 10, {20, 25, -90, -80}, 20, 5000);

    Eigen::VectorXd lat_bins(5), lon_bins(10);
    lat_bins << 20, 21, 22, 23, 24;
    lon_bins << -90, -89, -88, -87, -86, -85, -84, -83, -82, -81;

    REQUIRE((lat_bins - obj.get_latitude_bounds()).norm() < 10e-6);
    REQUIRE((lon_bins - obj.get_longitude_bounds()).norm() < 10e-6);

    ParticleList2D p_list;

    p_list.create_particle({21.5, -88.5});
    p_list.create_particle({24.5, -80.5});
    p_list.create_particle({20.5, -90});
    p_list.create_particle({25, -88.5});

    Eigen::MatrixXd particle_locations(4, 2);
    particle_locations << 21.5, -88.5,
                          24.5, -80.5,
                          20.5, -90.0,
                          25.0, -88.5;

    ParticleList2D p_list2;



    p_list2.create_particle({21.5, -88.5});
    p_list2.create_particle({24.5, -80.5});
    p_list2.create_particle({20.5, -90});
    p_list2.create_particle({23.5, -88.5});

    obj.set_observed_values(particle_locations);

    Eigen::MatrixXd bin(10, 5);
    bin << 0.25, 0, 0, 0, 0,
           0, 0.25, 0, 0, 0.25,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0.25;

    REQUIRE((bin - obj.get_observed_bins()).norm() < 10e-6);
    REQUIRE(obj.calculate_value(p_list) == Approx(0.0));
    REQUIRE(obj.calculate_value(p_list2) > 0.0);

}