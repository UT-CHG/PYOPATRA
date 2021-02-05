//
// Created by Georgia Stuart on 2/3/21.
//

#include <math.h>
#include "particle.h"

// Acceleration due to gravity in km/s^2
#define GRAVITY 0.0098

Particle::Particle()
    : latitude(0)
    , longitude(0)
    , diameter(0)
    , density(0)
    , depth(0)
{}

Particle::Particle(double latitute, double longitude, double diameter, double density, double depth)
    : latitude(latitute)
    , longitude(longitude)
    , diameter(diameter)
    , density(density)
    , depth(depth)
{}

void Particle::terminal_buoyancy_velocity() {
    //    Buoyancy method from Zheng and Yapa (2000)
    //    Diameter assumed to be in km
    if (diameter <= 0.001) {
        double Nd = (4.0 * current_node->water_density * (density - current_node->water_density) * GRAVITY * pow(diameter, 3)) / (3 * current_node->)

    } else if (diameter <= 0.015) {

    } else {

    }
}