//
// Created by Georgia Stuart on 2/3/21.
//

#include "particle.h"


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