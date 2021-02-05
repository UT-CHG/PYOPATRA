//
// Created by Georgia Stuart on 2/3/21.
//

#include <math.h>
#include <iostream>
#include <stdexcept>
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

double Particle::terminal_buoyancy_velocity() {
    //  Buoyancy method from Zheng and Yapa (2000)
    //  Diameter assumed to be in km
    if (diameter <= 0.001) {
        //  Archimedes number * 4/3
        double reynolds;
        double Nd = (4.0 * current_node->water_density * (density - current_node->water_density)
                * GRAVITY * pow(diameter, 3)) / (3 * pow(current_node->viscosity,2));

        if (Nd <= 73) {
            reynolds = Nd / 24.0 - 1.7569e-4 * pow(Nd, 2)
                    + 6.9252e-7 * pow(Nd, 3) - 2.3027e-10 * pow(Nd, 4);
        } else if (Nd <= 580) {
            double W = log10(Nd);
            reynolds = -1.7095 + 1.33438 * W - 0.11591 * pow(W, 2);
            reynolds = pow(10.0, reynolds);
        } else if (Nd <= 1.55e7) {
            double W = log10(Nd);
            reynolds = -1.81391 + 1.34671 * W - 0.12427 * pow(W, 2) + 0.006344 * pow(W, 3);
            reynolds = pow(10.0, reynolds);
        } else {
            throw std::runtime_error(std::string("Archimedes number is too large at ") + std::to_string(Nd) + std::string("."));
        }

        return (reynolds * current_node->viscosity) / (current_node->water_density * diameter);

    } else if (diameter <= 0.015) {
        double M = (GRAVITY * pow(current_node->viscosity, 4) * (density - current_node->water_density))
                / (pow(current_node->water_density, 2) * pow(interfacial_tension, 3));

    } else {

    }
}