//
// Created by Georgia Stuart on 2/3/21.
//

#include <math.h>
#include <iostream>
#include <stdexcept>
#include "particle.h"

// Acceleration due to gravity in km/s^2
#define GRAVITY 0.0098
// Water viscosity from Braida in kg/ms
#define WATER_VISCOSITY 0.0009

Particle::Particle()
    : latitude(0)
    , longitude(0)
    , diameter(0)
    , density(0)
    , depth(0)
{}

Particle::Particle(double latitute, double longitude, double diameter, double density, double depth, double interfacial_tension)
    : latitude(latitute)
    , longitude(longitude)
    , diameter(diameter)
    , density(density)
    , depth(depth)
    , interfacial_tension(interfacial_tension)
{}

double Particle::terminal_buoyancy_velocity() {
    //  Buoyancy method from Zheng and Yapa (2000)
    //  Diameter assumed to be in km
    if (diameter <= 0.001) {
        //  Archimedes number * 4/3
        double Nd = calculate_nd();
        double reynolds = calculate_reynolds(Nd);

        return (reynolds * current_node->viscosity) / (current_node->density * diameter);

    } else {
        double critical_diameter = calculate_critical_diameter();

    }
}

double Particle::calculate_nd() {
    return (4.0 * current_node->density * (density - current_node->density)
            * GRAVITY * pow(diameter, 3)) / (3 * pow(current_node->viscosity,2));
}

double Particle::calculate_reynolds(double Nd) {
    if (Nd <= 73) {
        return Nd / 24.0 - 1.7569e-4 * pow(Nd, 2)
                   + 6.9252e-7 * pow(Nd, 3) - 2.3027e-10 * pow(Nd, 4);
    } else if (Nd <= 580) {
        double W = log10(Nd);
        double log_reynolds = -1.7095 + 1.33438 * W - 0.11591 * pow(W, 2);
        return pow(10.0, log_reynolds);
    } else if (Nd <= 1.55e7) {
        double W = log10(Nd);
        double log_reynolds = -1.81391 + 1.34671 * W - 0.12427 * pow(W, 2) + 0.006344 * pow(W, 3);
        return pow(10.0, log_reynolds);
    } else {
        throw std::runtime_error(std::string("Archimedes number is too large at ") + std::to_string(Nd) + std::string("."));
    }
}

double Particle::calculate_critical_diameter() {
    // Calculate x1, y1
    double H = 59.3;
    double J = 0.94 * pow(H, 0.757);
    double M = calculate_morton_number();
    double terminal_velocity = current_node->viscosity / (current_node->density * diameter)
            * pow(M, -0.149) * (J - 0.857);

    double x1 = log10(diameter);
    double y1 = log10(terminal_velocity);

    // Calculate x2, y2
    double E0 = calculate_eotvos_number();
    H = 4.0/3.0 * E0 * pow(M, -0.149) * pow(current_node->viscosity / WATER_VISCOSITY, -0.14);
    J = H <= 59.3 ? 0.94 * pow(H, 0.757) : 3.42 * pow(H, 0.441);
    terminal_velocity = current_node->viscosity / (current_node->density * 0.015)
                        * pow(M, -0.149) * (J - 0.857);

    double x2 = log10(0.015);
    double y2 = log10(terminal_velocity);

    // Calculate a1, b1, a2, b2
    double a1 = 0.5;
    double b1 = log10(0.711 * pow(GRAVITY * (density - current_node->density) / current_node->density, 0.5));
    double a2 = (y2 - y1) / (x2 - x1);
    double b2 = y1 - a2 * x1;

    return pow(10.0, (b2 - b1) / (a1 - a2));
}

// Morton + Eotvos Characterize shape of bubbles in a continuous phase
double Particle::calculate_morton_number() {
    return (GRAVITY * pow(current_node->viscosity, 4) * (density - current_node->density))
           / (pow(current_node->density, 2) * pow(interfacial_tension, 3));
}

double Particle::calculate_eotvos_number() {
    return GRAVITY * (density - current_node->density) * pow(diameter, 2) / interfacial_tension;
}