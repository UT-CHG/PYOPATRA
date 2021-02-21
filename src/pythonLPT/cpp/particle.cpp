//
// Created by Georgia Stuart on 2/3/21.
//

#include <math.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include "particle.h"

// Acceleration due to gravity in m/s^2
#define GRAVITY (-9.8)

Particle::Particle()
    : diameter(0.0)
    , density(0.0)
    , interfacial_tension(0.0)
    , position()
    , current_node(nullptr)
{ }

Particle::Particle(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension)
    : diameter(diameter)
    , density(density)
    , interfacial_tension(interfacial_tension)
    , position(latitude, longitude, depth)
    , current_node(nullptr)
{ }

double Particle::terminal_buoyancy_velocity() {
    //  Buoyancy method from Zheng and Yapa (2000)
    //  Diameter assumed to be in meters
    if (diameter <= 0.001) {
        //  Archimedes number * 4/3
        double Nd = calculate_nd();
        double reynolds = calculate_reynolds(Nd);

        return (reynolds * current_node->viscosity) / (current_node->density * diameter);

    } else {
        double critical_diameter = calculate_critical_diameter();

        if (diameter <= critical_diameter) {
            double J;
            double Eo = calculate_eotvos_number();
            double M = calculate_morton_number();
            double H = 4.0 / 3.0 * Eo * pow(M, -0.149)
                    * pow(current_node->viscosity / calculate_water_viscosity(current_node->temperature), -0.14);

            if (2 < H && H <= 59.3) {
                J = 0.94 * pow(H, 0.757);
            } else if (H > 59.3) {
                J = 3.42 * pow(H, 0.441);
            } else {
                throw std::runtime_error(std::string("H >= 2 at ") + std::to_string(H) + std::string("."));
            }

            return current_node->viscosity / (current_node->density * diameter) * pow(M, -0.149) * (J - 0.857);
        } else {
            return 0.711 * sqrt(GRAVITY * diameter * (density - current_node->density) / current_node->density);
        }
    }
}

double Particle::calculate_nd() const {
    return (4.0 * current_node->density * (density - current_node->density)
            * GRAVITY * pow(diameter, 3)) / (3 * pow(current_node->viscosity,2));
}

/* static */ double Particle::calculate_reynolds(double Nd) {
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


double Particle::calculate_critical_diameter() const {
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
    H = 4.0/3.0 * E0 * pow(M, -0.149) * pow(current_node->viscosity / calculate_water_viscosity(current_node->temperature), -0.14);
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
double Particle::calculate_morton_number() const {
    return (GRAVITY * pow(current_node->viscosity, 4) * (density - current_node->density))
           / (pow(current_node->density, 2) * pow(interfacial_tension, 3));
}

double Particle::calculate_eotvos_number() const {
    return GRAVITY * (density - current_node->density) * pow(diameter, 2) / interfacial_tension;
}

// Water Viscosity
// From Huber et al 2009
// Temperature assumed to be Kelvin
// Density assumed to be kg m^-3

const double density_coefficients[6][7] = {
        {5.20094e-1,    2.22531e-1, -2.81378e-1,    1.61913e-1, -3.25372e-2,    0.0,        0.0},
        {8.50895e-2,    9.99115e-1, -9.06851e-1,    2.57399e-1, 0.0,            0.0,        0.0},
        {-1.08374,      1.88797,    -7.72479e-1,    0.0,        0.0,            0.0,        0.0},
        {-2.89555e-1,   1.26613,    -4.89827e-1,    0.0,        6.98452e-2,     0.0,        -4.35673e-3},
        {0.0,           0.0,        -2.57040e-1,    0.0,        0.0,            8.72102e-3, 0.0},
        {0.0,           1.20573e-1, 0.0,            0.0,        0.0,            0.0,        -5.93264e-4}
};

inline double dimensionless_temperature(double temperature) {
    return temperature / 647.096;
}

inline double dimensionless_density(double density) {
    return density / 322.0;
}

double dimensionless_viscosity_density_zero(double temperature) {
    double T = dimensionless_temperature(temperature);
    return (100.0 * sqrt(T)) / (1.67752 + 2.20462 / T + 0.6366564 / pow(T, 2) - 0.241605 / pow(T, 3));
}

/* static */ double Particle::calculate_water_viscosity(double temperature) {
    return
}