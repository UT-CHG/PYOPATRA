//
// Created by Georgia Stuart on 2/3/21.
//

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <iomanip>
#include "particle.h"
#include "mesh/mesh_base.h"

// Acceleration due to gravity in m/s^2
#define GRAVITY (-9.8)

Particle3D::Particle3D()
    : ParticleBase<3>()
    , diameter(0.0)
    , density(0.0)
    , interfacial_tension(0.0)
    , depth_index(0)
{ }

Particle3D::Particle3D(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension)
    : ParticleBase<3>()
    , diameter(diameter)
    , density(density)
    , interfacial_tension(interfacial_tension)
    , depth_index(0) // FIX THIS LATER!!!!
{
    location(0) = latitude;
    location(1) = longitude;
    location(2) = depth;
}

double Particle3D::terminal_buoyancy_velocity(double fluid_density, double fluid_viscosity, double water_viscosity) const {
    //  Buoyancy method from Zheng and Yapa (2000)
    //  Diameter assumed to be in meters

    if (diameter <= 0.001) {
        //  Archimedes number * 4/3
        double Nd = calculate_nd(fluid_density, fluid_viscosity);
        double reynolds = calculate_reynolds(Nd);

        return (reynolds * fluid_viscosity) / (fluid_density * diameter);

    } else {
        double critical_diameter = calculate_critical_diameter(fluid_density, fluid_viscosity, water_viscosity);

        if (diameter <= critical_diameter) {
            double J;
            double Eo = calculate_eotvos_number(diameter, fluid_density);
            double M = calculate_morton_number(fluid_density, fluid_viscosity);
            double H = 4.0 / 3.0 * Eo * pow(M, -0.149)
                    * pow(fluid_viscosity / water_viscosity, -0.14);

            if (2 < H && H <= 59.3) {
                J = 0.94 * pow(H, 0.757);
            } else if (H > 59.3) {
                J = 3.42 * pow(H, 0.441);
            } else {
                throw std::runtime_error(std::string("H >= 2 at ") + std::to_string(H) + std::string("."));
            }

            return fluid_viscosity / (fluid_density * diameter) * pow(M, -0.149) * (J - 0.857);
        } else {
            return 0.711 * sqrt(GRAVITY * diameter * (density - fluid_density) / fluid_density);
        }
    }
}

double Particle3D::calculate_nd(double fluid_density, double fluid_viscosity) const {
    return (4.0 * fluid_density * (density - fluid_density)
            * GRAVITY * pow(diameter, 3)) / (3 * pow(fluid_viscosity, 2));
}

/* static */ double Particle3D::calculate_reynolds(double Nd) {
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


double Particle3D::calculate_critical_diameter(double fluid_density, double fluid_viscosity, double water_viscosity) const {
    // Calculate x1, y1
    double H = 59.3;
    double J = 0.94 * pow(H, 0.757);
    double M = calculate_morton_number(fluid_density, fluid_viscosity);
    double temp_diameter = calculate_diameter_from_H(H, M, fluid_density, fluid_viscosity, water_viscosity);
    double terminal_velocity = fluid_viscosity / (fluid_density * temp_diameter)
                               * pow(M, -0.149) * (J - 0.857);

    double x1 = log10(temp_diameter);
    double y1 = log10(terminal_velocity);

    // Calculate x2, y2
    double E0 = calculate_eotvos_number(0.015, fluid_density);
    H = 4.0/3.0 * E0 * pow(M, -0.149) * pow(fluid_viscosity / water_viscosity, -0.14);
    J = H <= 59.3 ? 0.94 * pow(H, 0.757) : 3.42 * pow(H, 0.441);
    terminal_velocity = fluid_viscosity / (fluid_density * 0.015)
                        * pow(M, -0.149) * (J - 0.857);

    double x2 = log10(0.015);
    double y2 = log10(terminal_velocity);

    // Calculate a1, b1, a2, b2
    double a1 = 0.5;
    double b1 = log10(0.711 * pow(GRAVITY * (density - fluid_density) / fluid_density, 0.5));
    double a2 = (y2 - y1) / (x2 - x1);
    double b2 = y1 - a2 * x1;

    return pow(10.0, (b2 - b1) / (a1 - a2));
}

// Diameter from a give H, needed for computing critical diameter
double Particle3D::calculate_diameter_from_H(double H, double M, double fluid_density, double fluid_viscosity, double water_viscosity) const {
    double E0 = 3.0/4.0 * H * pow(M, 0.149) * pow(fluid_viscosity
                                                  / water_viscosity, 0.14);
    return sqrt(E0 * interfacial_tension / (GRAVITY * (density - fluid_density)));
}

// Morton + Eotvos Characterize shape of bubbles in a continuous phase
double Particle3D::calculate_morton_number(double fluid_density, double fluid_viscosity) const {
    return (GRAVITY * pow(fluid_viscosity, 4) * (density - fluid_density))
           / (pow(fluid_density, 2) * pow(interfacial_tension, 3));
}

double Particle3D::calculate_eotvos_number(double diameter_effective, double fluid_density) const {
    return GRAVITY * (density - fluid_density) * pow(diameter_effective, 2) / interfacial_tension;
}