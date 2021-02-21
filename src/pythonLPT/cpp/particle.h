//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef LPTCPP_PARTICLE_H
#define LPTCPP_PARTICLE_H

#include "mesh.h"

class Particle {
public:
    double latitude, longitude, diameter, density, depth, interfacial_tension;
    MeshNode *current_node;

    Particle();
    Particle(double latitute, double longitude, double diameter, double density, double depth, double interfacial_tension);
    double terminal_buoyancy_velocity();
    double calculate_nd() const;
    static double calculate_reynolds(double Nd);
    double calculate_critical_diameter() const;
    double calculate_morton_number() const;
    double calculate_eotvos_number() const;
};

#endif //LPTCPP_PARTICLE_H
