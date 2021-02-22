//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYOPATRA_PARTICLE_H
#define PYOPATRA_PARTICLE_H

#include "mesh.h"

class Particle {
public:
    double diameter, density, interfacial_tension;
    CoordinateD position;
    int depth_index;
    MeshNode *current_node;

    Particle();
    Particle(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension);
    double terminal_buoyancy_velocity();
    double calculate_nd() const;
    static double calculate_reynolds(double Nd);
    double calculate_critical_diameter() const;
    double calculate_morton_number() const;
    double calculate_eotvos_number() const;
};

#endif //PYOPATRA_PARTICLE_H
