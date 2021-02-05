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
    Particle(double latitute, double longitude, double diameter, double density, double depth);
    double terminal_buoyancy_velocity();
};

#endif //LPTCPP_PARTICLE_H
