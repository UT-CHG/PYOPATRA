//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYOPATRA_PARTICLE_H
#define PYOPATRA_PARTICLE_H

#include "mesh/mesh_element.h"
#include "illnode.h"

class Particle {
public:
    double diameter, density, interfacial_tension;
    CoordinateD position;
    int depth_index;
    MeshElement *current_mesh_node;

    // Enable Intrusive Linked List Structure
    ILLNode<Particle> node;

    Particle();
    Particle(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension);
    double terminal_buoyancy_velocity() const;
    double calculate_nd() const;
    static double calculate_reynolds(double Nd);
    double calculate_critical_diameter() const;
    double calculate_diameter_from_H(double H, double M) const;
    double calculate_morton_number() const;
    double calculate_eotvos_number(double diameter_effective) const;
};


#endif //PYOPATRA_PARTICLE_H
