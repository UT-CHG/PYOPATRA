//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYOPATRA_PARTICLE_H
#define PYOPATRA_PARTICLE_H

#include "mesh/mesh_element.h"
#include "illnode.h"

class MeshBase;

class Particle {
public:
    double diameter, density, interfacial_tension;
    Coordinate3D location;
    int depth_index;
    MeshElementCursor *current_mesh_node;

    // Enable Intrusive Linked List Structure
    ILLNode<Particle> node;

    Particle();
    Particle(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension);
    double terminal_buoyancy_velocity(double fluid_density, double fluid_viscosity, double water_viscosity) const;
    double calculate_nd(double fluid_density, double fluid_viscosity) const;
    static double calculate_reynolds(double Nd);
    double calculate_critical_diameter(double fluid_density, double fluid_viscosity, double water_viscosity) const;
    double calculate_diameter_from_H(double H, double M, double fluid_density, double fluid_viscosity, double water_viscosity) const;
    double calculate_morton_number(double fluid_density, double fluid_viscosity) const;
    double calculate_eotvos_number(double diameter_effective, double fluid_density) const;
};


#endif //PYOPATRA_PARTICLE_H
