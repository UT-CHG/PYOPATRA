//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYOPATRA_PARTICLE_H
#define PYOPATRA_PARTICLE_H

#include "mesh/mesh_element.h"
#include "illnode.h"

template <int dimension>
class ParticleBase {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

protected:
    Vector location;
    MeshElementCursor<dimension> *current_mesh_node;

    // Enable Intrusive Linked List Structure
    ILLNode<ParticleBase<dimension>> node;

public:
    ParticleBase() : location(Vector::Zero()), current_mesh_node(nullptr), node(this) {}
    const Vector get_location() const { return location; }
    [[nodiscard]] ILLNode<ParticleBase<dimension>>& get_node() { return node; }
};

template <int dimension>
class Particle: public ParticleBase<dimension> {};

class Particle3D: public ParticleBase<3> {
private:
    double diameter, density, interfacial_tension;
    int depth_index;
public:
    Particle3D();
    Particle3D(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension);

    double get_diameter() const { return diameter; }
    double get_density() const { return density; }
    double get_interfacial_tension() const { return interfacial_tension; }
    int get_depth_index() const { return depth_index; }

    double terminal_buoyancy_velocity(double fluid_density, double fluid_viscosity, double water_viscosity) const;
    double calculate_nd(double fluid_density, double fluid_viscosity) const;
    static double calculate_reynolds(double Nd);
    double calculate_critical_diameter(double fluid_density, double fluid_viscosity, double water_viscosity) const;
    double calculate_diameter_from_H(double H, double M, double fluid_density, double fluid_viscosity, double water_viscosity) const;
    double calculate_morton_number(double fluid_density, double fluid_viscosity) const;
    double calculate_eotvos_number(double diameter_effective, double fluid_density) const;
};

using Particle2D = Particle<2>;

#endif //PYOPATRA_PARTICLE_H
