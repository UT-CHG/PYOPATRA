//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYOPATRA_PARTICLE_H
#define PYOPATRA_PARTICLE_H

#include <Eigen/Dense>
#include "illnode.h"
#include "coordinate.h"

template <int dimension>
class ParticleBase {
public:
    using Vector = Eigen::Matrix<double, dimension, 1>;

protected:
    Vector location;
    size_t last_known_water_column_index;

    // Enable Intrusive Linked List Structure
    ILLNode<ParticleBase<dimension>> node;

public:
    ParticleBase() : location(Vector::Zero()), last_known_water_column_index(0), node(this) {}
    const Vector get_location() const { return location; }
    void set_location(const Vector& new_location) { location = new_location; }
    [[nodiscard]] ILLNode<ParticleBase<dimension>>& get_node() { return node; }
    void update_location(Vector& velocity, double time_delta) {
        location(1) += meters_to_longitude(velocity(1), location(0)) * 3600.0 * time_delta;
        location(0) += meters_to_latitude(velocity(0)) * 3600.0 * time_delta;
//        location += velocity * 6.0 / 185.0 * time_delta;
    }
    size_t get_last_known_water_column_index() { return last_known_water_column_index; }
    void set_water_column_index(size_t wc_index) { last_known_water_column_index = wc_index; }
    ParticleBase<dimension>* get_next() {
        if (node.next) {
            return node.next->owner;
        } else {
            return nullptr;
        }
    }
};

template <int dimension>
class Particle: public ParticleBase<dimension> {};

template<>
class Particle<3> :  public ParticleBase<3> {
private:
    double diameter, density, interfacial_tension;
    int depth_index;
public:
    Particle();
    Particle(double latitude, double longitude, double diameter, double density, double depth, double interfacial_tension);

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
using Particle3D = Particle<3>;

#endif //PYOPATRA_PARTICLE_H
