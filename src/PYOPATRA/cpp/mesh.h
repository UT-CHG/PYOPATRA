//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_H
#define PYOPATRA_MESH_H

#include <vector>
#include "coordinate.h"

class MeshNode {
public:
    int num_depth_layers;
    std::vector<double> density, temperature, water_viscosity, viscosity;
    CoordinateD velocity, location;

    MeshNode();
    explicit MeshNode(int num_depths);
    explicit MeshNode(double density, double temperature);
    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};

class Mesh {};

#endif //PYOPATRA_MESH_H
