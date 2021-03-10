//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef PYOPATRA_MESH_NODE_H
#define PYOPATRA_MESH_NODE_H

#include <vector>
#include <array>
#include "coordinate.h"

class MeshNode {
public:
    int num_depth_layers, mesh_index;
    std::vector<double> density, temperature, water_viscosity, viscosity;
    CoordinateD velocity, location;

    MeshNode();
    explicit MeshNode(int mesh_index);
    explicit MeshNode(int mesh_index, int num_depths);
    explicit MeshNode(int mesh_index, double density, double temperature);
    static double calculate_fluid_viscosity(double temperature, double density);
    static double calculate_pure_water_viscosity(double temperature);
};

class TriangularMeshNode: public MeshNode {
public:
    std::array<Location, 3> vertices;
};

#endif //PYOPATRA_MESH_NODE_H
