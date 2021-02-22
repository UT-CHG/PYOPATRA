//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef LPTCPP_MESH_H
#define LPTCPP_MESH_H

#include <vector>
#include "coordinate.h"

class MeshNode {
public:
    int num_depth_layers;
    std::vector<double> density, temperature, salinity, water_viscosity, viscosity;
    CoordinateD velocity, location;

    MeshNode();
    explicit MeshNode(int num_depths);
    explicit MeshNode(double density, double temperature, double salinity);
    static double calculate_water_viscosity(double temperature, double density);
    static double calculate_salt_water_viscosity(double temperature, double density, double salinity);
};

class Mesh {};

#endif //LPTCPP_MESH_H
