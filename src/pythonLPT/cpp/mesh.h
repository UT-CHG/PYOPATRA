//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef LPTCPP_MESH_H
#define LPTCPP_MESH_H

#include "coordinate.h"

class MeshNode {
public:
    double density, temperature, viscosity;
    CoodinateD velocity, location;

    MeshNode();
    MeshNode(double water_density, double temperature, double viscosity);
};

class Mesh {};

#endif //LPTCPP_MESH_H
