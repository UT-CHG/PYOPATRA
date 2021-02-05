//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef LPTCPP_MESH_H
#define LPTCPP_MESH_H

#include "coordinate.h"

class MeshNode {
public:
    double water_density, temperature, viscosity;
    CoodinateD velocity, location;
};

class Mesh {};

#endif //LPTCPP_MESH_H
