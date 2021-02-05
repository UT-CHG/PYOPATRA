//
// Created by Georgia Stuart on 2/4/21.
//

#include "mesh.h"

MeshNode::MeshNode() 
    : density(0.0)
    , temperature(0.0)
    , viscosity(0.0)
{}

MeshNode::MeshNode(double density, double temperature, double viscosity)
    : density(density)
    , temperature(temperature)
    , viscosity(viscosity)
{}