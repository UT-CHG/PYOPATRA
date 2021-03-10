//
// Created by Georgia on 3/10/21.
//

#ifndef PYOPATRA_MESH_BASE_H
#define PYOPATRA_MESH_BASE_H

#include "mesh_base.h"
#include "mesh_element.h"
#include "../particle.h"

class MeshBase {
public:
    virtual ~MeshBase() = default;
    virtual PolygonMeshElement* find_particle_location(Particle &particle) = 0;
};

#endif //PYOPATRA_MESH_BASE_H
