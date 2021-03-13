//
// Created by Georgia Stuart on 3/10/21.
//

#ifndef PYOPATRA_TRIANGULAR_MESH_H
#define PYOPATRA_TRIANGULAR_MESH_H

#include <vector>
#include <array>
#include <Eigen/Dense>

#include "mesh_base.h"



class TriangularMesh: public MeshBase {
//protected:
//    using AdjacencyElement = std::array<TriangularMeshElement*, 3>;
//    std::vector<AdjacencyElement> adjacency_list;
//
//public:
//    TriangularMesh();
//    explicit TriangularMesh(int num_elements);
//
//    TriangularMeshElement* find_particle_location(Particle &particle) override;
//    const std::vector<AdjacencyElement>& GetList() const { return adjacency_list; }
//    void setup_mesh(Eigen::Ref<Eigen::ArrayXXd> node_coordinates, Eigen::Ref<Eigen::ArrayXXi> element_nodes);
};

#endif //PYOPATRA_TRIANGULAR_MESH_H
