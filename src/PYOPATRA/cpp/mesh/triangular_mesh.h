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
protected:
    using AdjacencyElement = std::array<MeshElement*, 3>;
    std::vector<AdjacencyElement> adjacency_list;
    std::vector<>
public:
    TriangularMesh();
    explicit TriangularMesh(int num_elements);

    MeshElement* find_particle_location(Particle &particle) override;
    const std::vector<AdjacencyElement>& GetList() const;
    void setup_mesh(Eigen::Ref<Eigen::ArrayXXd> node_coordinates, Eigen::Ref<Eigen::ArrayXXi> element_nodes);
};

#endif //PYOPATRA_TRIANGULAR_MESH_H
