//
// Created by Georgia Stuart on 3/10/21.
//

#include "triangular_mesh.h"

TriangularMesh::TriangularMesh()
    : adjacency_list(0)
{}

TriangularMesh::TriangularMesh(int num_elements)
    : adjacency_list(num_elements)
{}