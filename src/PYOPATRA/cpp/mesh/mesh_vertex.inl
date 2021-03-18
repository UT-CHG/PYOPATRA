//
// Created by Georgia Stuart on 3/10/21.
//

#include <cmath>

template <int dimension>
MeshVertexBase<dimension>::MeshVertexBase()
    : location(Eigen::Matrix<double, dimension, 1>::Zero())
{}

template <int dimension>
MeshVertexBase<dimension>::MeshVertexBase(double latitude, double longitude)
{
    MeshVertexBase<dimension>();
    location(0) = latitude;
    location(1) = longitude;
}

template <int dimension>
bool operator==(const MeshVertexBase<dimension>& lhs, const MeshVertex<dimension>& rhs) {
    return (lhs.get_location() == rhs.get_location());
}
