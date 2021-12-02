//
// Created by Georgia Stuart on 3/10/21.
//

#include <cmath>

template <int dimension>
MeshVertexBase<dimension>::MeshVertexBase(int index, int num_timesteps, int num_wind_timesteps)
    : location(MeshVertexBase<dimension>::Vector::Zero())
    , velocity(index * num_timesteps)
    , diffusion_coefficient(index * num_timesteps)
    , wind(index * (num_wind_timesteps == 0 ? 2 : num_wind_timesteps))
{}

template <int dimension>
MeshVertexBase<dimension>::MeshVertexBase(double latitude, double longitude, int index, int num_timesteps, int num_wind_timesteps)
    : MeshVertexBase<dimension>(index, num_timesteps, num_wind_timesteps)
{
    location(0) = latitude;
    location(1) = longitude;
}

//template <int dimension>
//MeshVertexBase<dimension>::MeshVertexBase(double latitude, double longitude, Vector velocity)
//    : location(Eigen::Matrix<double, dimension, 1>::Zero())
//    , velocity(velocity)
//    , diffusion_coefficient(Eigen::Matrix<double, dimension, 1>::Zero())
//{
//    location(0) = latitude;
//    location(1) = longitude;
//}
//
//template <int dimension>
//MeshVertexBase<dimension>::MeshVertexBase(double latitude, double longitude, Vector velocity, Vector diffusion_coefficient)
//        : location(Eigen::Matrix<double, dimension, 1>::Zero())
//        , velocity(velocity)
//        , diffusion_coefficient(diffusion_coefficient)
//{
//    location(0) = latitude;
//    location(1) = longitude;
//}

template <int dimension>
bool operator==(const MeshVertexBase<dimension>& lhs, const MeshVertex<dimension>& rhs) {
    return (lhs.get_location() == rhs.get_location());
}
