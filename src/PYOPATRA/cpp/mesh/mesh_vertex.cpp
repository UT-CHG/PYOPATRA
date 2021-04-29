//
// Created by Georgia Stuart on 3/18/21.
//

#include "mesh_vertex.h"

#include <utility>

MeshVertex3D::MeshVertex(int num_timesteps)
        : MeshVertexBase<3>(num_timesteps)
        , density(num_timesteps, 0.0)
        , temperature(num_timesteps, 0.0)
        , water_viscosity(num_timesteps, 0.0)
        , viscosity(num_timesteps, 0.0)
{}


MeshVertex3D::MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_timesteps)
        : MeshVertex3D::MeshVertex(num_timesteps)
{
    location(0) = latitude;
    location(1) = longitude;
    location(2) = bathymetric_depth;
}

//MeshVertex3D::MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature, int num_timesteps)
//        : MeshVertexBase<3>(latitude, longitude, num_timesteps)
//        , density(density)
//        , temperature(temperature)
//        , water_viscosity(calculate_pure_water_viscosity(temperature))
//        , viscosity(calculate_fluid_viscosity(temperature, density))
//{
//    location(2) = bathymetric_depth;
//}

//MeshVertex3D::MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature, Vector velocity)
//        : MeshVertexBase<3>(latitude, longitude, std::move(velocity))
//        , density(density)
//        , temperature(temperature)
//        , water_viscosity(calculate_pure_water_viscosity(temperature))
//        , viscosity(calculate_fluid_viscosity(temperature, density))
//{
//    location(2) = bathymetric_depth;
//}
//
//MeshVertex3D::MeshVertex(double latitude, double longitude, double bathymetric_depth, double density, double temperature, Vector velocity, Vector diffusion_coefficient)
//        : MeshVertexBase<3>(latitude, longitude, std::move(velocity), std::move(diffusion_coefficient))
//        , density(density)
//        , temperature(temperature)
//        , water_viscosity(calculate_pure_water_viscosity(temperature))
//        , viscosity(calculate_fluid_viscosity(temperature, density))
//{
//    location(2) = bathymetric_depth;
//}

//void MeshVertex3D::set_temperature(double new_temperature) {
//    temperature = new_temperature;
//}
//
//void MeshVertex3D::set_density(double new_density) {
//    density = new_density;
//}

// Water Viscosity
// From Huber et al 2009
// Temperature assumed to be Kelvin
// Density assumed to be kg m^-3

const double density_coefficients[6][7] = {
        {5.20094e-1,    2.22531e-1, -2.81378e-1,    1.61913e-1, -3.25372e-2,    0.0,        0.0},
        {8.50895e-2,    9.99115e-1, -9.06851e-1,    2.57399e-1, 0.0,            0.0,        0.0},
        {-1.08374,      1.88797,    -7.72479e-1,    0.0,        0.0,            0.0,        0.0},
        {-2.89555e-1,   1.26613,    -4.89837e-1,    0.0,        6.98452e-2,     0.0,        -4.35673e-3},
        {0.0,           0.0,        -2.57040e-1,    0.0,        0.0,            8.72102e-3, 0.0},
        {0.0,           1.20573e-1, 0.0,            0.0,        0.0,            0.0,        -5.93264e-4}
};

inline double dimensionless_temperature(double temperature) {
    return temperature / 647.096;
}

inline double dimensionless_density(double density) {
    return density / 322.0;
}

double dimensionless_viscosity_density_zero(double temperature) {
    double T = dimensionless_temperature(temperature);
    return (100.0 * sqrt(T)) / (1.67752 + 2.20462 / T + 0.6366564 / pow(T, 2) - 0.241605 / pow(T, 3));
}

double dimensionless_viscosity_due_to_density(double temperature, double density) {
    double rho = dimensionless_density(density);
    double T = dimensionless_temperature(temperature);
    double viscosity = 0.0;

    for (int i = 0; i <= 5; i++) {
        double temp = 0.0;
        for (int j = 0; j <= 6; j++) {
            temp += density_coefficients[i][j] * pow(rho - 1.0, j);
        }
        viscosity += pow(1.0 / T - 1.0, i) * temp;
    }
    return exp(rho * viscosity);
}

/* static */ double MeshVertex3D::calculate_fluid_viscosity(double temperature, double density) {
    return 1.0e-6 * dimensionless_viscosity_density_zero(temperature) * dimensionless_viscosity_due_to_density(temperature, density);
}

// This equation is only recommended for the temperature range of 253.15 K to 383.15 K (-20 C to 110 C).
/* static */ double MeshVertex3D::calculate_pure_water_viscosity(double temperature) {
    double T = temperature / 300.0;
    return 1.0e-6 * (280.68 * pow(T, -1.9) + 511.45 * pow(T, -7.7) + 61.131 * pow(T, -19.6) + 0.45903 * pow(T, -40.0));
}