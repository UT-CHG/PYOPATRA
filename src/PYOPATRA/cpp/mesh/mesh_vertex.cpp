//
// Created by georgia on 3/10/21.
//

#include <cmath>

#include "mesh_vertex.h"


MeshVertex::MeshVertex()
    : location(0.0, 0.0, 0.0)
    , density(1, 1)
    , temperature(1, 1)
    , water_viscosity(1, 1)
    , viscosity(1, 1)
{}

MeshVertex::MeshVertex(double latitude, double longitude)
        : location(latitude, longitude, 0.0)
        , density(1, 1)
        , temperature(1, 1)
        , water_viscosity(1, 1)
        , viscosity(1, 1)
{}

MeshVertex::MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_depth_layers)
    : location(latitude, longitude, bathymetric_depth)
    , density(num_depth_layers, 1)
    , temperature(num_depth_layers, 1)
    , water_viscosity(num_depth_layers, 1)
    , viscosity(num_depth_layers, 1)
{}

MeshVertex::MeshVertex(double latitude, double longitude, double bathymetric_depth, int num_depth_layers, int num_time_steps)
        : location(latitude, longitude, bathymetric_depth)
        , density(num_depth_layers, num_time_steps)
        , temperature(num_depth_layers, num_time_steps)
        , water_viscosity(num_depth_layers, num_time_steps)
        , viscosity(num_depth_layers, num_time_steps)
{}

MeshVertex::MeshVertex(double latitude, double longitude, double bathymetric_depth, Eigen::Ref<const Eigen::ArrayXXd>& density, Eigen::Ref<const Eigen::ArrayXXd>& temperature)
    : location(latitude, longitude, bathymetric_depth)
    , density(density)
    , temperature(temperature)
{}

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

/* static */ double MeshVertex::calculate_fluid_viscosity(double temperature, double density) {
    return 1.0e-6 * dimensionless_viscosity_density_zero(temperature) * dimensionless_viscosity_due_to_density(temperature, density);
}

// This equation is only recommended for the temperature range of 253.15 K to 383.15 K (-20 C to 110 C).
/* static */ double MeshVertex::calculate_pure_water_viscosity(double temperature) {
    double T = temperature / 300.0;
    return 1.0e-6 * (280.68 * pow(T, -1.9) + 511.45 * pow(T, -7.7) + 61.131 * pow(T, -19.6) + 0.45903 * pow(T, -40.0));
}