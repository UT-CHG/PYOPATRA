//
// Created by Georgia Stuart on 7/26/21.
//

#ifndef PYOPATRA_OBJECTIVE_FUNCTIONS_H
#define PYOPATRA_OBJECTIVE_FUNCTIONS_H

#include "Eigen/Dense"
#include "../particle_list.h"

// At this time, only rectangular domains are implemented

template <int dimension>
class ObjectiveFunctionBase {
public:
    ObjectiveFunctionBase() = default;
    virtual ~ObjectiveFunctionBase() = default;

    virtual void set_observed_values(const ParticleList<dimension>& particles) = 0;
    virtual double calculate_value(const ParticleList<dimension>& particles) = 0;
};

template <int dimension>
class BinObjectiveFunctionBase : public ObjectiveFunctionBase<dimension> {
private:
    Eigen::VectorXd latitude_bounds, longitude_bounds;
    Eigen::MatrixXd bins, observed_bins;

public:
    BinObjectiveFunctionBase()
        : latitude_bounds(0)
        , longitude_bounds(0)
        , bins(0, 0)
        , observed_bins(0, 0) {}

    // bounds are {min latitude, max latitude, min longitude, max longitude}
    BinObjectiveFunctionBase(int num_bins_lat, int num_bins_lon, const Eigen::Vector4d& bounds)
        : latitude_bounds(num_bins_lat)
        , longitude_bounds(num_bins_lon)
        , bins(num_bins_lon, num_bins_lat)
        , observed_bins(num_bins_lon, num_bins_lat)
    {
        double lat_step = (bounds(1) - bounds(0)) / num_bins_lat;
        double lon_step = (bounds(3) - bounds(2)) / num_bins_lon;

        for (int i = 0; i < num_bins_lat; i++) {
            latitude_bounds[i] = bounds(0) + lat_step * i;
        }

        for (int i = 0; i < num_bins_lon; i++) {
            longitude_bounds[i] = bounds(2) + lon_step * i;
        }

        observed_bins.setZero();
        bins.setZero();
    }

    virtual double calculate_value(const ParticleList<dimension>& particles) {
        fill_bins(particles, bins);

        return calculate_value_impl();
    }

    void set_observed_values(const ParticleList<dimension>& particles) {
        fill_bins(particles, observed_bins);
    }

    Eigen::MatrixXd& get_observed_bins() { return observed_bins; }
    Eigen::MatrixXd& get_bins() { return bins; }
    Eigen::VectorXd& get_latitude_bounds() { return latitude_bounds; }
    Eigen::VectorXd& get_longitude_bounds() { return longitude_bounds; }

protected:
    virtual double calculate_value_impl() = 0;

    int get_1D_bin_coord(double position, const Eigen::VectorXd& bounds) {
        size_t len = bounds.size();
        size_t upper = len;
        size_t lower = 0;

        while (upper - lower > 1) {
            size_t half = (upper + lower) / 2;

            if (position >= bounds(half)) {
                lower = half;
            } else {
                upper = half;
            }
        }

        return lower;
    }

    void fill_bins(const ParticleList<dimension>& particles, Eigen::MatrixXd& bins_array) {
        bins_array.setZero();
        int lat_bin, lon_bin;
        auto current_particle = particles.get_head();

        while (current_particle) {
            lat_bin = get_1D_bin_coord(current_particle->get_location()(0), latitude_bounds);
            lon_bin = get_1D_bin_coord(current_particle->get_location()(1), longitude_bounds);

            bins_array(lon_bin, lat_bin) += 1;

            current_particle = current_particle->get_next();
        }

        bins_array /= particles.get_length();
    }
};

template <int dimension>
class SlicedWassersteinDistance : public BinObjectiveFunctionBase<dimension> {
private:
    using Parent = BinObjectiveFunctionBase<dimension>;

public:
    SlicedWassersteinDistance(int num_bins_lat, int num_bins_lon, const Eigen::Vector4d& bounds)
        : Parent(num_bins_lat, num_bins_lon, bounds)
    {}

private:
    double calculate_value_impl() {
        return 0.0;
    }
};

using SlicedWassersteinDistance2D = SlicedWassersteinDistance<2>;

#endif //PYOPATRA_OBJECTIVE_FUNCTIONS_H
