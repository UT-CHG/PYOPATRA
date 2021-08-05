//
// Created by Georgia Stuart on 7/26/21.
//

#ifndef PYOPATRA_OBJECTIVE_FUNCTIONS_H
#define PYOPATRA_OBJECTIVE_FUNCTIONS_H

#include "Eigen/Dense"
#include "mpi.h"

#include "../particle_list.h"
#include "../util.h"

// At this time, only rectangular domains are implemented

template <int dimension>
class ObjectiveFunctionBase {
public:
    ObjectiveFunctionBase() {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        ptr_wrapper = this;
    };
    virtual ~ObjectiveFunctionBase() = default;

//    virtual void set_observed_values(const ParticleList<dimension>& particles) = 0;
    virtual double calculate_value(const ParticleList<dimension>& particles) = 0;
    PointerWrapper<ObjectiveFunctionBase<dimension>> get_pointer_wrapper() { return ptr_wrapper; }
    void set_observed_values(const Eigen::Ref<Eigen::MatrixXd> particle_locations) {
        auto temp_list = new ParticleList<dimension>();

        for (int i = 0; i < particle_locations.rows(); i++) {
            temp_list->create_particle(particle_locations.row(i));
        }

        finish_observed_setup_impl(*temp_list);

        temp_list->delete_all_particles();
        delete temp_list;
    }
    int get_rank() { return rank; }

protected:
    int rank, world_size;
    virtual void finish_observed_setup_impl(ParticleList<dimension>& particles) = 0;
    PointerWrapper<ObjectiveFunctionBase<dimension>> ptr_wrapper;
};

template <int dimension>
class BinObjectiveFunctionBase : public ObjectiveFunctionBase<dimension> {
protected:
    Eigen::VectorXd latitude_bounds, longitude_bounds;
    Eigen::MatrixXd bins, observed_bins;

public:
    BinObjectiveFunctionBase()
        : ObjectiveFunctionBase<dimension>()
        , latitude_bounds(0)
        , longitude_bounds(0)
        , bins(0, 0)
        , observed_bins(0, 0)
    {}

    // bounds are {min latitude, max latitude, min longitude, max longitude}
    BinObjectiveFunctionBase(int num_bins_lat, int num_bins_lon, const Eigen::Vector4d& bounds)
        : ObjectiveFunctionBase<dimension>()
        , latitude_bounds(num_bins_lat)
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

    double calculate_value(const ParticleList<dimension>& particles) {
        fill_bins(particles, bins);
        return calculate_value_impl(particles);
    }

    Eigen::MatrixXd& get_observed_bins() { return observed_bins; }
    Eigen::MatrixXd& get_bins() { return bins; }
    Eigen::VectorXd& get_latitude_bounds() { return latitude_bounds; }
    Eigen::VectorXd& get_longitude_bounds() { return longitude_bounds; }

protected:
    virtual double calculate_value_impl(const ParticleList<dimension>& particles) = 0;

    void finish_observed_setup_impl(ParticleList<dimension>& particles) {
        fill_bins(particles, observed_bins);
        observed_bins /= particles.get_length();
    }

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
    }
};


// https://github.com/PythonOT/POT/blob/master/ot/sliced.py
// https://stats.stackexchange.com/questions/404775/calculate-earth-movers-distance-for-two-grayscale-images
// http://robotics.stanford.edu/~rubner/emd/default.htm
// https://stackoverflow.com/questions/5101004/python-code-for-earth-movers-distance
// https://github.com/scipy/scipy/blob/v1.7.0/scipy/stats/stats.py#L8245-L8319
// https://stackoverflow.com/questions/45740677/earth-mover-distance-between-numpy-1-d-histograms
// https://home.ttic.edu/~ssameer/Research/Papers/WEMD_CVPR08.pdf
// https://en.wikipedia.org/wiki/Earth_mover%27s_distance#Computing_the_EMD
template <int dimension>
class SlicedWassersteinDistance : public BinObjectiveFunctionBase<dimension> {
private:
    using Parent = BinObjectiveFunctionBase<dimension>;
    Eigen::VectorXd proj, sample_proj, observed_proj, emd_vec;
    Eigen::MatrixXd recv_bin;
    int num_proj;
    std::default_random_engine sw_generator;

public:
    SlicedWassersteinDistance(int num_bins_lat, int num_bins_lon, const Eigen::Vector4d& bounds, int num_proj, unsigned int seed)
        : Parent(num_bins_lat, num_bins_lon, bounds)
        , proj(num_bins_lat)
        , sample_proj(num_bins_lon)
        , observed_proj(num_bins_lon)
        , emd_vec(num_bins_lon)
        , recv_bin(num_bins_lon, num_bins_lat)
        , num_proj(num_proj)
        , sw_generator(seed)
    {
        this->num_proj /= Parent::ObjectiveFunctionBase::world_size;
    }

    void set_num_proj(int new_num_proj) { num_proj = new_num_proj; }

private:
    double calculate_value_impl(const ParticleList<dimension>& particles) {
        double sum = 0.0;
        int proc_num_particles = particles.get_length();
        int num_particles;

        MPI_Allreduce(&proc_num_particles, &num_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Parent::get_bins().data(), recv_bin.data(), recv_bin.rows() * recv_bin.cols(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        recv_bin /= num_particles;

        for (int i = 0; i < num_proj; i++) {
            proj = Eigen::VectorXd::NullaryExpr(proj.size(), [&]() { return unif(sw_generator); });
            proj.normalize();

            sample_proj = recv_bin * proj;
            observed_proj = Parent::get_observed_bins() * proj;

            sum += wasserstein_distance_1d(sample_proj, observed_proj);
        }

        double all_sum;

        MPI_Reduce(&sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (Parent::ObjectiveFunctionBase::rank == 0) {
            std::cout << "Proj: " << num_proj << ", world size: " << Parent::ObjectiveFunctionBase::world_size << std::endl;
            return all_sum / (num_proj * Parent::ObjectiveFunctionBase::world_size);
        } else {
            return 0.0;
        }
    }

    // https://en.wikipedia.org/wiki/Earth_mover%27s_distance#Computing_the_EMD
    double wasserstein_distance_1d(Eigen::VectorXd& dist1, Eigen::VectorXd& dist2) {
        emd_vec.setZero();

        for (int i = 0; i < emd_vec.size() - 1; i ++) {
            emd_vec(i + 1) = dist1(i) + emd_vec(i) - dist2(i);
        }

        return emd_vec.cwiseAbs().sum();
    }
};

using SlicedWassersteinDistance2D = SlicedWassersteinDistance<2>;
using ObjectiveFunction2DPtrWrapper = PointerWrapper<ObjectiveFunctionBase<2>>;

#endif //PYOPATRA_OBJECTIVE_FUNCTIONS_H
