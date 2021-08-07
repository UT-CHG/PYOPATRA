//
// Created by Georgia Stuart on 7/28/21.
//

#ifndef PYOPATRA_SOLVER_H
#define PYOPATRA_SOLVER_H

#include "Eigen/Dense"
#include "mpi.h"

#include "mesh/mesh.h"
#include "particle_list.h"
#include "inversion_tools/objective_functions.h"
#include "util.h"

template <int num_vertices_per_element, int dimension>
class Solver {
public:
    using SolverMesh = Mesh<num_vertices_per_element, dimension>;
    using SolverParticles = ParticleList<dimension>;
    using SolverObjFn = ObjectiveFunctionBase<dimension>;
    using Vector = Eigen::Matrix<double, dimension, 1>;

    Solver()
        : mesh(nullptr)
        , particles(nullptr)
        , obj_fn(nullptr)
        , current_time_step(0)
        , starting_time(measured_times(0))
        , time(starting_time)
        , measured_times(Eigen::VectorXd::Zero())
    {}

    Solver(PointerWrapper<SolverMesh>& mesh_ptr, PointerWrapper<SolverParticles>& particle_list_ptr, PointerWrapper<SolverObjFn>& obj_fn_ptr, Eigen::VectorXd measured_times)
        : mesh(mesh_ptr.get_pointer())
        , particles(particle_list_ptr.get_pointer())
        , obj_fn(obj_fn_ptr.get_pointer())
        , current_time_step(0)
        , starting_time(measured_times(0))
        , time(starting_time)
        , measured_times(measured_times)
    {
        update_particle_location_indices();
    }

    Solver(PointerWrapper<SolverMesh>& mesh_ptr, PointerWrapper<SolverParticles>& particle_list_ptr, Eigen::VectorXd measured_times)
            : mesh(mesh_ptr.get_pointer())
            , particles(particle_list_ptr.get_pointer())
            , obj_fn(nullptr)
            , current_time_step(0)
            , starting_time(measured_times(0))
            , time(starting_time)
            , measured_times(measured_times)
    {
        update_particle_location_indices();
    }

    double get_current_time() { return time; }

    void time_step(double time_delta) {
        update_particle_locations(time_delta);
    }

    void reset_solver() {
        time = starting_time;
        current_time_step = 0;
        particles->delete_all_particles();
    }

    void update_particle_locations(double time_delta) {
        auto current = particles->get_head();
        auto temp = current;
        size_t lower_bound = current_time_step;
        Vector velocity_update;

        time += time_delta;

        while (time >= measured_times[lower_bound + 1]) {
            lower_bound++;
        }

        current_time_step = lower_bound;

        while (current) {
            temp = current->get_next();
            mesh->get_water_columns()[current->get_last_known_water_column_index()]
                    .interpolate_velocity(mesh->get_elements(), mesh->get_vertices(), mesh->get_velocities_ptr(), mesh->get_diffusions_ptr(), current->get_location(), lower_bound, time_delta, time,
                                          measured_times[current_time_step], measured_times[current_time_step + 1], velocity_update);
            current->update_location(velocity_update, time_delta);
            update_particle_mesh_location(*current);
            current = temp;
        }
    }

    void update_particle_location_indices() {
        typename ParticleList<dimension>::ParticleN* current = particles->get_head();

        while (current) {
            update_particle_mesh_location(*current);
            current = current->get_next();
        }
    }

    void update_particle_mesh_location(typename ParticleList<dimension>::ParticleN& particle) {
        int new_col = mesh->locate_new_water_column(&mesh->get_water_columns()[particle.get_last_known_water_column_index()], particle.get_location());
        if (new_col >= 0) {
            particle.set_water_column_index(new_col);
        } else {
            // Particle has gone off the mesh, remove it
            particle.get_node().remove();
            particles->decrement_length();
            delete &particle;
            std::cout << "Removing a particle" << std::endl;
        }
    }

    double calculate_objective_value() {
        return obj_fn->calculate_value(*particles);
    }
private:
    SolverMesh *mesh;
    SolverParticles *particles;
    SolverObjFn *obj_fn;

    size_t current_time_step;
    double starting_time;
    double time;
    Eigen::VectorXd measured_times;
};

using TriangularMesh2DSolver = Solver<3, 2>;

#endif //PYOPATRA_SOLVER_H
