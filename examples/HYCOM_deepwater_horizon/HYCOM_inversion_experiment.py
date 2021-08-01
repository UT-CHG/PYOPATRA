import os
from datetime import date, timedelta
from configparser import ConfigParser
import h5py
import numpy as np

from PYOPATRA import *


if __name__ == '__main__':
    # Configuration settings
    # In hours
    time_delta = 1
    # Particle release per timedelta and process
    num_particles = 5
    # True particle release location
    true_particle_lon = -88.365997
    true_particle_lat = 28.736628
    # Starting particle release location
    particle_lon = -87.8
    particle_lat = 29.2
    prev_loc = np.array((particle_lat, particle_lon))
    proposed_loc = np.array((particle_lat, particle_lon))
    # Time elapsed
    total_days = 8 * 7
    total_time_steps = int(24 / time_delta * total_days) - 4
    # When to add particles (time steps, not hours)
    add_particles_time_step_interval = 3
    # How frequently to save particles (time steps, not hours)
    particle_save_interval = 3
    # Number of Particles at the end
    total_particles = (total_time_steps // add_particles_time_step_interval + 1) * num_particles
    # Frame interval
    frame_interval = 3
    # Number of samples to draw for MCMC process
    num_samples = 10000
    if rank == 0:
        samples = np.zeros((num_samples, 2))
        obj_values = np.zeros(num_samples)
    # Other MCMC specifiers
    precision_parameter = 0.1
    step_length = np.array((0.03, 0.03))

    times = ['000', '003', '006', '009', '012', '015', '018', '021']
    start_date = date(2010, 4, 20)

    file_prefix = os.path.dirname(os.path.realpath(__file__))
    hycom_files = []

    for days_since_start in range(total_days):
        date = start_date + timedelta(days=days_since_start)

        for time_index, time_str in enumerate(times):
            file = '{}/data/hycom_gomu_501_{}{:02d}{:02d}00_t{}.nc'.format(file_prefix, date.year, date.month, date.day, time_str)
            hycom_files.append(file)

    if rank == 0:
        print('Reading HYCOM files....')
    # Read HYCOM files
    hfp = HYCOMFileParser()
    hfp.read(hycom_files, diffusion_coefficient=10.0)

    if rank == 0:
        print('Setting up mesh...')
    # Set up 2D Triangular Mesh
    tm2d = TriangularMesh2D()
    tm2d.setup_mesh(hfp, 2)

    # Set up particles
    particles = ParticleList()

    # Set up objective function
    obj = SlicedWassersteinDistance(700, 1000, [hfp.latitude[0], hfp.latitude[-1], hfp.longitude[0], hfp.longitude[-1]], 5000, 3000)

    # Set up solver
    solver = Solver(hfp.times, tm2d, particles, obj)

    obs_particles = None
    # Set up objective function
    if rank == 0:
        with h5py.File("{}/data/observed_particles.hdf5".format(file_prefix), "r") as fp:
            obs_particles_temp = fp['particles'][:, :]

        obs_particles = obs_particles_temp[~np.all(obs_particles_temp == 0, axis=1)]

    obs_particles = comm.bcast(obs_particles, root=0)
    obj.set_observed_values(obs_particles)

    previous_log_likelihood = 0
    previous_obj_value = 0
    accepted = 0

    value = 0
    check = 0

    for sample in range(num_samples):
        if sample > 0:
            if rank == 0:
                print("Sample number {}".format(sample))
                proposed_loc = prev_loc +  np.random.randn(2) * step_length
            proposed_loc = comm.bcast(proposed_loc, root=0)

        # print('Time stepping...')
        current_num_particles = 0
        frame = 0
        # Time stepping
        for i in range(total_time_steps):
            # print('Time step {}'.format(i))

            # Inject more particles
            if i % add_particles_time_step_interval == 0:
                for j in range(num_particles):
                    particles.append_particle(proposed_loc[0], proposed_loc[1])
                    current_num_particles += 1

            solver.time_step(time_delta)

        obj_value = solver.calculate_objective_value()
        # obj_value = comm.bcast(obj_value, root=0)
        log_likelihood = -obj_value / (2 * precision_parameter**2)
        if rank == 0:
            print(obj_value, previous_obj_value, log_likelihood, previous_log_likelihood)

        if sample == 0:
            previous_log_likelihood = log_likelihood
            previous_obj_value = obj_value
        else:
            if rank == 0:
                value = np.exp(log_likelihood - previous_log_likelihood)
                check = np.random.random()

                print(value, check, value > check)

                if value > check:
                    previous_log_likelihood = log_likelihood
                    previous_obj_value = obj_value
                    prev_loc[:] = proposed_loc[:]
                    accepted += 1

        if rank == 0:
            samples[sample, :] = prev_loc[:]
            obj_values[sample] = previous_obj_value
        solver.reset_solver()
        comm.barrier()

    if rank == 0:
        print("Acceptance Ratio: {}".format(accepted / num_samples))

        with h5py.File("{}/data/mcmc_save_data_{}_samples.hdf5".format(file_prefix, num_samples), "w") as fp:
            fp.create_dataset('samples', data=samples)
            fp.create_dataset('objectives', data=obj_values)

