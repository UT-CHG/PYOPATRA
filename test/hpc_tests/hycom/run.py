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
    # Particle release per timedelta
    num_particles = 5
    # Particle release location
    particle_lon = -88.365997
    particle_lat = 28.736628
    release_loc = np.array((particle_lat, particle_lon))
    # Time elapsed
    total_days = 8 * 7
    total_time_steps = int(24 / time_delta * total_days) - 4
    # When to add particles (time steps, not hours)
    add_particles_time_step_interval = 3
    # How frequently to save particles (time steps, not hours)
    particle_save_interval = 3
    # Number of Particles at the end
    total_particles = (total_time_steps // add_particles_time_step_interval + 1) * (num_particles * size)
    # Frame interval
    frame_interval = 3

    # # Read in config file
    # config = ConfigParser()
    # config.read('hycom_setup.ini')

    # Set up HYCOM file paths
    times = ['000', '003', '006', '009', '012', '015', '018', '021']
    start_date = date(2010, 4, 20)

    file_prefix = os.path.dirname(os.path.realpath(__file__))
    hycom_files = []

    for days_since_start in range(total_days):
        date = start_date + timedelta(days=days_since_start)

        for time_index, time_str in enumerate(times):
            file = '{}/data/hycom_gomu_501_{}{:02d}{:02d}00_t{}.nc'.format(file_prefix, date.year, date.month, date.day, time_str)
            hycom_files.append(file)

    print('Reading HYCOM files....')
    # Read HYCOM files
    hfp = HYCOMFileParser()
    hfp.read(hycom_files, diffusion_coefficient=10.0)

    print('Setting up mesh...')
    # Set up 2D Triangular Mesh
    tm2d = TriangularMesh2D()
    tm2d.setup_mesh(hfp, 2)

    # Set up particles
    particles = ParticleList()

    # Set up solver
    solver = Solver(hfp.times, tm2d, particles)

    # Set up snapshot particle saving file
    if rank == 0:
        with h5py.File('{}/data/snapshots.hdf5'.format(file_prefix), 'w') as fp:
            fp.create_dataset('snapshots', (total_particles, 2, total_time_steps // frame_interval + 1))

    print('Time stepping...')
    current_num_particles = 0
    frame = 0
    # Time stepping
    for i in range(total_time_steps):
        print('Time step {}'.format(i))

        # Inject more particles
        if i % add_particles_time_step_interval == 0:
            for j in range(num_particles):
                # print('Appending particle at ({}, {})'.format(particle_lon, particle_lat))
                particles.append_particle(release_loc[0], release_loc[1])
                current_num_particles += 1

        # Time stepping
        solver.time_step(time_delta)

        if i % frame_interval == 0:
            print('Saving particles')
            particle_locations = particles.get_all_particle_locations()
            if rank == 0:
                with h5py.File('{}/data/snapshots.hdf5'.format(file_prefix), 'a') as fp:
                    fp['snapshots'][:particle_locations.shape[0], :, frame] = particle_locations
            frame += 1



    print('Saving particle locations...')
    particle_locations = particles.get_all_particle_locations()
    if rank == 0:
        with h5py.File('{}/data/observed_particles.hdf5'.format(file_prefix), 'w') as fp:
            fp.create_dataset('particles', particle_locations.shape, data=particle_locations)

