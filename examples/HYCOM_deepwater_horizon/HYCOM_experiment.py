import os
from datetime import date, timedelta
from configparser import ConfigParser

from PYOPATRA import *


if __name__ == '__main__':
    # Configuration settings
    # In hours
    timestep = 1
    # Particle release per timedelta
    num_particles = 10
    # Particle release location
    particle_lon = -88.365997
    particle_lat = 28.736628
    # Time elapsed
    total_days = 7
    total_time_steps = int(24 / timestep * total_days)
    # When to add particles (time steps, not hours)
    add_particles_time_step_interval = 1
    # How frequently to save particles (time steps, not hours)
    particle_save_interval = 3

    # Read in config file
    config = ConfigParser()
    config.read('hycom_setup.ini')

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

    print(hycom_files)

    # Read HYCOM files
    hfp = HYCOMFileParser()
    hfp.read(hycom_files, diffusion_coefficient=10.0)

    # Set up 2D Triangular Mesh
    tm2d = TriangularMesh2D()
    tm2d.setup_mesh(hfp, 2)

    # Set Up ParticleList
    particles = ParticleList()

    # Time stepping
    for i in range(total_time_steps):
        # Inject more particles
        if i % add_particles_time_step_interval == 0:
            for j in range(num_particles):
                particles.append_particle(particle_lat, particle_lon)

        # Time stepping
        particles.time_step(time_delta=timestep)

        # Save particle snapshot (not to disk)
        if i % particle_save_interval == 0:
            particles.snapshot()

    # Save particle snapshots to disk
    particles.save_hdf5('{}/data/hycom_particles.hdf5'.format(file_prefix))



