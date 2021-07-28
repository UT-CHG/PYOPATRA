import os
import numpy as np

from PYOPATRA import *


if __name__ == '__main__':
    file_prefix = os.path.dirname(os.path.realpath(__file__)) + '/hycom_data'
    hycom_files = [
        file_prefix + '/hycom_gomu_501_2010042000_t000.nc',
        file_prefix + '/hycom_gomu_501_2010042000_t003.nc',
        file_prefix + '/hycom_gomu_501_2010042000_t006.nc',
        file_prefix + '/hycom_gomu_501_2010042000_t009.nc',
        file_prefix + '/hycom_gomu_501_2010042000_t012.nc',
    ]

    hfp = HYCOMFileParser()
    hfp.read(hycom_files, diffusion_coefficient=10.0)

    tm2d = TriangularMesh2D()
    tm2d.setup_mesh(hfp, 2)

    particles = ParticleList()
    solver = Solver(tm2d, particles)

    particle_locations = np.array(
        [[10.5, -80.1],
         [3.14, -90.16],
         [5.5, -89.0]],
        order='F'
    )

    for row_index in range(particle_locations.shape[0]):
        particles.append_particle(particle_locations[row_index, 0], particle_locations[row_index, 1])