from configparser import ConfigParser
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

from .pyopatra_pybind import CppParticleList2D

class ParticleList:
    def __init__(self, particle_save_interval=None, dimensions=2):
        self.dimensions = dimensions
        self._cpp_particle_list = None

        if dimensions == 2:
            self._cpp_particle_list = CppParticleList2D()
        else:
            raise NotImplementedError('Only 2D currently implemented')

    def append_particle(self, latitude, longitude, depth=None):
        if self.dimensions == 2:
            self._cpp_particle_list.create_particle(np.array([latitude, longitude]))
        else:
            raise NotImplementedError('Only 2D currently implemented')

    def get_num_particles(self):
        return self._cpp_particle_list.get_length()

    def get_all_particle_locations(self):
        locations = comm.gather(self._cpp_particle_list.get_all_particle_locations(), root=0)

        if rank == 0:
            total_num_particles = 0
            for location in locations:
                total_num_particles += location.shape[0]

            snapshot = np.zeros((total_num_particles, 2))

            position = 0
            for location in locations:
                snapshot[position:position + location.shape[0], :] = location
                position += location.shape[0]

            return snapshot
        else:
            return None

    def get_all_particle_column_indices(self):
        locations = comm.gather(self._cpp_particle_list.get_all_particle_column_indices(), root=0)

        if rank == 0:
            total_num_particles = 0
            for location in locations:
                total_num_particles += location.shape[0]

            snapshot = np.zeros(total_num_particles, dtype=int)

            position = 0
            for location in locations:
                snapshot[position:position + location.shape[0]] = location
                position += location.shape[0]

            return snapshot
        else:
            return None

    def get_num_particles(self):
        return comm.reduce(self._cpp_particle_list.get_length(), root=0, op=MPI.SUM)

    def reset_particles(self):
        self._cpp_particle_list.reset_particles()

    def snapshot(self):
        pass

    def save_hdf5(self, filename):
        pass