from configparser import ConfigParser
import numpy as np

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

    def get_all_particle_locations(self):
        return self._cpp_particle_list.get_all_particle_locations()

    def get_all_particle_column_indices(self):
        return self._cpp_particle_list.get_all_particle_column_indices()

    def get_num_particles(self):
        return self._cpp_particle_list.get_length()

    def reset_particles(self):
        self._cpp_particle_list.reset_particles()

    def snapshot(self):
        pass

    def save_hdf5(self, filename):
        pass