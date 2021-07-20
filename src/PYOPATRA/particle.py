from configparser import ConfigParser
from .pyopatra_pybind import CppParticleList2D

class ParticleList:
    def __init__(self, particle_save_interval=None, dimensions=2):
        self.dimensions = dimensions
        self.cpp_particle_list = None

        if dimensions == 2:
            self.cpp_particle_list = CppParticleList2D()
        else:
            raise NotImplementedError('Only 2D currently implemented')

    # def time_step(self, time_delta=None):
    #     pass
    #
    # def append_particle(self, latitude, longitude, depth=None):
    #     if self.dimensions == 2:
    #         self.cpp_particle_list.create_particle([latitude, longitude])
    #     else:
    #         raise NotImplementedError('Only 2D currently implemented')
    #
    # def snapshot(self):
    #     pass
    #
    # def save_hdf5(self, filename):
    #     pass