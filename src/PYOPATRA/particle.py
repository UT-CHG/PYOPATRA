from configparser import ConfigParser

class ParticleList:
    def __init__(self, particle_save_interval=None):
        pass

    def time_step(self, time_delta=None):
        pass

    def append_particle(self, latitude, longitude):
        pass

    def snapshot(self):
        pass

    def save_hdf5(self, filename):
        pass