import numpy as np
from .pyopatra_pybind import MeshVertex2D

class Vertex:
    def __init__(self, latitude, longitude, velocity, dimensions):
        if dimensions == 2:
            self.vertex = MeshVertex2D(latitude, longitude, velocity)
        else:
            raise NotImplementedError('Dimensions other than 2 have not been implemented (yet)')

    def get_latitude(self):
        return self.vertex.get_latitude()

    def get_longitude(self):
        return self.vertex.get_longitude()

    def get_velocity(self):
        return self.vertex.get_velocity()
