import numpy as np
from .pyopatra_pybind import CppTriangularMesh2D


class TriangularMesh2D:
    def __init__(self):
        pass

    def setup_elements_and_adjacency_list(self, vertex_list, regular=None, depths_array=None):
        elements = []
        adjacency_list = []

        if regular is not None:
            if len(regular) == 2:
                for i in range(regular[0]):
                    for j in range(regular[1]):
                        pass



