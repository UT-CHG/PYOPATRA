import numpy as np

from PYOPATRA import FileParserBase, MeshVertex2D
from .pyopatra_pybind import CppTriangularMesh2D, TriangularMeshElement2D


class TriangularMesh2D:
    def __init__(self):
        self.vertex_list = None
        self._cpp_vertex_list = None

        self.adjacency_list = None
        self.element_list = None

        self.times = None
        self.regular_dimensions = None

    def setup_vertices(self, file_parser: FileParserBase):
        self.vertex_list = []
        self._cpp_vertex_list = []
        index = 0

        self.times = file_parser.times

        if file_parser.regular_dimensions is not None:
            self.regular_dimensions = file_parser.regular_dimensions
            for i in range(file_parser.regular_dimensions[0]):
                for j in range(file_parser.regular_dimensions[1]):
                    for time in range(len(self.times)):
                        self.vertex_list.append(MeshVertex2D(
                            file_parser.latitude[i],
                            file_parser.longitude[j],
                            file_parser.velocity[:, i * j, time],
                            file_parser.diffusion_coefficient[:, i * j, time],
                            len(file_parser.regular_dimensions)
                        ))

                        self._cpp_vertex_list.append(self.vertex_list[-1].vertex)
                        index += 1

    def setup_elements_and_adjacency_list(self, regular=None, depths_array=None):
        elements = []
        adjacency_list = []

        if self.vertex_list is None:
            raise AttributeError('vertex_list must be set prior to setting up elements and adjacency.')

        if regular is not None:
            index = 0

            if len(regular) == 2:
                for i in range(regular[0] - 1):
                    for j in range(regular[1] - 1):
                        for k in range(2):
                            for time in range(len(self.times)):
                                adj = []

                                # First Adjacency
                                # Top row, even triangles
                                if i == 0 and (index - time) % 2 == 0:
                                    adj.append(None)
                                # Other even triangles
                                elif (index - time) % 2 == 0:
                                    adj.append((((i - 1) * (regular[1] - 1) + j) * 2 + 1) * len(self.times) + time)
                                # Right side, odd triangles
                                elif j == regular[1] - 2:
                                    adj.append(None)
                                # Other odd triangles
                                else:
                                    adj.append(((index - time) // len(self.times) + 1) * len(self.times) + time)

                                # Second Adjacency
                                # Bottom row, odd triangles
                                if i == regular[0] - 2 and index % 2 == 1:
                                    adj.append(None)
                                # Other odd triangles
                                elif (index - time) % 2 == 1:
                                    adj.append((((i + 1) * (regular[1] - 1) + j) * 2) * len(self.times) + time)
                                # Even triangles
                                else:
                                    adj.append(((index - time) // len(self.times) + 1) * len(self.times) + time)

                                # Third Adjacency
                                # Left side, even triangles
                                if i == 0 and (index - time) % 2 == 0:
                                    adj.append(None)
                                # All other triangles
                                else:
                                    adj.append(((index - time) // len(self.times) - 1) * len(self.times) + time)

                                if (index - time) % 2 == 0:
                                    elements.append(TriangularMeshElement2D(
                                        self._cpp_vertex_list[(i * regular[1] + j) * len(self.times) + time],
                                        self._cpp_vertex_list[(i * regular[1] + j + 1) * len(self.times) + time],
                                        self._cpp_vertex_list[((i + 1) * regular[1] + j) * len(self.times) + time],
                                        index
                                    ))
                                else:
                                    elements.append(TriangularMeshElement2D(
                                        self._cpp_vertex_list[(i * regular[1] + j + 1) * len(self.times) + time],
                                        self._cpp_vertex_list[((i + 1) * regular[1] + j + 1) * len(self.times) + time],
                                        self._cpp_vertex_list[((i + 1) * regular[1] + j) * len(self.times) + time],
                                        index
                                    ))

                                adjacency_list.append(adj)
                                index += 1

        self.element_list = elements
        self.adjacency_list = adjacency_list






