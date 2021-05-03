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

        self._cpp_mesh = None

    def _setup_mesh(self, file_parser: FileParserBase):
        pass

    def _setup_vertices(self, file_parser: FileParserBase):
        self.vertex_list = []
        self._cpp_vertex_list = []
        index = 0

        self.times = file_parser.times

        if file_parser.regular_dimensions is not None:
            self.regular_dimensions = file_parser.regular_dimensions
            for i in range(file_parser.regular_dimensions[0]):
                for j in range(file_parser.regular_dimensions[1]):
                    self._cpp_mesh.get_vertices()[i * file_parser.regular_dimensions[1] + j]\
                        .set_location([file_parser.latitude[i], file_parser.longitude[j]])
                    for time in range(len(self.times)):
                        self._cpp_mesh.get_vertices()[i * file_parser.regular_dimensions[1] + j]\
                            .set_velocity(file_parser.velocity[:, i * self.regular_dimensions[1] + j, time], time)
                        self._cpp_mesh.get_vertices()[i * file_parser.regular_dimensions[1] + j] \
                            .set_diffusion_coefficient(file_parser.velocity[:, i * self.regular_dimensions[1] + j, time], time)
                    #     self.vertex_list.append(MeshVertex2D(
                    #         file_parser.latitude[i],
                    #         file_parser.longitude[j],
                    #         file_parser.velocity[:, i * self.regular_dimensions[1] + j, time],
                    #         file_parser.diffusion_coefficient[:, i * self.regular_dimensions[1] + j, time],
                    #         len(file_parser.regular_dimensions)
                    #     ))
                    #
                    #     self._cpp_vertex_list.append(self.vertex_list[-1].vertex)
                    #     index += 1

    def _setup_elements_and_adjacency_list(self, depths_array=None):
        elements = []
        adjacency_list = []
        water_columns = []

        if self.vertex_list is None:
            raise AttributeError('vertex_list must be set prior to setting up elements and adjacency.')

        if self.regular_dimensions is not None:
            index = 0

            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0] - 1):
                    for j in range(self.regular_dimensions[1] - 1):
                        for k in range(2):
                            for time in range(len(self.times)):
                                # adj = []
                                #
                                # # First Adjacency
                                # # Top row, even triangles
                                # if i == 0 and (index - time) % 2 == 0:
                                #     adj.append(None)
                                # # Other even triangles
                                # elif (index - time) % 2 == 0:
                                #     adj.append((((i - 1) * (self.regular_dimensions[1] - 1) + j) * 2 + 1) * len(self.times) + time)
                                # # Right side, odd triangles
                                # elif j == self.regular_dimensions[1] - 2:
                                #     adj.append(None)
                                # # Other odd triangles
                                # else:
                                #     adj.append(((index - time) // len(self.times) + 1) * len(self.times) + time)
                                #
                                # # Second Adjacency
                                # # Bottom row, odd triangles
                                # if i == self.regular_dimensions[0] - 2 and index % 2 == 1:
                                #     adj.append(None)
                                # # Other odd triangles
                                # elif (index - time) % 2 == 1:
                                #     adj.append((((i + 1) * (self.regular_dimensions[1] - 1) + j) * 2) * len(self.times) + time)
                                # # Even triangles
                                # else:
                                #     adj.append(((index - time) // len(self.times) + 1) * len(self.times) + time)
                                #
                                # # Third Adjacency
                                # # Left side, even triangles
                                # if i == 0 and (index - time) % 2 == 0:
                                #     adj.append(None)
                                # # All other triangles
                                # else:
                                #     adj.append(((index - time) // len(self.times) - 1) * len(self.times) + time)

                                if (index - time) % 2 == 0:
                                    elements.append(TriangularMeshElement2D(
                                        self._cpp_vertex_list[(i * self.regular_dimensions[1] + j) * len(self.times) + time],
                                        self._cpp_vertex_list[(i * self.regular_dimensions[1] + j + 1) * len(self.times) + time],
                                        self._cpp_vertex_list[((i + 1) * self.regular_dimensions[1] + j) * len(self.times) + time],
                                        index
                                    ))
                                else:
                                    elements.append(TriangularMeshElement2D(
                                        self._cpp_vertex_list[(i * self.regular_dimensions[1] + j + 1) * len(self.times) + time],
                                        self._cpp_vertex_list[((i + 1) * self.regular_dimensions[1] + j + 1) * len(self.times) + time],
                                        self._cpp_vertex_list[((i + 1) * self.regular_dimensions[1] + j) * len(self.times) + time],
                                        index
                                    ))

                                # adjacency_list.append(adj)
                                index += 1

        self.element_list = elements
        self.adjacency_list = adjacency_list

    def _setup_water_columns(self):
        if self.element_list is None or self.vertex_list is None:
            raise AttributeError('Vertex and elements must be set up prior to water columns.')
        
        if self.regular_dimensions is not None:
            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0] - 1):
                    for j in range(self.regular_dimensions[1] - 1):
                        for k in range(2):
                            pass

    def get_velocity_u(self):
        if self.regular_dimensions is not None:
            velocity = np.zeros((*self.regular_dimensions, len(self.times)))

            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0]):
                    for j in range(self.regular_dimensions[1]):
                        for t in range(len(self.times)):
                            velocity[i, j, t] = self.vertex_list[((i * self.regular_dimensions[1]) + j) * len(self.times) + t].get_velocity()[0]

            return velocity




