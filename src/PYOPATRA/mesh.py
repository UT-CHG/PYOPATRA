import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

from PYOPATRA import FileParserBase, MeshVertex2D
from .pyopatra_pybind import CppTriangularMesh2D, TriangularMeshElement2D


class MeshBase:
    def __init__(self):
        self.vertex_list = None
        self._cpp_vertex_list = None

        self.adjacency_list = None
        self.element_list = None

        self.times = None
        self.regular_dimensions = None

        self._cpp_mesh = None

        self.dimensions = None

    # Not efficient! For testing purposes only
    def get_vertex_locations(self):
        return self._cpp_mesh.get_vertex_locations()

    # Not efficient! For testing purposes only
    def get_velocities(self, time_index):
        return self._cpp_mesh.get_velocities(time_index)

    # def update_particle_locations(self):
    #     self._cpp_mesh.update_particle_locations()

    def get_water_column_pointer(self, index):
        return self._cpp_mesh.get_water_column_pointer(index)

    def get_water_column_adjacencies(self, index):
        return self._cpp_mesh.get_water_column_adjacencies(index)

    def check_water_column_adjacency(self, origin_index, destination_index, side):
        return self._cpp_mesh.check_water_column_adjacency(origin_index, destination_index, side)

    # def append_particle(self, location):
    #     if self.dimensions == 2:
    #         self._cpp_mesh.add_particle(location)
    #     else:
    #         raise NotImplementedError('Only 2D currently implemented')

    # Not efficient! Only use for testing purposes or infrequently
    # def get_all_particle_locations(self):
    #     return self._cpp_mesh.get_all_particle_locations()
    #
    # def save_particle_locations_hdf5(self, filename):
    #     pass
    #
    # def time_step(self, time_delta):
    #     self._cpp_mesh.time_step(time_delta)
    #
    # def reset_mesh(self):
    #     self._cpp_mesh.reset_mesh()
    #
    # def setup_objective_function(self, particle_locations, num_bins_lat_long=None, bounds=None, num_proj=None, seed=None):
    #     temp = np.array(particle_locations, order='F')
    #     self._cpp_mesh.create_sliced_wasserstein_distance(temp, num_bins_lat_long[0], num_bins_lat_long[1], np.array(bounds), num_proj, seed)
    #
    # def get_objective_value(self):
    #     return self._cpp_mesh.calculate_objective_function()


class TriangularMesh(MeshBase):
    def __init__(self):
        super().__init__()


class TriangularMesh2D(TriangularMesh):
    def __init__(self):
        super().__init__()

        self.dimensions = 2

    def setup_mesh(self, file_parser: FileParserBase, dimensions: int):
        self._setup_mesh(file_parser, dimensions)
        self._setup_vertices(file_parser)
        self._setup_elements_vertices()
        self._setup_water_columns(file_parser)

        comm.barrier()
        print('Finished Set Up')

    def _setup_mesh(self, file_parser: FileParserBase, dimensions: int):
        if dimensions == 2:
            # print(file_parser.num_elements, file_parser.num_vertices, file_parser.times)
            self._cpp_mesh = CppTriangularMesh2D(file_parser.num_elements, file_parser.num_vertices, file_parser.times)
        elif dimensions == 3:
            raise NotImplementedError('3D Meshes are not yet implemented.')
        else:
            raise NotImplementedError('Meshes must be 2D or 3D.')

    def _setup_vertices(self, file_parser: FileParserBase):
        self.vertex_list = []
        self._cpp_vertex_list = []
        index = 0

        self.times = file_parser.times

        num_time_steps = int(len(self.times) / size) + (1 if len(self.times) % size > rank else 0)

        ts_list = comm.gather(num_time_steps, root=0)
        if rank == 0:
            print(ts_list)
            ts_list_array = np.array(ts_list)
            print('Time steps:', len(file_parser.times), np.sum(ts_list_array))

        if rank == 0:
            for i, ts in enumerate(ts_list):
                if i == 0:
                    velocities = file_parser.velocity[:, :, :ts]
                    diffusions = file_parser.diffusion_coefficient[:, :, :ts]
                    start_time = 0
                else:
                    temp = file_parser.velocity[:, :, np.sum(ts_list_array[:i]):np.sum(ts_list_array[:i+1])]
                    comm.send(temp, dest=i, tag=0)
                    comm.send(np.sum(ts_list_array[:i]), dest=i, tag=1)
                    temp = file_parser.diffusion_coefficient[:, :, np.sum(ts_list_array[:i]):np.sum(ts_list_array[:i+1])]
                    comm.send(temp, dest=i, tag=2)
        else:
            velocities = comm.recv(source=0, tag=0)
            start_time = comm.recv(source=0, tag=1)
            diffusions = comm.recv(source=0, tag=2)
            print('Rank {} received velocity of shape {} and start time {}'.format(rank, velocities.shape, start_time))

        if file_parser.regular_dimensions is not None:
            self.regular_dimensions = file_parser.regular_dimensions
            for i in range(file_parser.regular_dimensions[0]):
                for j in range(file_parser.regular_dimensions[1]):
                    if rank == 0:
                        self._cpp_mesh.set_vertex_location(i * file_parser.regular_dimensions[1] + j,
                                                           [file_parser.latitude[i], file_parser.longitude[j]])
                    for time in range(velocities.shape[2]):
                        self._cpp_mesh.set_vertex_velocity(i * file_parser.regular_dimensions[1] + j, start_time + time,
                                                           velocities[:, i * self.regular_dimensions[1] + j, time])
                        self._cpp_mesh.set_vertex_diffusion(i * file_parser.regular_dimensions[1] + j, start_time + time,
                                                           diffusions[:, i * self.regular_dimensions[1] + j, time])

    def _setup_elements_vertices(self, depths_array=None):
        if self.regular_dimensions is not None:
            index = 0

            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0] - 1):
                    for j in range(self.regular_dimensions[1] - 1):
                        for k in range(2):
                            if k == 0:
                                self._cpp_mesh.set_element_vertex((i * (self.regular_dimensions[1] - 1) + j) * 2,
                                                                  0, 0,
                                                                  i * self.regular_dimensions[1] + j)

                                self._cpp_mesh.set_element_vertex((i * (self.regular_dimensions[1] - 1) + j) * 2,
                                                                  0, 1,
                                                                  i * self.regular_dimensions[1] + j + 1)

                                self._cpp_mesh.set_element_vertex((i * (self.regular_dimensions[1] - 1) + j) * 2,
                                                                  0, 2,
                                                                  (i + 1) * self.regular_dimensions[1] + j)

                            else:
                                self._cpp_mesh.set_element_vertex((i * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                  0, 0,
                                                                  i * self.regular_dimensions[1] + j + 1)

                                self._cpp_mesh.set_element_vertex((i * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                  0, 1,
                                                                  (i + 1) * self.regular_dimensions[1] + j + 1)

                                self._cpp_mesh.set_element_vertex((i * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                  0, 2,
                                                                  (i + 1) * self.regular_dimensions[1] + j)

    def _setup_water_columns(self, file_parser: FileParserBase):
        if self.regular_dimensions is not None:
            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0] - 1):
                    for j in range(self.regular_dimensions[1] - 1):
                        for k in range(2):
                            # First Adjacency
                            # Top row, even triangles
                            if k == 0 and i == 0:
                                # No adjacency, should stay a nullptr
                                pass

                            # Other even triangles
                            elif k == 0:
                                self._cpp_mesh.set_water_column_adjacency((i * (self.regular_dimensions[1] - 1) + j) * 2,
                                                                          ((i - 1) * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                          0)
                            # Right side, odd triangles
                            elif j == self.regular_dimensions[1] - 2:
                                # No adjancency, should stay a nullptr
                                pass

                            # Other odd triangles
                            else:
                                self._cpp_mesh.set_water_column_adjacency((i * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                          (i * (self.regular_dimensions[1] - 1) + j) * 2 + 2,
                                                                          0)

                            # Second Adjacency
                            # Bottom row, odd triangles
                            if i == self.regular_dimensions[0] - 2 and k == 1:
                                # No adjancency
                                pass

                            # Other odd triangles
                            elif k == 1:
                                self._cpp_mesh.set_water_column_adjacency((i * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                          ((i + 1) * (self.regular_dimensions[1] - 1) + j) * 2,
                                                                          1)

                            # Even triangles
                            else:
                                self._cpp_mesh.set_water_column_adjacency((i * (self.regular_dimensions[1] - 1) + j) * 2,
                                                                          (i * (self.regular_dimensions[1] - 1) + j) * 2 + 1,
                                                                          1)

                            # Third Adjacency
                            # Left side, even triangles
                            if j == 0 and k == 0:
                                # No adjacency
                                pass

                            # All other triangles
                            else:
                                self._cpp_mesh.set_water_column_adjacency((i * (self.regular_dimensions[1] - 1) + j) * 2 + k,
                                                                          (i * (self.regular_dimensions[1] - 1) + j) * 2 - 1 + k,
                                                                          2)

    def get_velocity_u(self):
        if self.regular_dimensions is not None:
            velocity = np.zeros((*self.regular_dimensions, len(self.times)))

            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0]):
                    for j in range(self.regular_dimensions[1]):
                        for t in range(len(self.times)):
                            velocity[i, j, t] = self.vertex_list[((i * self.regular_dimensions[1]) + j) * len(self.times) + t].get_velocity()[0]

            return velocity

