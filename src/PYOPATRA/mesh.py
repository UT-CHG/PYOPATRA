import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

sharedcomm = comm.Split_type(MPI.COMM_TYPE_SHARED, key=rank)
sharedsize = sharedcomm.Get_size()
sharedrank = sharedcomm.Get_rank()

from PYOPATRA import FileParserBase
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

        self.wind_coef = None

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
        if rank == 0:
            print('Finished Set Up')

    def _setup_mesh(self, file_parser: FileParserBase, dimensions: int):
        if dimensions == 2:
            self._cpp_mesh = CppTriangularMesh2D(file_parser.num_elements, file_parser.num_vertices, file_parser.times,
                                                 file_parser.wind_times, file_parser.wind_coef)
        elif dimensions == 3:
            raise NotImplementedError('3D Meshes are not yet implemented.')
        else:
            raise NotImplementedError('Meshes must be 2D or 3D.')

    def _setup_vertices(self, file_parser: FileParserBase):
        self.vertex_list = []
        self._cpp_vertex_list = []
        index = 0

        self.times = file_parser.times

        if file_parser.regular_dimensions is not None:
            self.regular_dimensions = file_parser.regular_dimensions

            num_time_steps = int(len(self.times) / sharedsize) + (1 if len(self.times) % sharedsize > sharedrank else 0)

            ts_list = sharedcomm.gather(num_time_steps, root=0)
            if sharedrank == 0:
                ts_list_array = np.array(ts_list)

                for i, ts in enumerate(ts_list):
                    if i == 0:
                        velocities = file_parser.velocity[:, :, :ts]
                        diffusions = file_parser.diffusion_coefficient[:, :, :ts]
                        start_time = 0
                    else:
                        temp = file_parser.velocity[:, :, np.sum(ts_list_array[:i]):np.sum(ts_list_array[:i+1])]
                        sharedcomm.send(temp, dest=i, tag=0)
                        sharedcomm.send(np.sum(ts_list_array[:i]), dest=i, tag=1)
                        temp = file_parser.diffusion_coefficient[:, :, np.sum(ts_list_array[:i]):np.sum(ts_list_array[:i+1])]
                        sharedcomm.send(temp, dest=i, tag=2)
            else:
                velocities = sharedcomm.recv(source=0, tag=0)
                start_time = sharedcomm.recv(source=0, tag=1)
                diffusions = sharedcomm.recv(source=0, tag=2)

            velocities = np.moveaxis(velocities, 0, 2).reshape((-1, 2))
            diffusions = np.moveaxis(diffusions, 0, 2).reshape((-1, 2))

            self._cpp_mesh.set_velocities(velocities)
            self._cpp_mesh.set_diffusions(diffusions)

            if sharedrank == 0:
                num_lat, num_lon = file_parser.regular_dimensions[:2]
                vertex_locations = np.zeros((num_lat, num_lon, 2))
                vertex_locations[:, :, 0] = file_parser.latitude.reshape((num_lat, 1))
                vertex_locations[:, :, 1] = file_parser.longitude.reshape((1, num_lon))
                self._cpp_mesh.set_vertex_locations(vertex_locations.reshape((-1, 2)))


            if sharedrank == 0:
                if file_parser.winds is None: return
                winds = file_parser.winds.reshape((-1, 2))
                self._cpp_mesh.set_winds(winds)

    def _setup_elements_vertices(self, depths_array=None):
        if self.regular_dimensions is not None:
            index = 0
            m, n = self.regular_dimensions[0] - 1, self.regular_dimensions[1] - 1
            vertices = np.zeros((m,n, 2, 3), dtype=np.int64)
            i = np.arange(m).reshape((m, 1))
            j = np.arange(n).reshape((1, n))
            # even elements
            vertices[:, :, 0, 0] = i * self.regular_dimensions[1] + j
            vertices[:, :, 0, 1] = i * self.regular_dimensions[1] + j + 1
            vertices[:, :, 0, 2] = (i + 1) * self.regular_dimensions[1] + j

            # odd ones
            vertices[:, :, 1, 0] = i * self.regular_dimensions[1] + j + 1
            vertices[:, :, 1, 1] = (i + 1) * self.regular_dimensions[1] + j + 1
            vertices[:, :, 1, 2] = (i + 1) * self.regular_dimensions[1] + j

            self._cpp_mesh.set_element_vertices(vertices.reshape((-1, 3)))

    def _setup_water_columns(self, file_parser: FileParserBase):
        if self.regular_dimensions is not None:
            if len(self.regular_dimensions) == 2:
                m, n = self.regular_dimensions[0] - 1, self.regular_dimensions[1] - 1
                inds = np.arange(m*n*2).reshape((m, n, 2))
                # First adjaency
                adjacencies = np.zeros((m, n, 2, 3), dtype=np.int64)
                adjacencies[:] = -1
                # even triangles
                adjacencies[1:, :, 0, 0] = inds[1:, :, 0] - n
                # odd ones
                adjacencies[:, :-1, 1, 0] = inds[:, :-1, 1] + 1

                # Second adjacency
                # even triangles
                adjacencies[:, :, 0, 1] = inds[:, :, 0] + 1
                # odd ones
                adjacencies[:-1, :, 1, 1] = inds[:-1, :, 1] + n
                
                # Third adjacency
                # even triangles
                adjacencies[:, 1:, 0, 2] = inds[:, 1:, 0] - 1
                # odd ones
                adjacencies[:, :, 1, 2] = inds[:, :, 1] - 1

                self._cpp_mesh.set_water_column_adjacencies(adjacencies.reshape((-1, 3)))

    def get_velocity_u(self):
        if self.regular_dimensions is not None:
            velocity = np.zeros((*self.regular_dimensions, len(self.times)))

            if len(self.regular_dimensions) == 2:
                for i in range(self.regular_dimensions[0]):
                    for j in range(self.regular_dimensions[1]):
                        for t in range(len(self.times)):
                            velocity[i, j, t] = self.vertex_list[((i * self.regular_dimensions[1]) + j) * len(self.times) + t].get_velocity()[0]

            return velocity

    def set_wind_coef(self, new_coef):
        if sharedrank == 0:
            self._cpp_mesh.set_wind_coef(new_coef)

