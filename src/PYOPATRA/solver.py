from .pyopatra_pybind import CppTriangularMesh2DSolver
from .particle import ParticleList
from .objective_functions import ObjectiveFunctionBase
from .mesh import MeshBase
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


class Solver:
    def __init__(self, times, mesh: MeshBase, particles: ParticleList, objective_function=None, wind_times=None):
        # print("times in solver {}".format(times))
        self.mesh = mesh
        self.particles = particles
        self.objective_function = objective_function

        if wind_times is None:
            wind_times = np.empty((0,))

        if objective_function is None:
            self._cpp_solver = CppTriangularMesh2DSolver(mesh._cpp_mesh.get_pointer_wrapper(),
                                                         particles._cpp_particle_list.get_pointer_wrapper(), times, wind_times)
        else:
            self._cpp_solver = CppTriangularMesh2DSolver(mesh._cpp_mesh.get_pointer_wrapper(),
                                                         particles._cpp_particle_list.get_pointer_wrapper(),
                                                         objective_function._cpp_obj_fn.get_pointer_wrapper(), times, wind_times)

    def time_step(self, time_delta):
        self._cpp_solver.time_step(time_delta)

    def reset_solver(self):
        self._cpp_solver.reset_solver()

    def update_particle_location_indices(self):
        self._cpp_solver.update_particle_location_indices()

    def calculate_objective_value(self):
        return self._cpp_solver.calculate_objective_value()

    def get_current_time(self):
        return self._cpp_solver.get_current_time()
