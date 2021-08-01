import numpy as np

from .pyopatra_pybind import CppSlicedWassersteinDistance2D
from .particle import ParticleList, comm, rank

class ObjectiveFunctionBase:
    def __init__(self):
        self._cpp_obj_fn = None

    def set_observed_values(self, observed_locations):
        self._cpp_obj_fn.set_observed_values(np.array(observed_locations, order='F'))


class SlicedWassersteinDistance(ObjectiveFunctionBase):
    def __init__(self, num_bins_lat, num_bins_lon, bounds, num_projections, rng_seed):
        super().__init__()
        self._cpp_obj_fn = CppSlicedWassersteinDistance2D(num_bins_lat, num_bins_lon, np.array(bounds), num_projections, rng_seed)
