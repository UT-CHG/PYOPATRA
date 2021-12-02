# import numpy as np
# from .pyopatra_pybind import CppMeshVertex2D
#
#
# class MeshVertex2D:
#     def __init__(self, latitude, longitude, velocity, diffusion_coefficient, dimensions):
#         if dimensions == 2:
#             self.vertex = CppMeshVertex2D(latitude, longitude, velocity, diffusion_coefficient)
#         else:
#             raise NotImplementedError('Dimensions other than 2 have not been implemented (yet)')
#
#     def get_location(self):
#         return self.vertex.get_location()
#
#     def get_latitude(self):
#         return self.vertex.get_latitude()
#
#     def get_longitude(self):
#         return self.vertex.get_longitude()
#
#     def get_velocity(self):
#         return self.vertex.get_velocity()
#
#     def get_diffusion_coefficient(self):
#         return self.vertex.get_diffusion_coefficient()
#
#     def set_location(self, new_location):
#         self.vertex.set_location(new_location)
#
#     def set_diffusion_coefficient(self, new_diffusion_coefficient):
#         self.vertex.set_diffusion_coefficient(new_diffusion_coefficient)
#
#     def set_velocity(self, new_velocity):
#         self.vertex.set_velocity(new_velocity)
