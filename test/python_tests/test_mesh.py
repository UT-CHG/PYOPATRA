import pytest
import os
import numpy as np

from PYOPATRA import *


class TestMesh:
    @pytest.fixture
    def triangular_2d_hycom_mesh(self):
        # file_prefix = os.path.dirname(os.path.realpath(__file__)) + '/hycom_data'
        # hycom_files = [
        #     file_prefix + '/hycom_gomu_501_2010042000_t000.nc',
        #     file_prefix + '/hycom_gomu_501_2010042000_t003.nc',
        #     file_prefix + '/hycom_gomu_501_2010042000_t006.nc',
        #     file_prefix + '/hycom_gomu_501_2010042000_t009.nc',
        #     file_prefix + '/hycom_gomu_501_2010042000_t012.nc',
        # ]
        #
        # hfp = HYCOMFileParser()
        # hfp.read(hycom_files, diffusion_coefficient=10.0)

        dummy_file = FileParserBase()
        dummy_file.vertices_per_polygon = 3
        dummy_file.times = np.array([3, 6, 9], dtype=int)
        dummy_file.regular_dimensions = (5, 4)
        dummy_file.num_vertices = 5 * 4
        dummy_file.num_elements = (5 - 1) * 2 * (4 - 1)
        dummy_file.latitude = np.array([1, 2, 3, 4, 5])
        dummy_file.longitude = np.array([11, 12, 13, 14])
        dummy_file.velocity = np.zeros((2, dummy_file.num_vertices, len(dummy_file.times)))
        dummy_file.velocity[0, :, 0] = np.array([
            [1.1, 1.2, 1.3, 1.4, 1.5],
            [2.1, 2.2, 2.3, 2.4, 2.5],
            [3.1, 3.2, 3.3, 3.4, 3.5],
            [4.1, 4.2, 4.3, 4.4, 4.5]
        ]).flatten()
        dummy_file.velocity[0, :, 1] = np.array([
            [2.1, 2.2, 2.3, 2.4, 2.5],
            [3.1, 3.2, 3.3, 3.4, 3.5],
            [4.1, 4.2, 4.3, 4.4, 4.5],
            [5.1, 5.2, 5.3, 5.4, 5.5]
        ]).flatten()
        dummy_file.velocity[0, :, 2] = np.array([
            [3.1, 3.2, 3.3, 3.4, 3.5],
            [4.1, 4.2, 4.3, 4.4, 4.5],
            [5.1, 5.2, 5.3, 5.4, 5.5],
            [6.1, 6.2, 6.3, 6.4, 6.5]
        ]).flatten()
        dummy_file.diffusion_coefficient = 10.0 * np.ones((2, dummy_file.num_vertices, len(dummy_file.times)))

        tm2d = TriangularMesh2D()
        tm2d.setup_mesh(dummy_file, 2)

        yield tm2d

    def test_particle_list_setup(self, triangular_2d_hycom_mesh):
        particle_locations = np.array(
            [[1.5, 11.5],
             [2.5, 11.5],
             [4.5, 13.5]],
            order='F'
        )

        for row_index in range(particle_locations.shape[0]):
            triangular_2d_hycom_mesh.append_particle(particle_locations[row_index, :])

        retrieved_particle_locations = triangular_2d_hycom_mesh.get_all_particle_locations()
        print(retrieved_particle_locations)
        assert np.linalg.norm(particle_locations - retrieved_particle_locations) < 10e-6

    def test_hycom_mesh_setup(self, triangular_2d_hycom_mesh):
        latitude = np.array([1, 2, 3, 4, 5])
        longitude = np.array([11, 12, 13, 14])

        X, Y = np.meshgrid(latitude, longitude)
        locations = np.zeros((len(X.flatten()), 2))
        locations[:, 0] = X.flatten(order='F')
        locations[:, 1] = Y.flatten(order='F')

        returned_locations = triangular_2d_hycom_mesh.get_vertex_locations()

        assert np.linalg.norm(locations - returned_locations) < 10e-6