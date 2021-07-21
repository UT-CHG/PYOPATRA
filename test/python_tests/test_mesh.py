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
        dummy_file.regular_dimensions = (4, 5)
        dummy_file.num_vertices = 5 * 4
        dummy_file.num_elements = (5 - 1) * 2 * (4 - 1)
        dummy_file.latitude = np.array([1, 2, 3, 4])
        dummy_file.longitude = np.array([11, 12, 13, 14, 15])
        dummy_file.velocity = np.zeros((2, dummy_file.num_vertices, len(dummy_file.times)), order='F')
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
        dummy_file.velocity[1, :, 2] = np.array([
            [1.1, 1.2, 1.3, 1.4, 1.5],
            [2.1, 2.2, 2.3, 2.4, 2.5],
            [3.1, 3.2, 3.3, 3.4, 3.5],
            [4.1, 4.2, 4.3, 4.4, 4.5]
        ]).flatten()
        dummy_file.velocity[1, :, 0] = np.array([
            [2.1, 2.2, 2.3, 2.4, 2.5],
            [3.1, 3.2, 3.3, 3.4, 3.5],
            [4.1, 4.2, 4.3, 4.4, 4.5],
            [5.1, 5.2, 5.3, 5.4, 5.5]
        ]).flatten()
        dummy_file.velocity[1, :, 1] = np.array([
            [3.1, 3.2, 3.3, 3.4, 3.5],
            [4.1, 4.2, 4.3, 4.4, 4.5],
            [5.1, 5.2, 5.3, 5.4, 5.5],
            [6.1, 6.2, 6.3, 6.4, 6.5]
        ]).flatten()

        dummy_file.diffusion_coefficient = 10.0 * np.ones((2, dummy_file.num_vertices, len(dummy_file.times)))

        tm2d = TriangularMesh2D()
        tm2d.setup_mesh(dummy_file, 2)

        yield tm2d, dummy_file

    def test_particle_list_setup(self, triangular_2d_hycom_mesh):
        mesh = triangular_2d_hycom_mesh[0]
        dummy_file = triangular_2d_hycom_mesh[1]

        particle_locations = np.array(
            [[1.5, 11.5],
             [2.5, 11.5],
             [4.5, 13.5]],
            order='F'
        )

        for row_index in range(particle_locations.shape[0]):
            mesh.append_particle(particle_locations[row_index, :])

        retrieved_particle_locations = mesh.get_all_particle_locations()
        print(retrieved_particle_locations)
        assert np.linalg.norm(particle_locations - retrieved_particle_locations) < 10e-6

    def test_hycom_mesh_setup(self, triangular_2d_hycom_mesh):
        mesh = triangular_2d_hycom_mesh[0]
        dummy_file = triangular_2d_hycom_mesh[1]

        latitude = dummy_file.latitude
        longitude = dummy_file.longitude

        X, Y = np.meshgrid(latitude, longitude)
        locations = np.zeros((len(X.flatten()), 2))
        locations[:, 0] = X.flatten(order='F')
        locations[:, 1] = Y.flatten(order='F')

        returned_locations = mesh.get_vertex_locations()
        returned_velocities = np.zeros((2, dummy_file.num_vertices, len(dummy_file.times)), order='F')

        for i in range(len(dummy_file.times)):
            temp = mesh.get_velocities(i)
            returned_velocities[0, :, i] = temp[:, 0]
            returned_velocities[1, :, i] = temp[:, 1]

        assert np.linalg.norm(locations - returned_locations) < 10e-6
        assert np.linalg.norm(dummy_file.velocity - returned_velocities) < 10e-6

        assert mesh._cpp_mesh.get_water_columns_size() == dummy_file.num_elements

        assert mesh.check_water_column_adjacency(0, 1, 1)
        assert mesh.check_water_column_adjacency(1, 0, 2)
        assert mesh.check_water_column_adjacency(1, 2, 0)
        assert mesh.check_water_column_adjacency(1, 8, 1)
        assert mesh.check_water_column_adjacency(2, 3, 1)
        assert mesh.check_water_column_adjacency(2, 1, 2)
        assert mesh.check_water_column_adjacency(10, 3, 0)
        assert mesh.check_water_column_adjacency(10, 11, 1)
        assert mesh.check_water_column_adjacency(10, 9, 2)
        assert mesh.check_water_column_adjacency(11, 12, 0)
        assert mesh.check_water_column_adjacency(11, 18, 1)
        assert mesh.check_water_column_adjacency(11, 10, 2)
        assert mesh.check_water_column_adjacency(23, 22, 2)
        assert mesh.check_water_column_adjacency(19, 20, 0)
        assert mesh.check_water_column_adjacency(19, 18, 2)






