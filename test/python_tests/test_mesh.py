import pytest
import os
import numpy as np

from PYOPATRA import *


class TestMesh:
    @pytest.fixture
    def triangular_2d_hycom_mesh(self):
        file_prefix = os.path.dirname(os.path.realpath(__file__)) + '/hycom_data'
        hycom_files = [
            file_prefix + '/hycom_gomu_501_2010042000_t000.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t003.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t006.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t009.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t012.nc',
        ]

        hfp = HYCOMFileParser()
        hfp.read(hycom_files, diffusion_coefficient=10.0)

        tm2d = TriangularMesh2D()
        tm2d.setup_vertices(hfp)
        tm2d.setup_elements_and_adjacency_list(regular=hfp.regular_dimensions)

        yield tm2d

    def test_hycom_mesh_setup(self, triangular_2d_hycom_mesh):
        assert len(triangular_2d_hycom_mesh.vertex_list) == 187186 * 5
        assert np.linalg.norm(triangular_2d_hycom_mesh.vertex_list[0].get_diffusion_coefficient() - np.array((10.0, 10.0))) < 10e-6
