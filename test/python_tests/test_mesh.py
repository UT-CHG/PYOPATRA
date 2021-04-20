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
        tm2d.setup_elements_and_adjacency_list()

        yield tm2d

    def test_hycom_mesh_setup(self, triangular_2d_hycom_mesh):
        assert len(triangular_2d_hycom_mesh.vertex_list) == 187186 * 5
        assert len(triangular_2d_hycom_mesh.element_list) == (346 - 1) * 2 * (541 - 1) * 5
        assert np.linalg.norm(triangular_2d_hycom_mesh.vertex_list[0].get_diffusion_coefficient() - np.array((10.0, 10.0))) < 10e-6
        assert triangular_2d_hycom_mesh.adjacency_list[0] == [None, 5, None]
        assert triangular_2d_hycom_mesh.adjacency_list[1] == [None, 6, None]
        assert triangular_2d_hycom_mesh.adjacency_list[5] == [10, (541 - 1) * 2 * 5, 0]

        assert triangular_2d_hycom_mesh.element_list[0].get_vertices()[0] is triangular_2d_hycom_mesh.vertex_list[0].vertex
        assert triangular_2d_hycom_mesh.element_list[0].get_vertices()[1] is triangular_2d_hycom_mesh.vertex_list[5].vertex
        assert triangular_2d_hycom_mesh.element_list[0].get_vertices()[2] is triangular_2d_hycom_mesh.vertex_list[541 * 5].vertex

        assert triangular_2d_hycom_mesh.element_list[1].get_vertices()[0] is triangular_2d_hycom_mesh.vertex_list[1].vertex
        assert triangular_2d_hycom_mesh.element_list[1].get_vertices()[1] is triangular_2d_hycom_mesh.vertex_list[6].vertex
        assert triangular_2d_hycom_mesh.element_list[1].get_vertices()[2] is triangular_2d_hycom_mesh.vertex_list[541 * 5 + 1].vertex

        assert triangular_2d_hycom_mesh.element_list[5].get_vertices()[0] is triangular_2d_hycom_mesh.vertex_list[5].vertex
        assert triangular_2d_hycom_mesh.element_list[5].get_vertices()[1] is triangular_2d_hycom_mesh.vertex_list[542 * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[5].get_vertices()[2] is triangular_2d_hycom_mesh.vertex_list[541 * 5].vertex

        assert triangular_2d_hycom_mesh.element_list[540 * 2 * 5].get_vertices()[0] is triangular_2d_hycom_mesh.vertex_list[(541 + 0) * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[540 * 2 * 5].get_vertices()[1] is triangular_2d_hycom_mesh.vertex_list[(541 + 1) * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[540 * 2 * 5].get_vertices()[2] is triangular_2d_hycom_mesh.vertex_list[(2 * 541 + 0) * 5].vertex

        assert triangular_2d_hycom_mesh.element_list[541 * 2 * 5].get_vertices()[0] is triangular_2d_hycom_mesh.vertex_list[(541 + 1) * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[541 * 2 * 5].get_vertices()[1] is triangular_2d_hycom_mesh.vertex_list[(541 + 2) * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[541 * 2 * 5].get_vertices()[2] is triangular_2d_hycom_mesh.vertex_list[(2 * 541 + 1) * 5].vertex

        assert triangular_2d_hycom_mesh.element_list[541 * 2 * 5 + 5].get_vertices()[0] is triangular_2d_hycom_mesh.vertex_list[(541 + 2) * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[541 * 2 * 5 + 5].get_vertices()[1] is triangular_2d_hycom_mesh.vertex_list[(2 * 541 + 2) * 5].vertex
        assert triangular_2d_hycom_mesh.element_list[541 * 2 * 5 + 5].get_vertices()[2] is triangular_2d_hycom_mesh.vertex_list[(2 * 541 + 1) * 5].vertex


