import pytest
from PYOPATRA import Vertex
import numpy as np


class TestMeshVertex:
    @pytest.fixture
    def vertex(self):
        yield Vertex(10, -90, np.array((3, 4)), 2)

    def test_get_latitude(self, vertex):
        assert vertex.get_latitude() == 10

    def test_get_longitude(self, vertex):
        assert vertex.get_longitude() == -90

    def test_location(self, vertex):
        assert np.linalg.norm(vertex.get_location() - np.array((10, -90))) < 10e-6
        vertex.set_location(np.array((15.5, -80)))
        assert np.linalg.norm(vertex.get_location() - np.array((15.5, -80))) < 10e-6

    def test_velocity(self, vertex):
        assert np.linalg.norm(vertex.get_velocity() - np.array((3, 4))) < 10e-6
        vertex.set_velocity(np.array((5, 10.5)))
        assert np.linalg.norm(vertex.get_velocity() - np.array((5, 10.5))) < 10e-6

    def test_diffusion_coefficient(self, vertex):
        assert np.linalg.norm(vertex.get_diffusion_coefficient() - np.array((0, 0))) < 10e-6
        vertex.set_diffusion_coefficient(np.array((10, 20)))
        assert np.linalg.norm(vertex.get_diffusion_coefficient() - np.array((10, 20))) < 10e-6

