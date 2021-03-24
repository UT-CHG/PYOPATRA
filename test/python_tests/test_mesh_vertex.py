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

    def test_get_velocity(self, vertex):
        assert np.linalg.norm(vertex.get_velocity() - np.array((3, 4))) < 10e-6;
