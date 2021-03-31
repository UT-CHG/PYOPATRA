#
#  Created by Georgia Stuart on 24 February 2021
#

import numpy as np
import netCDF4 as nc

class FileParserBase:
    """
    The base class for all file parsers.
    """
    def __init__(self):
        self.grid_name = None
        """Alpha-numeric grid identifier (name)"""
        self.num_elements = None
        """Number of grid elements"""
        self.num_vertices = None
        """Number of grid vertices"""
        self.vertices_per_polygon = None
        """Number of vertices that form each grid polygon"""
        self.vertices = None
        """Raw Vertex Information"""
        self.element_vertices = None
        """Defines the vertex indices that make up each element"""
        self.grid = None
        """Raw Grid Information"""
        self.boundary = None
        """Raw boundary information"""
        self.velocity = None
        """Raw velocity information"""
        self.temperature = None
        """Raw temperature information"""
        self.salinity = None
        """Raw salinity information"""
        self.density = None
        """Raw density information"""
        self.turbulence = None
        """Raw turbulence information"""
        self.times = None
        """Time spacing"""


class ADCIRCFileParser(FileParserBase):
    """
    File parsing for `ADCIRC <https://adcirc.org/>`_ input and output files.

    :param fort14: The ADCIRC `grid and boundary conditions file <https://adcirc.org/home/documentation/users-manual-v53/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14>`_
    :param fort44: The ADCIRC `3D Density, Temperature and/or Salinity at All Nodes in the Model Grid <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-density-temperature-and-or-salinity-nodes-model-grid-fort-44/>`_
    :param fort45: The ADCIRC `3D Velocity at All Nodes in the Model Grid <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-velocity-nodes-model-grid-fort-45/>`_
    :param fort46: The ADCIRC `3D Turbulence at All Nodes in the Model Grid <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-turbulence-nodes-model-grid-fort-46/>`_
    """
    def __init__(self, fort14=None, fort44=None, fort45=None, fort46=None):
        super().__init__()
        self.fort14 = fort14
        """ADCIRC `fort.14 <https://adcirc.org/home/documentation/users-manual-v53/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14>`_ file name and path"""
        self.fort44 = fort44
        """ADCIRC `fort.44 <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-density-temperature-and-or-salinity-nodes-model-grid-fort-44/>`_ file name and path"""
        self.fort45 = fort45
        """ADCIRC `fort.45 <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-velocity-nodes-model-grid-fort-45/>`_ file name and path"""
        self.fort46 = fort46
        """ADCIRC `fort.46 <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-turbulence-nodes-model-grid-fort-46/>`_ file name and path"""

    def read_grid_and_bc(self, file=None):
        """
        Reads the grid and boundary conditions from an `ADCIRC fort.14 file <https://adcirc.org/home/documentation/users-manual-v53/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14>`_

        :param file: The ADCIRC fort.14 file to be read. If None, uses the file specified in the constructor.
        """
        if not file:
            file = self.fort14

        if not file:
            raise FileNotFoundError('fort.14 file must be specified in either the ADCIRCFileParser constructor or in read_grid_and_bc function.')

        with open(file, 'r') as fp:
            self.grid_name = fp.readline().strip()
            num_elements_and_nodes = fp.readline().strip().split()
            self.num_elements = int(num_elements_and_nodes[0])
            self.num_vertices = int(num_elements_and_nodes[1])

            print('Reading {} grid with {} elements and {} vertices'.format(self.grid_name, self.num_elements, self.num_vertices))

            self.vertices = np.zeros((self.num_vertices, 3))

            for i in range(self.num_vertices):
                vertex = fp.readline().strip().split()
                index = int(vertex[0]) - 1
                coords = np.array(vertex[1:])
                self.vertices[index, :] = coords

            # We currently only support triangular ADCIRC meshes
            self.element_vertices = np.zeros((self.num_elements, 3), dtype='int')

            for i in range(self.num_elements):
                element = fp.readline().strip().split()
                index = int(element[0]) - 1
                vertices = np.array(element[2:], dtype='int') - 1
                self.element_vertices[index, :] = vertices

    def read_density_temperature_salinity(self, file=None):
        """
        Reads the density, temperature, and salinity from an `ADCIRC fort.44 file <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-density-temperature-and-or-salinity-nodes-model-grid-fort-44/>`_

        :param file: The ADCIRC fort.44 file to be read. If None, uses the file specified in the constructor.
        """
        pass

    def read_velocity(self, file=None):
        """
        Reads the velocity from an `ADCIRC fort.45 file <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-velocity-nodes-model-grid-fort-45/>`_

        :param file: The ADCIRC fort.45 file to be read. If None, uses the file specified in the constructor.
        """
        pass

    def read_turbulence(self, file=None):
        """
        Reads the turbulence from an `ADCIRC fort.46 file <https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/3d-turbulence-nodes-model-grid-fort-46/>`_

        :param file: The ADCIRC fort.46 file to be read. If None, uses the file specified in the constructor.
        """
        pass


class MITgcmFileParser(FileParserBase):
    pass


class HYCOMFileParser(FileParserBase):
    def __init__(self):
        super().__init__()

    def read(self, list_of_hycom_files, dimensions=2, triangulate=True):
        if triangulate:
            self.vertices_per_polygon = 3
        else:
            raise NotImplementedError('Only triangulated HYCOM data is implemented at this time.')

        self.times = np.zeros(len(list_of_hycom_files))

        with nc.Dataset(list_of_hycom_files[0]) as ds:
            self.num_vertices = ds['lat'].shape[0] * ds['lon'].shape[0]
            self.num_elements = (ds['lat'].shape[0] - 1) * 2 * ds['lon'].shape[0]

        if dimensions == 2:
            self.velocity = np.zeros((2, self.num_vertices, len(list_of_hycom_files)))
        else:
            raise NotImplementedError('Dimensions other than 2 have not been implemented.')

        for index, filename in enumerate(list_of_hycom_files):
            with nc.Dataset(filename) as ds:
                if dimensions == 2:
                    self.velocity[0, :, index] = ds['water_u'][0, 0, :, :].flatten()
                    self.velocity[1, :, index] = ds['water_v'][0, 0, :, :].flatten()
                    self.times[index] = ds['time'][0]
                else:
                    raise NotImplementedError('Dimensions other than 2 have not been implemented.')


class POMFileParser(FileParserBase):
    pass