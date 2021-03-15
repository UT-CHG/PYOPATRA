#
#  Created by Georgia Stuart on 24 February 2021
#


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
    pass


class POMFileParser(FileParserBase):
    pass