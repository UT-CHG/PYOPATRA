#
#  Created by Georgia Stuart on 24 February 2021
#

import numpy as np
import netCDF4 as nc
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
sharedcomm = comm.Split_type(MPI.COMM_TYPE_SHARED)
sharedrank = sharedcomm.Get_rank()
sharedmaster = rank if sharedrank == 0 else -1

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
        self.regular_dimensions = None
        """Dimensions for regular grid spacing (e.g., HYCOM)"""
        self.latitude = None
        """Latitude Numpy Array"""
        self.longitude = None
        """Longitude Numpy Array"""
        self.diffusion_coefficient = None
        """Diffusion Coefficient Matrix"""

        self.master_ranks = comm.gather(sharedmaster, root=0)
        # if rank == 0:
        #     print(self.master_ranks)


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

    def read(self, list_of_hycom_files, dimensions=2, triangulate=True, diffusion_coefficient=0.0):
        if triangulate:
            self.vertices_per_polygon = 3
        else:
            raise NotImplementedError('Only triangulated HYCOM data is implemented at this time.')

        self.times = np.zeros(len(list_of_hycom_files))

        if rank == 0:
            with nc.Dataset(list_of_hycom_files[0]) as ds:
                self.regular_dimensions = (ds['lat'].shape[0], ds['lon'].shape[0])
                self.num_vertices = ds['lat'].shape[0] * ds['lon'].shape[0]
                self.num_elements = (ds['lon'].shape[0] - 1) * 2 * (ds['lat'].shape[0] - 1)
                self.latitude = ds['lat'][:]
                self.longitude = ds['lon'][:]

        self.latitude = comm.bcast(self.latitude, root=0)
        self.longitude = comm.bcast(self.longitude, root=0)
        self.regular_dimensions = comm.bcast(self.regular_dimensions, root=0)
        self.num_vertices = comm.bcast(self.num_vertices, root=0)
        self.num_elements = comm.bcast(self.num_elements, root=0)

        if dimensions == 2:
            if rank == 0:
                self.velocity = np.zeros((2, self.num_vertices, len(list_of_hycom_files)))
                self.diffusion_coefficient = np.ones((2, self.num_vertices, len(list_of_hycom_files))) * diffusion_coefficient
        else:
            raise NotImplementedError('Dimensions other than 2 have not been implemented.')

        # TODO: Make diffusion coefficient more flexible

        if rank == 0:
            for index, filename in enumerate(list_of_hycom_files):
                with nc.Dataset(filename) as ds:
                    if dimensions == 2:
                        water_v = ds['water_v'][0, 0, :, :].flatten()
                        self.velocity[0, :, index] = water_v[:]
                        mv = self.velocity[0, :, index] == float(ds['water_v'].missing_value)
                        self.velocity[0, :, index][mv] = 0.0
                        self.diffusion_coefficient[0, :, index][mv] = 0.0

                        water_u = ds['water_u'][0, 0, :, :].flatten()
                        self.velocity[1, :, index] = water_u[:]
                        mv = self.velocity[1, :, index] == float(ds['water_u'].missing_value)
                        self.velocity[1, :, index][mv] = 0.0
                        self.diffusion_coefficient[0, :, index][mv] = 0.0

                        self.times[index] = ds['time'][0]
                    else:
                        raise NotImplementedError('Dimensions other than 2 have not been implemented.')

            for r in self.master_ranks:
                if r > 0:
                    comm.send(self.velocity, dest=r, tag=5)
                    comm.send(self.diffusion_coefficient, dest=r, tag=6)

        else:
            if sharedmaster > 0:
                self.velocity = comm.recv(source=0, tag=5)
                self.diffusion_coefficient = comm.recv(source=0, tag=6)

        self.times = comm.bcast(self.times, root=0)


class POMFileParser(FileParserBase):
    pass