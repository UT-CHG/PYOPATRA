import os
import numpy as np
import PYOPATRA as pta

class TestAdcircParsing:
    def test_read_grid_and_bc(self):
        afp = pta.ADCIRCFileParser(fort14=os.path.dirname(os.path.realpath(__file__)) + '/fort.14')
        afp.read_grid_and_bc()

        assert afp.grid_name == 'Pontchartrain Test'
        assert afp.num_elements == 620
        assert afp.num_vertices == 402

        assert np.allclose(afp.vertices[0, :], np.array((-89.4717160000, 30.0866470000, 3.0000000000)))
        assert np.allclose(afp.vertices[-1, :], np.array((-90.4216230000, 30.1109680000, 1.0691647500)))
        assert np.allclose(afp.element_vertices[0, :], np.array((0, 5, 6), dtype='int'))
        assert np.allclose(afp.element_vertices[-1, :], np.array((397, 398, 401), dtype='int'))

    # def test_read_density_temperature_salinity(self):
    #     assert False
    #
    # def test_read_velocity(self):
    #     assert False
    #
    # def test_read_turbulence(self):
    #     assert False

class TestHycomParsing:
    def test_read(self):
        file_prefix = os.path.dirname(os.path.realpath(__file__)) + '/hycom_data'
        hycom_files = [
            file_prefix + '/hycom_gomu_501_2010042000_t000.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t003.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t006.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t009.nc',
            file_prefix + '/hycom_gomu_501_2010042000_t012.nc',
        ]

        hfp = pta.HYCOMFileParser()
        hfp.read(hycom_files)
        print(hfp.velocity)
        assert hfp.velocity.shape == (2, 346 * 541, 5)
        assert hfp.times.shape == (5,)
        assert hfp.times[0] == 90288.0
        assert hfp.times[1] == 90288.0 + 3
        assert hfp.times[2] == 90288.0 + 6
        assert hfp.times[3] == 90288.0 + 9
        assert hfp.times[4] == 90288.0 + 12