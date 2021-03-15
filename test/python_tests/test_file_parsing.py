import os
import PYOPATRA as pta

class TestAdcircParsing:
    def test_read_grid_and_bc(self):
        afp = pta.ADCIRCFileParser(fort14=os.path.dirname(os.path.realpath(__file__)) + '/fort.14')
        afp.read_grid_and_bc()

        assert afp.grid_name == 'Pontchartrain Test'
        assert afp.num_elements == 620
        assert afp.num_vertices == 402

    # def test_read_density_temperature_salinity(self):
    #     assert False
    #
    # def test_read_velocity(self):
    #     assert False
    #
    # def test_read_turbulence(self):
    #     assert False
