import unittest
import PYOPATRA as pta


class ADCIRCFileTests(unittest.TestCase):
    def test_fort_14_reading(self):
        afp = pta.ADCIRCFileParser(fort14='fort.14')
        afp.read_grid_and_bc()

        self.assertEqual(afp.grid_name, 'Pontchartrain Test')
        self.assertEqual(afp.num_elements, 620)
        self.assertEqual(afp.num_vertices, 402)


if __name__ == '__main__':
    unittest.main()
