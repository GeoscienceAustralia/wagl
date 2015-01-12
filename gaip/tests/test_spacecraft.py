import unittest
from gaip import Spacecraft
import os


class SpacecraftTest(unittest.TestCase):

    def test_type(self):
        sat = Spacecraft.LS5
        self.assertTrue(isinstance(sat, Spacecraft))

    def test_name(self):
        # the 'name' property derives from being an enumeration member
        sat = Spacecraft.LS5
        self.assertTrue(isinstance(sat.name, str))
        self.assertEqual(sat.name, 'LS5')
        self.assertEqual(str(sat), 'Spacecraft.LS5')
        self.assertEqual(repr(sat), "Spacecraft.LS5 (Sensors: ['TM', 'MSS'])")

    def test_properties(self):
        sat = Spacecraft.LS8
        self.assertEqual(sat.omega, 0.001059)
        self.assertEqual(sat.tag, 'LS8')

    def test_iteration(self):
        names = ''
        for sat in Spacecraft:
            names = names + sat.name + ','

        self.assertEqual(names, 'LS5,LS7,LS8,')

    def test_access(self):
        craft_a = Spacecraft.LS7
        craft_b = Spacecraft['LS7']

        self.assertEqual(craft_a, craft_b)

    def no_test_immutable(self):
        """TODO: make Spacecraft members immutable"""
        sat = Spacecraft.LS8

        self.assertEqual(sat.altitude, 705000.0)
        try:
            sat.altitude = 710000.0
            self.fail("Spacecraft is immutable, should have raised a TypeError")
        except TypeError:
            pass  # expected a TypeError trying to change immutable Spacecraft object

if __name__ == '__main__':
    unittest.main()
