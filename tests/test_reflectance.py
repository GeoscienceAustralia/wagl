import unittest

from wagl.reflectance import sun_glint, sun_glint_correction


class TestReflectance(unittest.TestCase):

    def test_sunglint(self):
        self.assertAlmostEqual(sun_glint(5.0, 30.0, 200.0, 5.0), 0.03241, 4)

    def test_sunglintcorrection(self):
        self.assertAlmostEqual(sun_glint_correction(0.03, 0.8, 0.03), 0.006, 4)


if __name__ == '__main__':
    unittest.main()
