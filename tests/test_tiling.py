#!/usr/bin/env python

"""
Tests the tiling routines.
"""

from __future__ import absolute_import
import random
import unittest

import numpy

from wagl.tiling import generate_tiles


class TestGetTile3(unittest.TestCase):

    """Unit tests for generate_tiles function."""

    # pylint: disable=too-many-public-methods

    # Input data for tests. These are tuples (samples, lines, xtile, ytile)
    # used as arguments to the generate_tiles function during testing.

    normal_input = [(4000, 4000, 100, 100), (4001, 4001, 100, 100),
                    (4001, 4003, 97, 101)]

    boundary_input1 = [(4000, 1, 97, 101), (4000, 2, 101, 97),
                       (1, 4000, 100, 100), (2, 4000, 100, 100),
                       (1, 1, 97, 100), (2, 2, 100, 100),
                       (1, 2, 100, 100), (2, 1, 100, 100)]

    boundary_input2 = [(37, 41, 10, 1), (37, 41, 10, 2),
                       (40, 40, 1, 10), (40, 40, 2, 10),
                       (40, 40, 1, 1), (40, 40, 2, 2),
                       (40, 40, 1, 2), (40, 40, 1, 2)]

    exception_input1 = [(0, 4000, 97, 101), (4000, 0, 101, 97),
                        (-1, 4001, 100, 100), (4001, -1, 101, 97)]

    exception_input2 = [(4001, 4003, 0, 100), (4003, 4001, 100, 0)]
    # TODO have a test for negative tile sizes
    # exception_input2 = [(4001, 4003, 0, 100), (4003, 4001, 100, 0),
    #                     (4000, 4000, -1, 100), (4000, 4000, 100, -1)]

    def test_normal(self):
        """Test typical input:"""
        self.do_test(self.normal_input)

    def test_boundary_1(self):
        """Test small samples or lines:"""
        self.do_test(self.boundary_input1)

    def test_boundary_2(self):
        """Test small xtile or ytile:"""
        self.do_test(self.boundary_input2)

    def test_exception_1(self):
        """Test empty image:"""
        for (samples, lines, xtile, ytile) in self.exception_input1:
            tiles_list = list(generate_tiles(samples, lines, xtile, ytile))
            self.assertEqual(tiles_list, [], 'Expected an empty tile list.')

    def test_exception_2(self):
        """Test empty tiles:"""
        for (samples, lines, xtile, ytile) in self.exception_input2:
            self.assertRaises(ZeroDivisionError, generate_tiles,
                              samples, lines, xtile, ytile)

    def do_test(self, test_input):
        """Check sizes and coverage for a list of test input."""
        for (samples, lines, xtile, ytile) in test_input:
            tiles_list = list(generate_tiles(samples, lines, xtile, ytile))

            self.check_sizes(xtile, ytile, tiles_list)
            self.check_tiling(samples, lines, tiles_list)

    def check_sizes(self, xtile, ytile, tiles_list):
        """Check that the tiles in tiles_list have appropriate extents."""
        for tile in tiles_list:
            yse, xse = tile
            ystart, yend = yse
            xstart, xend = xse
            self.assertTrue(0 <= xstart < xend,
                            'Tile empty - xcoord: ' + repr(tile))
            self.assertTrue(0 <= ystart < yend,
                            'Tile empty - y coord: ' + repr(tile))
            self.assertLessEqual(xend - xstart, xtile,
                                 'Tile too big - x coord: ' + repr(tile))
            self.assertLessEqual(yend - ystart, ytile,
                                 'Tile too big - y coord: ' + repr(tile))

    def check_tiling(self, samples, lines, tiles_list):
        """Check the tiles in tiles_list for covarge and overlap."""

        # Each tile is given a random tag (between 1 and a maximum).
        # The test array is initialized to all zero then the tags are added
        # over each tiles extent. After this holes in the coverage should
        # have a value of zero, and overlaps should have a value that does
        # not corrispond to a vaild tag. The numpy unique function allows
        # this checking to be done quickly.

        (tag_array, tag_dict) = generate_tile_tags(tiles_list)

        test_array = generate_test_array(samples, lines, tiles_list, tag_array)

        # pylint: disable=unbalanced-tuple-unpacking
        (unique, unique_indices) = numpy.unique(test_array, return_index=True)
        # numpy.unique has a variable number of return values depending on
        # keyword flags. This confuses pylint.
        # pylint: enable=unbalanced-tuple-unpacking

        for (tag, flat_index) in zip(unique, unique_indices):
            index = numpy.unravel_index(flat_index, (lines, samples))
            self.assertGreater(tag, 0,
                               'Hole in coverage detected at ' + repr(index))
            self.assertIn(tag, tag_dict,
                          'Tile overlap detected at ' + repr(index))


def generate_tile_tags(tiles_list):
    """Generate a random tag for each tile.

    Returns a numpy array of the tags an a dictonary of tag->tile.

    """
    max_tag = 2**31 - 1
    tag_list = []
    tag_dict = {}
    for tile in tiles_list:
        tile_tag = random.randint(1, max_tag)
        tag_list.append(tile_tag)
        tag_dict[tile_tag] = tile
    tag_array = numpy.array(tag_list, dtype=numpy.uint32)
    return tag_array, tag_dict


def generate_test_array(samples, lines, tiles_list, tag_array):
    """Generate a test array of the same size as the arrays to be tiled.

    This array is initialized to all zero, then the tags are added
    to each tile's extent.

    """
    test_array = numpy.zeros((lines, samples), dtype=numpy.uint32)
    i = 0
    for tile in tiles_list:
        yse, xse = tile
        ystart, yend = yse
        xstart, xend = xse
        test_array[ystart:yend, xstart:xend] += tag_array[i]
        i += 1
    return test_array


def the_suite():
    """Returns a test suite of all the tests in this module."""

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestGetTile3)

    return suite


def run_the_tests():
    """Runs the tests defined in this module"""
    unittest.TextTestRunner(verbosity=2).run(the_suite())

if __name__ == '__main__':
    run_the_tests()
