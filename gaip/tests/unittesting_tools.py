#!/usr/bin/env python

import os
import sys
import unittest

from osgeo import gdal

class ParameterisedTestCase(unittest.TestCase):
    """ 
    TestCase classes that want to be parameterised should
    inherit from this class.

    Source for this code taken from:
    http://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases/

    Modified to suit our given parameters.
    """
    def __init__(self, methodName='runTest', reference_dir=None, test_dir=None,
                 decimal_precision=4, integer_precision=1):
        super(ParameterisedTestCase, self).__init__(methodName)

        # Reference and Testing directories
        self.reference_dir = reference_dir
        self.test_dir      = test_dir

        # Allowable numerical precision
        self.decimal_precision = decimal_precision
        self.integer_precision = integer_precision


    @staticmethod
    def parameterise(testcase_klass, reference_dir=None, test_dir=None,
                    decimal_precision=4, integer_precision=1):
        """
        Create a suite containing all tests taken from the given
        subclass, passing them the parameters 'reference_dir,
        test_dir, decimal_precision, integer_precision'.

        :param testcase_klass:
            A unittest.TestCase Class

        :param reference_dir:
            A full file pathname to the directory containing the
            reference data.

        :param test_dir:
            A full file pathname to the directory containing the
            test data.

        :param decimal_precision:
            The decimal precision to be used during array comparison.
            Default is 4, i.e. values must be correct up to 4 d.p. in
            order to Pass.

        :param integer_precision:
            The intger precision to be used during array comparison.
            Default is 1, i.e. values must be correct within 1 integer
            in order to Pass.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, reference_dir=reference_dir,
                          test_dir=test_dir,
                          decimal_precision=decimal_precision,
                          integer_precision=integer_precision))
        return suite

