#!/usr/bin/env python

import argparse
from collections import OrderedDict
import glob
import json
import os
from os.path import join as pjoin, abspath, basename
import unittest

import markdown
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import rasterio

from gaip import GriddedGeoBox
from gaip import write_img
from gaip.tests.unittesting_tools import ParameterisedTestCase
from idl_functions import histogram

"""
Notes:

What started out at as simple unittesting for the Lambertian-NBAR-TC products
has become more convoluted with the incorporation of auto report generation
for descisions to be made by the management board.
It might be more ideal to break this apart and have unittesting separated
from the auto report generation. JS 2015-07-24.
"""


class TestProductFileNames(ParameterisedTestCase):

    """
    Tests that the reference and the current datasets
    have matching files.
    """

    def test_lambert_files(self):
        """
        Check that the same number of reference and test files
        exist for the lambertian output.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'lambertian' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'lambertian' in f]

        # Change back to the original directory
        os.chdir(cwd)

        self.assertEqual(len(ref_files), len(test_files))


    def test_brdf_files(self):
        """
        Check that the same number of reference and test files
        exist for the brdf output.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'brdf' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'brdf' in f]

        # Change back to the original directory
        os.chdir(cwd)

        self.assertEqual(len(ref_files), len(test_files))


    def test_tc_files(self):
        """
        Check that the same number of reference and test files
        exist for the tc output.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'terrain' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'terrain' in f]

        # Change back to the original directory
        os.chdir(cwd)

        self.assertEqual(len(ref_files), len(test_files))


def test_compare_lmbrt_files(reference_directory, test_directory,
                             output_directory, percent):
    """
    Check that the lambertian outputs are roughly equal.
    """
    reflectance_ref_dir = pjoin(reference_directory, 'Reflectance_Outputs')
    reflectance_test_dir = pjoin(test_directory, 'Reflectance_Outputs')

    cwd = os.getcwd()

    # Get the reference files
    os.chdir(reflectance_ref_dir)
    files = glob.glob('*.bin')
    ref_files = [abspath(f) for f in files if 'lambertian' in f]

    # Get the test files
    os.chdir(reflectance_test_dir)
    files = glob.glob('*.bin')
    test_files = [abspath(f) for f in files if 'lambertian' in f]

    # Change back to the original directory
    os.chdir(cwd)

    ref_files.sort()
    test_files.sort()

    # Precision (Recycle the decmial_precision param)
    tolerance = 100 - percent

    # Setup lists to hold band names, min and max values
    diff_fnames = []
    min_diff = []
    max_diff = []

    # Open the first image and get the geobox
    with rasterio.open(ref_files[0], 'r') as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)

    # We need to calculate the difference images and determine the absolute
    # max and min in order to create a generic histogram that captures all
    # diff images on the same scale. This makes it easier to plot on the
    # same graph
    for i in range(len(ref_files)):
        ref_fname = ref_files[i]
        test_fname = test_files[i]
        fname = basename(test_fname)

        # Get the reference data
        with rasterio.open(ref_fname, 'r') as ds:
            ref_img = ds.read(1, masked=False)

        # Get the test data
        with rasterio.open(test_fname, 'r') as ds:
            test_img = ds.read(1, masked=False)

        # Calculate the difference image and get some basic stats
        diff = (ref_img - test_img).astype('int32')
        max_diff.append(diff.max())
        min_diff.append(diff.min())

        # Output the difference image to disk
        diff_out_fname = pjoin(output_directory,
                               fname.replace('.bin', '_diff.bin'))
        diff_fnames.append(diff_out_fname)
        write_img(diff, diff_out_fname, geobox=geobox)

        # Output the difference map
        out_fname = pjoin(output_directory,
                          fname.replace('.bin', '-Difference-Map.png'))
        plt.imshow(diff)
        title = fname.replace('.bin', '-Difference-Map')
        plt.title(title)
        plt.colorbar()
        plt.savefig(out_fname)
        plt.close()


    # Determine the min and max
    min_v = numpy.min(min_diff)
    max_v = numpy.max(max_diff)


    # Initialise dicts to hold the histogram and cumulative distributions
    hist_results = {}
    cdist_results = {}

    for diff_fname in diff_fnames:
        fname = basename(diff_fname)
        with rasterio.open(diff_fname, 'r') as ds:
            diff = ds.read(1, masked=False)

        h = histogram(diff, minv=min_v, maxv=max_v, omin='omin')
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float32')
        array_sz = diff.size
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]
        hist_results[fname] = hist
        cdist_results[fname] = cdf

        result = OrderedDict()
        result['Product'] = 'Lambertian'
        result['Filename'] = fname
        result['Min Difference'] = float(diff.min())
        result['Max Difference'] = float(diff.max())
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
            result['Acceptable Percent Difference?'] = True
        else:
            result['Acceptable Percent Difference?'] = False

        out_fname = pjoin(output_directory, fname.replace('.bin', '.json'))
        with open(out_fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)


    # Output the histograms
    out_fname = pjoin(output_directory, 'Lambertian-Histogram-Differences.png')
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for bname in hist_results:
        y = hist_results[bname]
        x = numpy.arange(y.shape[0]) + omin
        axis.plot(x, y, label=bname)
    plt.title('Difference Histograms for the Lambertian product')
    plt.xlabel('Pixel Difference Value')
    plt.ylabel('Count')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                     prop={'size': 10})
    plt.tight_layout()
    plt.savefig(out_fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    out_fname = pjoin(output_directory,
                      'Lambertian-Cumulative-Distribution-Differences.png')
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for bname in cdist_results:
        y = cdist_results[bname]
        x = numpy.arange(y.shape[0]) + omin
        axis.plot(x, y, label=bname)
    plt.title('Cumulative Distribution for Lambertian product')
    plt.xlabel('Pixel Difference Value')
    plt.ylabel('Percent %')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                     prop={'size': 10})
    plt.tight_layout()
    plt.savefig(out_fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def test_compare_brdf_files(reference_directory, test_directory,
                            output_directory, percent):
    """
    Check that the brdf outputs are roughly equal.
    """
    reflectance_ref_dir = pjoin(reference_directory, 'Reflectance_Outputs')
    reflectance_test_dir = pjoin(test_directory, 'Reflectance_Outputs')

    cwd = os.getcwd()

    # Get the reference files
    os.chdir(reflectance_ref_dir)
    files = glob.glob('*.bin')
    ref_files = [abspath(f) for f in files if 'brdf' in f]

    # Get the test files
    os.chdir(reflectance_test_dir)
    files = glob.glob('*.bin')
    test_files = [abspath(f) for f in files if 'brdf' in f]

    # Change back to the original directory
    os.chdir(cwd)

    ref_files.sort()
    test_files.sort()

    # Precision (Recycle the decmial_precision param)
    tolerance = 100 - percent

    # Setup lists to hold band names, min and max values
    diff_fnames = []
    min_diff = []
    max_diff = []

    # Open the first image and get the geobox
    with rasterio.open(ref_files[0], 'r') as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)

    # We need to calculate the difference images and determine the absolute
    # max and min in order to create a generic histogram that captures all
    # diff images on the same scale. This makes it easier to plot on the
    # same graph
    for i in range(len(ref_files)):
        ref_fname = ref_files[i]
        test_fname = test_files[i]
        fname = basename(test_fname)

        # Get the reference data
        with rasterio.open(ref_fname, 'r') as ds:
            ref_img = ds.read(1, masked=False)

        # Get the test data
        with rasterio.open(test_fname, 'r') as ds:
            test_img = ds.read(1, masked=False)

        # Calculate the difference image and get some basic stats
        diff = (ref_img - test_img).astype('int32')
        max_diff.append(diff.max())
        min_diff.append(diff.min())

        # Output the difference image to disk
        diff_out_fname = pjoin(output_directory,
                               fname.replace('.bin', '_diff.bin'))
        diff_fnames.append(diff_out_fname)
        write_img(diff, diff_out_fname, geobox=geobox)

        # Output the difference map
        out_fname = pjoin(output_directory,
                          fname.replace('.bin', '-Difference-Map.png'))
        plt.imshow(diff)
        title = fname.replace('.bin', '-Difference-Map')
        plt.title(title)
        plt.colorbar()
        plt.savefig(out_fname)
        plt.close()


    # Determine the min and max
    min_v = numpy.min(min_diff)
    max_v = numpy.max(max_diff)


    # Initialise dicts to hold the histogram and cumulative distributions
    hist_results = {}
    cdist_results = {}

    for diff_fname in diff_fnames:
        fname = basename(diff_fname)
        with rasterio.open(diff_fname, 'r') as ds:
            diff = ds.read(1, masked=False)

        h = histogram(diff, minv=min_v, maxv=max_v, omin='omin')
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float32')
        array_sz = diff.size
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]
        hist_results[fname] = hist
        cdist_results[fname] = cdf

        result = OrderedDict()
        result['Product'] = 'BRDF-Corrected'
        result['Filename'] = fname
        result['Min Difference'] = float(diff.min())
        result['Max Difference'] = float(diff.max())
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
            result['Acceptable Percent Difference?'] = True
        else:
            result['Acceptable Percent Difference?'] = False

        out_fname = pjoin(output_directory, fname.replace('.bin', '.json'))
        with open(out_fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)


    # Output the histograms
    out_fname = pjoin(output_directory,
                      'BRDF-Corrected-Histogram-Differences.png')
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for bname in hist_results:
        y = hist_results[bname]
        x = numpy.arange(y.shape[0]) + omin
        axis.plot(x, y, label=bname)
    plt.title('Difference Histograms for the BRDF-Corrected product')
    plt.xlabel('Pixel Difference Value')
    plt.ylabel('Count')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                     prop={'size': 10})
    plt.tight_layout()
    plt.savefig(out_fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    out_fname = pjoin(output_directory,
                      ('BRDF-Corrected-Cumulative-Distribution-'
                       'Differences.png'))
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for bname in cdist_results:
        y = cdist_results[bname]
        x = numpy.arange(y.shape[0]) + omin
        axis.plot(x, y, label=bname)
    plt.title('Cumulative Distribution for BRDF-Corrected product')
    plt.xlabel('Pixel Difference Value')
    plt.ylabel('Percent %')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                     prop={'size': 10})
    plt.tight_layout()
    plt.savefig(out_fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def test_compare_tc_files(reference_directory, test_directory,
                          output_directory, percent):
    """
    Check that the tc outputs are roughly equal.
    """
    reflectance_ref_dir = pjoin(reference_directory, 'Reflectance_Outputs')
    reflectance_test_dir = pjoin(test_directory, 'Reflectance_Outputs')

    cwd = os.getcwd()

    # Get the reference files
    os.chdir(reflectance_ref_dir)
    files = glob.glob('*.bin')
    ref_files = [abspath(f) for f in files if 'terrain' in f]

    # Get the test files
    os.chdir(reflectance_test_dir)
    files = glob.glob('*.bin')
    test_files = [abspath(f) for f in files if 'terrain' in f]

    # Change back to the original directory
    os.chdir(cwd)

    ref_files.sort()
    test_files.sort()

    # Precision (Recycle the decmial_precision param)
    tolerance = 100 - percent

    # Setup lists to hold band names, min and max values
    diff_fnames = []
    min_diff = []
    max_diff = []

    # Open the first image and get the geobox
    with rasterio.open(ref_files[0], 'r') as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)

    # We need to calculate the difference images and determine the absolute
    # max and min in order to create a generic histogram that captures all
    # diff images on the same scale. This makes it easier to plot on the
    # same graph
    for i in range(len(ref_files)):
        ref_fname = ref_files[i]
        test_fname = test_files[i]
        fname = basename(test_fname)

        # Get the reference data
        with rasterio.open(ref_fname, 'r') as ds:
            ref_img = ds.read(1, masked=False)

        # Get the test data
        with rasterio.open(test_fname, 'r') as ds:
            test_img = ds.read(1, masked=False)

        # Calculate the difference image and get some basic stats
        diff = (ref_img - test_img).astype('int32')
        max_diff.append(diff.max())
        min_diff.append(diff.min())

        # Output the difference image to disk
        diff_out_fname = pjoin(output_directory,
                               fname.replace('.bin', '_diff.bin'))
        diff_fnames.append(diff_out_fname)
        write_img(diff, diff_out_fname, geobox=geobox)

        # Output the difference map
        out_fname = pjoin(output_directory,
                          fname.replace('.bin', '-Difference-Map.png'))
        plt.imshow(diff)
        title = fname.replace('.bin', '-Difference-Map')
        plt.title(title)
        plt.colorbar()
        plt.savefig(out_fname)
        plt.close()


    # Determine the min and max
    min_v = numpy.min(min_diff)
    max_v = numpy.max(max_diff)


    # Initialise dicts to hold the histogram and cumulative distributions
    hist_results = {}
    cdist_results = {}

    for diff_fname in diff_fnames:
        fname = basename(diff_fname)
        with rasterio.open(diff_fname, 'r') as ds:
            diff = ds.read(1, masked=False)

        h = histogram(diff, minv=min_v, maxv=max_v, omin='omin')
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float32')
        array_sz = diff.size
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]
        hist_results[fname] = hist
        cdist_results[fname] = cdf

        result = OrderedDict()
        result['Product'] = 'Terrain-Corrected'
        result['Filename'] = fname
        result['Min Difference'] = float(diff.min())
        result['Max Difference'] = float(diff.max())
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
            result['Acceptable Percent Difference?'] = True
        else:
            result['Acceptable Percent Difference?'] = False

        out_fname = pjoin(output_directory, fname.replace('.bin', '.json'))
        with open(out_fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)


    # Output the histograms
    out_fname = pjoin(output_directory,
                      'Terrain-Corrected-Histogram-Differences.png')
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for bname in hist_results:
        y = hist_results[bname]
        x = numpy.arange(y.shape[0]) + omin
        axis.plot(x, y, label=bname)
    plt.title('Difference Histograms for the Terrain-Corrected product')
    plt.xlabel('Pixel Difference Value')
    plt.ylabel('Count')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                     prop={'size': 10})
    plt.tight_layout()
    plt.savefig(out_fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    out_fname = pjoin(output_directory,
                      ('Terrain-Corrected-Cumulative-'
                       'Distribution-Differences.png'))
    fig = plt.figure()
    axis = fig.add_subplot(111)
    for bname in cdist_results:
        y = cdist_results[bname]
        x = numpy.arange(y.shape[0]) + omin
        axis.plot(x, y, label=bname)
    plt.title('Cumulative Distribution for Terrain-Corrected product')
    plt.xlabel('Pixel Difference Value')
    plt.ylabel('Percent %')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                     prop={'size': 10})
    plt.tight_layout()
    plt.savefig(out_fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def produce_report(out_dir, reference_dir, test_dir):
    """
    Produce the markdown report document.
    """
    template = ("""# NBAR Comparison Report

## Reference Data:
*{ref_dir}*

## Test Data:
*{test_dir}*

## Lambertian Product

{lmbrt_data}

{lmbrt_diff_maps}

![lmbrt differences histogram](Lambertian-Histogram-Differences.png)

![lmbrt cumulative distribution](Lambertian-Cumulative-Distribution-Differences.png)

## BRDF Corrected Product

{brdf_data}

{brdf_diff_maps}

![brdf differences histogram](BRDF-Corrected-Histogram-Differences.png)

![brdf cumulative distribution](BRDF-Corrected-Cumulative-Distribution-Differences.png)

## BRDF & Terrain Corrected Product

{tc_data}

{tc_diff_maps}

![tc differences histogram](Terrain-Corrected-Histogram-Differences.png)

![tc cumulative distribution](Terrain-Corrected-Cumulative-Distribution-Differences.png)
""")

    # Get the JSON files
    cwd = os.getcwd()
    os.chdir(out_dir)
    json_files = glob.glob('*.json')

    # Get the lambertian json data
    lmbrt_json = [f for f in json_files if 'lambertian' in f]
    lmbrt_data = []
    for fname in lmbrt_json:
        with open(fname, 'r') as f:
            data = json.load(f, object_pairs_hook=OrderedDict)
        lmbrt_data.append(json.dumps(data, indent=4, separators=(',\n', ': ')))
    lmbrt_data_list = ['{}\n'.format(d) for d in lmbrt_data]
    lmbrt_data = ''.join(lmbrt_data_list)

    # Get the brdf json data
    brdf_json = [f for f in json_files if 'brdf' in f]
    brdf_data = []
    for fname in brdf_json:
        with open(fname, 'r') as f:
            data = json.load(f, object_pairs_hook=OrderedDict)
        brdf_data.append(json.dumps(data, indent=4, separators=(',\n', ': ')))
    brdf_data_list = ['{}\n'.format(d) for d in brdf_data]
    brdf_data = ''.join(brdf_data_list)

    # Get the terrain json data
    tc_json = [f for f in json_files if 'terrain' in f]
    tc_data = []
    for fname in tc_json:
        with open(fname, 'r') as f:
            data = json.load(f, object_pairs_hook=OrderedDict)
        tc_data.append(json.dumps(data, indent=4, separators=(',\n', ': ')))
    tc_data_list = ['{}\n'.format(d) for d in tc_data]
    tc_data = ''.join(tc_data_list)

    # Get the png difference maps
    png_files = glob.glob('*.png')

    # Get the lambertian difference maps
    lmbrt_diffs = [f for f in png_files if 'lambertian' in f]
    lmbrt_diffs = [f for f in lmbrt_diffs if 'Map' in f]
    lmbrt_diff_maps = []
    for i, m in enumerate(lmbrt_diffs):
        txt = '![lmbrt differences map {bn}]({mp})\n'
        lmbrt_diff_maps.append(txt.format(bn=i, mp=m))
    lmbrt_diffs = ''.join(lmbrt_diff_maps)

    # Get the brdf difference maps
    brdf_diffs = [f for f in png_files if 'brdf' in f]
    brdf_diffs = [f for f in brdf_diffs if 'Map' in f]
    brdf_diff_maps = []
    for i, m in enumerate(brdf_diffs):
        txt = '![brdf differences map {bn}]({mp})\n'
        brdf_diff_maps.append(txt.format(bn=i, mp=m))
    brdf_diffs = ''.join(brdf_diff_maps)

    # Get the terrain difference maps
    tc_diffs = [f for f in png_files if 'terrain' in f]
    tc_diffs = [f for f in tc_diffs if 'Map' in f]
    tc_diff_maps = []
    for i, m in enumerate(tc_diffs):
        txt = '![tc differences map {bn}]({mp})\n'
        tc_diff_maps.append(txt.format(bn=i, mp=m))
    tc_diffs = ''.join(tc_diff_maps)

    # Populate the markdown template
    md_template = template.format(lmbrt_data=lmbrt_data,
                                  lmbrt_diff_maps=lmbrt_diffs,
                                  brdf_data=brdf_data,
                                  brdf_diff_maps=brdf_diffs,
                                  tc_data=tc_data, tc_diff_maps=tc_diffs,
                                  ref_dir=reference_dir,
                                  test_dir=test_dir)

    # Output the nbar markdown
    out_fname = pjoin(out_dir, 'NBAR-Report.md')
    with open(out_fname, 'w') as outf:
        outf.write(md_template)

    # Convert to html and output
    md = markdown.Markdown()
    html = md.convert(md_template)

    out_fname = pjoin(out_dir, 'NBAR-Report.html')
    with open(out_fname, 'w') as outf:
        outf.write(html)

    os.chdir(cwd)


if __name__ == '__main__':

    desc = 'Test and compare the Lambertian, BRDF and TC outputs.'
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--reference_dir', required=True,
                        help='A directory path of a nbar reference output.')
    parser.add_argument('--test_dir', required=True,
                        help=('A directory path to the test output.'))
    parser.add_argument('--pct_tolerance', default=3, type=int,
                        help=('The tolerance for an acceptable percent '
                              'difference.'))
    parser.add_argument('--out_directory', required=True,
                        help='The output directory to contain the results.')

    parsed_args = parser.parse_args()

    reference_dir = parsed_args.reference_dir
    test_dir = parsed_args.test_dir
    tolerance = parsed_args.pct_tolerance
    out_dir = parsed_args.out_directory


    print "\nChecking that we have all the reference and test data files.\n"
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestProductFileNames,
                  reference_dir=reference_dir, test_dir=test_dir))
    unittest.TextTestRunner(verbosity=2).run(suite)


    # Run the difference tests
    print "\nTesting Lambertian\n"
    test_compare_lmbrt_files(reference_dir, test_dir, out_dir, tolerance)

    print "\nTesting BRDF\n"
    test_compare_brdf_files(reference_dir, test_dir, out_dir, tolerance)

    print "\nTesting Terrain\n"
    test_compare_tc_files(reference_dir, test_dir, out_dir, tolerance)

    # Produce the markdown and html docs
    produce_report(out_dir, reference_dir, test_dir)
