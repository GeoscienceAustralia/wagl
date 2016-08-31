#!/usr/bin/env python

import argparse
from collections import OrderedDict
import datetime
import glob
import json
import os
from os.path import join as pjoin, abspath, dirname
import unittest

import markdown
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import numpy.testing as npt
import rasterio

from gaip.tests.unittesting_tools import ParameterisedTestCaseFiles
from idl_functions import histogram


class TestFCProducts(ParameterisedTestCaseFiles):

    """
    Test and compare the FC product outputs.
    """

    def test_photosynthetic_veg_band(self):
        """
        Test and compare the differences between the
        reference and test photosynthetic vegetation
        output bands.
        """
        out_dir = self.output_directory

        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read(2, masked=False).astype('float')

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read(2, masked=False)

        # Precision
        tolerance = 100 - self.tolerance

        # Calculate the difference and get basic stats
        diff = (ref_img - test_img).astype('int')
        max_diff = diff.max()
        min_diff = diff.min()

        h = histogram(diff, minv=min_diff, maxv=max_diff, omin='omin')
        array_sz = ref_img.size
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float')
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]

        # Initialise a dict ready for output to JSON
        result = OrderedDict()
        result['Fractional Cover Index'] = 'Photsynthetic Vegetation'
        result['Min Difference'] = float(min_diff)
        result['Max Difference'] = float(max_diff)
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
           result['Acceptable Percent Difference?'] = True
        else:
           result['Acceptable Percent Difference?'] = False

        fname = pjoin(out_dir, 'Photosynthetic-Vegetation.json')
        with open(fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)

        # Output figures (difference map, histograms)
        fname = pjoin(out_dir, 'PV-Difference-Map.png')
        plt.imshow(diff)
        plt.title('Photosynthetic Vegetation Difference Map')
        plt.colorbar()
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'PV-Histogram-Differences.png')
        plt.plot(hist)
        plt.title('Histogram of Differences for Photosynthetic Vegetation')
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Count')
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'PV-Cumulative-Distribution-Differences.png')
        plt.plot(cdf)
        plt.title('Cumulative Distribution for Photosynthetic Vegetation')
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Percent %')
        plt.savefig(fname)
        plt.close()

        self.assertTrue(pct_no_diff > tolerance)


    def test_non_photosynthetic_veg_band(self):
        """
        Test and compare the differences between the
        reference and test non-photosynthetic vegetation
        output bands.
        """
        out_dir = self.output_directory

        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read(3, masked=False).astype('float')

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read(3, masked=False)

        # Precision
        tolerance = 100 - self.tolerance

        # Calculate the difference and get basic stats
        diff = (ref_img - test_img).astype('int')
        max_diff = diff.max()
        min_diff = diff.min()

        h = histogram(diff, minv=min_diff, maxv=max_diff, omin='omin')
        array_sz = ref_img.size
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float')
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]

        # Initialise a dict ready for output to JSON
        result = OrderedDict()
        result['Fractional Cover Index'] = 'Non Photsynthetic Vegetation'
        result['Min Difference'] = float(min_diff)
        result['Max Difference'] = float(max_diff)
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
           result['Acceptable Percent Difference?'] = True
        else:
           result['Acceptable Percent Difference?'] = False

        fname = pjoin(out_dir, 'Non-Photosynthetic-Vegetation.json')
        with open(fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)

        # Output figures (difference map, histograms)
        fname = pjoin(out_dir, 'NPV-Difference-Map.png')
        plt.imshow(diff)
        plt.title('Non Photosynthetic Vegetation Difference Map')
        plt.colorbar()
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'NPV-Histogram-Differences.png')
        plt.plot(hist)
        plt.title('Histogram of Differences for Non Photosynthetic Vegetation')
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Count')
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'NPV-Cumulative-Distribution-Differences.png')
        plt.plot(cdf)
        title = ('Cumulative Distribution of Differences for Non '
                 'Photsynthetic Vegetation')
        plt.title(title)
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Percent %')
        plt.savefig(fname)
        plt.close()

        self.assertTrue(pct_no_diff > tolerance)


    def test_bare_soil_band(self):
        """
        Test and compare the differences between the
        reference and test bare soil output bands.
        """
        out_dir = self.output_directory

        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read(1, masked=False).astype('float')

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read(1, masked=False)

        # Precision
        tolerance = 100 - self.tolerance

        # Calculate the difference and get basic stats
        diff = (ref_img - test_img).astype('int')
        max_diff = diff.max()
        min_diff = diff.min()

        h = histogram(diff, minv=min_diff, maxv=max_diff, omin='omin')
        array_sz = ref_img.size
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float')
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]

        # Initialise a dict ready for output to JSON
        result = OrderedDict()
        result['Fractional Cover Index'] = 'Bare Soil'
        result['Min Difference'] = float(min_diff)
        result['Max Difference'] = float(max_diff)
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
           result['Acceptable Percent Difference?'] = True
        else:
           result['Acceptable Percent Difference?'] = False

        fname = pjoin(out_dir, 'Bare-Soil.json')
        with open(fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)

        # Output figures (difference map, histograms)
        fname = pjoin(out_dir, 'BS-Difference-Map.png')
        plt.imshow(diff)
        plt.title('Difference Map')
        plt.colorbar()
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'BS-Histogram-Differences.png')
        plt.plot(hist)
        plt.title('Histogram of Differences for Bare Soil')
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Count')
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'BS-Cumulative-Distribution-Differences.png')
        plt.plot(cdf)
        plt.title('Cumulative Distribution of Differences for Bare Soil')
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Percent %')
        plt.savefig(fname)
        plt.close()

        self.assertTrue(pct_no_diff > tolerance)


    def test_unmixing_error_band(self):
        """
        Test and compare the differences between the
        reference and test unmixing error output bands.
        """
        out_dir = self.output_directory

        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read(4, masked=False).astype('float')

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read(4, masked=False)

        # Precision
        tolerance = 100 - self.tolerance

        # Calculate the difference and get basic stats
        diff = (ref_img - test_img).astype('int')
        max_diff = diff.max()
        min_diff = diff.min()

        h = histogram(diff, minv=min_diff, maxv=max_diff, omin='omin')
        array_sz = ref_img.size
        hist = h['histogram']
        omin = h['omin']
        cumu_h = numpy.cumsum(hist, dtype='float')
        cdf = (cumu_h / array_sz) * 100

        pct_no_diff = cdf[0 - omin]

        # Initialise a dict ready for output to JSON
        result = OrderedDict()
        result['Fractional Cover Index'] = 'Unmixing Error'
        result['Min Difference'] = float(min_diff)
        result['Max Difference'] = float(max_diff)
        result['Minimum Allowable Percent Difference'] = float(tolerance)
        result['Percent No Difference'] = float(pct_no_diff)

        if pct_no_diff > tolerance:
           result['Acceptable Percent Difference?'] = True
        else:
           result['Acceptable Percent Difference?'] = False

        fname = pjoin(out_dir, 'Umixing-Error.json')
        with open(fname, 'w') as outf:
            json.dump(result, outf, ensure_ascii=False, indent=4)

        # Output figures (difference map, histograms)
        fname = pjoin(out_dir, 'UE-Difference-Map.png')
        plt.imshow(diff)
        plt.title('Difference Map')
        plt.colorbar()
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'UE-Histogram-Differences.png')
        plt.plot(hist)
        plt.title('Histogram of Differences for Unmixing Error')
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Count')
        plt.savefig(fname)
        plt.close()

        fname = pjoin(out_dir, 'UE-Cumulative-Distribution-Differences.png')
        plt.plot(cdf)
        plt.xlabel('Pixel Difference Value')
        plt.ylabel('Percent %')
        plt.title('Cumulative Distribution of Differences for Unmixing Error')
        plt.savefig(fname)
        plt.close()

        self.assertTrue(pct_no_diff > tolerance)


def produce_report(out_dir):
    """
    Produce the markdown report document.
    """
    template = ("""# Fractional Cover Comparison Report

## Photosynthetic Vegetation

{pv_data}

![pv difference map](PV-Difference-Map.png)

![pv differences histogram](PV-Histogram-Differences.png)

![pv cumulative distribution](PV-Cumulative-Distribution-Differences.png)

## Non-Photosynthetic Vegetation

{npv_data}

![npv difference map](NPV-Difference-Map.png)

![npv differences histogram](NPV-Histogram-Differences.png)

![npv cumulative distribution](NPV-Cumulative-Distribution-Differences.png)

## Bare Soil

{bs_data}

![bs difference map](BS-Difference-Map.png)

![bs differences histogram](BS-Histogram-Differences.png)

![bs cumulative distribution](BS-Cumulative-Distribution-Differences.png)

## Unmixing Error

{ue_data}

![ue difference map](UE-Difference-Map.png)

![ue differences histogram](UE-Histogram-Differences.png)

![ue cumulative distribution](UE-Cumulative-Distribution-Differences.png)
""")

    # Get the PV json data
    fname = pjoin(out_dir, 'Photosynthetic-Vegetation.json')
    with open(fname, 'r') as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
    pv_data = json.dumps(data, indent=4, separators=(',\n', ': '))

    # Get the NPV json data
    fname = pjoin(out_dir, 'Non-Photosynthetic-Vegetation.json')
    with open(fname, 'r') as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
    npv_data = json.dumps(data, indent=4, separators=(',\n', ': '))

    # Get the BS json data
    fname = pjoin(out_dir, 'Bare-Soil.json')
    with open(fname, 'r') as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
    bs_data = json.dumps(data, indent=4, separators=(',\n', ': '))

    # Get the UE json data
    fname = pjoin(out_dir, 'Umixing-Error.json')
    with open(fname, 'r') as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
    ue_data = json.dumps(data, indent=4, separators=(',\n', ': '))

    md_template = template.format(pv_data=pv_data, npv_data=npv_data,
                                  bs_data=bs_data, ue_data=ue_data)

    out_fname = pjoin(out_dir, 'Fractional-Cover-Report.md')
    with open(out_fname, 'w') as outf:
        outf.write(md_template)

    md = markdown.Markdown()
    html = md.convert(md_template)

    out_fname = pjoin(out_dir, 'Fractional-Cover-Report.html')
    with open(out_fname, 'w') as outf:
        outf.write(html)


if __name__ == '__main__':

    desc = 'Test and compare the fractional cover outputs.'
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--reference_filename', required=True,
                        help='The filename for the referenc dataset.')
    parser.add_argument('--test_filename', required=True,
                        help='The filename for the test dataset.')
    parser.add_argument('--pct_tolerance', default=3, type=int,
                        help=('The tolerance for an acceptable percent '
                              'difference.'))
    parser.add_argument('--out_directory', required=True,
                        help='The output directory to contain the results.')

    parsed_args = parser.parse_args()

    reference_fname = parsed_args.reference_filename
    test_fname = parsed_args.test_filename
    tolerance = parsed_args.pct_tolerance
    out_dir = parsed_args.out_directory


    print "Testing the Fractional cover products."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCaseFiles.parameterise(TestFCProducts,
                  reference_fname=reference_fname, test_fname=test_fname,
                  tolerance=tolerance, output_directory=out_dir))
    unittest.TextTestRunner(verbosity=2).run(suite)

    # Produce the markdown and html docs
    produce_report(out_dir)
