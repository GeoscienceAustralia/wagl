#!/usr/bin/env python

# coding=utf-8
# A near literal port of the FMask matlab code to python using numpy/scipy and itk.
# - As a result of the literal port, original comments should still be intact - as is the original code structure.

# Transcription History
# 1.6.3 transcription Mitchell Wheeler
# 3.0 transcription Josh Sixsmith
#     + added skimage library
#     + better handling for imagery of different resolutions

import sys
import re
import gc
import time
import math
import logging
import datetime
import os.path
import argparse
import numpy
import numexpr
import scipy.stats
import scipy.signal
import scipy.ndimage.morphology
from osgeo import gdal
from skimage import morphology
from skimage import measure
from skimage import segmentation


# Sun earth distance look up table
sun_earth_distance = {1: 0.98331, 2: 0.98330, 3: 0.98330, 4: 0.98330, 5: 0.98330, 6: 0.98332, 7: 0.98333, 8: 0.98335, 9: 0.98338, 10: 0.98341, 11: 0.98345, 12: 0.98349, 13: 0.98354, 14: 0.98359, 15: 0.98365, 16: 0.98371, 17: 0.98378, 18: 0.98385, 19: 0.98393, 20: 0.98401, 21: 0.98410, 22: 0.98419, 23: 0.98428, 24: 0.98439, 25: 0.98449, 26: 0.98460, 27: 0.98472, 28: 0.98484, 29: 0.98496, 30: 0.98509, 31: 0.98523, 32: 0.98536, 33: 0.98551, 34: 0.98565, 35: 0.98580, 36: 0.98596, 37: 0.98612, 38: 0.98628, 39: 0.98645, 40: 0.98662, 41: 0.98680, 42: 0.98698, 43: 0.98717, 44: 0.98735, 45: 0.98755, 46: 0.98774, 47: 0.98794, 48: 0.98814, 49: 0.98835, 50: 0.98856, 51: 0.98877, 52: 0.98899, 53: 0.98921, 54: 0.98944, 55: 0.98966, 56: 0.98989, 57: 0.99012, 58: 0.99036, 59: 0.99060, 60: 0.99084, 61: 0.99108, 62: 0.99133, 63: 0.99158, 64: 0.99183, 65: 0.99208, 66: 0.99234, 67: 0.99260, 68: 0.99286, 69: 0.99312, 70: 0.99339, 71: 0.99365, 72: 0.99392, 73: 0.99419, 74: 0.99446, 75: 0.99474, 76: 0.99501, 77: 0.99529, 78: 0.99556, 79: 0.99584, 80: 0.99612, 81: 0.99640, 82: 0.99669, 83: 0.99697, 84: 0.99725, 85: 0.99754, 86: 0.99782, 87: 0.99811, 88: 0.99840, 89: 0.99868, 90: 0.99897, 91: 0.99926, 92: 0.99954, 93: 0.99983, 94: 1.00012, 95: 1.00041, 96: 1.00069, 97: 1.00098, 98: 1.00127, 99: 1.00155, 100: 1.00184, 101: 1.00212, 102: 1.00240, 103: 1.00269, 104: 1.00297, 105: 1.00325, 106: 1.00353, 107: 1.00381, 108: 1.00409, 109: 1.00437, 110: 1.00464, 111: 1.00492, 112: 1.00519, 113: 1.00546, 114: 1.00573, 115: 1.00600, 116: 1.00626, 117: 1.00653, 118: 1.00679, 119: 1.00705, 120: 1.00731, 121: 1.00756, 122: 1.00781, 123: 1.00806, 124: 1.00831, 125: 1.00856, 126: 1.00880, 127: 1.00904, 128: 1.00928, 129: 1.00952, 130: 1.00975, 131: 1.00998, 132: 1.01020, 133: 1.01043, 134: 1.01065, 135: 1.01087, 136: 1.01108, 137: 1.01129, 138: 1.01150, 139: 1.01170, 140: 1.01191, 141: 1.01210, 142: 1.01230, 143: 1.01249, 144: 1.01267, 145: 1.01286, 146: 1.01304, 147: 1.01321, 148: 1.01338, 149: 1.01355, 150: 1.01371, 151: 1.01387, 152: 1.01403, 153: 1.01418, 154: 1.01433, 155: 1.01447, 156: 1.01461, 157: 1.01475, 158: 1.01488, 159: 1.01500, 160: 1.01513, 161: 1.01524, 162: 1.01536, 163: 1.01547, 164: 1.01557, 165: 1.01567, 166: 1.01577, 167: 1.01586, 168: 1.01595, 169: 1.01603, 170: 1.01610, 171: 1.01618, 172: 1.01625, 173: 1.01631, 174: 1.01637, 175: 1.01642, 176: 1.01647, 177: 1.01652, 178: 1.01656, 179: 1.01659, 180: 1.01662, 181: 1.01665, 182: 1.01667, 183: 1.01668, 184: 1.01670, 185: 1.01670, 186: 1.01670,
                      187: 1.01670, 188: 1.01669, 189: 1.01668, 190: 1.01666, 191: 1.01664, 192: 1.01661, 193: 1.01658, 194: 1.01655, 195: 1.01650, 196: 1.01646, 197: 1.01641, 198: 1.01635, 199: 1.01629, 200: 1.01623, 201: 1.01616, 202: 1.01609, 203: 1.01601, 204: 1.01592, 205: 1.01584, 206: 1.01575, 207: 1.01565, 208: 1.01555, 209: 1.01544, 210: 1.01533, 211: 1.01522, 212: 1.01510, 213: 1.01497, 214: 1.01485, 215: 1.01471, 216: 1.01458, 217: 1.01444, 218: 1.01429, 219: 1.01414, 220: 1.01399, 221: 1.01383, 222: 1.01367, 223: 1.01351, 224: 1.01334, 225: 1.01317, 226: 1.01299, 227: 1.01281, 228: 1.01263, 229: 1.01244, 230: 1.01225, 231: 1.01205, 232: 1.01186, 233: 1.01165, 234: 1.01145, 235: 1.01124, 236: 1.01103, 237: 1.01081, 238: 1.01060, 239: 1.01037, 240: 1.01015, 241: 1.00992, 242: 1.00969, 243: 1.00946, 244: 1.00922, 245: 1.00898, 246: 1.00874, 247: 1.00850, 248: 1.00825, 249: 1.00800, 250: 1.00775, 251: 1.00750, 252: 1.00724, 253: 1.00698, 254: 1.00672, 255: 1.00646, 256: 1.00620, 257: 1.00593, 258: 1.00566, 259: 1.00539, 260: 1.00512, 261: 1.00485, 262: 1.00457, 263: 1.00430, 264: 1.00402, 265: 1.00374, 266: 1.00346, 267: 1.00318, 268: 1.00290, 269: 1.00262, 270: 1.00234, 271: 1.00205, 272: 1.00177, 273: 1.00148, 274: 1.00119, 275: 1.00091, 276: 1.00062, 277: 1.00033, 278: 1.00005, 279: 0.99976, 280: 0.99947, 281: 0.99918, 282: 0.99890, 283: 0.99861, 284: 0.99832, 285: 0.99804, 286: 0.99775, 287: 0.99747, 288: 0.99718, 289: 0.99690, 290: 0.99662, 291: 0.99634, 292: 0.99605, 293: 0.99577, 294: 0.99550, 295: 0.99522, 296: 0.99494, 297: 0.99467, 298: 0.99440, 299: 0.99412, 300: 0.99385, 301: 0.99359, 302: 0.99332, 303: 0.99306, 304: 0.99279, 305: 0.99253, 306: 0.99228, 307: 0.99202, 308: 0.99177, 309: 0.99152, 310: 0.99127, 311: 0.99102, 312: 0.99078, 313: 0.99054, 314: 0.99030, 315: 0.99007, 316: 0.98983, 317: 0.98961, 318: 0.98938, 319: 0.98916, 320: 0.98894, 321: 0.98872, 322: 0.98851, 323: 0.98830, 324: 0.98809, 325: 0.98789, 326: 0.98769, 327: 0.98750, 328: 0.98731, 329: 0.98712, 330: 0.98694, 331: 0.98676, 332: 0.98658, 333: 0.98641, 334: 0.98624, 335: 0.98608, 336: 0.98592, 337: 0.98577, 338: 0.98562, 339: 0.98547, 340: 0.98533, 341: 0.98519, 342: 0.98506, 343: 0.98493, 344: 0.98481, 345: 0.98469, 346: 0.98457, 347: 0.98446, 348: 0.98436, 349: 0.98426, 350: 0.98416, 351: 0.98407, 352: 0.98399, 353: 0.98391, 354: 0.98383, 355: 0.98376, 356: 0.98370, 357: 0.98363, 358: 0.98358, 359: 0.98353, 360: 0.98348, 361: 0.98344, 362: 0.98340, 363: 0.98337, 364: 0.98335, 365: 0.98333, 366: 0.98331}

# Replacement for original dir() function in this module.
# Renamed to avoid name collision with builtin.


def match_file(dir_path, pattern):

    # Ignore filenames starting with '.' character, e.g. rsync work files.

    for f in sorted([x for x in os.listdir(dir_path) if not x.startswith('.')]):
        if re.search(pattern, f):
            return os.path.join(dir_path, f)

    # No match -- complain loudly (or bail out??)

    logging.error('ERROR: %s.match_file("%s", "%s") found no match',
                  __file__.strip('.py'), dir_path, pattern
                  )
    return None


def im_info(filename):
    """
    A function to retrieve the geotransform and projection details using GDAL.
    The original MATLAB code handles it differently, we'll just implement something unique here
    and let GDAL automatically handle the read/write of the projection info.
    Less messy than the MATLAB code which assumes that everything is in UTM.
    Some products may come as Lat/Lon geographic projections.
    """

    img = gdal.Open(filename)
    geoT = img.GetGeoTransform()
    prj = img.GetProjection()
    return (geoT, prj)


def imread(filename, resample=False, samples=None, lines=None):
    img = gdal.Open(filename)
    band = img.GetRasterBand(1)
    if resample:
        driver = gdal.GetDriverByName('MEM')
        outds = driver.Create("", samples, lines, 1, band.DataType)
        outds.SetGeoTransform(img.GetGeoTransform())
        outds.SetProjection(img.GetProjection())
        gdal.ReprojectImage(img, outds)
        return outds.ReadAsArray()
    else:
        return band.ReadAsArray()


def imfill_pybuffer(img, ts):
    import itk

    # Convert img to ITK image
    image_type = itk.Image[itk.F, 2]

    # Have to transfer via a file for now... (until PyBuffer WrapITK extension
    # is fixed)
    inp = itk.PyBuffer[image_type].GetImageFromArray(img)

    # Run grayscale connected closing filter
    #filter = itk.GrayscaleConnectedClosingImageFilter[image_type, image_type].New()
    fltr = itk.GrayscaleFillholeImageFilter[image_type, image_type].New()
    fltr.SetInput(inp)
    fltr.Update()

    # Copy results from ITK image back to numpy array
    output_array = itk.PyBuffer[image_type].GetArrayFromImage(fltr.GetOutput())

    return output_array


def imfill(img, ts):
    import itk

    # Convert img to ITK image
    image_type = itk.Image[itk.SS, 2]

    # Have to transfer via a file for now... (until PyBuffer WrapITK extension
    # is fixed)
    logging.debug("serializing")
    scale_min = img.min()
    scale_max = img.max()
    c = gdal.GetDriverByName('GTiff').Create(
        "img_%s.tif" % ts, img.shape[1], img.shape[0], 1, gdal.GDT_UInt16)
    c.GetRasterBand(1).WriteArray(
        (((img - img.min()) / (img.max() - img.min())) * 32767.0).astype('int16'))
    c = None
    logging.debug("deserializing")
    c = itk.ImageFileReader[image_type].New()
    c.SetFileName("img_%s.tif" % ts)
    inp = c.GetOutput()
    #inp = itk.PyBuffer[image_type].GetImageFromArray(img)

    # Run grayscale connected closing filter
    logging.debug("processing")
    #fltr = itk.GrayscaleConnectedClosingImageFilter[image_type, image_type].New()
    fltr = itk.GrayscaleFillholeImageFilter[image_type, image_type].New()
    fltr.SetInput(inp)
    fltr.Update()

    # Save debug image (using ITK)
    logging.debug("serializing")
    test = itk.ImageFileWriter[image_type].New()
    test.SetFileName("holes_%s.tif" % ts)
    test.SetInput(fltr.GetOutput())
    test.Update()
    test = None

    logging.debug("deserializing")
    output_array = imread("holes_%s.tif" % ts)
    output_array = scale_min + \
        (output_array.astype('float32') / 32767.0) * (scale_max - scale_min)
    # Copy results from ITK image back to numpy array
    #output_array = itk.PyBuffer[image_type].GetArrayFromImage(fltr.GetOutput())

    os.remove("img_%s.tif" % ts)
    os.remove("holes_%s.tif" % ts)

    return output_array


def imfill_skimage(img):
    """
    Replicates the imfill function available within MATLAB.
    Based on the example provided in http://scikit-image.org/docs/dev/auto_examples/plot_holes_and_peaks.html#example-plot-holes-and-peaks-py.

    """

    seed = img.copy()

    # Define seed points and the start points for the erosion process.
    seed[1:-1, 1:-1] = img.max()

    # Define the mask; Probably unneeded.
    mask = img

    # Fill the holes
    filled = morphology.reconstruction(seed, mask, method='erosion')

    return filled


def lndhdrread(filename):
    """
    Load Landsat scene MTL file metadata.

    :param filename:
        A string containing the full path the scene's MTL file.
    """
    # Read in Landsat TM/ETM+ MTL header for Fmask
    # [Lmax,Lmin,Qcalmax,Qcalmin,ijdim_ref,ijdim_thm,reso_ref,reso_thm,ul,zen,azi,zc,Lnum,doy]=lndhdrread(filename)
    # Where:
    # Inputs:
    # filename='L*MTL.txt'
    # Ouenuts:
    # 1) Lmax = Max radiances
    # 2) Lmin = Min radiances
    # 3) Qcalmax = Max calibrated DNs
    # 4) Qcalmin = Min calibrated DNs
    # 5) ijdim_ref = [nrows,ncols] # dimension of optical bands
    # 6) ijdim_ref = [nrows,ncols] # dimension of thermal band
    # 7) reo_ref = 28/30 # resolution of optical bands
    # 8) reo_thm = 60/120 # resolution of thermal band
    # 9) ul = [upperleft_mapx upperleft_mapy]
    # 10) zen = solar zenith angle (degrees)
    # 11) azi = solar azimuth angle (degrees)
    # 12) zc = Zone Number
    # 13) Lnum = 4,5,or 7 Landsat sensor number
    # 14) doy = day of year (1,2,3,...,356)
    #
    ##
    # open and read hdr file
    data = {}
    fl = open(filename, 'r')
    file_lines = fl.readlines()
    for line in file_lines:
        values = line.split(' = ')
        if len(values) != 2:
            continue

        data[values[0].strip()] = values[1].strip().strip('"')

    fl.close()

    # Identify Landsat Number (Lnum = 4, 5 or 7)
    LID = data['SPACECRAFT_ID']
    Lnum = int(LID[len(LID) - 1])

    if ((Lnum >= 4) & (Lnum <= 7)):
        # LS8 variables only. The original MATLAB function returns all
        # variables.
        Refmax = None
        # LS8 variables only. The original MATLAB function returns all
        # variables.
        Refmin = None
        # Test for New/Old MTL file
        if not ('LANDSAT_SCENE_ID' in data.keys()):
            # Process Old version of MTL file

            # read in LMAX
            Lmax_B1 = numpy.float32(data['LMAX_BAND1'])
            Lmax_B2 = numpy.float32(data['LMAX_BAND2'])
            Lmax_B3 = numpy.float32(data['LMAX_BAND3'])
            Lmax_B4 = numpy.float32(data['LMAX_BAND4'])
            Lmax_B5 = numpy.float32(data['LMAX_BAND5'])
            if Lnum == 7:
                Lmax_B6 = numpy.float32(data['LMAX_BAND61'])
            else:
                Lmax_B6 = numpy.float32(data['LMAX_BAND6'])

            Lmax_B7 = numpy.float32(data['LMAX_BAND7'])
            Lmax = (
                Lmax_B1, Lmax_B2, Lmax_B3, Lmax_B4, Lmax_B5, Lmax_B6, Lmax_B7)

            # Read in LMIN
            Lmin_B1 = numpy.float32(data['LMIN_BAND1'])
            Lmin_B2 = numpy.float32(data['LMIN_BAND2'])
            Lmin_B3 = numpy.float32(data['LMIN_BAND3'])
            Lmin_B4 = numpy.float32(data['LMIN_BAND4'])
            Lmin_B5 = numpy.float32(data['LMIN_BAND5'])
            if Lnum == 7:
                Lmin_B6 = numpy.float32(data['LMIN_BAND61'])
            else:
                Lmin_B6 = numpy.float32(data['LMIN_BAND6'])

            Lmin_B7 = numpy.float32(data['LMIN_BAND7'])
            Lmin = (
                Lmin_B1, Lmin_B2, Lmin_B3, Lmin_B4, Lmin_B5, Lmin_B6, Lmin_B7)

            # Read in QCALMAX
            Qcalmax_B1 = numpy.float32(data['QCALMAX_BAND1'])
            Qcalmax_B2 = numpy.float32(data['QCALMAX_BAND2'])
            Qcalmax_B3 = numpy.float32(data['QCALMAX_BAND3'])
            Qcalmax_B4 = numpy.float32(data['QCALMAX_BAND4'])
            Qcalmax_B5 = numpy.float32(data['QCALMAX_BAND5'])
            if Lnum == 7:
                Qcalmax_B6 = numpy.float32(data['QCALMAX_BAND61'])
            else:
                Qcalmax_B6 = numpy.float32(data['QCALMAX_BAND6'])

            Qcalmax_B7 = numpy.float32(data['QCALMAX_BAND7'])
            Qcalmax = (Qcalmax_B1, Qcalmax_B2, Qcalmax_B3,
                       Qcalmax_B4, Qcalmax_B5, Qcalmax_B6, Qcalmax_B7)

            # Read in QCALMIN
            Qcalmin_B1 = numpy.float32(data['QCALMIN_BAND1'])
            Qcalmin_B2 = numpy.float32(data['QCALMIN_BAND2'])
            Qcalmin_B3 = numpy.float32(data['QCALMIN_BAND3'])
            Qcalmin_B4 = numpy.float32(data['QCALMIN_BAND4'])
            Qcalmin_B5 = numpy.float32(data['QCALMIN_BAND5'])
            if Lnum == 7:
                Qcalmin_B6 = numpy.float32(data['QCALMIN_BAND61'])
            else:
                Qcalmin_B6 = numpy.float32(data['QCALMIN_BAND6'])

            Qcalmin_B7 = numpy.float32(data['QCALMIN_BAND7'])
            Qcalmin = (Qcalmin_B1, Qcalmin_B2, Qcalmin_B3,
                       Qcalmin_B4, Qcalmin_B5, Qcalmin_B6, Qcalmin_B7)

            # Read in nrows & ncols of optical bands
            Sample_ref = int(data['PRODUCT_SAMPLES_REF'])
            Line_ref = int(data['PRODUCT_LINES_REF'])
            # record ijdimension of optical bands
            ijdim_ref = (Line_ref, Sample_ref)

            Sample_thm = int(data['PRODUCT_SAMPLES_THM'])
            Line_thm = int(data['PRODUCT_LINES_THM'])
            # record thermal band dimensions (i,j)
            ijdim_thm = (Line_thm, Sample_thm)

            # Read in resolution of optical and thermal bands
            reso_ref = numpy.float32(data['GRID_CELL_SIZE_REF'])
            reso_thm = numpy.float32(data['GRID_CELL_SIZE_THM'])

            # Read in UTM Zone Number
            zc = numpy.float32(data['ZONE_NUMBER'])
            # Read in Solar Azimuth & Elevation angle (degrees)
            azi = numpy.float32(data['SUN_AZIMUTH'])
            zen = 90 - numpy.float32(data['SUN_ELEVATION'])
            # Read in upperleft mapx,y
            ulx = numpy.float32(data['PRODUCT_UL_CORNER_MAPX'])
            uly = numpy.float32(data['PRODUCT_UL_CORNER_MAPY'])
            ul = (ulx, uly)
            # Read in date of year
            char_doy = data['DATEHOUR_CONTACT_PERIOD']
            doy = int(char_doy[2:5])
        else:
            # Process New version of MTL file

            # read in LMAX
            Lmax_B1 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_1'])
            Lmax_B2 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_2'])
            Lmax_B3 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_3'])
            Lmax_B4 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_4'])
            Lmax_B5 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_5'])
            if Lnum == 7:
                Lmax_B6 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_6_VCID_1'])
            else:
                Lmax_B6 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_6'])

            Lmax_B7 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_7'])
            Lmax = (
                Lmax_B1, Lmax_B2, Lmax_B3, Lmax_B4, Lmax_B5, Lmax_B6, Lmax_B7)

            # Read in LMIN
            Lmin_B1 = numpy.float32(data['RADIANCE_MINIMUM_BAND_1'])
            Lmin_B2 = numpy.float32(data['RADIANCE_MINIMUM_BAND_2'])
            Lmin_B3 = numpy.float32(data['RADIANCE_MINIMUM_BAND_3'])
            Lmin_B4 = numpy.float32(data['RADIANCE_MINIMUM_BAND_4'])
            Lmin_B5 = numpy.float32(data['RADIANCE_MINIMUM_BAND_5'])
            if Lnum == 7:
                Lmin_B6 = numpy.float32(data['RADIANCE_MINIMUM_BAND_6_VCID_1'])
            else:
                Lmin_B6 = numpy.float32(data['RADIANCE_MINIMUM_BAND_6'])

            Lmin_B7 = numpy.float32(data['RADIANCE_MINIMUM_BAND_7'])
            Lmin = (
                Lmin_B1, Lmin_B2, Lmin_B3, Lmin_B4, Lmin_B5, Lmin_B6, Lmin_B7)

            # Read in QCALMAX
            Qcalmax_B1 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_1'])
            Qcalmax_B2 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_2'])
            Qcalmax_B3 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_3'])
            Qcalmax_B4 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_4'])
            Qcalmax_B5 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_5'])
            if Lnum == 7:
                Qcalmax_B6 = numpy.float32(
                    data['QUANTIZE_CAL_MAX_BAND_6_VCID_1'])
            else:
                Qcalmax_B6 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_6'])

            Qcalmax_B7 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_7'])
            Qcalmax = (Qcalmax_B1, Qcalmax_B2, Qcalmax_B3,
                       Qcalmax_B4, Qcalmax_B5, Qcalmax_B6, Qcalmax_B7)

            # Read in QCALMIN
            Qcalmin_B1 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_1'])
            Qcalmin_B2 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_2'])
            Qcalmin_B3 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_3'])
            Qcalmin_B4 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_4'])
            Qcalmin_B5 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_5'])
            if Lnum == 7:
                Qcalmin_B6 = numpy.float32(
                    data['QUANTIZE_CAL_MIN_BAND_6_VCID_1'])
            else:
                Qcalmin_B6 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_6'])

            Qcalmin_B7 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_7'])
            Qcalmin = (Qcalmin_B1, Qcalmin_B2, Qcalmin_B3,
                       Qcalmin_B4, Qcalmin_B5, Qcalmin_B6, Qcalmin_B7)

            # Read in nrows & ncols of optical bands
            Sample_ref = int(data['REFLECTIVE_SAMPLES'])
            Line_ref = int(data['REFLECTIVE_LINES'])
            # record ijdimension of optical bands
            ijdim_ref = (Line_ref, Sample_ref)

            Sample_thm = int(data['THERMAL_SAMPLES'])
            Line_thm = int(data['THERMAL_LINES'])
            # record thermal band dimensions (i,j)
            ijdim_thm = (Line_thm, Sample_thm)

            # Read in resolution of optical and thermal bands
            reso_ref = numpy.float32(data['GRID_CELL_SIZE_REFLECTIVE'])
            reso_thm = numpy.float32(data['GRID_CELL_SIZE_THERMAL'])

            # Read in UTM Zone Number
            zc = numpy.float32(data['UTM_ZONE'])
            # Read in Solar Azimuth & Elevation angle (degrees)
            azi = numpy.float32(data['SUN_AZIMUTH'])
            zen = 90 - numpy.float32(data['SUN_ELEVATION'])
            # Read in upperleft mapx,y
            ulx = numpy.float32(data['CORNER_UL_PROJECTION_X_PRODUCT'])
            uly = numpy.float32(data['CORNER_UL_PROJECTION_Y_PRODUCT'])
            ul = (ulx, uly)
            # Read in date of year
            char_doy = data['LANDSAT_SCENE_ID']
            doy = int(char_doy[14:16])

    elif (Lnum == 8):
            # Retrieve LS8 info
        Lmax_B2 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_2'])
        Lmax_B3 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_3'])
        Lmax_B4 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_4'])
        Lmax_B5 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_5'])
        Lmax_B6 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_6'])
        Lmax_B7 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_7'])
        Lmax_B9 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_9'])
        Lmax_B10 = numpy.float32(data['RADIANCE_MAXIMUM_BAND_10'])

        Lmax = (Lmax_B2, Lmax_B3, Lmax_B4, Lmax_B5,
                Lmax_B6, Lmax_B7, Lmax_B9, Lmax_B10)

        # Read in LMIN
        Lmin_B2 = numpy.float32(data['RADIANCE_MINIMUM_BAND_2'])
        Lmin_B3 = numpy.float32(data['RADIANCE_MINIMUM_BAND_3'])
        Lmin_B4 = numpy.float32(data['RADIANCE_MINIMUM_BAND_4'])
        Lmin_B5 = numpy.float32(data['RADIANCE_MINIMUM_BAND_5'])
        Lmin_B6 = numpy.float32(data['RADIANCE_MINIMUM_BAND_6'])
        Lmin_B7 = numpy.float32(data['RADIANCE_MINIMUM_BAND_7'])
        Lmin_B9 = numpy.float32(data['RADIANCE_MINIMUM_BAND_9'])
        Lmin_B10 = numpy.float32(data['RADIANCE_MINIMUM_BAND_10'])

        Lmin = (Lmin_B2, Lmin_B3, Lmin_B4, Lmin_B5,
                Lmin_B6, Lmin_B7, Lmin_B9, Lmin_B10)

        # Read in QCALMAX
        Qcalmax_B2 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_2'])
        Qcalmax_B3 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_3'])
        Qcalmax_B4 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_4'])
        Qcalmax_B5 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_5'])
        Qcalmax_B6 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_6'])
        Qcalmax_B7 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_7'])
        Qcalmax_B9 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_9'])
        Qcalmax_B10 = numpy.float32(data['QUANTIZE_CAL_MAX_BAND_10'])

        Qcalmax = (Qcalmax_B2, Qcalmax_B3, Qcalmax_B4, Qcalmax_B5,
                   Qcalmax_B6, Qcalmax_B7, Qcalmax_B9, Qcalmax_B10)

        # Read in QCALMIN
        Qcalmin_B2 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_2'])
        Qcalmin_B3 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_3'])
        Qcalmin_B4 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_4'])
        Qcalmin_B5 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_5'])
        Qcalmin_B6 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_6'])
        Qcalmin_B7 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_7'])
        Qcalmin_B9 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_9'])
        Qcalmin_B10 = numpy.float32(data['QUANTIZE_CAL_MIN_BAND_10'])

        Qcalmin = (Qcalmin_B2, Qcalmin_B3, Qcalmin_B4, Qcalmin_B5,
                   Qcalmin_B6, Qcalmin_B7, Qcalmin_B9, Qcalmin_B10)

        # Read in Refmax
        Refmax_B2 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_2'])
        Refmax_B3 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_3'])
        Refmax_B4 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_4'])
        Refmax_B5 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_5'])
        Refmax_B6 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_6'])
        Refmax_B7 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_7'])
        Refmax_B9 = numpy.float32(data['REFLECTANCE_MAXIMUM_BAND_9'])

        Refmax = (Refmax_B2, Refmax_B3, Refmax_B4,
                  Refmax_B5, Refmax_B6, Refmax_B7, Refmax_B9)

        # Read in Refmin
        Refmin_B2 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_2'])
        Refmin_B3 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_3'])
        Refmin_B4 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_4'])
        Refmin_B5 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_5'])
        Refmin_B6 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_6'])
        Refmin_B7 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_7'])
        Refmin_B9 = numpy.float32(data['REFLECTANCE_MINIMUM_BAND_9'])

        Refmin = (Refmin_B2, Refmin_B3, Refmin_B4,
                  Refmin_B5, Refmin_B6, Refmin_B7, Refmin_B9)

        # Read in nrows & ncols of optical bands
        Sample_ref = int(data['REFLECTIVE_SAMPLES'])
        Line_ref = int(data['REFLECTIVE_LINES'])
        # record ijdimension of optical bands
        ijdim_ref = (Line_ref, Sample_ref)

        Sample_thm = int(data['THERMAL_SAMPLES'])
        Line_thm = int(data['THERMAL_LINES'])
        # record thermal band dimensions (i,j)
        ijdim_thm = (Line_thm, Sample_thm)

        # Read in resolution of optical and thermal bands
        reso_ref = numpy.float32(data['GRID_CELL_SIZE_REFLECTIVE'])
        reso_thm = numpy.float32(data['GRID_CELL_SIZE_THERMAL'])

        # Read in UTM Zone Number
        zc = numpy.float32(data['UTM_ZONE'])
        # Read in Solar Azimuth & Elevation angle (degrees)
        azi = numpy.float32(data['SUN_AZIMUTH'])
        zen = 90 - numpy.float32(data['SUN_ELEVATION'])
        # Read in upperleft mapx,y
        ulx = numpy.float32(data['CORNER_UL_PROJECTION_X_PRODUCT'])
        uly = numpy.float32(data['CORNER_UL_PROJECTION_Y_PRODUCT'])
        ul = (ulx, uly)
        # Read in date of year
        char_doy = data['LANDSAT_SCENE_ID']
        doy = int(char_doy[13:16])
    else:
        raise Exception('This sensor is not Landsat 4, 5, 7, or 8!')

    if ((doy < 1) or (doy > 366)):
        raise ValueError(
            'Invalid Day of Year metadata value - expected (1,366) got %s' % doy)

    # The new version returns Lmax,Lmin,Qcalmax,Qcalmin,Refmax,Refmin,ijdim_ref,ijdim_thm,reso_ref,reso_thm,ul,zen,azi,zc,Lnum,doy
    # return
    # (Lmax,Lmin,Qcalmax,Qcalmin,ijdim_ref,ijdim_thm,reso_ref,reso_thm,ul,zen,azi,zc,Lnum,doy)
    return (Lmax, Lmin, Qcalmax, Qcalmin, Refmax, Refmin, ijdim_ref, ijdim_thm, reso_ref, reso_thm, ul, zen, azi, zc, Lnum, doy)


def nd2toarbt(filename, images=None):
    """
    Load metadata from MTL file & calculate reflectance values for scene bands.

    :param filename:
        A string containing the file path of the MTL file for the landsat scene.

    :param images:
        A numpy.ndarray of pre-calculated reflectance values for each landsat band, to be used instead of calculating our own.
    """
    Lmax, Lmin, Qcalmax, Qcalmin, Refmax, Refmin, ijdim_ref, ijdim_thm, reso_ref, reso_thm, ul, zen, azi, zc, Lnum, doy = lndhdrread(
        filename)

    base = os.path.dirname(filename)

    # LPGS Upper left corner alignment (see Landsat handbook for detail)
    # Changed from (ul[0]-15,ul[1]+15), GA products are 25m, this should also
    # allow for other resolutions as well
    ul = (ul[0] - float(reso_ref) / 2, ul[1] + float(reso_ref) / 2)
    resolu = (reso_ref, reso_ref)

    if ((Lnum >= 4) & (Lnum <= 7)):

        # Band6
        if Lnum == 7:
            n_B6 = match_file(base, '.*B61.*')
        else:
            n_B6 = match_file(base, '.*B6.*')

        # Check that the thermal band resolution matches the reflectance bands.
        ref_lines, ref_samples = ijdim_ref
        thm_lines, thm_samples = ijdim_thm
        if ((thm_lines != ref_lines) | (thm_samples != ref_samples)):
            im_B6 = imread(
                n_B6, resample=True, samples=ref_samples, lines=ref_lines).astype(numpy.float32)
        else:
            im_B6 = imread(n_B6).astype(numpy.float32)

        # convert Band6 from radiance to BT
        # fprintf('From Band 6 Radiance to Brightness Temperature\n')
        # see G. Chander et al. RSE 113 (2009) 893-903
        K1_L4 = 671.62
        K2_L4 = 1284.30
        K1_L5 = 607.76
        K2_L5 = 1260.56
        K1_L7 = 666.09
        K2_L7 = 1282.71

        if Lnum == 7:
            K1 = K1_L7
            K2 = K2_L7
        elif Lnum == 5:
            K1 = K1_L5
            K2 = K2_L5
        elif Lnum == 4:
            K1 = K1_L4
            K2 = K2_L4

        if images != None:
            im_B1 = images[0, :, :].astype(numpy.float32)
            im_B2 = images[1, :, :].astype(numpy.float32)
            im_B3 = images[2, :, :].astype(numpy.float32)
            im_B4 = images[3, :, :].astype(numpy.float32)
            im_B5 = images[4, :, :].astype(numpy.float32)
            im_B7 = images[6, :, :].astype(numpy.float32)
            del images

            # find pixels that are saturated in the visible bands
            B1Satu = im_B1 == 255.0
            B2Satu = im_B2 == 255.0
            B3Satu = im_B3 == 255.0

            # only processing pixesl where all bands have values (id_mssing)
            id_missing = numexpr.evaluate(
                "(im_B1 == 0.0) | (im_B2 == 0.0) | (im_B3 == 0.0) | (im_B4 == 0.0) |(im_B5 == 0.0) | (im_B6 == 0.0) | (im_B7 == 0.0)")

        else:
            # Band1
            n_B1 = match_file(base, '.*B1.*')
            im_B1 = imread(n_B1).astype(numpy.float32)
            # Band2
            n_B2 = match_file(base, '.*B2.*')
            im_B2 = imread(n_B2).astype(numpy.float32)
            # Band3
            n_B3 = match_file(base, '.*B3.*')
            im_B3 = imread(n_B3).astype(numpy.float32)
            # Band4
            n_B4 = match_file(base, '.*B4.*')
            im_B4 = imread(n_B4).astype(numpy.float32)
            # Band5
            n_B5 = match_file(base, '.*B5.*')
            im_B5 = imread(n_B5).astype(numpy.float32)
            # Band7
            n_B7 = match_file(base, '.*B7.*')
            im_B7 = imread(n_B7).astype(numpy.float32)

            # Retrieve the projection and geotransform info from the blue band
            # (B1 LS 4,5,7)
            geoT, prj = im_info(n_B1)

            # find pixels that are saturated in the visible bands
            B1Satu = im_B1 == 255.0
            B2Satu = im_B2 == 255.0
            B3Satu = im_B3 == 255.0

            # only processing pixesl where all bands have values (id_mssing)
            id_missing = numexpr.evaluate(
                "(im_B1 == 0.0) | (im_B2 == 0.0) | (im_B3 == 0.0) | (im_B4 == 0.0) | (im_B5 == 0.0) | (im_B6 == 0.0) | (im_B7 == 0.0)")

            # ND to radiance first
            im_B1 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B1 - Qmi) + Lmi", {
                                     'Lma': Lmax[0], 'Lmi': Lmin[0], 'Qma': Qcalmax[0], 'Qmi': Qcalmin[0]}, locals())
            im_B2 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B2 - Qmi) + Lmi", {
                                     'Lma': Lmax[1], 'Lmi': Lmin[1], 'Qma': Qcalmax[1], 'Qmi': Qcalmin[1]}, locals())
            im_B3 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B3 - Qmi) + Lmi", {
                                     'Lma': Lmax[2], 'Lmi': Lmin[2], 'Qma': Qcalmax[2], 'Qmi': Qcalmin[2]}, locals())
            im_B4 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B4 - Qmi) + Lmi", {
                                     'Lma': Lmax[3], 'Lmi': Lmin[3], 'Qma': Qcalmax[3], 'Qmi': Qcalmin[3]}, locals())
            im_B5 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B5 - Qmi) + Lmi", {
                                     'Lma': Lmax[4], 'Lmi': Lmin[4], 'Qma': Qcalmax[4], 'Qmi': Qcalmin[4]}, locals())
            im_B6 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B6 - Qmi) + Lmi", {
                                     'Lma': Lmax[5], 'Lmi': Lmin[5], 'Qma': Qcalmax[5], 'Qmi': Qcalmin[5]}, locals())
            im_B7 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B7 - Qmi) + Lmi", {
                                     'Lma': Lmax[6], 'Lmi': Lmin[6], 'Qma': Qcalmax[6], 'Qmi': Qcalmin[6]}, locals())

            # radiance to TOA reflectances
            # fprintf('From Radiances to TOA ref\n')
            #  # Solar Spectral Irradiances from LEDAPS
            #  esun_L7=[1969.000, 1840.000, 1551.000, 1044.000, 225.700, -1.0, 82.07]
            #  esun_L5=[1957.0, 1826.0, 1554.0, 1036.0, 215.0, -1.0, 80.67]
            #  esun_L4=[1957.0, 1825.0, 1557.0, 1033.0, 214.9, -1.0, 80.72]

            # see G. Chander et al. RSE 113 (2009) 893-903
            esun_L7 = [
                1997.000, 1812.000, 1533.000, 1039.000, 230.800, -1.0, 84.90]
            esun_L5 = [1983.0, 1796.0, 1536.0, 1031.0, 220.0, -1.0, 83.44]
            esun_L4 = [1983.0, 1795.0, 1539.0, 1028.0, 219.8, -1.0, 83.49]

            if Lnum == 7:
                ESUN = esun_L7
            elif Lnum == 5:
                ESUN = esun_L5
            elif Lnum == 4:
                ESUN = esun_L4

            #  # Interpolate earth-sun distance with day of year from LEDAPS
            #  dsun_table_doy = [1,15,32,46,60,74,91,106,121,135,152,166,182,196,213,227,242,258,274,288,305,319,335,349,366]
            #  dsun_table_dis=  [0.9832,0.9836,0.9853,0.9878,0.9909,0.9945,0.9993,1.0033,1.0076,1.0109,1.0140,1.0158,1.0167,
            #  1.0165,1.0149,1.0128,1.0092,1.0057,1.0011,0.9972,0.9925,0.9892,0.9860,0.9843,0.9833]
            #
            #  for i=1:length(dsun_table_doy)-1
            #   if doy >=dsun_table_doy(i) and doy <=dsun_table_doy(i+1)
            #      break
            #  end
            #  end
            #
            #  dsun_doy=dsun_table_dis(i)+
            #  (dsun_table_dis(i+1)-dsun_table_dis(i))*(doy-dsun_table_doy(i))/(dsun_table_doy(i+1)-dsun_table_doy(i))

            # earth-sun distance see G. Chander et al. RSE 113 (2009) 893-903
            dsun_doy = sun_earth_distance[doy]

            # compute TOA reflectances
            # converted from degrees to radiance
            s_zen = math.radians(zen)
            stack = {
                'a': numpy.float32(10000.0 * math.pi),
                'b': numpy.float32(dsun_doy * dsun_doy),
                'c': numpy.float32(math.cos(s_zen))
            }

            im_B1 = numexpr.evaluate(
                "a * im_B1 * b / (sun * c)", dict(stack.items() + {'sun': numpy.float32(ESUN[0])}.items()), locals())
            im_B2 = numexpr.evaluate(
                "a * im_B2 * b / (sun * c)", dict(stack.items() + {'sun': numpy.float32(ESUN[1])}.items()), locals())
            im_B3 = numexpr.evaluate(
                "a * im_B3 * b / (sun * c)", dict(stack.items() + {'sun': numpy.float32(ESUN[2])}.items()), locals())
            im_B4 = numexpr.evaluate(
                "a * im_B4 * b / (sun * c)", dict(stack.items() + {'sun': numpy.float32(ESUN[3])}.items()), locals())
            im_B5 = numexpr.evaluate(
                "a * im_B5 * b / (sun * c)", dict(stack.items() + {'sun': numpy.float32(ESUN[4])}.items()), locals())
            im_B7 = numexpr.evaluate(
                "a * im_B7 * b / (sun * c)", dict(stack.items() + {'sun': numpy.float32(ESUN[6])}.items()), locals())

        # convert from Kelvin to Celcius with 0.01 scale_facor
        im_B6 = numexpr.evaluate("a * ((K2 / log((K1 / im_B6) + one)) - b)", {'a': numpy.float32(
            100), 'b': numpy.float32(273.15), 'one': numpy.float32(1.0)}, locals())

        # get data ready for Fmask
        im_B1[id_missing] = -9999
        im_B2[id_missing] = -9999
        im_B3[id_missing] = -9999
        im_B4[id_missing] = -9999
        im_B5[id_missing] = -9999
        im_B6[id_missing] = -9999
        im_B7[id_missing] = -9999
        del id_missing

        images = numpy.array(
            [im_B1, im_B2, im_B3, im_B4, im_B5, im_B7], 'float32')
        del im_B1, im_B2, im_B3, im_B4, im_B5, im_B7

        # We'll modify the return argument for the Python implementation
        # (geoT,prj) are added to the list
        return [im_B6, images, ijdim_ref, ul, zen, azi, zc, B1Satu, B2Satu, B3Satu, resolu, geoT, prj]
    elif (Lnum == 8):
        n_B10 = match_file(base, '.*B10.*')
        # Check that the thermal band resolution matches the reflectance bands.
        ref_lines, ref_samples = ijdim_ref
        thm_lines, thm_samples = ijdim_thm

        if ((thm_lines != ref_lines) | (thm_samples != ref_samples)):
            im_B10 = imread(
                n_B10, resample=True, samples=ref_samples, lines=ref_lines).astype(numpy.float32)
        else:
            im_B10 = imread(n_B10).astype(numpy.float32)

        # Band2
        n_B2 = match_file(base, '.*B2.*')
        im_B2 = imread(n_B2).astype(numpy.float32)
        # Band3
        n_B3 = match_file(base, '.*B3.*')
        im_B3 = imread(n_B3).astype(numpy.float32)
        # Band4
        n_B4 = match_file(base, '.*B4.*')
        im_B4 = imread(n_B4).astype(numpy.float32)
        # Band5
        n_B5 = match_file(base, '.*B5.*')
        im_B5 = imread(n_B5).astype(numpy.float32)
        # Band6
        n_B6 = match_file(base, '.*B6.*')
        im_B6 = imread(n_B6).astype(numpy.float32)
        # Band7
        n_B7 = match_file(base, '.*B7.*')
        im_B7 = imread(n_B7).astype(numpy.float32)
        # Band9
        n_B9 = match_file(base, '.*B9.*')
        im_B9 = imread(n_B9).astype(numpy.float32)

        # Retrieve the projection and geotransform info from the blue band (B2
        # in LS8)
        geoT, prj = im_info(n_B2)

        # only processing pixesl where all bands have values (id_mssing)
        id_missing = numexpr.evaluate(
            "(im_B2 == 0.0) | (im_B3 == 0.0) | (im_B4 == 0.0) | (im_B5 == 0.0) | (im_B6 == 0.0) | (im_B7 == 0.0) | (im_B9 == 0.0) | (im_B10 == 0.0)")

        # find pixels that are saturated in the visible bands
        B1Satu = im_B2 == 65535.0
        B2Satu = im_B3 == 65535.0
        B3Satu = im_B4 == 65535.0

        # DN to TOA reflectance with 0.0001 scale_factor
        # This formulae is similar to that used for LS 4,5,7. But is different to that given by
        # https://landsat.usgs.gov/Landsat8_Using_Product.php : Noted JS
        # 2013/11/28
        print 'From DNs to TOA ref & BT\n'
        im_B2 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B2 - Qmi) + Rmi", {
                                 'Rma': Refmax[0], 'Rmi': Refmin[0], 'Qma': Qcalmax[0], 'Qmi': Qcalmin[0]}, locals())
        im_B3 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B3 - Qmi) + Rmi", {
                                 'Rma': Refmax[1], 'Rmi': Refmin[1], 'Qma': Qcalmax[1], 'Qmi': Qcalmin[1]}, locals())
        im_B4 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B4 - Qmi) + Rmi", {
                                 'Rma': Refmax[2], 'Rmi': Refmin[2], 'Qma': Qcalmax[2], 'Qmi': Qcalmin[2]}, locals())
        im_B5 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B5 - Qmi) + Rmi", {
                                 'Rma': Refmax[3], 'Rmi': Refmin[3], 'Qma': Qcalmax[3], 'Qmi': Qcalmin[3]}, locals())
        im_B6 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B6 - Qmi) + Rmi", {
                                 'Rma': Refmax[4], 'Rmi': Refmin[4], 'Qma': Qcalmax[4], 'Qmi': Qcalmin[4]}, locals())
        im_B7 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B7 - Qmi) + Rmi", {
                                 'Rma': Refmax[5], 'Rmi': Refmin[5], 'Qma': Qcalmax[5], 'Qmi': Qcalmin[5]}, locals())
        im_B9 = numexpr.evaluate("((Rma - Rmi) / (Qma - Qmi)) * (im_B9 - Qmi) + Rmi", {
                                 'Rma': Refmax[6], 'Rmi': Refmin[6], 'Qma': Qcalmax[6], 'Qmi': Qcalmin[6]}, locals())
        im_B10 = numexpr.evaluate("((Lma - Lmi) / (Qma - Qmi)) * (im_B10 - Qmi) + Lmi", {
                                  'Lma': Lmax[7], 'Lmi': Lmin[7], 'Qma': Qcalmax[7], 'Qmi': Qcalmin[7]}, locals())

        s_zen = numpy.deg2rad(zen)
        im_B2 = numexpr.evaluate("10000 * im_B2 / cos(s_zen)")
        im_B3 = numexpr.evaluate("10000 * im_B3 / cos(s_zen)")
        im_B4 = numexpr.evaluate("10000 * im_B4 / cos(s_zen)")
        im_B5 = numexpr.evaluate("10000 * im_B5 / cos(s_zen)")
        im_B6 = numexpr.evaluate("10000 * im_B6 / cos(s_zen)")
        im_B7 = numexpr.evaluate("10000 * im_B7 / cos(s_zen)")
        im_B9 = numexpr.evaluate("10000 * im_B9 / cos(s_zen)")

        # convert Band6 from radiance to BT
        # fprintf('From Band 6 Radiance to Brightness Temperature\n');
        K1_B10 = numpy.float32(774.89)
        K2_B10 = numpy.float32(1321.08)
        one = numpy.float32(1)

        im_B10 = numexpr.evaluate("K2_B10 / log((K1_B10 / im_B10) + one)")

        # convert from Kelvin to Celcius with 0.01 scale_factor
        K = numpy.float32(273.15)
        im_B10 = numexpr.evaluate("100 * (im_B10 - K)")

        # get data ready for Fmask
        im_B2[id_missing] = -9999
        im_B3[id_missing] = -9999
        im_B4[id_missing] = -9999
        im_B5[id_missing] = -9999
        im_B6[id_missing] = -9999
        im_B7[id_missing] = -9999
        im_B9[id_missing] = -9999
        im_B10[id_missing] = -9999
        del id_missing

        images = numpy.array(
            [im_B2, im_B3, im_B4, im_B5, im_B6, im_B7, im_B9], 'float32')
        del im_B2, im_B3, im_B4, im_B5, im_B6, im_B7, im_B9

        # We'll modify the return argument for the Python implementation
        # (geoT,prj) are added to the list
        return [im_B10, images, ijdim_ref, ul, zen, azi, zc, B1Satu, B2Satu, B3Satu, resolu, geoT, prj]

    else:
        raise Exception('This sensor is not Landsat 4, 5, 7, or 8!')


def plcloud(filename, cldprob=22.5, num_Lst=None, images=None, shadow_prob=False, mask=None, aux_data={}):
    """
    Calculates a cloud mask for a landsat 5/7 scene.

    :param filename:
        A string containing the file path of the landsat scene MTL file.

    :param cldprob:
        The cloud probability for the scene (defaults to 22.5%).

    :param num_Lst:
        The Landsat satellite number.

    :param images:
        A numpy.ndarray of pre-loaded, scaled, and corrected bands.

    :param log_filename:
        A string containing the file path of the output log produced by FMask.

    :param shadow_prob:
        A flag indicating if the shadow probability should be calculated or not (required by FMask cloud shadow). Type Bool.

    :return:
        Tuple (zen,azi,ptm, temperature band (celcius*100),t_templ,t_temph, water mask, snow mask, cloud mask , shadow probability,dim,ul,resolu,zc).
    """
    start_time = time.time()

    Temp, data, dim, ul, zen, azi, zc, satu_B1, satu_B2, satu_B3, resolu, geoT, prj = nd2toarbt(
        filename, images)

    if num_Lst < 8:  # Landsat 4~7
        Thin_prob = 0  # there is no contribution from the new bands
    else:
        Thin_prob = numexpr.evaluate(
            "cirrus / 400", {'cirrus': data[-1]}, locals())

    Cloud = numpy.zeros(dim, 'uint8')  # cloud mask
    Snow = numpy.zeros(dim, 'uint8')  # Snow mask
    WT = numpy.zeros(dim, 'uint8')  # Water msk

    # process only the overlap area
    if mask == None:
        mask = Temp > -9999
    else:
        mask = mask.astype('bool')

    Shadow = numpy.zeros(dim, 'uint8')  # shadow mask

    data1 = data[0, :, :]
    data2 = data[1, :, :]
    data3 = data[2, :, :]
    data4 = data[3, :, :]
    data5 = data[4, :, :]
    data6 = data[5, :, :]

    NDVI = numexpr.evaluate("(data4 - data3) / (data4 + data3)")
    NDSI = numexpr.evaluate("(data2 - data5) / (data2 + data5)")

    NDVI[numexpr.evaluate("(data4 + data3) == 0")] = 0.01
    NDSI[numexpr.evaluate("(data2 + data5) == 0")] = 0.01

    # saturation in the three visible bands
    satu_Bv = numexpr.evaluate("(satu_B1 | satu_B2 | satu_B3)")
    del satu_B1
    # Basic cloud test
    idplcd = numexpr.evaluate(
        "(NDSI < 0.8) & (NDVI < 0.8) & (data6 > 300) & (Temp < 2700)")

    # Snow test
    # It takes every snow pixels including snow pixel under thin clouds or icy
    # clouds
    Snow[numexpr.evaluate(
        "(NDSI > 0.15) & (Temp < 1000) & (data4 > 1100) & (data2 > 1000)")] = 1
    #Snow[mask == 0] = 255
    # Water test
    # Zhe's water test (works over thin cloud)
    WT[numexpr.evaluate(
        "((NDVI < 0.01) & (data4 < 1100)) | ((NDVI < 0.1) & (NDVI > 0) & (data4 < 500))")] = 1
    WT[mask == 0] = 255
    # ################################################ Whiteness test
    # visible bands flatness (sum(abs)/mean < 0.6 => brigt and dark cloud )
    visimean = numexpr.evaluate("(data1 + data2 + data3) / 3 ")
    whiteness = numexpr.evaluate(
        "(abs(data1 - visimean) + abs(data2 - visimean)+ abs(data3 - visimean)) / visimean")
    del visimean

    # update idplcd
    whiteness[satu_Bv] = 0  # If one visible is saturated whiteness == 0
    idplcd &= whiteness < 0.7

    # Haze test
    HOT = numexpr.evaluate("data1 - 0.5 * data3 - 800")  # Haze test
    idplcd &= numexpr.evaluate("(HOT > 0) | satu_Bv")
    del HOT  # need to find thick warm cloud

    # Ratio4/5>0.75 cloud test
    idplcd &= numexpr.evaluate("(data4 / data5) > 0.75")

    # Cirrus tests from Landsat 8
    idplcd |= numexpr.evaluate("Thin_prob > 0.25")

    ####################################constants##########################
    l_pt = 0.175  # low percent
    h_pt = 1 - l_pt  # high percent
    # (temperature & snow test )
    # test whether use thermal or not
    idclr = numexpr.evaluate("(idplcd == False) & (mask == 1)")
    ptm = 100 * idclr.sum() / mask.sum()  # percent of del pixel
    idlnd = numexpr.evaluate("idclr & (WT == False)")
    idwt = numexpr.evaluate("idclr & (WT == True)")  # &data(:,:,6)<=300;
    lndptm = 100 * idlnd.sum() / mask.sum()

    logging.debug('idlnd.sum(): %s', idlnd.sum())
    logging.debug('lndptm: %s', lndptm)
    sys.stdout.flush()

    if ptm <= 0.1:  # no thermal test => meanless for snow detection (0~1)
        # fprintf('Clear pixel NOT exist in this scene (del prct =
        # #.2f)\n',ptm)
        Cloud[idplcd] = 1  # all cld

        # mask out the non-contiguous pixels
        Cloud[~(mask)] = 0

        # # improving by majority filtering
        # Cloud=bwmorph(Cloud,'majority')# exclude <5/9
        #Cloud = scipy.signal.convolve2d(Cloud,numpy.ones((3,3),Cloud.dtype.name))[1:-1,1:-1]
        #Cloud = (Cloud > 4).astype('uint8')
        # Applying twice, makes a cleaner result
        #Cloud = scipy.signal.convolve2d(Cloud,numpy.ones((3,3),Cloud.dtype.name))[1:-1,1:-1]
        #Cloud = Cloud > 4

        Shadow[Cloud == 0] = 1
        Temp = -1
        t_templ = -1
        t_temph = -1
    else:
        # fprintf('Clear pixel EXIST in this scene (del prct = #.2f)\n',ptm)
        # (temperature test )
        if lndptm >= 0.1:
            F_temp = Temp[idlnd]  # get land temperature
            #       fprintf('Land temperature\n')
        else:
            F_temp = Temp[idclr]  # get del temperature
            #        fprintf('Clear temperature\n')

        # Get cloud prob over water
        # temperature test (over water)
        # F_wtemp = Temp[numexpr.evaluate("(WT == 1) & (data6 <= 300)")] # get
        # del water temperature
        F_wtemp = Temp[idwt]
        if len(F_wtemp) == 0:
            t_wtemp = 0
        else:
            t_wtemp = scipy.stats.scoreatpercentile(F_wtemp, 100 * h_pt)
        wTemp_prob = numexpr.evaluate('(t_wtemp - Temp) / 400')
        wTemp_prob[numexpr.evaluate('wTemp_prob < 0')] = 0

        # Brightness test (over water)
        t_bright = 1100
        Brightness_prob = data5 / t_bright
        Brightness_prob[Brightness_prob > 1] = 1
        Brightness_prob[Brightness_prob < 0] = 0

        # Final prob mask (water)
        # cloud over water probability
        wfinal_prob = numexpr.evaluate(
            '100 * wTemp_prob * Brightness_prob + 100 * Thin_prob')
        wclr_max = scipy.stats.scoreatpercentile(
            wfinal_prob[idwt], 100 * h_pt) + cldprob  # dynamic threshold (land)
        # wclr_max=50;% fixed threshold (water)

        # release memory
        del wTemp_prob
        del Brightness_prob

        # Temperature test
        t_buffer = 4 * 100
        if len(F_temp) != 0:
            # 0.175 percentile background temperature (low)
            t_templ = scipy.stats.scoreatpercentile(F_temp, 100 * l_pt)
            # 0.825 percentile background temperature (high)
            t_temph = scipy.stats.scoreatpercentile(F_temp, 100 * h_pt)
        else:
            t_templ = 0
            t_temph = 0

        t_tempL = t_templ - t_buffer
        t_tempH = t_temph + t_buffer
        Temp_l = t_tempH - t_tempL
        Temp_prob = (t_tempH - Temp) / Temp_l
        # Temperature can have prob > 1
        Temp_prob[Temp_prob < 0] = 0
        # Temp_prob(Temp_prob > 1) = 1

        NDSI[numexpr.evaluate('satu_B2 & (NDSI < 0)')] = 0
        NDVI[numexpr.evaluate('satu_B3 & (NDVI > 0)')] = 0

        Vari_prob = 1 - \
            numpy.maximum(
                numpy.maximum(numpy.absolute(NDSI), numpy.absolute(NDVI)), whiteness)

        # release memory
        logging.debug("FMASK releasing memory")
        del satu_B2
        del satu_B3
        del NDSI
        del NDVI
        del whiteness

        # Final prob mask (land)
        final_prob = 100 * Temp_prob * Vari_prob + 100 * \
            Thin_prob  # cloud over land probability
        clr_max = scipy.stats.scoreatpercentile(
            final_prob[idlnd], 100 * h_pt) + cldprob  # dynamic threshold (land)

        # release memory
        del Vari_prob
        del Temp_prob
        del Thin_prob

        logging.debug('cldprob: %s', cldprob)
        logging.debug('clr_max: %s', clr_max)
        logging.debug('t_templ: %s', t_templ)
        sys.stdout.flush()

        # fprintf('pcloud probability threshold (land) = .2f#\n',clr_max)
        # cloud over land : (idplcd & (final_prob > clr_max) & (WT == 0))
        # thin cloud over water : (idplcd & (wfinal_prob > wclr_max) & (WT == 1))
        # high prob cloud (land) : (final_prob > 99.0) & (WT == 0)
        # extremly cold cloud : (Temp < t_templ - 3500)
        id_final_cld = numexpr.evaluate(
            '(idplcd & (final_prob > clr_max) & (WT == 0)) | (idplcd & (wfinal_prob > wclr_max) & (WT == 1)) | (Temp < t_templ - 3500)')

        # Star with potential cloud mask
        # # potential cloud mask
        Cloud[id_final_cld] = 1

        # release memory
        del final_prob
        del wfinal_prob
        del id_final_cld

        # Start with potential cloud shadow mask
        if shadow_prob:
            # band 4 flood fill
            nir = data4.astype('float32')
            # estimating background (land) Band 4 ref
            backg_B4 = scipy.stats.scoreatpercentile(nir[idlnd], 100.0 * l_pt)
            nir[mask == 0] = backg_B4
            # fill in regional minimum Band 4 ref
            nir = imfill_skimage(nir)
            nir = nir - data4

            # band 5 flood fill
            swir = data5
            # estimating background (land) Band 4 ref
            backg_B5 = scipy.stats.scoreatpercentile(swir[idlnd], 100.0 * l_pt)
            swir[mask == 0] = backg_B5
            # fill in regional minimum Band 5 ref
            swir = imfill_skimage(swir)
            swir = swir - data5

            # compute shadow probability
            shadow_prob = numpy.minimum(nir, swir)
            # release remory
            del nir
            del swir

            Shadow[shadow_prob > 200] = 1
            # release remory
            del shadow_prob

        # Cloud[idplcd==True]=1 # all cld

        #*************************************************************************************#
        #*************************************************************************************#
        #*************************************************************************************#
        #************The following code may be removed as the new code has changed************#
        # Mask the non-contiguous pixels
        #Cloud[~(mask)] = 0

        # improving by majority filtering
        # ERROR: not aware of a similar filter in scipy (though one may very well exist)
        # Doing convolution of all surrounding pixels & filtering those > 4, which should in theory have the same result
        #Cloud = scipy.signal.convolve2d(Cloud, numpy.ones((3,3), Cloud.dtype.name))[1:-1,1:-1]
        #Cloud = (Cloud > 4).astype('uint8')
        # Applying twice, makes a cleaner result
        #Cloud = scipy.signal.convolve2d(Cloud, numpy.ones((3,3), Cloud.dtype.name))[1:-1,1:-1]
        # Cloud=(Cloud>4).astype('uint8')
        # 3rd, still some single pixels at tile edges
        # Cloud=scipy.signal.convolve2d(Cloud,numpy.ones((3,3),Cloud.dtype.name))[1:-1,1:-1]
        #Cloud = Cloud > 4
        #************The above code may be removed as the new code has changed****************#
        #*************************************************************************************#
        #*************************************************************************************#
        #*************************************************************************************#

    del data
    images = None
    gc.collect()

    # refine Water mask - Zhe's water mask (no confusion water/cloud)
    WT[numexpr.evaluate("(WT == 1) & (Cloud == 0)")] = 1
    # bwmorph changed Cloud to Binary
    Cloud = Cloud.astype('uint8')
    Cloud[mask == 0] = 255
    Shadow[mask == 0] = 255
    processing_time = time.time() - start_time

    gc.collect()
    cloud_mask = (Cloud == 1) & mask

    if ptm > 0.1:
        cloud_temp = Temp[cloud_mask]

        logging.debug("FMASK snow percent: %f" %
                      ((float(Snow[mask].sum()) / float(mask.sum())) * 100.0))
        aux_data['fmask_snow_percent'] = float(
            Snow[mask].sum()) * 100.0 / float(mask.sum())
        logging.debug("FMASK cloud mean: %f C" %
                      (numpy.mean(cloud_temp) / 100.0))
        aux_data['fmask_cloud_mean_temp_deg_C'] = numpy.mean(
            cloud_temp) / 100.0

        if numpy.sum(cloud_mask) > 0:
            cloud_stddev = numpy.std(
                cloud_temp, dtype='float64', ddof=1) / 100.0
            pct_upper = numpy.percentile(cloud_temp, 97.5) / 100.0
            pct_lower = numpy.percentile(cloud_temp, 83.5) / 100.0
            pct_upper_max = numpy.percentile(cloud_temp, 98.75) / 100.0

            logging.debug("FMASK Standard Deviation: %f C" % cloud_stddev)
            aux_data['FMASK_std_dev_degC'] = cloud_stddev
            logging.debug("FMASK 97.5 percentile: %f C" % pct_upper)
            aux_data['FMASK_97.5_percentile'] = pct_upper
            logging.debug("FMASK 83.5 percentile: %f C" % pct_lower)
            aux_data['FMASK_83.5_percentile'] = pct_lower
            logging.debug("FMASK 98.75 percentile: %f C" % pct_upper_max)
            aux_data['FMASK_98.75_percentile'] = pct_upper_max

    cloud_skew = 0.0  # TODO

    logging.debug("FMASK Final Cloud Layer Percent: %f" %
                  ((float(Cloud[cloud_mask].sum()) / float(mask.sum())) * 100.0))
    aux_data['FMASK_cloud_layer_percent'] = (
        float(Cloud[cloud_mask].sum()) / float(mask.sum())) * 100.0
    aux_data['FMASK_processing_time'] = processing_time

    # We'll modify the return argument for the Python implementation
    # (geoT,prj) are added to the list
    return (zen, azi, ptm, Temp, t_templ, t_temph, WT, Snow, Cloud, Shadow, dim, ul, resolu, zc, geoT, prj)


def fcssm(Sun_zen, Sun_azi, ptm, Temp, t_templ, t_temph, Water, Snow, plcim, plsim, ijDim, resolu, ZC, cldpix, sdpix, snpix):
    """
    NEW:
    fcssm(dir_im,Sun_zen,Sun_azi,ptm,Temp,...
        t_templ,t_temph,Water,Snow,plcim,plsim,ijDim,jiUL,resolu,ZC,cldpix,sdpix,snpix)
    ORIGINAL:
    fcssm_1_6sav(dir_im,Sun_zen,Sun_azi,ptm,Temp,...
        t_templ,t_temph,Water,Snow,plcim,plsim,ijDim,jiUL,resolu,ZC,cldpix,sdpix)
    """
    """
    Calculates the cloud shadow mask for a scene, given solar geometry information, the thermal band for the scene & a cloud mask.

    :param Sun_zen:
        Solar Elevation angle (degrees).

    :param Sun_azi:
        Solar Azimuth angle (degrees).

    :param ptm:
        Percentage of deleted pixels.

    :param Temp:
        A numpy.ndarray containing the temperature band for the landsat scene (Celcius*100).
    :param t_templ:
        0.175 percentile background temperature (low).

    :param t_temph:
        0.825 percentile background temperature (high).

    :param Water:
        A numpy.ndarray of type Bool containing the water mask calculated by FMask.

    :param Snow:
        A numpy.ndarray of type Bool containing the snow mask calculated by FMask.

    :param plcim:
        A numpy.ndarray of type Bool containing the cloud mask calculated by FMask.

    :param plsim:
        A numpy.ndarray of type Bool containing the cloud shadow mask calculated by FMask.

    :param ijDim:
        A tuple containing the resolution of the scene bands (height, width).

    :param resolu:
        A tuple (number, numpber).

    :param ZC:
        The UTM Zone Number of the scene.

    :param cldpix:
        A number for the cloud mask dilation (in pixels).

    :param sdpix:
        A number for the cloud shadow mask dilation (in pixels)
    """
    # Function for Cloud, cloud Shadow, and Snow Masking 1.6.3sav
    # History of revisions:
    # cloud shadow do not have to overlap with potential cloud shadow layer (Zhe Zhu 04/24/2011)
    # exclude small cloud object <= 25 pixels (zhe Zhu 3/07/2011)
    # dilate shadow again (3 pixels as default) (Zhe Zhu 12/23/2010);
    # similarity < 0.95 (Zhe Zhu 11/06/2010)
    # boosts data by >5/9 (Zhe Zhu 12/08/2009)
    # use temperature to narrow iteration height (Zhe Zhu 12/09/2009)
    # fixed bug for height (Zhe Zhu 12/09/2009)
    # cloud DEM by thermal in cloud and shadow match (Zhe Zhu 1/03/2009)

    # solar elevation angle
    Sun_ele = 90.0 - Sun_zen
    sun_ele_rad = math.radians(Sun_ele)
    # solar azimuth anngle
    Sun_tazi = Sun_azi - 90.0
    sun_tazi_rad = math.radians(Sun_tazi)
    # assume resolu.x=resolu.y
    sub_size = resolu[0]
    win_height = ijDim[0]
    win_width = ijDim[1]

    # potential cloud & shadow layer
    cloud_test = numpy.zeros(ijDim, 'uint8')
    shadow_test = numpy.zeros(ijDim, 'uint8')
    # matched cloud & shadow layer
    shadow_cal = numpy.zeros(ijDim, 'uint8')
    cloud_cal = numpy.zeros(ijDim, 'uint8')
    # cloud_height=zeros(ijDim)# cloud relative height (m)
    # boundary layer
    boundary_test = numpy.zeros(ijDim, 'uint8')
    # final cloud, shadow and snow mask
    cs_final = numpy.zeros(ijDim, 'uint8')

    # get potential mask values
    shadow_test[plsim == 1] = 1  # plshadow layer
    del plsim  # empty memory

    boundary_test[plcim < 255] = 1  # boundary layer
    cloud_test[plcim == 1] = 1  # plcloud layer
    del plcim  # empty memory

    # revised percent of cloud on the scene after plcloud
    revised_ptm = numpy.sum(cloud_test) / numpy.sum(boundary_test)
    # no t test  => more than 98 # clouds and partly cloud over land
    # => no match => rest are definite shadows

    # cloud covers more than 90# of the scene
    # => no match => rest are definite shadows
    # fprintf('Cloud and cloud shadow matching ...\n')

    if ptm <= 0.1 or revised_ptm >= 0.90:
        #     fprintf('No Shadow Match due to too much cloud (>90 percent)\n')
        cloud_cal[cloud_test == True] = 1
        shadow_cal[cloud_test == False] = 1
        similar_num = -1
        #   height_num=-1

    else:
        #     fprintf('Shadow Match in processing\n')

        # define constants
        Tsimilar = 0.30
        Tbuffer = 0.95  # threshold for matching buffering
        max_similar = 0.95  # max similarity threshold
        num_cldoj = 3  # minimum matched cloud object (pixels)
        num_pix = 3  # number of inward pixes (90m) for cloud base temperature

        # enviromental lapse rate 6.5 degrees/km
        # dry adiabatic lapse rate 9.8 degrees/km
        rate_elapse = 6.5  # degrees/km
        rate_dlapse = 9.8  # degrees/km

        #     fprintf('Set cloud similarity = #.3f\n',Tsimilar)
        #     fprintf('Set matching buffer = #.3f\n',Tbuffer)
        #     fprintf('Shadow match for cloud object >= #d pixels\n',num_cldoj)

        i_step = 2 * sub_size * math.tan(sun_ele_rad)  # move 2 pixel at a time

        # get moving direction
        (rows, cols) = numpy.nonzero(boundary_test)
        (y_ul, num) = (rows.min(), rows.argmin())
        x_ul = cols[num]

        (y_lr, num) = (rows.max(), rows.argmax())
        x_lr = cols[num]

        (x_ll, num) = (cols.min(), cols.argmin())
        y_ll = rows[num]

        (x_ur, num) = (cols.max(), cols.argmax())
        y_ur = rows[num]

        # get view angle geometry
        (A, B, C, omiga_par, omiga_per) = viewgeo(float(x_ul), float(y_ul), float(
            x_ur), float(y_ur), float(x_ll), float(y_ll), float(x_lr), float(y_lr))

        # Segmentate each cloud
        #     fprintf('Cloud segmentation & matching\n')
        (segm_cloud_init, segm_cloud_init_features) = scipy.ndimage.measurements.label(
            cloud_test, scipy.ndimage.morphology.generate_binary_structure(2, 2))

        # filter out cloud object < than num_cldoj pixels
        morphology.remove_small_objects(
            segm_cloud_init, num_cldoj, in_place=True)
        segm_cloud, fw, inv = segmentation.relabel_from_one(segm_cloud_init)
        num = numpy.max(segm_cloud)

        # NOTE: properties is deprecated as of version 0.9 and all properties are computed. Currently using version 0.8.2. If this version or a later versionproves too slow, I'll implement another method. JS 16/12/2013
        # The properties is taking approx 3min, I can cut that down to just a
        # few seconds using another method, but will leave for the time being.
        s = measure.regionprops(
            segm_cloud,  properties=['Area', 'Coordinates'])

        # Use iteration to get the optimal move distance
        # Calulate the moving cloud shadow

        # height_num=zeros(num) # cloud relative height (m)
        similar_num = numpy.zeros(num)  # cloud shadow match similarity (m)

        # Newer method of looping through the cloud objects/segments JS
        # 16/12/2013
        for cloud_type in s:

            cld_area = cloud_type['Area']
            cld_label = cloud_type['Label']

            num_pixels = cld_area

            # Have re-formatted the arrays to be Python style memory ordering (2,num_pixels) JS 16/12/2013
            # moving cloud xys
            XY_type = numpy.zeros((2, num_pixels), dtype='uint32')

            # record the max threshold moving cloud xys
            tmp_XY_type = numpy.zeros((2, num_pixels), dtype='uint32')

            # corrected for view angle xys
            # Leave as float for the time being
            tmp_xys = numpy.zeros((2, num_pixels))

            # record this original ids
            orin_cid = (
                cloud_type['Coordinates'][:, 0], cloud_type['Coordinates'][:, 1])

            # Temperature of the cloud object
            temp_obj = Temp[orin_cid]

            # assume object is round r_obj is radium of object
            r_obj = math.sqrt(cld_area / math.pi)

            # number of inward pixes for correct temperature
            #        num_pix=8
            pct_obj = math.pow(r_obj - num_pix, 2) / math.pow(r_obj, 2)
            # pct of edge pixel should be less than 1
            pct_obj = numpy.minimum(pct_obj, 1)
            t_obj = scipy.stats.mstats.mquantiles(temp_obj, pct_obj)

            # put the edge of the cloud the same value as t_obj
            temp_obj[temp_obj > t_obj] = t_obj

            # wet adiabatic lapse rate 6.5 degrees/km
            # dry adiabatic lapse rate 9.8 degrees/km
            #        rate_wlapse=6.5# degrees/km
            #        rate_dlapse=9.8# degrees/km

            Max_cl_height = 12000  # Max cloud base height (m)
            Min_cl_height = 200  # Min cloud base height (m)

            # refine cloud height range (m)
            Min_cl_height = max(
                Min_cl_height, 10 * (t_templ - 400 - t_obj) / rate_dlapse)
            Max_cl_height = min(Max_cl_height, 10 * (t_temph + 400 - t_obj))

            # initialize height and similarity info
            record_h = 0.0
            record_thresh = 0.0

            # iterate in height (m)
            for base_h in numpy.arange(Min_cl_height, Max_cl_height, i_step):
                # Get the true postion of the cloud
                # calculate cloud DEM with initial base height
                h = (10 * (t_obj - temp_obj) / rate_elapse + base_h)
                tmp_xys[1, :], tmp_xys[0, :] = mat_truecloud(orin_cid[1], orin_cid[
                                                             0], h, A, B, C, omiga_par, omiga_per)  # Function is returned as (x_new,y_new)

                # shadow moved distance (pixel)
                # i_xy=h*cos(sun_tazi_rad)/(sub_size*math.tan(sun_ele_rad))
                i_xy = h / (sub_size * math.tan(sun_ele_rad))

                if Sun_azi < 180:
                    XY_type[1, :] = numpy.round(
                        tmp_xys[1, :] - i_xy * math.cos(sun_tazi_rad))  # X is for j,1
                    XY_type[0, :] = numpy.round(
                        tmp_xys[0, :] - i_xy * math.sin(sun_tazi_rad))  # Y is for i,0
                else:
                    XY_type[1, :] = numpy.round(
                        tmp_xys[1, :] + i_xy * math.cos(sun_tazi_rad))  # X is for j,1
                    XY_type[0, :] = numpy.round(
                        tmp_xys[0, :] + i_xy * math.sin(sun_tazi_rad))  # Y is for i,0

                tmp_j = XY_type[1, :]  # col
                tmp_i = XY_type[0, :]  # row

                # the id that is out of the image
                out_id = (tmp_i < 0) | (tmp_i >= win_height) | (
                    tmp_j < 0) | (tmp_j >= win_width)
                out_all = numpy.sum(out_id)

                tmp_ii = tmp_i[out_id == 0]
                tmp_jj = tmp_j[out_id == 0]

                tmp_id = [tmp_ii, tmp_jj]

                # the id that is matched (exclude original cloud)
                match_id = numexpr.evaluate("(b_test == 0) | ((seg != label) & ((cld_test > 0) | (shad_test == 1)))", {'b_test': boundary_test[
                                            tmp_id], 'seg': segm_cloud[tmp_id], 'label': cld_label, 'cld_test': cloud_test[tmp_id], 'shad_test': shadow_test[tmp_id]})
                matched_all = numpy.sum(match_id) + out_all

                # the id that is the total pixel (exclude original cloud)
                total_id = segm_cloud[tmp_id] != cld_label
                total_all = numpy.sum(total_id) + out_all

                thresh_match = numpy.float32(matched_all) / total_all
                if (thresh_match >= (Tbuffer * record_thresh)) and (base_h < (Max_cl_height - i_step)) and (record_thresh < 0.95):
                    if thresh_match > record_thresh:
                        record_thresh = thresh_match
                        record_h = h

                elif record_thresh > Tsimilar:
                    # -1 to account for the zero based index used by Python (MATLAB is 1 one based).
                    similar_num[cld_label - 1] = record_thresh
                    i_vir = record_h / (sub_size * math.tan(sun_ele_rad))
                    # height_num=record_h

                    if Sun_azi < 180:
                        tmp_XY_type[1, :] = numpy.round(
                            tmp_xys[1, :] - i_vir * math.cos(sun_tazi_rad))  # X is for col j,2
                        tmp_XY_type[0, :] = numpy.round(
                            tmp_xys[0, :] - i_vir * math.sin(sun_tazi_rad))  # Y is for row i,1
                    else:
                        tmp_XY_type[1, :] = numpy.round(
                            tmp_xys[1, :] + i_vir * math.cos(sun_tazi_rad))  # X is for col j,2
                        tmp_XY_type[0, :] = numpy.round(
                            tmp_xys[0, :] + i_vir * math.sin(sun_tazi_rad))  # Y is for row i,1

                    tmp_scol = tmp_XY_type[1, :]
                    tmp_srow = tmp_XY_type[0, :]

                    # put data within range
                    tmp_srow[tmp_srow < 0] = 0
                    tmp_srow[tmp_srow >= win_height] = win_height - 1
                    tmp_scol[tmp_scol < 0] = 0
                    tmp_scol[tmp_scol >= win_width] = win_width - 1

                    # sub2ind(ijDim,tmp_srow,tmp_scol)
                    tmp_sid = [tmp_srow, tmp_scol]
                    # give shadow_cal=1
                    shadow_cal[tmp_sid] = 1
                    # record matched cloud
                    # cloud_cal(orin_cid)=1
                    # cloud_height[orin_cid]=record_h
                    break
                else:
                    record_thresh = 0.0

        # # dilate each cloud and shadow object by 3 and 6 pixel outward in 8 connect directions
        #    cldpix=3 # number of pixels to be dilated for cloud
        #    sdpix=3 # number of pixels to be dilated for shadow
        # fprintf('Dilate #d pixels for cloud & #d pixels for shadow
        # objects\n',cldpix,sdpix)

        # NOTE: An alternative for a structuring element would be to use the iterations parameter of binary_dilation()
        # The number of iterations is equal to the number of dilations if using
        # a 3x3 structuring element. JS 16/12/2013
        SEc = 2 * cldpix + 1
        SEc = numpy.ones((SEc, SEc), 'uint8')
        SEs = 2 * sdpix + 1
        SEs = numpy.ones((SEs, SEs), 'uint8')
        SEsn = 2 * snpix + 1
        SEsn = numpy.ones((SEsn, SEsn), 'uint8')

        # dilate shadow first
        # NOTE: The original transcription returned the inverse, i.e.
        # cloud_shadow = 0 rather than 1. We'll try inverting it outside this
        # function in order to preserve the original return values of Fmask
        shadow_cal = scipy.ndimage.morphology.binary_dilation(
            shadow_cal, structure=SEs)

        #     # find shadow within plshadow
        #     shadow_cal(shadow_test~=1)=0
        #     # dilate shadow again with the more accurate cloud shadow
        #     shadow_cal=imdilate(shadow_cal,SEs)

        segm_cloud_tmp = numexpr.evaluate("segm_cloud != 0")
        # NOTE: The original transcription returned the inverse, i.e. cloud = 0
        # rather than 1. We'll try inverting it outside this function in order
        # to preserve the original return values of Fmask
        cloud_cal = scipy.ndimage.morphology.binary_dilation(
            segm_cloud_tmp, structure=SEc)

        Snow = scipy.ndimage.morphology.binary_dilation(Snow, structure=SEsn)

    cs_final[Water == 1] = 1
    # mask from plcloud
    # step 1 snow or unknow
    cs_final[Snow == 1] = 3  # snow
    # step 2 shadow above snow and everyting
    cs_final[shadow_cal == 1] = 2  # shadow
    # step 3 cloud above all
    cs_final[cloud_cal == 1] = 4  # cloud
    cs_final[boundary_test == 0] = 255

    # record cloud and cloud shadow percent
    tmpcs = ((cs_final == 1) | (cs_final == 3)).astype('uint8')
    cspt = 100.0 * (tmpcs.sum() / boundary_test.sum())

    return (similar_num, cspt, shadow_cal, cs_final)

# viewgeo function


def viewgeo(x_ul, y_ul, x_ur, y_ur, x_ll, y_ll, x_lr, y_lr):
    # imput "x",j
    # imput "y",i
    # imput cloud height "h"

    x_u = (x_ul + x_ur) / 2
    x_l = (x_ll + x_lr) / 2
    y_u = (y_ul + y_ur) / 2
    y_l = (y_ll + y_lr) / 2

    # get k of the upper left and right points
    K_ulr = (y_ul - y_ur) / (x_ul - x_ur)
    # get k of the lower left and right points
    K_llr = (y_ll - y_lr) / (x_ll - x_lr)
    K_aver = (K_ulr + K_llr) / 2
    omiga_par = math.atan(K_aver)  # get the angle of parallel lines k (in pi)

    # AX(j)+BY(i)+C=0
    A = y_u - y_l
    B = x_l - x_u
    C = y_l * x_u - x_l * y_u

    # get the angle which is perpendicular to the trace line
    omiga_per = math.atan(B / A)
    return (A, B, C, omiga_par, omiga_per)

# mat_truecloud function


def mat_truecloud(x, y, h, A, B, C, omiga_par, omiga_per):
    # imput "x",j col
    # imput "y",i row
    # imput cloud height "h"
    H = 705000  # average Landsat 7 height (m)
    # from the cetral perpendicular (unit: pixel)
    dist = (A * x + B * y + C) / math.sqrt(A * A + B * B)
    dist_par = dist / math.cos(omiga_per - omiga_par)
    dist_move = dist_par * h / H  # cloud move distance (m)
    delt_x = dist_move * math.cos(omiga_par)
    delt_y = dist_move * math.sin(omiga_par)

    x_new = x + delt_x  # new x, j
    y_new = y + delt_y  # new y, i

    return (x_new, y_new)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Computes the Fmask algorithm. Cloud, cloud shadow and Fmask combined (contains thecloud, cloud shadow and snow masks in a single array) are output to disk.')
    parser.add_argument(
        '--mtl', required=True, help='The full file path to the Landsat MTL file.')
    parser.add_argument('--cldprob', type=float, default=22.5,
                        help='The cloud probability for the scene. Default is 22.5 percent.')
    parser.add_argument('--cldpix', type=int, default=3,
                        help='The number of pixels to be dilated for the cloud mask. Default is 3.')
    parser.add_argument('--sdpix', type=int, default=3,
                        help='The number of pixels to be dilated for the cloud shadow mask. Default is 3.')
    parser.add_argument('--snpix', type=int, default=3,
                        help='The number of pixels to be dilated for the snow mask. Default is 3.')
    parser.add_argument('--outdir', required=True,
                        help='The full file path of the output directory that will contain the Fmask results.')

    parsed_args = parser.parse_args()
    mtl = parsed_args.mtl
    cldprob = parsed_args.cldprob
    cldpix = parsed_args.cldpix
    sdpix = parsed_args.sdpix
    snpix = parsed_args.snpix
    outdir = parsed_args.outdir

    # Check that the MTL file exists
    assert os.path.exists(mtl), "Invalid filename: %s" % mtl

    # Check that the output directory exists
    assert os.path.exists(outdir), "Directory doesn't exist: %s" % outdir

    # Create the output filenames
    log_fname = os.path.join(outdir, 'FMASK_LOGFILE.txt')
    cloud_fname = os.path.join(outdir, 'fmask_cloud')
    cloud_shadow_fname = os.path.join(outdir, 'fmask_cloud_shadow')
    fmask_fname = os.path.join(outdir, 'fmask')

    # Open the MTL file.
    # The original MATLAB code opens the file twice to retrieve the Landsat number.
    # It would be better to open it once and restructure the function
    # parameters.
    data = {}
    fl = open(mtl, 'r')
    file_lines = fl.readlines()
    for line in file_lines:
        values = line.split(' = ')
        if len(values) != 2:
            continue

        data[values[0].strip()] = values[1].strip().strip('"')

    fl.close()

    # Identify Landsat Number (Lnum = 4, 5 or 7)
    LID = data['SPACECRAFT_ID']
    Lnum = int(LID[len(LID) - 1])

    st = datetime.datetime.now()
    zen, azi, ptm, Temp, t_templ, t_temph, WT, Snow, Cloud, Shadow, dim, ul, resolu, zc, geoT, prj = plcloud(
        mtl, cldprob, num_Lst=Lnum, shadow_prob=True, log_filename=log_fname)
    et = datetime.datetime.now()
    print 'time taken for plcloud function: ', et - st
    st = datetime.datetime.now()
    similar_num, cspt, shadow_cal, cs_final = fcssm(
        zen, azi, ptm, Temp, t_templ, t_temph, WT, Snow, Cloud, Shadow, dim, resolu, zc, cldpix, sdpix, snpix)
    et = datetime.datetime.now()
    print 'time taken for fcssm function: ', et - st

    c = gdal.GetDriverByName('ENVI').Create(
        cloud_fname, Cloud.shape[1], Cloud.shape[0], 1, gdal.GDT_Byte)
    c.SetGeoTransform(geoT)
    c.SetProjection(prj)
    c.GetRasterBand(1).WriteArray(Cloud * 255)
    c = None

    c = gdal.GetDriverByName('ENVI').Create(
        cloud_shadow_fname, shadow_cal.shape[1], shadow_cal.shape[0], 1, gdal.GDT_Byte)
    c.SetGeoTransform(geoT)
    c.SetProjection(prj)
    c.GetRasterBand(1).WriteArray(shadow_cal * 255)
    c = None

    c = gdal.GetDriverByName('ENVI').Create(
        fmask_fname, cs_final.shape[1], cs_final.shape[0], 1, gdal.GDT_Byte)
    c.SetGeoTransform(geoT)
    c.SetProjection(prj)
    c.GetRasterBand(1).WriteArray(cs_final)
    c = None

    # TODO: Save water/snow masks?
