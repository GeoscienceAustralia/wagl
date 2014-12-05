import commands
import datetime
import os
import re

from gaip import constants
from gaip import BRDFLoader
from gaip import get_brdf_dirs_modis
from gaip import get_brdf_dirs_pre_modis
from gaip import read_subset
from gaip import write_img


def findFile(fileList, bandWL, factor):
    for file in fileList:
        if file.find(bandWL) != -1 and file.find(factor) != -1:
            return file
    return None


def get_brdf_data(acquisition, brdf_primary_path, brdf_secondary_path,
                  work_path):
    """
    Calculates the mean BRDF value for each band wavelength of your
    sensor, for each BRDF factor ['geo', 'iso', 'vol'] that covers
    your image extents.

    :param acquisition:
        An instance of an acquisitions object.

    :param brdf_primary_path:
        A string containing the full file system path to your directory
        containing the source BRDF files.  The BRDF directories are
        assumed to be yyyy.mm.dd naming convention.

    :param brdf_secondary_path:
        A string containing the full file system path to your directory
        containing the Jupp-Li backup BRDF data.  To be used for
        pre-MODIS and potentially post-MODIS acquisitions.

    :param work_path:
        A string containing the full file system path to your NBAR
        working directory. Intermediate BRDF files will be saved to
        work_path/brdf_intermediates/.

    :return:
        A dictionary with tuple (band, factor) as the keys. Each key
        represents the band of your satllite/sensor and brdf factor.
        Each key contains a dictionary with the following keys:
        data_source -> BRDF
        data_file -> File system path to the location of the selected
            BRDF wavelength and factor combination.
        value -> The mean BRDF value covering your image extents.
    """
    # Retrieve the satellite and sensor for the acquisition
    satellite = acquisition.spacecraft_id
    sensor = acquisition.sensor_id

    # Get the required BRDF LUT & factors list
    nbar_constants = constants.NBARConstants(satellite, sensor)

    brdf_lut = nbar_constants.getBRDFlut()
    brdf_factors = nbar_constants.getBRDFfactors()

    # Compute the geobox
    geobox = acquisition.gridded_geo_box()

    # Get the date of acquisition
    dt = acquisition.scene_center_datetime.date()

    # Get the boundary extents of the image
    # Each is a co-ordinate pair of (x, y)
    UL_Lon = geobox.ul_lonlat[0]
    UL_Lat = geobox.ul_lonlat[1]
    UR_Lon = geobox.ur_lonlat[0]
    UR_Lat = geobox.ur_lonlat[1]
    LR_Lon = geobox.lr_lonlat[0]
    LR_Lat = geobox.lr_lonlat[1]
    LL_Lon = geobox.ll_lonlat[0]
    LL_Lat = geobox.ll_lonlat[1]

    # Use maximal axis-aligned extents for BRDF mean value calculation.
    # Note that latitude min-max logic is valid for the Southern
    # hemisphere only.
    nw = (min(UL_Lon, LL_Lon), max(UL_Lat, UR_Lat))
    se = (max(LR_Lon, UR_Lon), min(LL_Lat, LR_Lat))

    # Compare the scene date and MODIS BRDF start date to select the
    # BRDF data root directory.
    # Scene dates outside the range of the CSIRO mosaic data
    # (currently 2000-02-18 through 2013-01-09) should use the pre-MODIS,
    # Jupp-Li BRDF.
    brdf_dir_list = sorted(os.listdir(brdf_primary_path))
    brdf_dir_range = [brdf_dir_list[0], brdf_dir_list[-1]]
    brdf_range = [datetime.date(*[int(x) for x in y.split('.')])
                  for y in brdf_dir_range]

    use_JuppLi_brdf = (dt < brdf_range[0] or dt > brdf_range[1])

    if use_JuppLi_brdf:
        brdf_base_dir = brdf_secondary_path
        brdf_dirs = get_brdf_dirs_pre_modis(brdf_base_dir, dt)
    else:
        brdf_base_dir = brdf_primary_path
        brdf_dirs = get_brdf_dirs_modis(brdf_base_dir, dt)

    # The following hdfList code was resurrected from the old SVN repo. JS
    # get all HDF files in the input dir
    dbDir = os.path.join(brdf_base_dir, brdf_dirs)
    three_tup = os.walk(dbDir)
    hdfList = []
    for (hdfHome, dirlist, filelist) in three_tup:
        for file in filelist:
            if file.endswith(".hdf.gz") or file.endswith(".hdf"):
                hdfList.append(file)

    # Initialise the brdf dictionary to store the results
    brdf_dict = {}

    # Create a BRDF directory in the work path to store the intermediate
    # files such as format conversion and subsets.
    brdf_out_path = os.path.join(work_path, 'brdf_intermediates')
    if not os.path.exists(brdf_out_path):
        os.makedirs(brdf_out_path)

    # Loop over each defined band and each BRDF factor
    for band in brdf_lut.keys():
        bandwl = brdf_lut[band]  # Band wavelength
        for factor in brdf_factors:
            hdfFileName = findFile(hdfList, bandwl, factor)

            hdfFile = os.path.join(hdfHome, hdfFileName)

            # Test if the file exists and has correct permissions
            try:
                with open(hdfFile, 'rb') as f:
                    pass
            except IOError as e:
                print "Unable to open file %s" % hdfFile

            # Unzip if we need to
            if hdfFile.endswith(".hdf.gz"):
                hdf_file = os.path.join(
                    work_path,
                    re.sub(".hdf.gz", ".hdf",
                           os.path.basename(hdfFile)))
                gunzipCmd = "gunzip -c %s > %s" % (hdfFile, hdf_file)
                (status, msg) = commands.getstatusoutput(gunzipCmd)
                assert status == 0, "gunzip failed: %s" % msg
            else:
                hdf_file = hdfFile

            # the following now converts the file format and outputs a subset.
            # this should proove useful for debugging and testing.

            # Load the file
            brdf_object = BRDFLoader(hdf_file, UL=nw, LR=se)

            # setup the output filename
            out_fname = '_'.join(['Band', str(band), bandwl, factor])
            out_fname = os.path.join(brdf_out_path, out_fname)

            # Convert the file format
            brdf_object.convert_format(out_fname)

            print 'x'

            # Read the subset and geotransform that corresponds to the subset
            subset, geobox_subset = read_subset(out_fname, (UL_Lat, UL_Lat),
                                                (UR_Lon, UR_Lat),
                                                (LR_Lon, LR_Lat),
                                                (LL_Lon, LL_Lat))

            print 'x'

            # The brdf_object has the scale and offsets so calculate the mean
            # through the brdf_object
            brdf_mean_value = brdf_object.get_mean(subset)

            print 'x'

            # Output the brdf subset
            out_fname_subset = out_fname + '_subset'
            write_img(subset, out_fname_subset, geobox=geobox_subset)

            # Remove temporary unzipped file
            if hdf_file.find(work_path) == 0:
                os.remove(hdf_file)

            # Add the brdf filename and mean value to brdf_dict
            brdf_dict[(band, factor)] = {'data_source': 'BRDF',
                                         'data_file': hdfFile,
                                         'value': brdf_mean_value}

    return brdf_dict
