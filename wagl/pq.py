#!/usr/bin/env python

"""
Contains the pixel quality workflow as well as a few utilities.
"""

from __future__ import absolute_import, print_function
import os
from os.path import join as pjoin, dirname
from glob import glob
import logging
import numpy
import h5py
import rasterio

from wagl.acca_cloud_masking import calc_acca_cloud_mask
from wagl.acquisition import acquisitions
from wagl.cloud_shadow_masking import cloud_shadow
from wagl.constants import BandType, DatasetName
from wagl.constants import PQAConstants, PQbits
from wagl.constants import ArdProducts as AP
from wagl.contiguity_masking import set_contiguity_bit
from wagl.fmask_cloud_masking_wrapper import fmask_cloud_mask
from wagl.hdf5 import dataset_compression_kwargs, write_h5_image
from wagl.land_sea_masking import set_land_sea_bit
from wagl.metadata import create_pq_yaml
from wagl.saturation_masking import set_saturation_bits
from wagl.temperature import get_landsat_temperature


def can_pq(level1, acq_parser_hint=None):
    """
    A simple test to check if we can process a scene through the
    pq pipeline.

    :param level1:
        An `str` containing the file path name to the directory
        containing the level-1 data.

    :return:
        True if the scene can be processed through PQ, else False.
    """
    supported = ['LANDSAT_5', 'LANDSAT_7', 'LANDSAT_8']
    acq = acquisitions(level1, acq_parser_hint).get_acquisitions()[0]
    return acq.platform_id in supported


class PQAResult(object):
    """
    Represents the PQA result
    """

    def __init__(self, shape, aGriddedGeoBox, dtype=numpy.uint16, aux_data=None):
        """
        Constructor

        Arguments:
            :shape:
                the shape of the numpy array holding the data
            :aGriddedGeoBox:
                an instance of class GriddedGeoBox providing the
                spatial location, scale and coordinate refernced system
                for this PQAResult
            :dtype:
                the datatype of the array
            :aux_data:
                a dictionary hold auxillary data associated with
                the PQAResult object. These may represent metadata
                elements that could be written to the final output
                file

        """
        assert shape is not None

        self.test_set = set()
        self.array = numpy.zeros(shape, dtype=dtype)
        self.dtype = dtype
        self.bitcount = self.array.itemsize * 8
        self.aux_data = aux_data or {}
        self.geobox = aGriddedGeoBox

    def set_mask(self, mask, bit_index, unset_bits=False):
        """Takes a boolean mask array and sets the bit in the result array.
        """
        assert mask.shape == self.array.shape, \
            "Mask shape %s does not match result array %s" % (mask.shape, self.array.shape)
        assert mask.dtype == numpy.bool, 'Mask must be of type bool'
        assert 0 <= bit_index < self.bitcount, 'Invalid bit index'
        assert bit_index not in self.test_set, 'Bit %d already set' % bit_index
        self.test_set.add(bit_index)

        c = sum(sum(mask))
        logging.debug('Setting result for bit %d, masking %d pixels', bit_index, c)
        numpy.bitwise_or(self.array, (mask << bit_index).astype(self.dtype),
                         self.array) # Set any 1 bits
        if unset_bits:
            numpy.bitwise_and(self.array,
                              ~(~mask << bit_index).astype(self.dtype),
                              self.array) # Clear any 0 bits

    def get_mask(self, bit_index):
        """
        Return boolean mask for specified bit index
        """
        assert 0 <= bit_index < self.bitcount, 'Invalid bit index'
        assert bit_index in self.test_set, 'Test %d not run' % bit_index
        return (self.array & (1 << bit_index)) > 0

    def add_to_aux_data(self, new_data=None):
        """
        Add the elements in the supplied dictionary to this objects
        aux_data property
        """
        self.aux_data.update(new_data or {})

    def save_as_tiff(self, path, crs=None):
        """
        Save the PQ result and attribute information in a GeoTiff.
        """
        os.makedirs(dirname(path))
        height, width = self.array.shape
        kwargs = {'driver': 'GTiff',
                  'width': width,
                  'height': height,
                  'count': 1,
                  'crs': self.geobox.crs.ExportToWkt(),
                  'transform': self.geobox.transform,
                  'dtype': 'uint16'}
        with rasterio.open(path, mode='w', **kwargs) as ds:
            ds.write_band(1, self.array)
            ds.update_tags(1, **self.aux_data)

    def save_as_h5_dataset(self, out_group, acq, product, compression):
        """
        Save the PQ result and attribute information in a HDF5
        `IMAGE` Class dataset.
        """
        kwargs = dataset_compression_kwargs(compression=compression,
                                            chunks=acq.tile_size)
        attrs = self.aux_data.copy()
        attrs['crs_wkt'] = self.geobox.crs.ExportToWkt()
        attrs['geotransform'] = self.geobox.transform.to_gdal()
        dname = DatasetName.PQ_FMT.value.format(product=product.value)
        write_h5_image(self.array, dname, out_group, attrs=attrs, **kwargs)

    @property
    def test_list(self):
        """Returns a sorted list of all bit indices which have been set
        """
        return sorted(self.test_set)

    @property
    def test_string(self):
        """Returns a string showing all bit indices which have been set
        """
        bit_list = ['0'] * self.bitcount
        for test_index in self.test_set:
            bit_list[test_index] = '1' # Show bits as big-endian
            # bit_list[15 - test_index] = '1' # Show bits as little-endian

        return ''.join(bit_list)


def _run_pq(level1, out_fname, scene_group, land_sea_path, compression, acq_parser_hint=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(out_fname) as fid:
        grp = fid[scene_group]
        run_pq(level1, grp, land_sea_path, grp, compression, AP.NBAR, acq_parser_hint)
        run_pq(level1, grp, land_sea_path, grp, compression, AP.NBART, acq_parser_hint)

def run_pq(level1, input_group, land_sea_path, out_group, compression='lzf',
           product=AP.NBAR, acq_parser_hint=None):
    """
    Runs the PQ workflow and saves the result in the same file as
    given by the `standardised_data_fname` parameter.

    :param level1:
        A `str` containing the file path name to the directory
        containing the level-1 data.

    :param input_group:
        The root HDF5 `Group` object containing the surface
        reflectance data that can be accessible via the enum
        specifier `constants.DatasetName.REFLECTANCE_FMT`.

    :param land_sea_path:
        A `str` containing the file path name to the directory
        containing the land/sea rasters.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.PQ_FMT

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param product:
        An enum representing the product type to use as input into the
        PQ algorithm. Default is `ArdProduct.nbar`.

    :return:
        None; the pixel quality result is stored in the same file
        as given by the `standardised_data_fname` parameter.
    """
    container = acquisitions(level1, acq_parser_hint)
    acqs = container.get_acquisitions()
    geo_box = acqs[0].gridded_geo_box()

    # filter out unwanted acquisitions
    acqs = [acq for acq in acqs if acq.band_type != BandType.PANCHROMATIC]

    platform_id = acqs[0].platform_id
    sensor = acqs[0].sensor_id

    # constants to be use for this PQA computation
    pq_const = PQAConstants(sensor)

    # get the band names based on acquisition description
    spectral_bands = []
    band_descriptions = ["Blue",
                         "Green",
                         "Red",
                         "NIR",
                         "SWIR 1",
                         "SWIR 2"]
    for band_desc in band_descriptions:
        spectral_bands.append([a.band_name for a in acqs if
                               a.desc == band_desc][0])


    # track the bits that have been set (tests that have been run)
    tests_run = {'band_1_saturated': False,
                 'band_2_saturated': False,
                 'band_3_saturated': False,
                 'band_4_saturated': False,
                 'band_5_saturated': False,
                 'band_6_saturated': False,
                 'band_7_saturated': False,
                 'contiguity': False,
                 'land_obs': False,
                 'cloud_acca': False,
                 'cloud_fmask': False,
                 'cloud_shadow_acca': False,
                 'cloud_shadow_fmask': False}

    # the PQAResult object for this run
    pqa_result = PQAResult(geo_box.shape, geo_box)

    # Saturation
    bits_set = set_saturation_bits(acqs, pq_const, pqa_result)
    for bit in bits_set:
        tests_run[PQbits(bit).name] = True

    # contiguity
    set_contiguity_bit(acqs, platform_id, pq_const, pqa_result)
    tests_run['contiguity'] = True

    # land/sea
    ancillary = set_land_sea_bit(geo_box, pq_const, pqa_result, land_sea_path)
    tests_run['land_obs'] = True

    contiguity_mask = (pqa_result.array & (1 << pq_const.contiguity)) > 0

    # fmask cloud mask
    if pq_const.run_cloud:
        aux_data = {}   # for collecting result metadata

        # TODO: pass in scene metadata via acquisitions or mtl reader
        mtl = glob(pjoin(level1, '*/*_MTL.txt'))[0]
        mask = fmask_cloud_mask(mtl, null_mask=contiguity_mask,
                                sat_tag=platform_id, aux_data=aux_data)

        # set the result
        pqa_result.set_mask(mask, pq_const.fmask)
        pqa_result.add_to_aux_data(aux_data)

        tests_run['cloud_fmask'] = True
    else:
        logging.warning('FMASK Not Run! %s sensor not configured for the '
                        'FMASK algorithm.', sensor)

    temperature = get_landsat_temperature(acqs, pq_const)

    # read NBAR data
    fmt = DatasetName.REFLECTANCE_FMT.value
    dname = fmt.format(product=product.value, band_name=spectral_bands[0])
    blue_dataset = input_group[dname]
    dname = fmt.format(product=product.value, band_name=spectral_bands[1])
    green_dataset = input_group[dname]
    dname = fmt.format(product=product.value, band_name=spectral_bands[2])
    red_dataset = input_group[dname]
    dname = fmt.format(product=product.value, band_name=spectral_bands[3])
    nir_dataset = input_group[dname]
    dname = fmt.format(product=product.value, band_name=spectral_bands[4])
    swir1_dataset = input_group[dname]
    dname = fmt.format(product=product.value, band_name=spectral_bands[5])
    swir2_dataset = input_group[dname]

    # acca cloud mask
    if pq_const.run_cloud:
        aux_data = {}   # for collecting result metadata
        mask = calc_acca_cloud_mask(blue_dataset, green_dataset,
                                    red_dataset, nir_dataset,
                                    swir1_dataset, swir2_dataset,
                                    temperature, pq_const,
                                    contiguity_mask, aux_data)

        # set the result
        pqa_result.set_mask(mask, pq_const.acca)
        pqa_result.add_to_aux_data(aux_data)

        tests_run['cloud_acca'] = True
    else:
        logging.warning('ACCA Not Run! %s sensor not configured for the '
                        'ACCA algorithm.', sensor)


    # parameters for cloud shadow masks
    land_sea_mask = pqa_result.get_mask(pq_const.land_sea)
    temperature = get_landsat_temperature(acqs, pq_const)

    # cloud shadow using the cloud mask generated by ACCA
    if pq_const.run_cloud_shadow:
        aux_data = {}   # for collecting result metadata

        cloud_mask = pqa_result.get_mask(pq_const.acca)
        sun_az_deg = acqs[0].sun_azimuth
        sun_elev_deg = acqs[0].sun_elevation

        mask = cloud_shadow(blue_dataset, green_dataset, red_dataset,
                            nir_dataset, swir1_dataset, swir2_dataset,
                            temperature, cloud_mask, geo_box, sun_az_deg,
                            sun_elev_deg, pq_const,
                            land_sea_mask=land_sea_mask,
                            contiguity_mask=contiguity_mask,
                            cloud_algorithm='ACCA',
                            growregion=True, aux_data=aux_data)

        pqa_result.set_mask(mask, pq_const.acca_shadow)
        pqa_result.add_to_aux_data(aux_data)

        tests_run['cloud_shadow_acca'] = True
    else: # OLI/TIRS only
        logging.warning('Cloud Shadow Algorithm Not Run! %s sensor not '
                        'configured for the cloud shadow '
                        'algorithm.', sensor)

    # cloud shadow using the cloud mask generated by FMASK
    if pq_const.run_cloud_shadow: # TM/ETM/OLI_TIRS
        aux_data = {}   # for collecting result metadata

        cloud_mask = pqa_result.get_mask(pq_const.fmask)
        sun_az_deg = acqs[0].sun_azimuth
        sun_elev_deg = acqs[0].sun_elevation

        mask = cloud_shadow(blue_dataset, green_dataset, red_dataset,
                            nir_dataset, swir1_dataset, swir2_dataset,
                            temperature, cloud_mask, geo_box, sun_az_deg,
                            sun_elev_deg, pq_const,
                            land_sea_mask=land_sea_mask,
                            contiguity_mask=contiguity_mask,
                            cloud_algorithm='FMASK',
                            growregion=True, aux_data=aux_data)

        pqa_result.set_mask(mask, pq_const.fmask_shadow)
        pqa_result.add_to_aux_data(aux_data)

        tests_run['cloud_shadow_fmask'] = True
    else: # OLI/TIRS only
        logging.warning('Cloud Shadow Algorithm Not Run! %s sensor not '
                        'configured for the cloud shadow '
                        'algorithm.', sensor)

    # write the pq result as an accompanying dataset to the standardised data
    pqa_result.save_as_h5_dataset(out_group, acqs[0], product, compression)
    # TODO: move metadata yaml creation outside of this func
    create_pq_yaml(acqs[0], ancillary, tests_run, out_group)
