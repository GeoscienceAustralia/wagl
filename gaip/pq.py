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
import gaip
from gaip.acca_cloud_masking import calc_acca_cloud_mask
from gaip.acquisition import acquisitions
from gaip.cloud_shadow_masking import cloud_shadow
from gaip.constants import BandType, DatasetName, NBARConstants
from gaip.constants import PQAConstants, PQbits
from gaip.contiguity_masking import set_contiguity_bit
from gaip.fmask_cloud_masking_wrapper import fmask_cloud_mask
from gaip.hdf5 import dataset_compression_kwargs, write_h5_image, write_scalar
from gaip.land_sea_masking import set_land_sea_bit
from gaip.metadata import create_pq_yaml
from gaip.saturation_masking import set_saturation_bits
from gaip.temperature import get_landsat_temperature


def can_pq(level1):
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
    acq = acquisitions(level1).get_acquisitions()[0]
    return acq.spacecraft_id in supported


class PQAResult(object):
    """
    Represents the PQA result
    """

    def __init__(self, shape, aGriddedGeoBox, dtype=numpy.uint16, aux_data={}):
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
        self.aux_data = aux_data
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
        logging.debug('Setting result for bit %d, masking %d pixels' % (bit_index, c))
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

    def add_to_aux_data(self, new_data={}):
        """
        Add the elements in the supplied dictionary to this objects
        aux_data property
        """
        self.aux_data.update(new_data)

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

    def save_as_h5_dataset(self, out_fname, compression):
        """
        Save the PQ result and attribute information in a HDF5
        `IMAGE` Class dataset.
        """
        with h5py.File(out_fname) as fid:
            chunks = (1, self.geobox.x_size())
            kwargs = dataset_compression_kwargs(compression=compression,
                                                chunks=chunks)
            attrs = self.aux_data.copy()
            attrs['crs_wkt'] = self.geobox.crs.ExportToWkt()
            attrs['geotransform'] = self.geobox.transform.to_gdal()
            dname = DatasetName.pixel_quality.value
            write_h5_image(self.array, dname, fid, attrs=attrs, **kwargs)
       
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


def run_pq(level1, standardised_data_fname, land_sea_path, compression='lzf'):
    """
    Runs the PQ workflow and saves the result in the same file as
    given by the `standardised_data_fname` parameter.

    :param level1:
        A `str` containing the file path name to the directory
        containing the level-1 data.

    :param standardised_data_fname:
        A `str` containing the file path name to the file containing
        the standardised data (surface reflectance and surface
        brightness temperature).

    :param land_sea_path:
        A `str` containing the file path name to the directory
        containing the land/sea rasters.

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :return:
        None; the pixel quality result is stored in the same file
        as given by the `standardised_data_fname' parameter.
    """
    container = acquisitions(level1)
    acqs = container.get_acquisitions()
    geo_box = acqs[0].gridded_geo_box()

    # filter out unwanted acquisitions
    acqs = [acq for acq in acqs if acq.band_type != BandType.Panchromatic]

    spacecraft_id = acqs[0].spacecraft_id
    sensor = acqs[0].sensor_id

    # constants to be use for this PQA computation 
    pq_const = PQAConstants(sensor)
    nbar_const = NBARConstants(spacecraft_id, sensor)
    avail_bands = nbar_const.get_nbar_lut()

    # TODO: better method of band number access
    spectral_bands = avail_bands[1:] if pq_const.oli_tirs else avail_bands

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
    set_contiguity_bit(acqs, spacecraft_id, pq_const, pqa_result)
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
                                sat_tag=spacecraft_id, aux_data=aux_data)

        # set the result
        pqa_result.set_mask(mask, pq_const.fmask)
        pqa_result.add_to_aux_data(aux_data)

        tests_run['cloud_fmask'] = True
    else:
        logging.warning(('FMASK Not Run! {} sensor not configured for the '
                         'FMASK algorithm.').format(sensor))
    
    temperature = get_landsat_temperature(acqs, pq_const)

    # read NBAR data
    dname_fmt = DatasetName.reflectance_fmt.value
    with h5py.File(standardised_data_fname, 'r') as fid:
        dname = dname_fmt.format(product='brdf', band=spectral_bands[0])
        blue_dataset = fid[dname]
        dname = dname_fmt.format(product='brdf', band=spectral_bands[1])
        green_dataset = fid[dname]
        dname = dname_fmt.format(product='brdf', band=spectral_bands[2])
        red_dataset = fid[dname]
        dname = dname_fmt.format(product='brdf', band=spectral_bands[3])
        nir_dataset = fid[dname]
        dname = dname_fmt.format(product='brdf', band=spectral_bands[4])
        swir1_dataset = fid[dname]
        dname = dname_fmt.format(product='brdf', band=spectral_bands[5])
        swir2_dataset = fid[dname]

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
            logging.warning(('ACCA Not Run! {} sensor not configured for the '
                             'ACCA algorithm.').format(sensor))


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
            logging.warning(('Cloud Shadow Algorithm Not Run! {} sensor not '
                             'configured for the cloud shadow '
                             'algorithm.').format(sensor))

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
            logging.warning(('Cloud Shadow Algorithm Not Run! {} sensor not '
                             'configured for the cloud shadow '
                             'algorithm.').format(sensor))

    # write the pq result as an accompanying dataset to the standardised data
    pqa_result.save_as_h5_dataset(standardised_data_fname, compression)

    with h5py.File(standardised_data_fname) as out_fid:
        create_pq_yaml(acqs[0], ancillary, tests_run, out_fid)
