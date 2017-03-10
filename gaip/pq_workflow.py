#!/bin/env python
# 
# Runs Pixel Quality workflow against specified input data
#
#

from __future__ import absolute_import
from datetime import datetime as dt
import os
from os.path import join as pjoin
from os.path import basename
import gc
from glob import glob
import logging
import subprocess
import yaml
from enum import Enum
import numpy
import h5py
import luigi
import gaip
from gaip.calculate_reflectance import DATASET_NAME_FMT
from gaip.nbar_workflow import TerrainCorrection


class PQbits(Enum):
    band_1_saturated = 0
    band_2_saturated = 1
    band_3_saturated = 2
    band_4_saturated = 3
    band_5_saturated = 4
    band_6_saturated = 5
    band_7_saturated = 6
    contiguity = 7
    land_obs = 8
    cloud_acca = 9
    cloud_fmask = 10
    cloud_shadow_acca = 11
    cloud_shadow_fmask = 12


class PixelQualityTask(luigi.Task):

    l1t_path = luigi.Parameter()
    work_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()
    land_sea_path = luigi.Parameter()
    compression = luigi.Parameter(default='lzf', significant=False)

    def output(self):
        return luigi.LocalTarget(pjoin(self.work_root, 'pq.h5'))

    def requires(self):
        return TerrainCorrection(self.level1, self.work_root, self.granule,
                                 self.group) 

    def run(self):
        # read L1T data
        logging.debug("Getting L1T acquisition data")
        l1t_acqs = gaip.acquisitions(self.l1t_path).get_acquisitions(group=self.group)

        # get the selected acquisitions and assciated band data and 
        # GriddedGeoBox. The latter provides the spatial context for the
        # band data
        l1t_acqs = [acq for acq in l1t_acqs if acq.band_type != gaip.PAN]
        l1t_data, geo_box = gaip.stack_data(l1t_acqs)

        spacecraft_id = l1t_acqs[0].spacecraft_id
        sensor = l1t_acqs[0].sensor_id

        # constants to be use for this PQA computation 
        pq_const = gaip.PQAConstants(sensor)
        nbar_const = gaip.constants.NBARConstants(spacecraft_id, sensor)
        avail_bands = nbar_const.get_nbar_lut()

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
        pqaResult = gaip.PQAResult(l1t_data[0].shape, geo_box)

        # Saturation
        bits_set = gaip.set_saturation_bits(l1t_data, pq_const, pqaResult)
        for bit in bits_set:
            tests_run[PQbits(bit).name] = True

        # contiguity
        gaip.set_contiguity_bit(l1t_data, spacecraft_id, pq_const, pqaResult)
        tests_run['contiguity'] = True

        # land/sea
        md = gaip.set_land_sea_bit(geo_box, pq_const, pqaResult,
                                   self.land_sea_path)
        tests_run['land_obs'] = True

        # get temperature data from thermal band in prepartion for cloud detection
        kelvin_band = gaip.get_landsat_temperature(l1t_data, l1t_acqs,
                                                   pq_const)
  
        contiguity_mask = (pqaResult.array & (1 << pq_const.contiguity)) > 0

        # we are done with L1T data so recover the memory
        del l1t_data
        gc.collect()
  
        # fmask cloud mask
        if pq_const.run_cloud:
            mask = None
            aux_data = {}   # for collecting result metadata
            
            # TODO: pass in scene metadata via Dale's new MTL reader
            mtl = glob(os.path.join(self.l1t_path, '*/*_MTL.txt'))[0]
            mask = gaip.fmask_cloud_mask(mtl, null_mask=contiguity_mask,
                                         sat_tag=spacecraft_id,
                                         aux_data=aux_data)

            # set the result
            pqaResult.set_mask(mask, pq_const.fmask)
            pqaResult.add_to_aux_data(aux_data)

            tests_run['cloud_fmask'] = True
        else:
            logging.warning(('FMASK Not Run! {} sensor not configured for the '
                             'FMASK algorithm.').format(sensor))
        
        # read NBAR data
        nbar_fid = h5py.File(self.input(), 'r')
        nbar_data = numpy.zeros(geo_box.shape, dtype='int16')
        for i, bn in enumerate(avail_bands):
            dataset_name = DATASET_NAME_FMT.format(product='brdf', band=bn)
            nbar_fid[dataset_name].read_direct(nbar_data[i])

        # acca cloud mask
        if pq_const.run_cloud:
            mask = None
            aux_data = {}   # for collecting result metadata
            if pq_const.oli_tirs:
                mask = gaip.calc_acca_cloud_mask(nbar_data[1:, :, :],
                                                 kelvin_band, pq_const,
                                                 contiguity_mask, aux_data)
            else: # TM or ETM
                mask = gaip.calc_acca_cloud_mask(nbar_data, kelvin_band,
                                                 pq_const, contiguity_mask,
                                                 aux_data)

            # set the result
            pqaResult.set_mask(mask, pq_const.acca)
            pqaResult.add_to_aux_data(aux_data)

            tests_run['cloud_acca'] = True
        else:
            logging.warning(('ACCA Not Run! {} sensor not configured for the '
                             'ACCA algorithm.').format(sensor))


        # parameters for cloud shadow masks
        land_sea_mask = pqaResult.get_mask(pq_const.land_sea)

        # acca cloud shadow
        if pq_const.run_cloud_shadow: # TM/ETM/OLI_TIRS
            mask = None
            aux_data = {}   # for collecting result metadata

            cloud_mask = pqaResult.get_mask(pq_const.acca)
            sun_az_deg = l1t_acqs[0].sun_azimuth
            sun_elev_deg = l1t_acqs[0].sun_elevation
            if pq_const.oli_tirs:
                mask = gaip.cloud_shadow(nbar_data[1 :, :, :], kelvin_band,
                                         cloud_mask, geo_box, sun_az_deg,
                                         sun_elev_deg, pq_const,
                                         land_sea_mask=land_sea_mask,
                                         contiguity_mask=contiguity_mask,
                                         cloud_algorithm='ACCA',
                                         growregion=True, aux_data=aux_data)
            else: # TM or ETM
                mask = gaip.cloud_shadow(nbar_data, kelvin_band, cloud_mask,
                                         geo_box, sun_az_deg, sun_elev_deg,
                                         pq_const, land_sea_mask=land_sea_mask,
                                         contiguity_mask=contiguity_mask,
                                         cloud_algorithm='ACCA',
                                         growregion=True, aux_data=aux_data)

            pqaResult.set_mask(mask, pq_const.acca_shadow)
            pqaResult.add_to_aux_data(aux_data)

            tests_run['cloud_shadow_acca'] = True
        else: # OLI/TIRS only
            logging.warning(('Cloud Shadow Algorithm Not Run! {} sensor not '
                             'configured for the cloud shadow '
                             'algorithm.').format(sensor))

        # FMASK cloud shadow
        if pq_const.run_cloud_shadow: # TM/ETM/OLI_TIRS
            mask = None
            aux_data = {}   # for collecting result metadata

            cloud_mask = pqaResult.get_mask(pq_const.fmask)
            sun_az_deg = l1t_acqs[0].sun_azimuth
            sun_elev_deg = l1t_acqs[0].sun_elevation
            if pq_const.oli_tirs:
                mask = gaip.cloud_shadow(nbar_data[1:, :, :], kelvin_band,
                                         cloud_mask, geo_box, sun_az_deg,
                                         sun_elev_deg, pq_const,
                                         land_sea_mask=land_sea_mask,
                                         contiguity_mask=contiguity_mask,
                                         cloud_algorithm='FMASK',
                                         growregion=True, aux_data=aux_data)
            else: # TM or ETM
                mask = gaip.cloud_shadow(nbar_data, kelvin_band, cloud_mask,
                                         geo_box, sun_az_deg, sun_elev_deg,
                                         pq_const, land_sea_mask=land_sea_mask,
                                         contiguity_mask=contiguity_mask,
                                         cloud_algorithm='FMASK',
                                         growregion=True, aux_data=aux_data)

            pqaResult.set_mask(mask, pq_const.fmask_shadow)
            pqaResult.add_to_aux_data(aux_data)

            tests_run['cloud_shadow_fmask'] = True
        else: # OLI/TIRS only
            logging.warning(('Cloud Shadow Algorithm Not Run! {} sensor not '
                             'configured for the cloud shadow '
                             'algorithm.').format(sensor))

        # metadata
        system_info = {}
        proc = subprocess.Popen(['uname', '-a'], stdout=subprocess.PIPE)
        system_info['node'] = proc.stdout.read()
        system_info['time_processed'] = dt.utcnow()

        source_info = {}
        source_info['source_l1t'] = self.l1t_path
        source_info['source_nbar'] = self.work_root

        algorithm = {}
        algorithm['software_version'] = gaip.get_version()
        algorithm['software_repository'] = ('https://github.com/'
                                            'GeoscienceAustralia/'
                                            'ga-neo-landsat-processor.git')
        algorithm['pq_doi'] = 'http://dx.doi.org/10.1109/IGARSS.2013.6723746'
        
        metadata = {}
        metadata['system_information'] = system_info
        metadata['source_data'] = source_info
        metadata['algorithm_information'] = algorithm
        metadata['ancillary'] = md
        metadata['tests_run'] = tests_run

        with open(pjoin(self.work_root, "pq_metadata.yaml"), 'w') as src:
            yaml.dump(metadata, src, default_flow_style=False)

        # write PQA file as output
        with self.output().temporary_path() as out_fname:
            pqaResult.save_as_h5_dataset(out_fname, self.compression)


class PQ(luigi.WrapperTask):

    """Kicks off NBAR tasks for each level1 entry."""

    level1_csv = luigi.Parameter()
    output_directory = luigi.Parameter()
    work_extension = luigi.Parameter(default='.gaip-work', significant=False)

    def requires(self):
        with open(self.level1_csv) as src:
            level1_scenes = src.readlines()

        for scene in level1_scenes:
            work_name = basename(scene) + self.work_extension
            work_root = pjoin(self.output_directory, work_name)
            container = gaip.acquisitions(scene)
            for granule in container.granules:
                for group in container.groups:
                    yield PixelQualityTask(scene, work_root, granule, group)

        
if __name__ == '__main__':
    luigi.run()
