#!/bin/env python
# 
# Runs Pixel Quality workflow against specified input data
#
#

import argparse
import luigi.contrib.mpi as mpi
import luigi
import logging
import os
import sys
import gc
import re
import gaip


from glob import glob


L1T_PATTERN = '(?P<spacecraft_id>LS\d)_(?P<sensor_id>\w+)_(?P<product_type>\w+)' \
    '_(?P<product_id>P\d+)_GA(?P<product_code>.*)-(?P<station_id>\d+)_' \
    '(?P<wrs_path>\d+)_(?P<wrs_row>\d+)_(?P<acquisition_date>\d{8})'
PAT = re.compile(L1T_PATTERN)


def nbar_name_from_l1t(l1t_fname):
    """
    Return an NBAR file name given a L1T file name or None if 
    invalid L1T name.
    """
    m = PAT.match(l1t_fname)
    if m :
        sensor_id =  m.group('sensor_id')
        sensor_id = sensor_id.replace('OLITIRS', 'OLI_TIRS')
        return "%s_%s_NBAR_P54_GANBAR01-%s_%s_%s_%s" % ( \
            m.group('spacecraft_id'),
            sensor_id,
            m.group('station_id'),
            m.group('wrs_path'),
            m.group('wrs_row'),
            m.group('acquisition_date'))

def pqa_name_from_l1t(l1t_fname):
    """
    Return a PQA file name given a L1T file name or None if 
    invalid L1T name.
    """
    m = PAT.match(l1t_fname)
    if m :
        sensor_id =  m.group('sensor_id')
        sensor_id = sensor_id.replace('OLITIRS', 'OLI_TIRS')
        return "%s_%s_PQ_P55_GAPQ01-%s_%s_%s_%s" % ( \
            m.group('spacecraft_id'),
            sensor_id,
            m.group('station_id'),
            m.group('wrs_path'),
            m.group('wrs_row'),
            m.group('acquisition_date'))


class PixelQualityTask(luigi.Task):

    l1t_path = luigi.Parameter()
    nbar_path = luigi.Parameter()
    land_sea_path = luigi.Parameter()
    pq_path = luigi.Parameter()

    def output(self):
        return PQDataset(self.pq_path)

    def requires(self):
        return NBARTask(self.nbar_path) 

    def run(self):
        logging.info("In PixelQualityTask.run method, L1T=%s NBAR=%s, output=%s" %\
           (self.l1t_path, self.input().nbar_path, self.output().path))

        # read L1T data
        logging.debug("Getting L1T acquisition data")
        l1t_acqs = gaip.acquisitions(self.l1t_path)

        # get the selected acquisitions and assciated band data and 
        # GriddedGeoBox. The latter provides the spatial context for the
        # band data

        l1t_acqs, l1t_data, geo_box = gaip.stack_data(l1t_acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN))
        logging.debug("l1t_data shape=%s, geo_box= %s" % \
            (str(l1t_data.shape), str(geo_box)))


        spacecraft_id = l1t_acqs[0].spacecraft_id
        sensor = l1t_acqs[0].sensor_id
        logging.debug("Satellite is %s, sensor is %s" % (spacecraft_id, sensor ))

        # constants to be use for this PQA computation 

        logging.debug("setting constants for sensor=%s" % (sensor, ))
        pq_const = gaip.PQAConstants(sensor)

        # the PQAResult object for this run

        pqaResult = gaip.PQAResult(l1t_data[0].shape, geo_box)

        # Saturation

        logging.debug("setting saturation bits")
        gaip.setSaturationBits(l1t_data, pq_const, pqaResult)
        logging.debug("done setting saturation bits")

        # contiguity

        logging.debug("setting contiguity bit")
        gaip.setContiguityBit(l1t_data, spacecraft_id, pq_const, pqaResult)
        logging.debug("done setting contiguity bit")

        # land/sea

        logging.debug("setting land/sea bit")
 #       affine = geo_box.affine
        gaip.setLandSeaBit(geo_box, pq_const, pqaResult, self.land_sea_path)
        logging.debug("done setting land/sea bit")

        # get temperature data from thermal band in prepartion for cloud detection

        logging.debug("calculating kelvin band")
        kelvin_band = gaip.get_landsat_temperature(l1t_data, l1t_acqs, pq_const)
  
        contiguity_mask = (pqaResult.array & (1 << pq_const.contiguity)) > 0

        # we are done with L1T data so recover the memory
   
        del l1t_data
        gc.collect()
  
        # fmask cloud mask

        logging.debug("calculating fmask cloud mask")
        if pq_const.run_cloud:
            mask = None
            aux_data = {}   # for collecting result metadata
            
            # TODO: pass in scene metadata via Dale's new MTL reader
            mtl = glob(os.path.join(self.l1t_path, \
                'scene01/*_MTL.txt'))[0] # Crude but effective
            mask = gaip.FMaskCloudMask(mtl, null_mask=contiguity_mask, sat_tag=spacecraft_id, \
                aux_data=aux_data)

            # set the result
            pqaResult.set_mask(mask, pq_const.fmask)
            pqaResult.add_to_aux_data(aux_data)
        else:
            logging.warning('FMASK Not Run! %s sensor not configured for the FMASK algorithm.' \
               % (sensor, ))
        
        logging.debug("done calculating fmask cloud mask")

        # read NBAR data
        logging.debug("Getting NBAR acquisition data")
        nbar_acqs = gaip.acquisitions(self.nbar_path)

        # get the selected acquisitions and associated band data and 
        # GriddedGeoBox. The latter provides the spatial context for the
        # band data

        nbar_acqs, nbar_data, geo_box = gaip.stack_data(nbar_acqs, \
            filter=(lambda acq: True))
        logging.debug("nbar_data shape=%s, geo_box= %s" % \
            (str(nbar_data.shape), str(geo_box)))

        # acca cloud mask

        logging.debug("calculating acca cloud mask")
        if pq_const.run_cloud:
            mask = None
            aux_data = {}   # for collecting result metadata
            if pq_const.oli_tirs:
                mask = gaip.calc_acca_cloud_mask(nbar_data[1:,:,:], kelvin_band, pq_const, \
                    contiguity_mask, aux_data)
            else: # TM or ETM
                mask = gaip.calc_acca_cloud_mask(nbar_data, kelvin_band, pq_const, \
                    contiguity_mask, aux_data)

            # set the result
            pqaResult.set_mask(mask, pq_const.acca)
            pqaResult.add_to_aux_data(aux_data)
        else:
            logging.warning('ACCA Not Run! %s sensor not configured for the ACCA algorithm.' \
               % (sensor, ))


        # parameters for cloud shadow masks

        land_sea_mask = pqaResult.get_mask(pq_const.land_sea)

        # acca cloud shadow

        logging.debug("calculating ACCA cloud shadow mask")
        if pq_const.run_cloud_shadow: # TM/ETM/OLI_TIRS
            mask = None
            aux_data = {}   # for collecting result metadata

            cloud_mask = pqaResult.get_mask(pq_const.acca)
            sun_az_deg = l1t_acqs[0].sun_azimuth
            sun_elev_deg = l1t_acqs[0].sun_elevation
            if pq_const.oli_tirs:
                mask = gaip.Cloud_Shadow(nbar_data[1:,:,:], kelvin_band, cloud_mask, geo_box, \
                    sun_az_deg, sun_elev_deg, pq_const, \
                    land_sea_mask=land_sea_mask, contiguity_mask=contiguity_mask, \
                    cloud_algorithm='ACCA', growregion=True, aux_data=aux_data)
            else: # TM or ETM
                mask = gaip.Cloud_Shadow(nbar_data, kelvin_band, cloud_mask, geo_box, \
                    sun_az_deg, sun_elev_deg, pq_const, \
                    land_sea_mask=land_sea_mask, contiguity_mask=contiguity_mask, \
                    cloud_algorithm='ACCA', growregion=True, aux_data=aux_data)

            pqaResult.set_mask(mask, pq_const.acca_shadow)
            pqaResult.add_to_aux_data(aux_data)

        else: # OLI/TIRS only
            logger.warning('Cloud Shadow Algorithm Not Run! %s sensor not configured for the cloud shadow algorithm.' \
                % (sensor, ))

        logging.debug("done calculating ACCA cloud shadow mask")

        # FMASK cloud shadow

        logging.debug("calculating FMASK cloud shadow mask")
        if pq_const.run_cloud_shadow: # TM/ETM/OLI_TIRS
            mask = None
            aux_data = {}   # for collecting result metadata

            cloud_mask = pqaResult.get_mask(pq_const.fmask)
            sun_az_deg = l1t_acqs[0].sun_azimuth
            sun_elev_deg = l1t_acqs[0].sun_elevation
            if pq_const.oli_tirs:
                mask = gaip.Cloud_Shadow(nbar_data[1:,:,:], kelvin_band, cloud_mask, geo_box, \
                    sun_az_deg, sun_elev_deg, pq_const, \
                    land_sea_mask=land_sea_mask, contiguity_mask=contiguity_mask,
                    cloud_algorithm='FMASK', growregion=True, aux_data=aux_data)
            else: # TM or ETM
                mask = gaip.Cloud_Shadow(nbar_data, kelvin_band, cloud_mask, geo_box, \
                    sun_az_deg, sun_elev_deg, pq_const, \
                    land_sea_mask=land_sea_mask, contiguity_mask=contiguity_mask,
                    cloud_algorithm='FMASK', growregion=True, aux_data=aux_data)

            pqaResult.set_mask(mask, pq_const.fmask_shadow)
            pqaResult.add_to_aux_data(aux_data)

        else: # OLI/TIRS only
            logger.warning('Cloud Shadow Algorithm Not Run! %s sensor not configured for the cloud shadow algorithm.' \
                % (sensor, ))

        logging.debug("done calculating FMASK cloud shadow mask")

        # write PQA file as output

        logging.debug("saving PQA result GeoTiff")
        pqa_output_path = os.path.join(self.pq_path, "pqa.tif")
        pqaResult.save_as_tiff(pqa_output_path)
        logging.debug("done saving PQA result GeoTiff")
            

class PQDataset(luigi.Target):

    def __init__(self, path):
        self.path = path

    def exists(self):
        return os.path.exists(self.path)

class NBARTask(luigi.ExternalTask):
    nbar_path = luigi.Parameter()

    def output(self):
        return NBARdataset(self.nbar_path)


class NBARdataset(luigi.Target):

    def __init__(self, nbar_path):
        self.nbar_path = nbar_path
        self.acquisitions = gaip.acquisitions(self.nbar_path)

    def exists(self):
        return os.path.exists(self.nbar_path)

def is_valid_directory(parser, arg):
    """Used by argparse"""
    if not os.path.exists(arg):
        parser.error("%s does not exist" % (arg, ))
    else:
        return arg

def nbar_name_from(l1t_fname):
    """
    Return an NBAR file name given a L1T file name
    """

if __name__ == '__main__':

    # command line arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("--l1t_path", help="path to directory containing OTH datasets", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--nbar_path", help="path to directory containing NBAR datasets", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--land_sea_path", help="path to directory containing Land/Sea datasets", \
        default='/g/data/v10/eoancillarydata/Land_Sea_Rasters', \
        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--out_path", help="path to directory where PQA dataset is to be written", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--log_path", help="path to directory where where log files will be written", \
        default='.', type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--debug", help="selects more detail logging (default is INFO)", \
        default=False, action='store_true')

    args = parser.parse_args()

    # setup logging
    
    logfile = "%s/run_pq_%s_%d.log" % (args.log_path, os.uname()[1], os.getpid())
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(filename=logfile, level=logging_level, \
        format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info("run_pq started")


    logging.info('l1t_path=%s' % (args.l1t_path, ))
    logging.info('nbar_path=%s' % (args.nbar_path, ))
    logging.info('land_sea_path=%s' % (args.land_sea_path, ))
    logging.info('out_path=%s' % (args.out_path, ))
    logging.info('log_path=%s' % (args.log_path, ))

    # create the task list based on L1T files to processa

    tasks = []
    for l1t_file in [f for f in os.listdir(args.l1t_path) if '_OTH_' in f]:
        l1t_dataset_path = os.path.join(args.l1t_path, l1t_file)
        nbar_dataset_path = os.path.join(args.nbar_path, nbar_name_from_l1t(l1t_file))
        pqa_dataset_path = os.path.join(args.out_path, pqa_name_from_l1t(l1t_file))

        tasks.append( \
            PixelQualityTask( \
                l1t_dataset_path,
                nbar_dataset_path,
                args.land_sea_path,
                pqa_dataset_path \
            ) \
        )
        
        print nbar_dataset_path

    mpi.run(tasks)

