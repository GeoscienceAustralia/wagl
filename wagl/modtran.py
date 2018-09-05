#!/usr/bin/env python

"""
MODTRAN drivers
---------------

"""

from __future__ import absolute_import, print_function
import os
from os.path import join as pjoin, exists, dirname
import subprocess
import glob
from posixpath import join as ppjoin
import numpy
import h5py
import pandas as pd
import json
import fnmatch
from wagl.constants import Workflow, BandType, DatasetName, GroupName, Albedos
from wagl.constants import POINT_FMT, ALBEDO_FMT, POINT_ALBEDO_FMT
from wagl.constants import AtmosphericCoefficients as AC
from wagl.hdf5 import write_dataframe, read_h5_table, create_external_link
from wagl.hdf5 import VLEN_STRING, write_scalar, H5CompressionFilter
import wagl.modtran_profile_json as mpjson


class JsonEncoder(json.JSONEncoder):
    """
    A wrapper class to address the issue of json encoding error
    This class handles  the json serializing error for numpy 
    datatype: 'float32' and numpy arrays
    """
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            raise OSError('Json Encoding Error')

def prepare_modtran(acquisitions, coordinate, albedos, basedir, modtran_exe):
    """
    Prepares the working directory for a MODTRAN execution.
    """

    
    data_dir = pjoin('/g/data/v10/private/modules/MODTRAN/MODTRAN-6.0.1/MODTRAN6.0','DATA')

    if not exists(data_dir):
        raise OSError('Cannot find MODTRAN')

    point_dir = pjoin(basedir, POINT_FMT.format(p=coordinate))

    for albedo in albedos:
        if albedo == Albedos.ALBEDO_TH:
            band_type = BandType.THERMAL
        else:
            band_type = BandType.REFLECTIVE

        acq = [acq for acq in acquisitions if acq.band_type == band_type][0]

        modtran_work = pjoin(point_dir, ALBEDO_FMT.format(a=albedo.value))

        if not exists(modtran_work):
            os.makedirs(modtran_work)

        symlink_dir = pjoin(modtran_work, 'DATA')

        if exists(symlink_dir):
            os.unlink(symlink_dir)

        os.symlink(data_dir, symlink_dir)

        out_fname = pjoin(modtran_work, acq.spectral_filter_file)  # filter file name needs to be put in DATA directory
        print(out_fname)
        response = acq.spectral_response(as_list=True)

        with open(out_fname, 'wb') as src:
            src.writelines(response)


def _format_json(acquisitions, satellite_solar_angles_fname,
                longitude_latitude_fname, ancillary_fname, out_fname, workflow):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_angles_fname, 'r') as sat_sol_fid,\
        h5py.File(longitude_latitude_fname, 'r') as lon_lat_fid,\
        h5py.File(ancillary_fname, 'r') as anc_fid,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = anc_fid[GroupName.ANCILLARY_GROUP.value]
        grp2 = sat_sol_fid[GroupName.SAT_SOL_GROUP.value]
        grp3 = lon_lat_fid[GroupName.LON_LAT_GROUP.value]
        json_data, _ = format_json(acquisitions, grp1, grp2, grp3, workflow, fid)

    return json_data


def format_json(acquisitions, ancillary_group, satellite_solar_group,
               lon_lat_group, workflow, out_group):
    """
    Creates json files for the albedo (0) and thermal
    """
    # angles data
    sat_view = satellite_solar_group[DatasetName.SATELLITE_VIEW.value]
    sat_azi = satellite_solar_group[DatasetName.SATELLITE_AZIMUTH.value]
    longitude = lon_lat_group[DatasetName.LON.value]
    latitude = lon_lat_group[DatasetName.LAT.value]

    # retrieve the averaged ancillary if available
    anc_grp = ancillary_group.get(GroupName.ANCILLARY_AVG_GROUP.value)
    if anc_grp is None:
        anc_grp = ancillary_group

    # ancillary data
    coordinator = ancillary_group[DatasetName.COORDINATOR.value]
    aerosol = anc_grp[DatasetName.AEROSOL.value][()]
    water_vapour = anc_grp[DatasetName.WATER_VAPOUR.value][()]
    ozone = anc_grp[DatasetName.OZONE.value][()]
    elevation = anc_grp[DatasetName.ELEVATION.value][()]

    npoints = coordinator.shape[0]
    view = numpy.zeros(npoints, dtype='float32')
    azi = numpy.zeros(npoints, dtype='float32')
    lat = numpy.zeros(npoints, dtype='float64')
    lon = numpy.zeros(npoints, dtype='float64')

    for i in range(npoints):
        yidx = coordinator['row_index'][i]
        xidx = coordinator['col_index'][i]
        view[i] = sat_view[yidx, xidx]
        azi[i] = sat_azi[yidx, xidx]
        lat[i] = latitude[yidx, xidx]
        lon[i] = longitude[yidx, xidx]

    view_corrected = 180 - view
    azi_corrected = azi + 180
    rlon = 360 - lon

    # check if in western hemisphere
    idx = rlon >= 360
    rlon[idx] -= 360

    idx = (180 - view_corrected) < 0.1
    view_corrected[idx] = 180
    azi_corrected[idx] = 0

    idx = azi_corrected > 360
    azi_corrected[idx] -= 360

    # get the modtran profiles to use based on the centre latitude
    _, centre_lat = acquisitions[0].gridded_geo_box().centre_lonlat

    if out_group is None:
        out_group = h5py.File('atmospheric-inputs.h5', 'w')

    if GroupName.ATMOSPHERIC_INPUTS_GRP.value not in out_group:
        out_group.create_group(GroupName.ATMOSPHERIC_INPUTS_GRP.value)

    group = out_group[GroupName.ATMOSPHERIC_INPUTS_GRP.value]
    iso_time = acquisitions[0].acquisition_datetime.isoformat()
    group.attrs['acquisition-datetime'] = iso_time

    json_data = {}
    # setup the json files required by MODTRAN
    if workflow == Workflow.STANDARD or workflow == Workflow.NBAR:
        acqs = [a for a in acquisitions if a.band_type == BandType.REFLECTIVE]

        for p in range(npoints):

            for alb in Workflow.NBAR.albedos:
                name = 'POINT-%i-ALBEDO-%s' %(p,str(alb.value))
                water = water_vapour
                ozone = ozone
                filter_function = acqs[0].spectral_filter_file
                visibility = -aerosol
                elevation = elevation
                sat_height = acquisitions[0].altitude / 1000.0
                sat_view = view_corrected[p]
                doy =  acquisitions[0].julian_day()
                binary = False
                input_data = {'water':water,
                              'ozone':ozone,
                              'filter_function':filter_function,
                              'visibility':visibility,
                              'elevation': elevation,
                              'sat_height':sat_height,
                              'sat_view':sat_view,
                              'doy':doy,
                              'description':'Input file for MODTRAN',
                              'file_format':'json'}
                
                albedo = float(alb.value)
                lats = lat[p]
                lons = rlon[p]
                time = acquisitions[0].decimal_hour()
                sat_azimuth = azi_corrected[p]
                input_data['albedo'] = albedo
                input_data['lat'] = lats
                input_data['lon'] = lons
                input_data['time'] = time
                input_data['sat_azimuth'] = sat_azimuth

                if centre_lat < -23.0:
                    data = mpjson.midlat_summer_albedo(name, water, ozone, visibility, doy, lats, lons, time,
                                                       sat_azimuth, elevation, sat_height, sat_view, albedo,
                                                       filter_function,binary)
                else:
                    data = mpjson.tropical_albedo(name, water, ozone, visibility, doy, lats, lons, time, sat_azimuth,
                                                  sat_height, elevation, sat_view, albedo, filter_function, binary)

                json_data[(p, alb)] = data

                data = json.dumps(data, cls=JsonEncoder, indent=4)
                dname = ppjoin(POINT_FMT.format(p=p),
                               ALBEDO_FMT.format(a=alb.value),
                               DatasetName.MODTRAN_INPUT.value)

                write_scalar(numpy.string_(data), dname, group, input_data)

    # create json for sbt if it has been collected
    if ancillary_group.attrs.get('sbt-ancillary'):
        dname = ppjoin(POINT_FMT, DatasetName.ATMOSPHERIC_PROFILE.value)
        acqs = [a for a in acquisitions if a.band_type == BandType.THERMAL]

        for p in range(npoints):
            name = 'POINT-%i-ALBEDO-TH' %p
            #atmospheric_profile = []
            atmos_profile = read_h5_table(ancillary_group, dname.format(p=p))
            n_layers = atmos_profile.shape[0] + 6
            elevation = atmos_profile.iloc[0]['GeoPotential_Height']
            
            prof_alt,prof_pres,prof_temp,prof_water = [],[],[],[]
            
            for i, row in atmos_profile.iterrows():
                prof_alt.append(row['GeoPotential_Height'])
                prof_pres.append(row['Pressure'])
                prof_temp.append(row['Temperature'])
                prof_water.append(row['Relative_Humidity'])
                
                
            ozone = ozone
            filter_function = acqs[0].spectral_filter_file
            visibility = -aerosol
            gpheight= elevation
            n =  n_layers
            sat_height = acquisitions[0].altitude / 1000.0
            sat_view = view_corrected[p]
            binary = False
            input_data = {'ozone':ozone,
                          'filter_function':filter_function,
                          'visibility':visibility,
                          'gpheight':elevation,
                          'n': n,
                          'sat_height':sat_height,
                          'sat_view':sat_view,
                          'prof_alt':prof_alt,
                          'prof_pres':prof_pres,
                          'prof_temp':prof_temp,
                          'prof_water':prof_water,
                          'description':'Input File for MODTRAN',
                          'file_format':'json'}


            data = mpjson.thermal_transmittance(name, ozone, n, prof_alt, prof_pres, prof_temp, prof_water,
                                                visibility, sat_height, gpheight, sat_view, filter_function, binary)
            json_data[(p, Albedos.ALBEDO_TH)] = data

            data = json.dumps(data, cls=JsonEncoder, indent=4)
            out_dname = ppjoin(POINT_FMT.format(p=p),
                               ALBEDO_FMT.format(a=Albedos.ALBEDO_TH.value),
                               DatasetName.MODTRAN_INPUT.value)
            write_scalar(numpy.string_(data), out_dname, group, input_data)

    # attach location info to each point Group
    for p in range(npoints):
        lonlat = (coordinator['longitude'][p], coordinator['latitude'][p])
        group[POINT_FMT.format(p=p)].attrs['lonlat'] = lonlat

    return json_data, out_group


def _run_modtran(acquisitions, modtran_exe, basedir, point, albedos, workflow,
                 npoints, atmospheric_inputs_fname, out_fname,
                 compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(atmospheric_inputs_fname, 'r') as atmos_fid,\
        h5py.File(out_fname, 'w') as fid:

        atmos_grp = atmos_fid[GroupName.ATMOSPHERIC_INPUTS_GRP.value]
        run_modtran(acquisitions, atmos_grp, workflow, npoints, point, albedos,
                    modtran_exe, basedir, fid, compression, filter_opts)


def run_modtran(acquisitions, atmospherics_group, workflow, npoints, point,
                albedos, modtran_exe, basedir, out_group,
                compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    Run MODTRAN and channel results.
    """
    lonlat = atmospherics_group[POINT_FMT.format(p=point)].attrs['lonlat']

    # determine the output group/file
    if out_group is None:
        fid = h5py.File('atmospheric-results.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    # initial attributes
    base_attrs = {'Point': point,
                  'lonlat': lonlat,
                  'datetime': acquisitions[0].acquisition_datetime}

    base_path = ppjoin(GroupName.ATMOSPHERIC_RESULTS_GRP.value,
                       POINT_FMT.format(p=point))

    # what atmospheric calculations have been run and how many points
    group_name = GroupName.ATMOSPHERIC_RESULTS_GRP.value
    if group_name not in fid:
        fid.create_group(group_name)

    fid[group_name].attrs['npoints'] = npoints
    applied = workflow == Workflow.STANDARD or workflow == Workflow.NBAR
    fid[group_name].attrs['nbar_atmospherics'] = applied
    applied = workflow == Workflow.STANDARD or workflow == Workflow.SBT
    fid[group_name].attrs['sbt_atmospherics'] = applied

    acqs = acquisitions
    for albedo in albedos:
        base_attrs['Albedo'] = albedo.value
        workpath = pjoin(basedir, POINT_FMT.format(p=point),
                         ALBEDO_FMT.format(a=albedo.value))
        
        json_mod_infile  = pjoin(workpath, ''.join([POINT_ALBEDO_FMT.format(p=point,a=albedo.value), '.json']))
        
        data_dir = pjoin(workpath, 'DATA')

        group_path = ppjoin(base_path, ALBEDO_FMT.format(a=albedo.value))

        subprocess.check_call([modtran_exe, json_mod_infile, data_dir],cwd=workpath)


        chn_fname = glob.glob(pjoin(workpath, '*.chn'))[0]    # need to change how .chn file is read
        tp6_fname = glob.glob(pjoin(workpath,'*.tp6'))[0]

        
        if albedo == Albedos.ALBEDO_TH:
            acq = [acq for acq in acqs if acq.band_type == BandType.THERMAL][0]
            
            channel_data = read_modtran_channel(chn_fname,tp6_fname, acq, albedo)

            attrs = base_attrs.copy()
            dataset_name = DatasetName.UPWARD_RADIATION_CHANNEL.value
            attrs['description'] = ('Upward radiation channel output from '
                                    'MODTRAN')
            dset_name = ppjoin(group_path, dataset_name)
            write_dataframe(channel_data[0], dset_name, fid, compression,
                            attrs=attrs, filter_opts=filter_opts)

            # downward radiation
            attrs = base_attrs.copy()
            dataset_name = DatasetName.DOWNWARD_RADIATION_CHANNEL.value
            attrs['description'] = ('Downward radiation channel output from '
                                    'MODTRAN')
            dset_name = ppjoin(group_path, dataset_name)
            write_dataframe(channel_data[1], dset_name, fid, compression,
                            attrs=attrs, filter_opts=filter_opts)
        else:
            acq = [acq for acq in acqs if
                   acq.band_type == BandType.REFLECTIVE][0]

            channel_data = read_modtran_channel(chn_fname,tp6_fname, acq, albedo) # needs changing to read json

            attrs = base_attrs.copy()
            dataset_name = DatasetName.CHANNEL.value
            attrs['description'] = 'Channel output from MODTRAN'
            dset_name = ppjoin(group_path, dataset_name)
            write_dataframe(channel_data[0], dset_name, fid, compression,
                            attrs=attrs, filter_opts=filter_opts)

            # solar zenith angle at surface
            attrs = base_attrs.copy()
            dataset_name = DatasetName.SOLAR_ZENITH_CHANNEL.value
            attrs['description'] = 'Solar zenith angle at different atmosphere levels'
            dset_name = ppjoin(group_path,dataset_name)
            write_dataframe(channel_data[1],dset_name,fid,compression,attrs=attrs,filter_opts=filter_opts)
            

    # metadata for a given point
    alb_vals = [alb.value for alb in workflow.albedos]
    fid[base_path].attrs['lonlat'] = lonlat
    fid[base_path].attrs['datetime'] = acqs[0].acquisition_datetime.isoformat()
    fid[base_path].attrs.create('albedos', data=alb_vals, dtype=VLEN_STRING)

    if out_group is None:
        return fid

def _calculate_coefficients(atmosheric_results_fname, out_fname,
                            compression=H5CompressionFilter.LZF,
                            filter_opts=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(atmosheric_results_fname, 'r') as atmos_fid,\
        h5py.File(out_fname, 'w') as fid:

        results_group = atmos_fid[GroupName.ATMOSPHERIC_RESULTS_GRP.value]
        calculate_coefficients(results_group, fid, compression, filter_opts)


def calculate_coefficients(atmospheric_results_group, out_group,
                           compression=H5CompressionFilter.LZF,
                           filter_opts=None):
    """
    Calculate the atmospheric coefficients from the MODTRAN output
    and used in the BRDF and atmospheric correction.
    Coefficients are computed for each band for each each coordinate
    for each atmospheric coefficient. The atmospheric coefficients can be
    found in `Workflow.STANDARD.atmos_coefficients`.

    :param atmospheric_results_group:
        The root HDF5 `Group` that contains the atmospheric results
        from each MODTRAN run.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The datasets will be formatted to the HDF5 TABLE specification
        and the dataset names will be as follows:

        * DatasetName.NBAR_COEFFICIENTS (if Workflow.STANDARD or Workflow.NBAR)
        * DatasetName.SBT_COEFFICIENTS (if Workflow.STANDARD or Workflow.SBT)

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    nbar_coefficients = pd.DataFrame()
    sbt_coefficients = pd.DataFrame()

    channel_data = upward = downward = None

    # Initialise the output group/file
    if out_group is None:
        fid = h5py.File('atmospheric-coefficients.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    res = atmospheric_results_group
    npoints = res.attrs['npoints']
    nbar_atmos = res.attrs['nbar_atmospherics']
    sbt_atmos = res.attrs['sbt_atmospherics']

    for point in range(npoints):
        point_grp = res[POINT_FMT.format(p=point)]
        lonlat = point_grp.attrs['lonlat']
        timestamp = pd.to_datetime(point_grp.attrs['datetime'])
        grp_path = ppjoin(POINT_FMT.format(p=point), ALBEDO_FMT)
        
        if nbar_atmos:

            channel_path = ppjoin(grp_path.format(a=Albedos.ALBEDO_0.value),
                                  DatasetName.CHANNEL.value)
            channel_data = read_h5_table(res, channel_path)
            
            
            channel_solar_angle_path = ppjoin(grp_path.format(a=Albedos.ALBEDO_0.value),
                                  DatasetName.SOLAR_ZENITH_CHANNEL.value)
            
            channel_solar_angle = read_h5_table(res, channel_solar_angle_path)

            
        if sbt_atmos:
            dname = ppjoin(grp_path.format(a=Albedos.ALBEDO_TH.value),
                           DatasetName.UPWARD_RADIATION_CHANNEL.value)
            
            upward = read_h5_table(res, dname)
            
            dname = ppjoin(grp_path.format(a=Albedos.ALBEDO_TH.value),
                           DatasetName.DOWNWARD_RADIATION_CHANNEL.value)
            downward = read_h5_table(res, dname)

        
        kwargs = {'channel_data': channel_data,
                  'solar_zenith_angle':channel_solar_angle,
                  'upward_radiation': upward,
                  'downward_radiation': downward}

        result = coefficients(**kwargs)

        # insert some datetime/geospatial fields
        if result[0] is not None:
            result[0].insert(0, 'POINT', point)
            result[0].insert(1, 'LONGITUDE', lonlat[0])
            result[0].insert(2, 'LATITUDE', lonlat[1])
            result[0].insert(3, 'DATETIME', timestamp)
            nbar_coefficients = nbar_coefficients.append(result[0])

        if result[1] is not None:
            result[1].insert(0, 'POINT', point)
            result[1].insert(1, 'LONGITUDE', lonlat[0])
            result[1].insert(2, 'LATITUDE', lonlat[1])
            result[1].insert(3, 'DATETIME', pd.to_datetime(timestamp))
            sbt_coefficients = sbt_coefficients.append(result[1])

    nbar_coefficients.reset_index(inplace=True)
    sbt_coefficients.reset_index(inplace=True)

    attrs = {'npoints': npoints}
    description = "Coefficients derived from the VNIR solar irradiation."
    attrs['description'] = description
    dname = DatasetName.NBAR_COEFFICIENTS.value

    if GroupName.COEFFICIENTS_GROUP.value not in fid:
        fid.create_group(GroupName.COEFFICIENTS_GROUP.value)

    group = fid[GroupName.COEFFICIENTS_GROUP.value]
    
    if nbar_atmos:
        write_dataframe(nbar_coefficients, dname, group, compression,
                        attrs=attrs, filter_opts=filter_opts)

    description = "Coefficients derived from the THERMAL solar irradiation."
    attrs['description'] = description
    dname = DatasetName.SBT_COEFFICIENTS.value

    if sbt_atmos:
        write_dataframe(sbt_coefficients, dname, group, compression,
                        attrs=attrs, filter_opts=filter_opts)

    if out_group is None:
        return fid


def coefficients(channel_data=None,solar_zenith_angle=None,upward_radiation=None,downward_radiation=None):
    """
    Calculate the coefficients for a given point.
    Calculate the atmospheric coefficients from the MODTRAN output
    and used in the BRDF and atmospheric correction.
    Coefficients are computed for each band.
    The atmospheric coefficients can be found in
    `Workflow.STANDARD.atmos_coefficients`.


    :param channel_data:
        A `pandas.DataFrame` containing the channel data for that
        point, and structured as returned by the
        `read_modtran_channel` function.
        Only used for NBAR calculations.
        
    :param solar_zenith_angle:
        A 'pandas.DataFrame' containing the value for solar zenith angle 
        for atmosphere layers(altitudes)
   
    :param upward_radiation:
        A `pandas.DataFrame` containing the upward radiation data for
        that point, and structured as returned by the
        `read_modtran_channel` function.
        Only used for SBT calculations.

    :param downward_radiation:
        A `pandas.DataFrame` containing the downward radiation data for
        that point, and structured as returned by the
        `read_modtran_channel` function.
        Only used for SBT calculations.

    :return:
        A `tuple` (nbar_coefficients, sbt_coefficients) whereby each
        item is a `pandas.DataFrame` containing the atmospheric
        coefficients for each band.
        If `accumulation_albedo_0` is None, then the first item in
        the returned `tuple` will be None.
        If `upward_radiation` is None, then the second item in the
        returned `tuple` will be None.
    """
    

    nbar = sbt = None
    
    if channel_data is not None:

        # calculate transmittance using channel data
        columns = [v.value for v in Workflow.NBAR.atmos_coefficients]
        nbar = pd.DataFrame(columns=columns, index=channel_data.index)
        cszen0 = numpy.cos(numpy.radians(float(solar_zenith_angle['solar_zenith'][-1:])))
        csnsrf = numpy.cos(numpy.radians(float(solar_zenith_angle['solar_zenith'][0])))
        tv = channel_data['24']
        ts = ((channel_data['19'] / channel_data['18']) * csnsrf/numpy.pi) / tv
        Ts = channel_data['21'] / tv
        Tv = tv + channel_data['22']/Ts
        E0_cozen = (channel_data['18'] * numpy.pi/channel_data['8']) * (cszen0/csnsrf)
        A_prime = E0_cozen * Ts * Tv/numpy.pi

        # nbar coefficients
        nbar[AC.FS.value] = ts / Ts
        nbar[AC.FV.value] = tv / Tv
        nbar[AC.A.value] = A_prime * 10000000
        nbar[AC.B.value] = channel_data['4'] * 10000000
        nbar[AC.S.value] = channel_data['23']
        nbar[AC.DIR.value] = E0_cozen * ts * 10000000
        nbar[AC.DIF.value] = (Ts - ts) * E0_cozen * 10000000
        nbar[AC.TS.value] = ts
        
    if upward_radiation is not None:
        columns = [v.value for v in Workflow.SBT.atmos_coefficients]
        columns.extend(['TRANSMITTANCE-DOWN']) # Currently not required
        sbt = pd.DataFrame(columns=columns, index=upward_radiation.index)

        # sbt coefficients
        sbt[AC.PATH_UP.value] = upward_radiation['4'] * 10000000
        sbt[AC.TRANSMITTANCE_UP.value] = upward_radiation['15']
        sbt[AC.PATH_DOWN.value] = downward_radiation['4'] * 10000000
        sbt['TRANSMITTANCE-DOWN'] = downward_radiation['15']

    return nbar, sbt

def read_spectral_response(fname, as_list=False, spectral_range=None):

    """
    Read the spectral response function text file used during
    MODTRAN processing.

    :param fname:
        A `str` containing the full file path name, or an opened
        `file` buffer.

    :param as_list:
        A `bool` indicating whether or not to return the spectral
        response data as a list instead of a `pd.DataFrame`.
        Default is `False` which returns a `pd.DataFrame`.

    :param spectral_range:
        A `list` or `generator` of the [start, stop, step] for the
        spectral range to be used in defining the spectral response.
        Default is [2600, 349, -1].

    :return:
        A `pd.DataFrame` containing the spectral response
        function.
    """
    if isinstance(fname, str):
        with open(fname, 'r') as src:
            lines = src.readlines()
    else:
        lines = fname.readlines()

    if as_list:
        return lines

    lines = [line.strip().decode('utf-8') for line in lines]

    # find the starting locations of each band description label
    ids = []
    for i, val in enumerate(lines):
        if 'B' in val:
            ids.append(i)

    # get the spectral response data up to band n-1
    response = {}
    for i, idx in enumerate(ids[0:-1]):
        data = numpy.array([l.split('  ') for l in lines[idx+1:ids[i+1]]],
                           dtype='float')
        df = pd.DataFrame({'band_name': lines[idx],
                           'wavelength': data[:, 0],
                           'response': data[:, 1]})
        response[lines[idx]] = df

    # get spectral response data for band n
    idx = ids[-1]
    data = numpy.array([l.split('  ') for l in lines[idx+1:]], dtype='float')
    df = pd.DataFrame({'band_name': lines[idx],
                       'wavelength': data[:, 0],
                       'response': data[:, 1]})
    response[lines[idx]] = df

    if spectral_range is None:
        wavelengths = range(2600, 349, -1)
    else:
        wavelengths = list(spectral_range)

    for band in response:
        base_df = pd.DataFrame({'wavelength': wavelengths,
                                'response': 0.0,
                                'band_name': band},
                               index=wavelengths)
        df = response[band]
        base_df.ix[df['wavelength'], 'response'] = df['response'].values

        response[band] = base_df

    spectral_response = pd.concat(response, names=['band_name', 'wavelength'])
    spectral_response.drop(['band_name', 'wavelength'], inplace=True, axis=1)

    return spectral_response

def _get_solar_angles(tp6_fname):
    """
    Read a MODTRAN output '*.tp6' ascii file to extract
    solar zenith angles
    
    :param tp6_fname:
        A 'str' containing the full file pathname of the tp6 
        data file 
        
    :return: 
        solar zenith angle at surface and top of atmosphere 
    
    """
    
    f = open(tp6_fname,'r')
    cnt = 0
    cnt_start = 0

    for line in f.readlines(): 
        cnt +=1 
        if fnmatch.fnmatch(line,'*SINGLE SCATTER SOLAR PATH GEOMETRY TABLE FOR MULTIPLE SCATTERING*'):
            cnt_start = cnt
            continue

    f.close()
    f = open(tp6_fname,'r')
    
    lines = list(f.readlines()[cnt_start+4: cnt_start + 40])  # this constant is (40-4 = num of atm 
                                                              # layer)
    data = numpy.vstack([ln.rstrip().split() for ln in lines]).T
      
    solar_zenith = data[3] # solar zenith angle at all atmosphere layers
    
    return solar_zenith

    
def read_modtran_channel(fname, tp6_fname, acquisition, albedo):
    """
    Read a MODTRAN output `*.chn` ascii file.
    
    :param fname:
        A `str` containing the full file pathname of the channel
        data file.
    
    :param tp6_fname:
        A 'str' containing the full file pathname of the tp6 
        data file 


    :param acquisition:
        An instance of an acquisition object.

    :param albedo:
        An albedo identifier from either Workflow.nbar.albedos or
        Workflow.SBT.albedos

    :return:
        A `pandas.DataFrame` containing the channel data, and index
        by the `band_name`.
    """
    response = acquisition.spectral_response()
    nbands = response.index.get_level_values('band_name').unique().shape[0]

   
    if albedo == Albedos.ALBEDO_TH:
        
        upward_radiation = pd.read_csv(fname, skiprows=5, header=None,
                                       delim_whitespace=True, nrows=nbands)
        
        downward_radiation = pd.read_csv(fname, skiprows=10+nbands,
                                         header=None, delim_whitespace=True,
                                         nrows=nbands)
        
        upward_radiation['band_name'] = upward_radiation[17]
        downward_radiation['band_name'] = downward_radiation[17]
        upward_radiation.drop(17, inplace=True, axis=1)
        downward_radiation.drop(17, inplace=True, axis=1)
        upward_radiation.set_index('band_name', inplace=True)
        downward_radiation.set_index('band_name', inplace=True)
        upward_radiation.columns = upward_radiation.columns.astype(str)
        downward_radiation.columns = downward_radiation.columns.astype(str)

        return upward_radiation, downward_radiation
    
    
    #get solar zenith angle at all layers from *.tp6 file 
    solar_zenith = _get_solar_angles(tp6_fname)
   
    df_sz_angle = pd.DataFrame()
   
    df_sz_angle['solar_zenith'] = solar_zenith
    
    
    chn_data = pd.read_csv(fname, skiprows=5, header=None, nrows=nbands,
                           delim_whitespace=True)

    chn_data['band_name'] = chn_data[26] 
    chn_data.drop(26,inplace=True,axis=1)
    
    chn_data.set_index('band_name',inplace=True)
    chn_data.columns = chn_data.columns.astype(str)

    return chn_data,df_sz_angle


def link_atmospheric_results(input_targets, out_fname, npoints, workflow):
    """
    Uses h5py's ExternalLink to combine the atmospheric results into
    a single file.

    :param input_targets:
        A `list` of luigi.LocalTargets.

    :param out_fname:
        A `str` containing the output filename.

    :param npoints:
        An `int` containing the number of points (vertices) used for
        evaluating the atmospheric conditions.

    :param workflow:
        An Enum given by wagl.constants.Workflow.

    :return:
        None. Results from each file in `input_targets` are linked
        into the output file.
    """
    base_group_name = GroupName.ATMOSPHERIC_RESULTS_GRP.value
    nbar_atmospherics = False
    sbt_atmospherics = False
    attributes = []
    for fname in input_targets:
        with h5py.File(fname.path, 'r') as fid:
            points = list(fid[base_group_name].keys())

            # copy across several attributes on the POINT Group
            # as the linking we do here links direct to a dataset
            # which will create the required parent Groups
            for point in points:
                group = ppjoin(base_group_name, point)
                attributes.append((point,
                                   fid[group].attrs['lonlat'],
                                   fid[group].attrs['datetime'],
                                   fid[group].attrs['albedos']))
                
        for point in points:
            for albedo in workflow.albedos:
                if albedo == Albedos.ALBEDO_TH:
                    datasets = [DatasetName.UPWARD_RADIATION_CHANNEL.value,
                                DatasetName.DOWNWARD_RADIATION_CHANNEL.value]
                    sbt_atmospherics = True
                else:
                    datasets = [DatasetName.CHANNEL.value,
                                DatasetName.SOLAR_ZENITH_CHANNEL.value]
                    nbar_atmospherics = True

                grp_path = ppjoin(base_group_name, point,
                                  ALBEDO_FMT.format(a=albedo.value))

                for dset in datasets:
                    dname = ppjoin(grp_path, dset)
                    create_external_link(fname.path, dname, out_fname, dname)

    with h5py.File(out_fname) as fid:
        group = fid[GroupName.ATMOSPHERIC_RESULTS_GRP.value]
        group.attrs['npoints'] = npoints
        group.attrs['nbar_atmospherics'] = nbar_atmospherics
        group.attrs['sbt_atmospherics'] = sbt_atmospherics

        # assign the lonlat attribute for each POINT Group
        for point, lonlat, date_time, albedos in attributes:
            group[point].attrs['lonlat'] = lonlat
            group[point].attrs['datetime'] = date_time
            group[point].attrs.create('albedos', data=albedos,
                                      dtype=VLEN_STRING)
