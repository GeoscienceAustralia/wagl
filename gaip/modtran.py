"""
MODTRAN drivers
---------------

"""
from __future__ import absolute_import, print_function
import os
from os.path import join as pjoin, exists, dirname, basename, splitext
from posixpath import join as ppjoin
import subprocess
import glob

import numpy
from scipy.io import FortranFile
import h5py
import pandas as pd

from gaip.constants import Model, POINT_FMT, ALBEDO_FMT, POINT_ALBEDO_FMT
from gaip.hdf5 import write_dataframe, read_table, create_external_link
from gaip.modtran_profiles import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO
from gaip.modtran_profiles import MIDLAT_SUMMER_TRANSMITTANCE, SBT_FORMAT
from gaip.modtran_profiles import TROPICAL_TRANSMITTANCE, THERMAL_TRANSMITTANCE


def create_modtran_dirs(coords, albedos, modtran_root, modtran_exe_root,
                        workpath_format, input_format):
    """Create all modtran subdirectories. and input files."""

    if not exists(modtran_root):
        os.makedirs(modtran_root)

    data_dir = pjoin(modtran_exe_root, 'DATA')
    if not exists(data_dir):
        raise OSError('Cannot find MODTRAN')

    for coord in coords:
        for albedo in albedos:
            modtran_work = workpath_format.format(coord=coord, albedo=albedo)
            modtran_work = pjoin(modtran_root, modtran_work)
            mod5root_in = input_format.format(coord=coord, albedo=albedo)
            mod5root_in = pjoin(modtran_root, mod5root_in)

            if not exists(modtran_work):
                os.makedirs(modtran_work)

            with open(mod5root_in, 'w') as outfile:
                outfile.write(coord + '_alb_' + albedo + '\n')

            symlink_dir = pjoin(modtran_work, 'DATA')
            if exists(symlink_dir):
                os.unlink(symlink_dir)

            os.symlink(data_dir, symlink_dir)


def prepare_modtran(acquisition, coordinate, albedo, modtran_work,
                    modtran_exe):
    """
    Prepares the working directory for a MODTRAN execution.
    """
    data_dir = pjoin(dirname(modtran_exe), 'DATA')
    if not exists(data_dir):
        raise OSError('Cannot find MODTRAN')

    out_fname = pjoin(modtran_work, 'mod5root.in')

    with open(out_fname, 'w') as src:
        src.write(POINT_ALBEDO_FMT.format(p=coordinate, a=albedo) + '\n')

    symlink_dir = pjoin(modtran_work, 'DATA')
    if exists(symlink_dir):
        os.unlink(symlink_dir)

    os.symlink(data_dir, symlink_dir)

    # TODO: write the spectral response function
    out_fname = pjoin(modtran_work, acquisition.spectral_filter_file)
    response = acquisition.spectral_response(as_list=True)
    with open(out_fname, 'wb') as src:
        src.writelines(response)


def _format_tp5(acquisition, satellite_solar_angles_fname,
                longitude_latitude_fname, ancillary_fname, out_fname,
                nbar_tp5=True):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_angles_fname, 'r') as sat_sol,\
        h5py.File(longitude_latitude_fname, 'r') as lon_lat_ds,\
        h5py.File(ancillary_fname, 'r') as anc_ds,\
        h5py.File(out_fname, 'w') as fid:

        # angles data
        view_dset = sat_sol['satellite-view']
        azi_dset = sat_sol['satellite-azimuth']
        lon_dset = lon_lat_ds['longitude']
        lat_dset = lon_lat_ds['latitude']

        # ancillary data
        coord_dset = anc_ds['coordinator']
        aerosol = anc_ds['aerosol'][()]
        water_vapour = anc_ds['water-vapour'][()]
        ozone = anc_ds['ozone'][()]
        elevation = anc_ds['elevation'][()]

        if anc_ds.attrs.get('sbt-ancillary'):
            sbt_ancillary = {}
            dname = ppjoin(POINT_FMT, 'atmospheric-profile')
            for i in range(coord_dset.shape[0]):
                sbt_ancillary[i] = read_table(anc_ds, dname.format(p=i))
        else:
            sbt_ancillary = None

        tp5_data, metadata = format_tp5(acquisition, coord_dset, view_dset,
                                        azi_dset, lat_dset, lon_dset, ozone,
                                        water_vapour, aerosol, elevation,
                                        coord_dset.shape[0], sbt_ancillary,
                                        nbar_tp5)

        group = fid.create_group('modtran-inputs')
        iso_time = acquisition.scene_centre_date.isoformat()
        group.attrs['acquisition-datetime'] = iso_time

        for key in metadata:
            dname = ppjoin(POINT_FMT.format(p=key[0]),
                           ALBEDO_FMT.format(a=key[1]), 'tp5_data')
            str_data = numpy.string_(tp5_data[key])
            dset = group.create_dataset(dname, data=str_data)
            for k in metadata[key]:
                dset.attrs[k] = metadata[key][k]

        # attach some meaningful location information to the point groups
        lon = coord_dset['longitude']
        lat = coord_dset['latitude']
        for i in range(coord_dset.shape[0]):
            group[POINT_FMT.format(p=i)].attrs['lonlat'] = (lon[i], lat[i])

    return tp5_data


def format_tp5(acquisition, coordinator, view_dataset, azi_dataset,
               lat_dataset, lon_dataset, ozone, vapour, aerosol, elevation,
               npoints, sbt_ancillary=None, nbar_tp5=True):
    """
    Creates str formatted tp5 files for the albedo (0, 1) and
    transmittance (t).
    """
    geobox = acquisition.gridded_geo_box()
    filter_file = acquisition.spectral_filter_file
    cdate = acquisition.scene_centre_date
    doy = int(cdate.strftime('%j'))
    altitude = acquisition.altitude / 1000.0  # in km
    dechour = acquisition.decimal_hour

    view = numpy.zeros(npoints, dtype='float32')
    azi = numpy.zeros(npoints, dtype='float32')
    lat = numpy.zeros(npoints, dtype='float64')
    lon = numpy.zeros(npoints, dtype='float64')

    for i in range(npoints):
        yidx = coordinator['row_index'][i]
        xidx = coordinator['col_index'][i]
        idx = (slice(yidx, yidx + 1), slice(xidx, xidx + 1))
        view[i] = view_dataset[idx][0, 0]
        azi[i] = azi_dataset[idx][0, 0]
        lat[i] = lat_dataset[idx][0, 0]
        lon[i] = lon_dataset[idx][0, 0]

    view_cor = 180 - view
    azi_cor = azi + 180
    rlon = 360 - lon
    
    # check if in western hemisphere
    wh = rlon >= 360
    rlon[wh] -= 360
    
    wh = (180 - view_cor) < 0.1
    view_cor[wh] = 180
    azi_cor[wh] = 0
    
    wh = azi_cor > 360
    azi_cor[wh] -= 360

    # get the modtran profiles to use based on the centre latitude 
    _, centre_lat = geobox.centre_lonlat
    if centre_lat < -23.0:
        albedo_profile = MIDLAT_SUMMER_ALBEDO
        trans_profile = MIDLAT_SUMMER_TRANSMITTANCE
    else:
        albedo_profile = TROPICAL_ALBEDO
        trans_profile = TROPICAL_TRANSMITTANCE

    # we'll only cater for MODTRAN to output binary form
    binary = 'T'

    tp5_data = {}
    metadata = {}

    # write the tp5 files required for input into MODTRAN
    if nbar_tp5:
        for i in range(npoints):
            for alb in Model.nbar.albedos:
                input_data = {'water': vapour,
                              'ozone': ozone,
                              'filter_function': filter_file,
                              'visibility': -aerosol,
                              'elevation': elevation,
                              'sat_height': altitude,
                              'sat_view': view_cor[i],
                              'doy': doy,
                              'binary': binary}
                if alb == Model.nbar.albedos[2]:
                    input_data['albedo'] = 0.0
                    input_data['sat_view_offset'] = 180.0-view_cor[i]
                    data = trans_profile.format(**input_data)
                else:
                    input_data['albedo'] = float(alb)
                    input_data['lat'] = lat[i]
                    input_data['lon'] = rlon[i]
                    input_data['time'] = dechour
                    input_data['sat_azimuth'] = azi_cor[i]
                    data = albedo_profile.format(**input_data)

                tp5_data[(i, alb)] = data
                metadata[(i, alb)] = input_data

    # tp5 for sbt; the current logic for NBAR uses 9 coordinator points
    # and sbt uses 25 coordinator points
    # as such points [0, 9) in nbar will not be the same [0, 9) points in
    # the sbt coordinator
    # hopefully the science side of the algorithm will be re-engineered
    # so as to ensure a consistant logic between the two products

    if sbt_ancillary is not None:
        for p in range(npoints):
            atmospheric_profile = []
            atmos_profile = sbt_ancillary[p]
            n_layers = atmos_profile.shape[0] + 6
            elevation = atmos_profile.iloc[0]['GeoPotential_Height']
            for i, row in atmos_profile.iterrows():
                input_data = {'gpheight': row['GeoPotential_Height'],
                              'pressure': row['Pressure'],
                              'airtemp': row['Temperature'],
                              'humidity': row['Relative_Humidity'],
                              'zero': 0.0}
                atmospheric_profile.append(SBT_FORMAT.format(**input_data))
            
            input_data = {'ozone': ozone,
                          'filter_function': filter_file,
                          'visibility': -aerosol,
                          'gpheight': elevation,
                          'n': n_layers,
                          'sat_height': altitude,
                          'sat_view': view_cor[p],
                          'binary': binary,
                          'data_array': ''.join(atmospheric_profile)}

            data = THERMAL_TRANSMITTANCE.format(**input_data)
            tp5_data[(p, Model.sbt.albedos[0])] = data
            metadata[(p, Model.sbt.albedos[0])] = input_data

    return tp5_data, metadata


def _run_modtran(acquisition, modtran_exe, workpath, point, albedo,
                 atmospheric_inputs_fname, out_fname, compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(atmospheric_inputs_fname, 'r') as fid:
        grp_path = ppjoin('modtran-inputs', POINT_FMT.format(p=point))
        lonlat = fid[grp_path].attrs['lonlat']

    rfid = run_modtran(acquisition, modtran_exe, workpath, point, albedo,
                       lonlat, out_fname, compression)

    rfid.close()
    return
        

def run_modtran(acquisition, modtran_exe, workpath, point, albedo, lonlat=None,
                out_fname=None, compression='lzf'):
    """
    Run MODTRAN and return the flux and channel results.
    """
    subprocess.check_call([modtran_exe], cwd=workpath)

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('modtran-results.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    fid.attrs['point'] = point
    fid.attrs['albedo'] = albedo

    # base group pathname
    group_path = ppjoin(POINT_FMT.format(p=point), ALBEDO_FMT.format(a=albedo))

    if lonlat is None:
        lonlat = (numpy.nan, numpy.nan)

    # initial attributes
    attrs = {'Point': point, 'Albedo': albedo, 'lonlat': lonlat}

    if albedo != Model.sbt.albedos[0]:
        flux_fname = glob.glob(pjoin(workpath, '*_b.flx'))[0]
        flux_data, altitudes = read_modtran_flux(flux_fname)

        # ouput the flux data
        dset_name = ppjoin(group_path, 'flux')
        attrs['Description'] = 'Flux output from MODTRAN'
        write_dataframe(flux_data, dset_name, fid, attrs=attrs)

        # output the altitude data
        attrs['Description'] = 'Altitudes output from MODTRAN'
        attrs['altitude levels'] = altitudes.shape[0]
        attrs['units'] = 'km'
        dset_name = ppjoin(group_path, 'altitudes')
        write_dataframe(altitudes, dset_name, fid, attrs=attrs)

        # accumulate the solar irradiance
        transmittance = True if albedo == Model.nbar.albedos[2] else False
        response = acquisition.spectral_response()
        accumulated = calculate_solar_radiation(flux_data, response,
                                                altitudes.shape[0],
                                                transmittance)

        dset_name = ppjoin(group_path, 'solar-irradiance')
        description = ("Accumulated solar irradiation for point {} "
                       "and albedo {}.")
        attrs['Description'] = description.format(point, albedo),
        write_dataframe(accumulated, dset_name, fid, compression, attrs=attrs)

    chn_fname = glob.glob(pjoin(workpath, '*.chn'))[0]
    channel_data = read_modtran_channel(chn_fname, acquisition, albedo)

    if albedo == Model.sbt.albedos[0]:
        # upward radiation
        attrs['Description'] = 'Upward radiation channel output from MODTRAN'
        dset_name = ppjoin(group_path, 'upward-radiation-channel')
        write_dataframe(channel_data[0], dset_name, fid, attrs=attrs)

        # downward radiation
        attrs['Description'] = 'Downward radiation channel output from MODTRAN'
        dset_name = ppjoin(group_path, 'downward-radiation-channel')
        write_dataframe(channel_data[1], dset_name, fid, attrs=attrs)
    else:
        # output the channel data
        attrs['Description'] = 'Channel output from MODTRAN'
        dset_name = ppjoin(group_path, 'channel')
        write_dataframe(channel_data, dset_name, fid, attrs=attrs)


    # meaningful location description
    fid[POINT_FMT.format(p=point)].attrs['lonlat'] = lonlat

    return fid


def _calculate_coefficients(atmospheric_fname, out_fname, compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    nbar_albedos = Model.nbar.albedos
    with h5py.File(atmospheric_fname, 'r') as fid:
        nbar_atmos = fid.attrs['nbar atmospherics']
        sbt_atmos = fid.attrs['sbt atmospherics']

        # initialise dicts to hold the data for each point
        accumulation_albedo_0 = {} if nbar_atmos else None
        accumulation_albedo_1 = {} if nbar_atmos else None
        accumulation_albedo_t = {} if nbar_atmos else None
        channel_data = {} if nbar_atmos else None
        upward = {} if sbt_atmos else None
        downward = {} if sbt_atmos else None

        for point in range(fid.attrs['npoints']):
            grp_path = ppjoin(POINT_FMT.format(p=point), ALBEDO_FMT)
            if nbar_atmos:
                albedo_0_path = ppjoin(grp_path.format(a=nbar_albedos[0]),
                                       'solar-irradiance')
                albedo_1_path = ppjoin(grp_path.format(a=nbar_albedos[1]),
                                       'solar-irradiance')
                albedo_t_path = ppjoin(grp_path.format(a=nbar_albedos[2]),
                                       'solar-irradiance')
                channel_path = ppjoin(grp_path.format(a=nbar_albedos[0]),
                                      'channel')

                accumulation_albedo_0[point] = read_table(fid, albedo_0_path)
                accumulation_albedo_1[point] = read_table(fid, albedo_1_path)
                accumulation_albedo_t[point] = read_table(fid, albedo_t_path)
                channel_data[point] = read_table(fid, channel_path)
            if sbt_atmos:
                dname = ppjoin(grp_path.format(a=Model.sbt.albedos[0]),
                               'upward-radiation-channel')
                upward[point] = read_table(fid, dname)
                dname = ppjoin(grp_path.format(a=Model.sbt.albedos[0]),
                               'downward-radiation-channel')
                downward[point] = read_table(fid, dname)
        
        kwargs = {'accumulation_albedo_0': accumulation_albedo_0,
                  'accumulation_albedo_1': accumulation_albedo_1,
                  'accumulation_albedo_t': accumulation_albedo_t,
                  'channel_data': channel_data,
                  'upward': upward,
                  'downward': downward,
                  'npoints': fid.attrs['npoints'],
                  'out_fname': out_fname,
                  'compression': compression}
                  
        rfid = calculate_coefficients(**kwargs)

    rfid.close()
    return
        

def calculate_coefficients(accumulation_albedo_0=None,
                           accumulation_albedo_1=None,
                           accumulation_albedo_t=None, channel_data=None,
                           upward=None, downward=None, npoints=None,
                           out_fname=None, compression='lzf'):
    """
    Calculate the atmospheric coefficients from the MODTRAN output
    and used in the BRDF and atmospheric correction.
    Coefficients are computed for each band for each each coordinate
    for each factor.  The factors are:
    ['fs', 'fv', 'a', 'b', 's', 'dir', 'dif', 'ts'].

    :param accumulation_albedo_0:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the solar accumulated
        irradiance (for albedo 0) and structured as returned by the
        `calculate_solar_radiation` function.
        Only used for NBAR calculations.

    :param accumulation_albedo_1:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the solar accumulated
        irradiance (for albedo 1) and structured as returned by the
        `calculate_solar_radiation` function.
        Only used for NBAR calculations.

    :param accumulation_albedo_t:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the solar accumulated
        irradiance (for albeod t; transmittance) and structured as
        returned by the `calculate_solar_radiation` function.
        Only used for NBAR calculations.

    :param channel_data:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the channel data for that
        point, and structured as returned by the
        `read_modtran_channel` function.
        Only used for NBAR calculations.

    :param upward:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the upward radiation data for
        that point, and structured as returned by the
        `read_modtran_channel` function.
        Only used for SBT calculations.

    :param downward:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the downward radiation data for
        that point, and structured as returned by the
        `read_modtran_channel` function.
        Only used for SBT calculations.

    :param npoints:
        An integer containing the number of location points over
        which MODTRAN was run.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.

        2 datasets formatted to the HDF5 TABLE specification
        named:

        * nbar-coefficients (if accumulation_albedo_0 is not None)
        * sbt-coefficients (if upward is not None)

    :notes:
        If accumulation_albedo_0 is None, then nbar-coefficients will
        not be calculated.
        If upward is None, then sbt-coefficients will not be
        calculated.
    """
    # Initialise the output file
    if out_fname is None:
        fid = h5py.File('coefficients.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    attrs = {'Number of atmospheric points': npoints}

    if accumulation_albedo_0 is not None:
        result = pd.DataFrame()
        for point in range(npoints):
            # MODTRAN channel output .chn file (albedo 0)
            data1 = channel_data[point]

            # solar radiation file (albedo 0)
            data3 = accumulation_albedo_0[point]

            # solar radiation file (albedo 1)
            data4 = accumulation_albedo_1[point]

            # solar radiation file (transmittance mode)
            data5 = accumulation_albedo_t[point]

            # calculate
            diff_0 = data3['diffuse'] * 10000000.0
            diff_1 = data4['diffuse'] * 10000000.0
            dir_0 = data3['direct'] * 10000000.0
            dir_1 = data4['direct'] * 10000000.0
            dir_t = data5['direct']
            dir0_top = data3['direct_top'] * 10000000.0
            dirt_top = data5['direct_top']
            tv_total = data5['transmittance']
            ts_total = (diff_0 + dir_0) / dir0_top
            ts_dir = dir_0 / dir0_top
            tv_dir = dir_t / dirt_top

            # TODO: better descriptive names
            columns = ['point',
                       'fs',
                       'fv',
                       'a',
                       'b',
                       's',
                       'dir',
                       'dif',
                       'ts']
            df = pd.DataFrame(columns=columns, index=data1.index)

            df['point'] = point
            df['fs'] = ts_dir / ts_total
            df['fv'] = tv_dir / tv_total
            df['a'] = (diff_0 + dir_0) / numpy.pi * tv_total
            df['b'] = data1['3'] * 10000000
            df['s'] = 1 - (diff_0 + dir_0) / (diff_1 + dir_1)
            df['dir'] = dir_0
            df['dif'] = diff_0
            df['ts'] = ts_dir

            result = result.append(df)

        result.reset_index(inplace=True)

        attrs['Description'] = ("Coefficients derived from the VNIR "
                                "accumulated solar irradiation.")


        write_dataframe(result, 'nbar-coefficients', fid, compression,
                        attrs=attrs)

    if upward is not None:
        upward_radiation = pd.concat(upward, names=['point'])
        downward_radiation = pd.concat(downward, names=['point'])

        result = pd.DataFrame(index=upward_radiation.index)
        result['path_up'] = upward_radiation['3'] * 10000000
        result['transmittance_up'] = upward_radiation['14']
        result['path_down'] = downward_radiation['3'] * 10000000
        result['transmittance_down'] = downward_radiation['14']
        result.reset_index(inplace=True)

        attrs = {}
        attrs['Description'] = ("Coefficients derived from the THERMAL "
                                "solar irradiation.")
        attrs['Number of atmospheric points'] = npoints

        write_dataframe(result, 'sbt-coefficients', fid, compression,
                        attrs=attrs)

    return fid


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
        df = pd.DataFrame({'band_id': lines[idx],
                           'wavelength': data[:, 0],
                           'response': data[:, 1]})
        response[lines[idx]] = df

    # get spectral response data for band n
    idx = ids[-1]
    data = numpy.array([l.split('  ') for l in lines[idx+1:]], dtype='float')
    df = pd.DataFrame({'band_id': lines[idx],
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
                                'band_id': band},
                               index=wavelengths)
        df = response[band]
        base_df.ix[df['wavelength'], 'response'] = df['response'].values

        response[band] = base_df

    spectral_response = pd.concat(response, names=['band_id', 'wavelength'])
    spectral_response.drop(['band_id', 'wavelength'], inplace=True, axis=1)

    return spectral_response


def read_modtran_flux(fname):
    """
    Read a MODTRAN output `*_b.flx` binary file.

    :param fname:
        A `str` containing the full file pathname of the flux
        data file.

    :return:
        Two `pandas.DataFrame's`. The first contains the spectral flux
        table data, and the second is contains the atmospheric height
        levels in km.
    """
    # define a datatype for the hdr info
    hdr_dtype = numpy.dtype([('record_length', 'int32'),
                             ('spectral_unit', 'S1'),
                             ('relabs', 'S1'),
                             ('linefeed', 'S1'),
                             ('mlflx', 'int32'),
                             ('iv1', 'float32'),
                             ('band_width', 'float32'),
                             ('fwhm', 'float32'),
                             ('ifwhm', 'float32')])

    # datatype for the dataframe containing the flux data
    flux_dtype = numpy.dtype([('upward_diffuse', 'float64'),
                              ('downward_diffuse', 'float64'),
                              ('direct_solar', 'float64')])

    with open(fname, 'rb') as src:
        # read the hdr record
        hdr_data = numpy.fromfile(src, hdr_dtype, count=1)

        # maximum flux levels at a spectral grid point
        levels = hdr_data['mlflx'][0] + 1

        # define a datatype to read a record containing flux data
        dtype = numpy.dtype([('wavelength', 'float64'),
                             ('flux_data', 'float64', (levels, 3))])

        # read the rest of the hdr which contains the altitude data
        altitude = numpy.fromfile(src, dtype='float32', count=levels)

        # read the record length end value
        _ = numpy.fromfile(src, 'int32', count=1)

        # initialise the FORTRAN read
        ffile = FortranFile(src)

        # read data from 2600 down to 350
        flux = {}
        wavelength_steps = range(2600, 349, -1)
        for wv in wavelength_steps:
            data = ffile.read_record(dtype)
            df = pd.DataFrame(numpy.zeros(levels, dtype=flux_dtype))
            #df['wavelength'] = data['wavelength'][0]
            df['upward_diffuse'] = data['flux_data'].squeeze()[:, 0]
            df['downward_diffuse'] = data['flux_data'].squeeze()[:, 1]
            df['direct_solar'] = data['flux_data'].squeeze()[:, 2]
            flux[wv] = df

    # concatenate into a single table
    flux_data = pd.concat(flux, names=['wavelength', 'level'])

    # setup a dataframe for the altitude
    altitude = pd.DataFrame({'altitude': altitude})
    altitude.index.name = 'layer'

    return flux_data, altitude


def read_modtran_channel(fname, acquisition, albedo):
    """
    Read a MODTRAN output `*.chn` ascii file.

    :param fname:
        A `str` containing the full file pathname of the channel
        data file.

    :param acquisition:
        An instance of an acquisition object.

    :param albedo:
        An albedo identifier from either Model.nbar.albedos or
        Model.sbt.albedos

    :return:
        A `pandas.DataFrame` containing the channel data, and index
        by the `band_id`.
    """
    if albedo == Model.sbt.albedos[0]:
        response = acquisition.spectral_response()
        nbands = response.index.get_level_values('band_id').unique().shape[0]
        upward_radiation = pd.read_csv(fname, skiprows=5, header=None,
                                       delim_whitespace=True, nrows=nbands)
        downward_radiation = pd.read_csv(fname, skiprows=10+nbands,
                                         header=None, delim_whitespace=True,
                                         nrows=nbands)
        upward_radiation['band_id'] = (upward_radiation[16] + ' ' + 
                                       upward_radiation[17].astype(str))
        downward_radiation['band_id'] = (downward_radiation[16] + ' ' +
                                         downward_radiation[17].astype(str))
        upward_radiation.drop([16, 17], inplace=True, axis=1)
        downward_radiation.drop([16, 17], inplace=True, axis=1)
        upward_radiation.set_index('band_id', inplace=True)
        downward_radiation.set_index('band_id', inplace=True)
        upward_radiation.columns = upward_radiation.columns.astype(str)
        downward_radiation.columns = downward_radiation.columns.astype(str)

        return upward_radiation, downward_radiation
    else:
        chn_data = pd.read_csv(fname, skiprows=5, header=None,
                               delim_whitespace=True)
        chn_data['band_id'] = chn_data[20] + ' ' + chn_data[21].astype(str)
        chn_data.drop([20, 21], inplace=True, axis=1)
        chn_data.set_index('band_id', inplace=True)
        chn_data.columns = chn_data.columns.astype(str)

        return chn_data


def calculate_solar_radiation(flux_data, spectral_response, levels=36,
                              transmittance=False):
    """
    Retreive the flux data from the MODTRAN output `*.flx`, and
    calculate the solar radiation.

    The solar radiation will be calculated for each of the bands
    contained within the spectral response dataset.

    :param flux_data:
        A `pandas.DataFrame` structured as if read from the
        `read_modtran_flux` function.

    :param spectral_response:
        A `pandas.DataFrame` containing the spectral response
        and structured as if read from the `read_spectral_response`
        function.

    :param levels:
        The number of atmospheric levels. Default is 36.

    :param transmittance:
        If set to `True`, then calculate the solar radiation in
        transmittance mode. Default is to calculate from albedo.

    :return:
        A `pandas.DataFrame` containing the solar radiation
        accumulation.
    """
    # index location of the top atmospheric level
    idx = levels - 1

    # group via the available bands
    band_index = spectral_response.index.get_level_values('band_id')
    groups = spectral_response.groupby(band_index)

    # output dataframe
    # later on the process can be refined to only evaluate the bands
    # we wish to process
    if transmittance:
        columns = ['diffuse',
                   'direct',
                   'diffuse_top',
                   'direct_top',
                   'transmittance']
    else:
        columns = ['diffuse', 'direct', 'direct_top']
    df = pd.DataFrame(columns=columns, index=groups.groups.keys(),
                      dtype='float64')

    # start wavelength and end wavelength, eg 2600 & 350 respectively
    st_wl = spectral_response.index[0][1]
    ed_wl = spectral_response.index[-1][1]

    # indices for flux at bottom and top of atmosphere layers
    wv_idx = range(st_wl - 1, ed_wl - 1, -1)
    wv_idx2 = [(i, 0) for i in wv_idx]
    wv_idx3 = [(i, idx) for i in wv_idx]

    # loop over each band and get the solar radiation
    for band, grp in groups:
        # df.ix[band, 'band_id'] = band

        # downward diffuse at bottom of atmospheric levels
        diffuse_bottom = (grp.ix[band, st_wl]['response'] *
                          flux_data.ix[st_wl, 'downward_diffuse'][0] +
                          grp.ix[band, ed_wl]['response'] *
                          flux_data.ix[ed_wl, 'downward_diffuse'][0]) / 2

        # direct solar at bottom of atmospheric levels
        direct_bottom = (grp.ix[band, st_wl]['response'] *
                         flux_data.ix[st_wl, 'direct_solar'][0] +
                         grp.ix[band, ed_wl]['response'] *
                         flux_data.ix[ed_wl, 'direct_solar'][0]) / 2

        # direct solar at top of atmospheric levels
        direct_top = (grp.ix[band, st_wl]['response'] *
                      flux_data.ix[st_wl, 'direct_solar'][idx] +
                      grp.ix[band, ed_wl]['response'] *
                      flux_data.ix[ed_wl, 'direct_solar'][idx]) / 2

        response_sum = (grp.ix[band, st_wl]['response'] +
                        grp.ix[band, ed_wl]['response']) / 2

        # Fuqin's code now loops over each wavelength, in -1 decrements
        # we'll use indices rather than a loop
        response_subs = grp.ix[band].ix[wv_idx]['response'].values
        flux_data_subs = flux_data.ix[wv_idx2]

        response_sum = response_sum + response_subs.sum()

        df.ix[band, 'diffuse'] = (((flux_data_subs['downward_diffuse'].values *
                                    response_subs).sum() + diffuse_bottom) /
                                  response_sum)
        df.ix[band, 'direct'] = (((flux_data_subs['direct_solar'].values *
                                   response_subs).sum() + direct_bottom) /
                                 response_sum)

        # direct solar at top of atmospheric levels
        flux_data_subs = flux_data.ix[wv_idx3]
        df.ix[band, 'direct_top'] = (((flux_data_subs['direct_solar'].values *
                                       response_subs).sum() + direct_top) /
                                     response_sum)

        if transmittance:
            # downward diffuse at top of atmospheric levels
            diffuse_top = (grp.ix[band, st_wl]['response'] *
                           flux_data.ix[st_wl, 'downward_diffuse'][idx] +
                           grp.ix[band, ed_wl]['response'] *
                           flux_data.ix[ed_wl, 'downward_diffuse'][idx]) / 2

            edo_top = flux_data_subs['downward_diffuse'].values
            df.ix[band, 'diffuse_top'] = ((edo_top * response_subs).sum() +
                                          diffuse_top) / response_sum
            t_result = ((df.ix[band, 'diffuse'] + df.ix[band, 'direct']) /
                        (df.ix[band, 'diffuse_top'] +
                         df.ix[band, 'direct_top']))
            df.ix[band, 'transmittance'] = t_result

    df.sort_index(inplace=True)

    return df


def create_solar_irradiance_tables(fname, out_fname, compression='lzf'):
    """
    Writes the accumulated solar irradiance table into a HDF5 file.
    The file is opened in 'a' mode, allowing multiple tables to be
    added. If a table already exists within the file it is removed.
    """
    dset_name = splitext(basename(fname))[0]
    with h5py.File(out_fname, 'a') as fid1:
        # check for a previous run
        if dset_name in fid1:
            del fid1[dset_name]

        with h5py.File(fname, 'r') as fid2:
            fid2.copy(dset_name, fid1)

    return


def link_atmospheric_results(input_targets, out_fname, npoints):
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

    :return:
        None. Results from each file in `input_targets` are linked
        into the output file.
    """
    nbar_atmospherics = False
    sbt_atmospherics = False
    for fname in input_targets:
        with h5py.File(fname.path, 'r') as fid:
            point = fid.attrs['point']
            albedo = fid.attrs['albedo']

        if albedo == Model.sbt.albedos[0]:
            datasets = ['upward-radiation-channel',
                        'downward-radiation-channel']
            sbt_atmospherics = True
        else:
            datasets = ['flux',
                        'altitudes',
                        'solar-irradiance',
                        'channel']
            nbar_atmospherics = True

        grp_path = ppjoin(POINT_FMT.format(p=point),
                          ALBEDO_FMT.format(a=albedo))

        for dset in datasets:
            dname = ppjoin(grp_path, dset)
            create_external_link(fname.path, dname, out_fname, dname)

    with h5py.File(out_fname) as fid:
        fid.attrs['npoints'] = npoints
        fid.attrs['nbar atmospherics'] = nbar_atmospherics
        fid.attrs['sbt atmospherics'] = sbt_atmospherics

    return
