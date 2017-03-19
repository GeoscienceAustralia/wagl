"""
MODTRAN drivers
---------------

"""
from __future__ import absolute_import, print_function, unicode_literals
import os
from os.path import join as pjoin, exists, dirname, basename, splitext
from posixpath import join as ppjoin
import subprocess
import glob

import numpy
from scipy.io import FortranFile
import h5py
import pandas as pd

from gaip.hdf5 import write_dataframe
from gaip.hdf5 import read_table
from gaip.modtran_profiles import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO
from gaip.modtran_profiles import MIDLAT_SUMMER_TRANSMITTANCE
from gaip.modtran_profiles import TROPICAL_TRANSMITTANCE


POINT_FMT = 'point-{p}'
ALBEDO_FMT = 'albedo-{a}'
POINT_ALBEDO_FMT = ''.join([POINT_FMT, '-', ALBEDO_FMT])


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
    with open(out_fname, 'w') as src:
        src.writelines(response)


def _format_tp5(acquisition, satellite_solar_angles_fname,
                longitude_fname, latitude_fname, ancillary_fname, out_fname,
                npoints, albedos):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_angles_fname, 'r') as sat_sol,\
        h5py.File(longitude_fname, 'r') as lon_ds,\
        h5py.File(latitude_fname, 'r') as lat_ds,\
        h5py.File(ancillary_fname, 'r') as anc_ds,\
        h5py.File(out_fname, 'w') as fid:

        # angles data
        view_dset = sat_sol['satellite-view']
        azi_dset = sat_sol['satellite-azimuth']
        coord_dset = sat_sol['coordinator']
        lon_dset = lon_ds['longitude']
        lat_dset = lat_ds['latitude']

        # ancillary data
        aerosol = anc_ds['aerosol'][()]
        water_vapour = anc_ds['water-vapour'][()]
        ozone = anc_ds['ozone'][()]
        elevation = anc_ds['elevation'][()]

        tp5_data, metadata = format_tp5(acquisition, coord_dset, view_dset,
                                        azi_dset, lat_dset, lon_dset, ozone,
                                        water_vapour, aerosol, elevation,
                                        npoints, albedos)

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

    return tp5_data


def format_tp5(acquisition, coordinator, view_dataset, azi_dataset,
               lat_dataset, lon_dataset, ozone, vapour, aerosol, elevation,
               npoints, albedos):
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
    for i in range(npoints):
        for alb in albedos:
            input_data = {'water': vapour,
                          'ozone': ozone,
                          'filter_function': filter_file,
                          'visibility': -aerosol,
                          'elevation': elevation,
                          'sat_height': altitude,
                          'sat_view': view_cor[i],
                          'doy': doy,
                          'binary': binary}
            if alb == 't':
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

    return tp5_data, metadata


def _run_modtran(modtran_exe, workpath, point, albedo, out_fname,
                 compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    rfid = run_modtran(modtran_exe, workpath, point, albedo, out_fname,
                       compression)

    rfid.close()
    return
        

def run_modtran(modtran_exe, workpath, point, albedo, out_fname=None,
                compression='lzf'):
    """
    Run MODTRAN and return the flux and channel results.
    """
    subprocess.check_call([modtran_exe], cwd=workpath)
    flux_fname = glob.glob(pjoin(workpath, '*_b.flx'))[0]
    chn_fname = glob.glob(pjoin(workpath, '*.chn'))[0]

    flux_data, altitudes = read_modtran_flux(flux_fname)

    # read the channel file
    channel_data = read_modtran_channel(chn_fname)

    # initial attributes
    attrs = {'Point': point, 'Albedo': albedo}

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('modtran-results.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    # base group pathname
    group_path = ppjoin(POINT_FMT.format(p=point), ALBEDO_FMT.format(a=albedo))

    # ouput the flux data
    dset_name = ppjoin(group_path, 'flux')
    attrs['Description'] = 'Flux output from MODTRAN'
    write_dataframe(flux_data, dset_name, fid, attrs=attrs)

    # output the channel data
    attrs['Description'] = 'Channel output from MODTRAN'
    dset_name = ppjoin(group_path, 'channel')
    write_dataframe(channel_data, dset_name, fid, attrs=attrs)

    # output the altitude data
    attrs['Description'] = 'Altitudes output from MODTRAN'
    attrs['altitude levels'] = altitudes.shape[0]
    attrs['units'] = 'km'
    dset_name = ppjoin(group_path, 'altitudes')
    write_dataframe(altitudes, dset_name, fid, attrs=attrs)

    return fid


def _calculate_coefficients(accumulated_fname, npoints, out_fname,
                            compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(accumulated_fname, 'r') as fid:

        # initialise dicts to hold the data for each point
        accumulation_albedo_0 = {}
        accumulation_albedo_1 = {}
        accumulation_albedo_t = {}
        channel_data = {}

        for point in range(npoints):
            grp_path = ppjoin(POINT_FMT, ALBEDO_FMT)
            albedo_0_path = ppjoin(grp_path.format(p=point, a='0'),
                                   'solar-irradiance')
            albedo_1_path = ppjoin(grp_path.format(p=point, a='1'),
                                   'solar-irradiance')
            albedo_t_path = ppjoin(grp_path.format(p=point, a='t'),
                                   'solar-irradiance')
            channel_path = ppjoin(grp_path.format(p=point, a='0'), 'channel')

            accumulation_albedo_0[point] = read_table(fid, albedo_0_path)
            accumulation_albedo_1[point] = read_table(fid, albedo_1_path)
            accumulation_albedo_t[point] = read_table(fid, albedo_t_path)
            channel_data[point] = read_table(fid, channel_path)
        
        rfid = calculate_coefficients(accumulation_albedo_0,
                                      accumulation_albedo_1,
                                      accumulation_albedo_t, channel_data,
                                      npoints, out_fname, compression)

    rfid.close()
    return
        

def calculate_coefficients(accumulation_albedo_0, accumulation_albedo_1,
                           accumulation_albedo_t, channel_data, npoints,
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

    :param accumulation_albedo_1:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the solar accumulated
        irradiance (for albedo 1) and structured as returned by the
        `calculate_solar_radiation` function.

    :param accumulation_albedo_t:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the solar accumulated
        irradiance (for albeod t; transmittance) and structured as
        returned by the `calculate_solar_radiation` function.

    :param channel_data:
        A `dict` containing range(npoints) as the keys, and the values
        a `pandas.DataFrame` containing the channel data for that
        point, and structured as returned by the
        `read_modtran_channel` function.

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

        Containing a dataset formatted to the HDF5 TABLE specification
        named:

        * coefficients
    """
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

    attrs = {}
    attrs['Description'] = ("Coefficients derived from the "
                            "accumulated solar irradiation.")
    attrs['Number of atmospheric points'] = npoints

    # Initialise the output file
    if out_fname is None:
        fid = h5py.File('coefficients.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    write_dataframe(result, 'coefficients', fid, compression,
                    attrs=attrs)

    return fid


def read_spectral_response(fname, as_list=False):
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

    :return:
        A `pd.DataFrame` containing the spectral response
        function.
    """
    if isinstance(fname, file):
        lines = fname.readlines()
    else:
        # open the text file
        with open(fname, 'r') as src:
            lines = src.readlines()

    if as_list:
        return lines

    lines = [line.strip() for line in lines]

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

    wavelengths = range(2600, 349, -1)
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


def read_modtran_channel(fname):
    """
    Read a MODTRAN output `*.chn` ascii file.

    :param fname:
        A `str` containing the full file pathname of the channel
        data file.

    :return:
        A `pandas.DataFrame` containing the channel data, and index
        by the `band_id`.
    """
    chn_data = pd.read_csv(fname, skiprows=5, header=None,
                           delim_whitespace=True)
    chn_data['band_id'] = chn_data[20] + ' ' + chn_data[21].astype(str)
    chn_data.drop([20, 21], inplace=True, axis=1)
    chn_data.set_index('band_id', inplace=True)
    chn_data.columns = chn_data.columns.astype(str)

    return chn_data


def _calculate_solar_radiation(acquisition, flux_fnames, out_fname,
                               compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    description = ("Accumulated solar irradiation for point {} "
                   "and albedo {}.")

    spectral_response = acquisition.spectral_response()

    with h5py.File(out_fname, 'w') as fid:
        for key in flux_fnames:
            flux_fname = flux_fnames[key].path
            point, albedo = key
            group_path = ppjoin(POINT_FMT.format(p=point),
                                ALBEDO_FMT.format(a=albedo))
            transmittance = True if albedo == 't' else False
            flux_dataset_name = ppjoin(group_path, 'flux')
            atmos_dataset_name = ppjoin(group_path, 'altitudes')
            channel_dname = ppjoin(group_path, 'channel')

            # link in the channel data for later easy access
            fid[channel_dname] = h5py.ExternalLink(flux_fname, channel_dname)

            # retrieve the flux data and the number of atmospheric levels
            with h5py.File(flux_fname, 'r') as fid2:
                flux_data = read_table(fid2, flux_dataset_name)
                levels = fid2[atmos_dataset_name].attrs['altitude levels']

            # accumulate solar irradiance
            df = calculate_solar_radiation(flux_data, spectral_response,
                                           levels, transmittance)

            # output
            dataset_name = ppjoin(group_path, 'solar-irradiance')
            attrs = {'Description': description.format(point, albedo),
                     'Point': point,
                     'Albedo': albedo}
            write_dataframe(df, dataset_name, fid, compression, attrs=attrs)

    return


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
