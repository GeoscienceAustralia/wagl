"""
MODTRAN drivers
---------------

"""
import os
from os.path import join as pjoin, exists, abspath, dirname, basename, splitext
from posixpath import join as ppjoin
import subprocess

import numpy
from scipy.io import FortranFile
import h5py
import pandas
import rasterio
import gaip
from gaip import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO
from gaip import MIDLAT_SUMMER_TRANSMITTANCE, TROPICAL_TRANSMITTANCE
from gaip import attach_image_attributes
from gaip import dataset_compression_kwargs
from gaip import write_h5_image
from gaip import write_dataframe


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


def prepare_modtran(coordinate, albedo, modtran_work, modtran_exe):
    """
    Prepares the working directory for a MODTRAN execution.
    """
    data_dir = pjoin(dirname(modtran_exe), 'DATA')
    if not exists(data_dir):
        raise OSError('Cannot find MODTRAN')

    out_fname = pjoin(modtran_work, 'mod5root.in')

    with open(out_fname, 'w') as src:
        src.write(coordinate + '_alb_' + albedo + '\n')

    symlink_dir = pjoin(modtran_work, 'DATA')
    if exists(symlink_dir):
        os.unlink(symlink_dir)

    os.symlink(data_dir, symlink_dir)


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
        dataset_fmt = '{point}/alb_{albedo}/tp5_data'

        for key in metadata:
            dname = dataset_fmt.format(point=key[0], albedo=key[1])
            str_data = numpy.string_(tp5_data[key])
            dset = group.create_dataset(dname, data=str_data)
            for k in metadata[key]:
                dset.attrs[k] = metadata[key][k]

    return tp5_data


def format_tp5(acquisition, coordinator, view_dataset, azi_dataset,
               lat_dataset, lon_dataset, ozone, vapour, aerosol, elevation,
               npoints, albedos):
    """Creates str formatted tp5 files for the albedo (0, 1) and transmittance (t)."""
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

    for i in range(1, npoints + 1):
        yidx = coordinator['row_index'][i]
        xidx = coordinator['col_index'][i]
        idx = (slice(yidx -1, yidx), slice(xidx -1, xidx))
        view[i-1] = view_dataset[idx][0, 0]
        azi[i-1] = azi_dataset[idx][0, 0]
        lat[i-1] = lat_dataset[idx][0, 0]
        lon[i-1] = lon_dataset[idx][0, 0]

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


def run_modtran(modtran_exe, workpath):
    """Run MODTRAN."""
    subprocess.check_call([modtran_exe], cwd=workpath)


def _calculate_coefficients(npoints, chn_input_fmt, accumulated_fname, mod_root,
                            out_fname, compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    # TODO: read the accumlated data
    # with h5py.File(accumulated_fname, 'r') as fid:
        
    res1, res2 = calculate_coefficients(npoints, chn_input_fmt,
                                        accumulated_fname, mod_root)
    attrs1 = {}
    attrs1['Description'] = ("Coefficients derived from the "
                             "accumulated solar irradiation.")
    attrs2 = {}
    attrs2['Description'] = ("Coefficients derived from the "
                             "accumulated solar irradiation, and formatted "
                             "into a 4x4 grid per MODTRAN factor.")
    attrs2['Grid Layout'] = ("TL, TM, ML, MM\n"
                             "TM, TR, MM, MR\n"
                             "ML, MM, BL, BM\n"
                             "MM, MR, BM, BR")

    with h5py.File(out_fname, 'w') as fid:
        write_dataframe(res1, 'coefficients_format_1', fid, compression,
                        attrs=attrs1)
        write_dataframe(res2, 'coefficients_format_2', fid, compression,
                        attrs=attrs2)

    return
        

def calculate_coefficients(npoints, chn_input_fmt, dir_input_fmt, mod_root):
    """
    Calculate the atmospheric coefficients from the MODTRAN output
    and used in the BRDF and atmospheric correction.
    Coefficients are computed for each band for each each coordinate
    for each factor.  The factors are:
    ['fs', 'fv', 'a', 'b', 's', 'dir', 'dif', 'ts'].

    :param acqs:
        A `list` of acquisitions.

    :param coords:
        A `list` of `string` coordinates indicating the location
        within an array, eg.
        ["TL", "TM", "TR", "ML", "MM", "MR", "BL", "BM", "BR"]

    :param chn_input_fmt:
        A `string` format for the MODTRAN *.chn output file.
        eg '{coord}/alb_{albedo}/{coord}_alb_{albedo}.chn'.

    :param dir_input_fmt:
        A `string` format for the MODTRAN *.dir output file.
        eg '{coord}_alb_{albedo}.dir'.

    :param mod_root:
        A `string` containing the full file pathname to the MODTRAN
        output directory.

    :return:
        A `pandas.DataFrame` object containing the coefficients of each
        factor for each band.
    """
    result = {}
    for point in range(npoints):
        # MODTRAN output .chn file (albedo 0)
        fname1 = pjoin(mod_root, chn_input_fmt.format(point=point, albedo=0))

        # **********UNUSED**********
        # MODTRAN output .chn file (albedo 1)
        # fname2 = pjoin(cwd, chn_input_fmt.format(coord=coord, albedo=1))

        # solar radiation file (albedo 0)
        fname3 = pjoin(mod_root, dir_input_fmt.format(point=point, albedo=0))

        # solar radiation file (albedo 1)
        fname4 = pjoin(mod_root, dir_input_fmt.format(point=point, albedo=1))

        # solar radiation file (transmittance mode)
        fname5 = pjoin(mod_root, dir_input_fmt.format(point=point, albedo='t'))

        # read the data
        # **********UNUSED**********
        # data2 = pandas.read_csv(fname2, skiprows=5, header=None,
        #                         delim_whitespace=True)

        # solar radiation (albedo 0, albedo 1, transmittance mode)
        data3 = pandas.read_csv(fname3, header=0, sep='\t')
        data4 = pandas.read_csv(fname4, header=0, sep='\t')
        data5 = pandas.read_csv(fname5, header=0, sep='\t')

        # set the index to be the band name
        # we didn't write the index out previously as we'll try to keep
        # the same format so Fuqin can use it within her code
        data3.set_index('band', inplace=True, drop=False)
        data4.set_index('band', inplace=True, drop=False)
        data5.set_index('band', inplace=True, drop=False)

        # MODTRAN output .chn file (albedo 0)
        data1 = pandas.read_csv(fname1, skiprows=5, header=None,
                                delim_whitespace=True, nrows=data3.shape[0])

        fmt = 'BAND {}'
        band_idx = [fmt.format(val) for key, val in data1[21].iteritems()]
        data1['band'] = band_idx
        data1.set_index('band', inplace=True, drop=False)

        # calculate
        diff_0 = data3['diffuse'] * 10000000.0
        diff_1 = data4['diffuse'] * 10000000.0
        dir_0 = data3['direct'] * 10000000.0
        dir_1 = data4['direct'] * 10000000.0
        dir_t = data5['direct']
        dir0_top = data3['directtop'] * 10000000.0
        dirt_top = data5['directtop']
        tv_total = data5['transmittance']
        ts_total = (diff_0 + dir_0) / dir0_top
        ts_dir = dir_0 / dir0_top
        tv_dir = dir_t / dirt_top

        columns = ['band',
                   'fs',
                   'fv',
                   'a',
                   'b',
                   's',
                   'dir',
                   'dif',
                   'ts']
        df = pandas.DataFrame(columns=columns, index=band_idx)

        df['band'] = data1[21]
        df['fs'] = ts_dir / ts_total
        df['fv'] = tv_dir / tv_total
        df['a'] = (diff_0 + dir_0) / numpy.pi * tv_total
        df['b'] = data1[3] * 10000000
        df['s'] = 1 - (diff_0 + dir_0) / (diff_1 + dir_1)
        df['dir'] = dir_0
        df['dif'] = diff_0
        df['ts'] = ts_dir

        result[point] = df

    # set the band numbers as a searchable index
    for key in result:
        result[key].set_index('band', inplace=True, drop=False)

    # combine all results into a single pandas.DataFrame
    df = pandas.concat(result)
    groups = df.groupby('band')

    # reformat the derived coefficients into another format layout
    # specifically a 4x4 layout, for each factor, for each band
    """
    TL, TM, ML, MM
    TM, TR, MM, MR
    ML, MM, BL, BM
    MM, MR, BM, BR
    """

    factors = columns[1:]

    result2 = {}
    for bn, grp in groups:
        for factor in factors:
            s1 = [grp.ix[('TL', bn)][factor],
                  grp.ix[('TM', bn)][factor],
                  grp.ix[('ML', bn)][factor],
                  grp.ix[('MM', bn)][factor]]
            s2 = [grp.ix[('TM', bn)][factor],
                  grp.ix[('TR', bn)][factor],
                  grp.ix[('MM', bn)][factor],
                  grp.ix[('MR', bn)][factor]]
            s3 = [grp.ix[('ML', bn)][factor],
                  grp.ix[('MM', bn)][factor],
                  grp.ix[('BL', bn)][factor],
                  grp.ix[('BM', bn)][factor]]
            s4 = [grp.ix[('MM', bn)][factor],
                  grp.ix[('MR', bn)][factor],
                  grp.ix[('BM', bn)][factor],
                  grp.ix[('BR', bn)][factor]]

            sdata = {'s1': s1, 's2': s2, 's3': s3, 's4': s4}
            df_reformat = pandas.DataFrame(sdata)

            result2[(bn, factor)] = df_reformat

    # re-create basic single index tables
    result = pandas.concat(result, names=['coordinate', 'band'])
    result.reset_index(level='coordinate')
    result.reset_index(inplace=True)

    result2 = pandas.concat(result2, names=['band_number', 'factor'])
    result2.reset_index(level=['band_number', 'factor'], inplace=True)
    result2.reset_index(inplace=True)

    return result, result2


def _bilinear_interpolate(acq, factor, sat_sol_angles_fname,
                          coefficients_fname, out_fname, compression):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    band = acq.band_num
    geobox = acq.gridded_geo_box()

    with h5py.File(sat_sol_angles_fname, 'r') as sat_sol,\
        h5py.File(coefficients_fname, 'r') as coef:

        coord_dset = sat_sol['coordinator']
        centre_dset = sat_sol['centreline']
        box_dset = sat_sol['boxline']

        coef_dset = coef['coefficients_format_2']
        wh = ((coef_dset['band_number'] == band) &
              (coef_dset['factor'] == factor))
        coefficients = pandas.DataFrame(coef_dset[wh])

        data = bilinear_interpolate(acq, coord_dset, box_dset, centre_dset,
                                    coefficients)

    with h5py.File(out_fname, 'w') as fid:
        dset_name = splitext(basename(out_fname))[0]
        kwargs = dataset_compression_kwargs(compression=compression,
                                            chunks=(1, geobox.x_size()))
        no_data = -999
        kwargs['fillvalue'] = no_data
        attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
                 'geotransform': geobox.affine.to_gdal(),
                 'no_data_value': no_data}
        desc = ("Contains the bi-linearly interpolated result of factor {}"
                "for band {} from sensor {}.")
        attrs['Description'] = desc.format(factor, band, acq.satellite_name)
        write_h5_image(data, dset_name, fid, attrs, **kwargs)

    return


def bilinear_interpolate(acq, coordinator_dataset, boxline_dataset,
                         centreline_dataset, coefficients):
    """Perform bilinear interpolation."""
    geobox = acq.gridded_geo_box()
    cols, rows = geobox.get_shape_xy()

    coord = numpy.zeros((9, 2), dtype='int')
    coord[:, 0] = coordinator_dataset['row_index'][:]
    coord[:, 1] = coordinator_dataset['col_index'][:]
    centre = boxline_dataset['bisection'][:]
    start = boxline_dataset['start'][:]
    end = boxline_dataset['end'][:]

    # get the individual atmospheric components
    s1 = coefficients['s1'][:]
    s2 = coefficients['s2'][:]
    s3 = coefficients['s3'][:]
    s4 = coefficients['s4'][:]

    res = numpy.zeros((rows, cols), dtype='float32')
    gaip.bilinear(cols, rows, coord, s1, s2, s3, s4, start, end,
                  centre, res.transpose())

    return res


def read_spectral_response(fname):
    """
    Read the spectral response function text file used during
    MODTRAN processing.

    :param fname:
        A `str` containing the full file path name.

    :return:
        A `pandas.DataFrame` containing the spectral response
        function.
    """
    # open the text file
    with open(fname, 'r') as src:
        lines = src.readlines()

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
        df = pandas.DataFrame({'band_description': lines[idx],
                               'wavelength': data[:, 0],
                               'response': data[:, 1]})
        response[lines[idx]] = df

    # get spectral response data for band n
    idx = ids[-1]
    data = numpy.array([l.split('  ') for l in lines[idx+1:]], dtype='float')
    df = pandas.DataFrame({'band_description': lines[idx],
                           'wavelength': data[:, 0],
                           'response': data[:, 1]})
    response[lines[idx]] = df

    wavelengths = range(2600, 349, -1)
    for band in response:
        base_df = pandas.DataFrame({'wavelength': wavelengths,
                                    'response': 0.0,
                                    'band_description': band},
                                   index=wavelengths)
        df = response[band]
        base_df.ix[df['wavelength'], 'response'] = df['response'].values

        response[band] = base_df

    spectral_response = pandas.concat(response)

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
    flux_dtype = numpy.dtype([('wavelength', 'float64'),
                              ('upward_diffuse', 'float64'),
                              ('downwar_diffuse', 'float64'),
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
            df = pandas.DataFrame(numpy.zeros(levels, dtype=flux_dtype))
            df['wavelength'] = data['wavelength'][0]
            df['upward_difuse'] = data['flux_data'].squeeze()[:, 0]
            df['downward_diffuse'] = data['flux_data'].squeeze()[:, 1]
            df['direct_solar'] = data['flux_data'].squeeze()[:, 2]
            flux[wv] = df

    # concatenate into a single table
    flux_data = pandas.concat(flux)

    return flux_data, altitude


def _calculate_solar_radiation(flux_fname, response_fname, transmittance,
                               out_fname, compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    df = calculate_solar_radiation(flux_fname, response_fname, transmittance)
    dset_name = splitext(basename(out_fname))[0]
    coordinate, _, albedo = dset_name.split('-')

    with h5py.File(out_fname, 'w') as fid:

        description = ("Accumulated solar irradiation for coordinate {} "
                       "and albedo {}.")
        attrs = {'Description': description.format(coordinate, albedo),
                 'Coordinate': coordinate,
                 'Albedo': albedo}
        write_dataframe(df, dset_name, fid, compression, attrs=attrs)

    return


def calculate_solar_radiation(flux_fname, response_fname, transmittance=False):
    """
    Retreive the flux data from the MODTRAN output `*.flx`, and
    calculate the solar radiation.

    The solar radiation will be calculated for each of the bands
    contained within the spectral response dataset.

    :param flux_fname:
        A `str` containing the full file path name of the flux dataset
        output by MODTRAN.

    :param response_fname:
        A `str` containing the full file path name of the spectral
        response dataset.

    :param transmittance:
        If set to `True`, then calculate the solar radiation in
        transmittance mode. Default is to calculate from albedo.

    :return:
        A `pandas.DataFrame` containing the solar radiation
        accumulation.
    """
    # get the flux and spectral response datasets
    flux_data, altitudes = read_modtran_flux(flux_fname)
    response = read_spectral_response(response_fname)

    # No. of atmospheric levels & index location of the top atmospheric level
    levels = altitudes.shape[0]
    idx = levels - 1

    # group via the available bands
    groups = response.groupby('band_description')

    # output dataframe
    # later on the process can be refined to only evaluate the bands
    # we wish to process
    if transmittance:
        columns = ['band',
                   'diffuse',
                   'direct',
                   'diffusetop',
                   'directtop',
                   'transmittance']
    else:
        columns = ['band',
                   'diffuse',
                   'direct',
                   'directtop']
    df = pandas.DataFrame(columns=columns, index=groups.groups.keys())

    # indices for flux at bottom and top of atmosphere layers
    wv_idx = range(2599, 350, -1)
    wv_idx2 = [(i, 0) for i in wv_idx]
    wv_idx3 = [(i, idx) for i in wv_idx]

    # loop over each band and get the solar radiation
    for band, grp in groups:
        df.ix[band, 'band'] = band

        # downward diffuse at bottom of atmospheric levels
        diffuse_bottom = (grp.ix[band, 2600]['response'] *
                          flux_data.ix[2600, 'downward_diffuse'][0] +
                          grp.ix[band, 350]['response'] *
                          flux_data.ix[350, 'downward_diffuse'][0]) / 2

        # direct solar at bottom of atmospheric levels
        direct_bottom = (grp.ix[band, 2600]['response'] *
                         flux_data.ix[2600, 'direct_solar'][0] +
                         grp.ix[band, 350]['response'] *
                         flux_data.ix[350, 'direct_solar'][0]) / 2

        # direct solar at top of atmospheric levels
        direct_top = (grp.ix[band, 2600]['response'] *
                      flux_data.ix[2600, 'direct_solar'][idx] +
                      grp.ix[band, 350]['response'] *
                      flux_data.ix[350, 'direct_solar'][idx]) / 2

        response_sum = (grp.ix[band, 2600]['response'] +
                        grp.ix[band, 350]['response']) / 2

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
        df.ix[band, 'directtop'] = (((flux_data_subs['direct_solar'].values *
                                      response_subs).sum() + direct_top) /
                                    response_sum)

        if transmittance:
            # downward diffuse at top of atmospheric levels
            diffuse_top = (grp.ix[band, 2600]['response'] *
                           flux_data.ix[2600, 'downward_diffuse'][idx] +
                           grp.ix[band, 350]['response'] *
                           flux_data.ix[350, 'downward_diffuse'][idx]) / 2

            edo_top = flux_data_subs['downward_diffuse'].values
            df.ix[band, 'diffusetop'] = ((edo_top * response_subs).sum() +
                                         diffuse_top) / response_sum
            t_result = ((df.ix[band, 'diffuse'] + df.ix[band, 'direct']) /
                        (df.ix[band, 'diffusetop'] + df.ix[band, 'directtop']))
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


def link_bilinear_data(data, out_fname):
    """
    Links the individual bilinearly interpolated results into a
    single file for easier access.
    """
    for key in data:
        # band, factor = key
        fname = data[key]
        base_dname = splitext(basename(fname))[0]

        # do we need two group levels?
        # dset_name = ppjoin(band, factor, base_dname)

        with h5py.File(out_fname, 'w') as fid:
            fid[base_dname] = h5py.ExternalLink(fname, base_dname)

    return
