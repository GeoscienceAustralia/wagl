"""
MODTRAN drivers
---------------

"""
import os
from os.path import join as pjoin, exists, abspath, dirname
import subprocess

import numpy
from scipy.io import FortranFile
import pandas
import rasterio
import gaip
from gaip import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO
from gaip import MIDLAT_SUMMER_TRANSMITTANCE, TROPICAL_TRANSMITTANCE

BIN_DIR = abspath(pjoin(dirname(__file__), '..', 'bin'))


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


def create_satellite_filter_file(acquisitions, satfilter_path, target):
    """Generate satellite filter input file."""
    refbands = [a for a in acquisitions if a.band_type == gaip.REF]
    filterfile = acquisitions[0].spectral_filter_file
    filterpath = os.path.join(satfilter_path, filterfile)

    with open(target, 'w') as outfile:
        outfile.write("%i\n" % len(refbands))
        outfile.write("%s\n" % filterpath)

    return target


# TODO: once validated, this can function can be deprecated
def write_modtran_input(acquisitions, modtran_input_file, ozone, vapour,
                        aerosol, elevation):
    """Generate modtran input file."""
    acq = acquisitions[0]
    filter_file = acq.spectral_filter_file
    cdate = acq.scene_centre_date
    altitude = acq.altitude / 1000.0  # in km
    dechour = acq.decimal_hour

    with open(modtran_input_file, 'w') as outfile:
        outfile.write("%f\n" % ozone)
        outfile.write("%f\n" % vapour)
        outfile.write("DATA/%s\n" % filter_file)
        outfile.write("-%f\n" % aerosol)
        outfile.write("%f\n" % elevation)
        outfile.write("Annotation, %s\n" % cdate.strftime('%Y-%m-%d'))
        outfile.write("%d\n" % altitude)
        outfile.write("%d\n" % int(cdate.strftime('%j')))
        outfile.write("%f\n" % dechour)


# TODO: once validated, this can function can be deprecated
# as we can write direct to the tp5 template
def write_modtran_inputs(acquisition, coordinator, view_fname, azi_fname,
                         lat_fname, lon_fname, ozone, vapour, aerosol,
                         elevation, coords, albedos, out_fname_fmt):
    filter_file = acquisition.spectral_filter_file
    cdate = acquisition.scene_centre_date
    altitude = acquisition.altitude / 1000.0  # in km
    dechour = acquisition.decimal_hour
    coord = pandas.read_csv(coordinator, header=None, sep=r'\s+\s+',
                            engine='python', names=['row', 'col'])

    with rasterio.open(view_fname) as view_ds,\
        rasterio.open(azi_fname) as azi_ds,\
        rasterio.open(lat_fname) as lat_ds,\
        rasterio.open(lon_fname) as lon_ds:

        npoints = len(coords)
        view = numpy.zeros(npoints, dtype='float32')
        azi = numpy.zeros(npoints, dtype='float32')
        lat = numpy.zeros(npoints, dtype='float64')
        lon = numpy.zeros(npoints, dtype='float64')

        for i in range(1, npoints + 1):
            yidx = coord['row'][i]
            xidx = coord['col'][i]
            idx = ((yidx -1, yidx), (xidx -1, xidx))
            view[i-1] = view_ds.read(1, window=idx)[0, 0]
            azi[i-1] = azi_ds.read(1, window=idx)[0, 0]
            lat[i-1] = lat_ds.read(1, window=idx)[0, 0]
            lon[i-1] = lon_ds.read(1, window=idx)[0, 0]

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

    for i, p in enumerate(coords):
        for alb in albedos:
            out_fname = out_fname_fmt.format(coord=p, albedo=alb)
            with open(out_fname, 'w') as src:
                src.write("{:.8f}\n".format(float(alb)))
                src.write("{:.14f}\n".format(ozone))
                src.write("{:.14f}\n".format(vapour))
                src.write("DATA/{}\n".format(filter_file))
                src.write("-{:.14f}\n".format(aerosol))
                src.write("{:.14f}\n".format(elevation))
                src.write("Annotation, {}\n".format(cdate.strftime('%Y-%m-%d')))
                src.write("{:.14f}\n".format(altitude))
                src.write("{:f}\n".format(view_cor[i]))
                src.write("{:d}\n".format(int(cdate.strftime('%j'))))
                src.write("{:.14f}\n".format(lat[i]))
                src.write("{:.14f}\n".format(rlon[i]))
                src.write("{:.14f}\n".format(dechour))
                src.write("{:f}\n".format(azi_cor[i]))


def write_tp5(acquisition, coordinator, view_fname, azi_fname,
              lat_fname, lon_fname, ozone, vapour, aerosol, elevation,
              coords, albedos, out_fname_fmt):
    """Writes the tp5 files for the albedo (0, 1) and transmittance (t)."""
    geobox = acquisition.gridded_geo_box()
    filter_file = acquisition.spectral_filter_file
    cdate = acquisition.scene_centre_date
    doy = int(cdate.strftime('%j'))
    altitude = acquisition.altitude / 1000.0  # in km
    dechour = acquisition.decimal_hour
    coord = pandas.read_csv(coordinator, header=None, sep=r'\s+\s+',
                            engine='python', names=['row', 'col'])

    with rasterio.open(view_fname) as view_ds,\
        rasterio.open(azi_fname) as azi_ds,\
        rasterio.open(lat_fname) as lat_ds,\
        rasterio.open(lon_fname) as lon_ds:

        npoints = len(coords)
        view = numpy.zeros(npoints, dtype='float32')
        azi = numpy.zeros(npoints, dtype='float32')
        lat = numpy.zeros(npoints, dtype='float64')
        lon = numpy.zeros(npoints, dtype='float64')

        for i in range(1, npoints + 1):
            yidx = coord['row'][i]
            xidx = coord['col'][i]
            idx = ((yidx -1, yidx), (xidx -1, xidx))
            view[i-1] = view_ds.read(1, window=idx)[0, 0]
            azi[i-1] = azi_ds.read(1, window=idx)[0, 0]
            lat[i-1] = lat_ds.read(1, window=idx)[0, 0]
            lon[i-1] = lon_ds.read(1, window=idx)[0, 0]

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
    centre_lon, centre_lat = geobox.centre_lonlat
    if centre_lat < -23.0:
        albedo_profile = MIDLAT_SUMMER_ALBEDO
        trans_profile = MIDLAT_SUMMER_TRANSMITTANCE
    else:
        albedo_profile = TROPICAL_ALBEDO
        trans_profile = TROPICAL_TRANSMITTANCE

    # we'll only cater for MODTRAN to output binary form
    binary = 'T'

    # write the tp5 files required for input into MODTRAN
    for i, p in enumerate(coords):
        for alb in albedos:
            out_fname = out_fname_fmt.format(coord=p, albedo=alb)
            if alb == 't':
                data = trans_profile.format(albedo=0.0,
                                            water=vapour,
                                            ozone=ozone,
                                            filter_function=filter_file,
                                            visibility=-aerosol,
                                            elevation=elevation,
                                            sat_height=altitude,
                                            sat_view=view_cor[i],
                                            doy=doy,
                                            sat_view_offset=180.0-view_cor[i],
                                            binary=binary)
            else:
                data = albedo_profile.format(albedo=float(alb),
                                             water=vapour,
                                             ozone=ozone,
                                             filter_function=filter_file,
                                             visibility=-aerosol,
                                             elevation=elevation,
                                             sat_height=altitude,
                                             sat_view=view_cor[i],
                                             doy=doy,
                                             lat=lat[i],
                                             lon=rlon[i],
                                             time=dechour,
                                             sat_azimuth=azi_cor[i],
                                             binary=binary)
            with open(out_fname, 'w') as src:
                src.write(data)


def run_modtran(modtran_exe, workpath):
    """Run MODTRAN."""
    subprocess.check_call([modtran_exe], cwd=workpath)


def calculate_coefficients(coords, chn_input_fmt, dir_input_fmt,
                           output_fmt, output_reformat, cwd):
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

    :param output_fmt:
        A `string` format for the output filename.
        eg '{coord}_alb.txt'. If set to `None`, then a `dictionary`
        with the `coords` as the keys will be returned.

    :param output_reformat:
        A `string` format for the reformatted output file.
        eg '{factor}_out_b_{band}.txt'. If set to `None`, then a
        `dictionary` with the combination of (band, factor) as the
        keys will be returned.

    :param cwd:
        A `string` containing the full file pathname to the MODTRAN
        output directory.

    :return:
        A `tuple` of two `dictionaries`. The first containing the
        calculated coefficients with the `coords` as the keys.
        The second containing the reformated coefficients and
        a `tuple` of (band, factor) as the keys.
    """

    result = {}
    for coord in coords:
        # MODTRAN output .chn file (albedo 0)
        fname1 = pjoin(cwd, chn_input_fmt.format(coord=coord, albedo=0))

        # **********UNUSED**********
        # MODTRAN output .chn file (albedo 1)
        # fname2 = pjoin(cwd, chn_input_fmt.format(coord=coord, albedo=1))

        # solar radiation file (albedo 0)
        fname3 = pjoin(cwd, dir_input_fmt.format(coord=coord, albedo=0))

        # solar radiation file (albedo 1)
        fname4 = pjoin(cwd, dir_input_fmt.format(coord=coord, albedo=1))

        # solar radiation file (transmittance mode)
        fname5 = pjoin(cwd, dir_input_fmt.format(coord=coord, albedo='t'))

        # output file
        if output_fmt is not None:
            out_fname = pjoin(cwd, output_fmt.format(coord=coord))

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

        # output to disk; tab delimited
        if output_fmt is not None:
            df.to_csv(out_fname, sep='\t', index=False)

        result[coord] = df

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

            if output_reformat is not None:
                out_fname = pjoin(cwd, output_reformat.format(factor=factor,
                                                              band=bn))
                df_reformat.to_csv(out_fname, index=False, header=False,
                                   sep='\t')

            result2[(bn, factor)] = df_reformat

    return result, result2


def bilinear_interpolate(acqs, factors, coordinator, boxline, centreline,
                         input_fmt, output_fmt, workpath):
    """Perform bilinear interpolation."""

    bands = [a.band_num for a in acqs]
    geobox = gaip.gridded_geo_box(acqs[0])
    cols, rows = geobox.get_shape_xy()

    # dataframes for the coords, scene centreline, boxline
    coords = pandas.read_csv(coordinator, header=None, sep=r'\s+\s+',
                             engine='python', skiprows=1, names=['row', 'col'])
    cent = pandas.read_csv(centreline, skiprows=2, header=None, sep=r'\s+\s+',
                           engine='python',
                           names=['line', 'centre', 'npoints', 'lat', 'lon'])
    box = pandas.read_csv(boxline, header=None, sep=r'\s+\s+', engine='python',
                          names=['line', 'cstart', 'cend'])

    coord = numpy.zeros((9, 2), dtype='int')
    coord[:, 0] = coords.row.values
    coord[:, 1] = coords.col.values
    centre = cent.centre.values
    start = box.cstart.values
    end = box.cend.values

    # Initialise the dict to store the locations of the bilinear outputs
    bilinear_outputs = {}

    for band in bands:
        for factor in factors:
            fname = output_fmt.format(factor=factor, band=band)
            fname = pjoin(workpath, fname)
            atmospheric_fname = input_fmt.format(factor=factor, band=band)
            bilinear_outputs[(band, factor)] = fname

            # atmospheric paramaters
            atmos = pandas.read_csv(atmospheric_fname, header=None,
                                    sep='\t', engine='python',
                                    names=['s1', 's2', 's3', 's4'])

            # get the individual atmospheric components
            s1 = atmos.s1.values
            s2 = atmos.s2.values
            s3 = atmos.s3.values
            s4 = atmos.s4.values

            res = numpy.zeros((rows, cols), dtype='float32')
            gaip.bilinear(cols, rows, coord, s1, s2, s3, s4, start, end,
                          centre, res.transpose())

            # Output the result to disk
            gaip.write_img(res, fname, fmt='GTiff', geobox=geobox, nodata=-999,
                           compress='deflate', options={'zlevel': 1})

    return bilinear_outputs


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
