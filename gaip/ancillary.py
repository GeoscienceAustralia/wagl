"""
Ancillary dataset retrieval and storage
"""

import logging
import subprocess
import re
import os
from os.path import join as pjoin, splitext, exists, abspath, dirname
from os.path import pardir, basename
from posixpath import join as ppjoin
import datetime
import glob

import h5py
import pandas
from geopandas import GeoSeries
import numpy
import rasterio
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely import wkt
import gaip
from gaip import attach_attributes
from gaip import write_scalar
from gaip import write_dataframe
from gaip import read_meatadata_tags
from gaip import read_table


log = logging.getLogger()


ECWMF_LEVELS = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200,
                225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
                775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]


def get_4d_idx(day):
    """
    A small utility function for indexing into a 4D dataset
    represented as a 3D dataset.
    [month, level, y, x], where level contains 37 levels, and day
    contains 28, 29, 30 or 31 days.
    """
    start = 1 + 37 * (day - 1)
    stop = start + 37
    return range(start, stop, 1)


def kelvin_2_celcius(kelvin):
    """
    A small utility function for converting degrees Kelvin to
    degrees Celcius.
    """
    return kelvin - 273.15


def relative_humdity(surface_temp, dewpoint_temp, kelvin=True):
    """
    Calculates relative humidity given a surface temperature and
    dewpoint temperature.
    """
    if kelvin:
        surf_t = kelvin_2_celcius(surface_temp)
        dew_t = kelvin_2_celcius(dewpoint_temp)
    else:
        surf_t = surface_temp
        dew_t = dewpoint_temp

    rh = 100 * ((112.0 - 0.1 * surf_t + dew_t) / (112.0 + 0.9 * surf_t)) **8

    return rh


def collect_thermal_ancillary(acquisition, dewpoint_path, temperature_2m_path,
                              surface_pressure_path, geopotential_path,
                              temperature_path, relative_humidity_path,
                              invariant_fname, out_fname=None,
                              compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('ecwmf-ancillary.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    dt = acquisition.scene_center_datetime
    geobox = acquisition.gridded_geo_box()
    lonlat = geobox.centre_lonlat

    # get data located at the surface
    dew = ecwmf_dewpoint_temperature(dewpoint_path, lonlat, dt)
    t2m = ecwmf_temperature_2metre(temperature_2m_path, lonlat, dt)
    sfc_prs = ecwmf_surface_pressure(surface_pressure_path, lonlat, dt)
    sfc_hgt = ecwmf_elevation(invariant_fname, lonlat)
    rh = relative_humdity(t2m[0], dew[0])

    # output the scalar data along with the attrs
    write_scalar(dew[0], 'dewpoint-temperature-2metre', fid, dew[1])
    write_scalar(t2m[0], 'temperature-2metre', fid, t2m[1])
    write_scalar(sfc_prs[0], 'surface-pressure', fid, sfc_prs[1])
    write_scalar(sfc_hgt[0], 'surface-geopotential-height', fid, sfc_hgt[1])
    attrs = {'Description': 'Relative Humidity calculated at the surface'}
    write_scalar(rh, 'surface-relative-humidity', fid, attrs)

    # get the data from each of the pressure levels (1 -> 1000 ISBL)
    gph_df, gph_md = ecwmf_geo_potential(geopotential_path, lonlat, dt)
    tmp_df, tmp_md = ecwmf_temperature(temperature_path, lonlat, dt)
    rh_df, rh_md = ecwmf_relative_humidity(relative_humidity_path, lonlat, dt)

    write_dataframe(gph_df, 'geopotential', fid, compression, attrs=gph_md)
    write_dataframe(tmp_df, 'temperature', fid, compression, attrs=tmp_md)
    write_dataframe(rh_df, 'relative-humidity', fid, compression, attrs=rh_md)

    # combine the surface and higher pressure layers into a single array
    cols = ['GeoPotential_Height', 'Pressure', 'Temperature',
            'Relative_Humidity']
    df = pandas.DataFrame(columns=cols, index=range(rh_df.shape[0]+1))
    df['GeoPotential_Height'].iloc[1:] = gph_df['GeoPotential_Height'].copy()
    df['Pressure'].iloc[1:] = ECWMF_LEVELS
    df['Temperature'].iloc[1:] = tmp_df['Temperature'].copy()
    df['Relative_Humidity'].iloc[1:] = rh_df['Relative_Humidity'].copy()

    # insert the surface level
    df['GeoPotential_Height'].iloc[0] = sfc_hgt[0]
    df['Pressure'].iloc[0] = sfc_prs[0]
    df['Temperature'].iloc[0] = t2m[0]
    df['Relative_Humidity'].iloc[0] = rh

    description = ('Combined Surface and Pressure Layer data retrieved from '
                   'the ECWMF catalogue.')
    attrs = {'Description': description,
             'Date used for querying ECWMF': dt}
    write_dataframe(df, 'atmospheric-profile', fid, compression, attrs=attrs)

    return fid


def collect_ancillary_data(acquisition, aerosol_path, water_vapour_path,
                           ozone_path, dem_path, brdf_path,
                           brdf_premodis_path, out_fname=None,
                           compression='lzf', work_path=''):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    def _format_brdf_attrs(factor):
        """Converts BRDF shortnames to longnames"""
        scattering_names = {'iso': 'isometric',
                            'vol': 'volumetric',
                            'geo': 'geometric'}

        attrs = {}
        description = ("Bidirectional Reflectance Distribution Function "
                       "for the {} scattering fraction.")
        attrs['Description'] = description.format(scattering_names[factor])

        return attrs

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('ancillary.h5', driver='core', backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    dt = acquisition.scene_center_datetime
    geobox = acquisition.gridded_geo_box()

    aerosol = get_aerosol_data(acquisition, aerosol_path)
    write_scalar(aerosol[0], 'aerosol', fid, aerosol[1])

    wv = get_water_vapour(acquisition, water_vapour_path)
    write_scalar(wv[0], 'water-vapour', fid, wv[1])

    ozone = get_ozone_data(acquisition, geobox.centre_lonlat, dt)
    write_scalar(ozone[0], 'ozone', fid, ozone[1])

    elev = get_elevation_data(geobox.centre_lonlat, dem_path)
    write_scalar(elev[0], 'elevation', fid, elev[1])

    # brdf
    group = fid.create_group('brdf-image-datasets')
    data = gaip.get_brdf_data(acquisition, brdf_path, brdf_premodis_path,
                              group, compression=compression,
                              work_path=work_path)
    dname_format = "BRDF-Band-{band}-{factor}"
    for key in data:
        band, factor = key
        attrs = _format_brdf_attrs(factor)
        for k in attrs:
            data[key][k] = attrs[k]
        dname = dname_format.format(band=band, factor=factor)
        fid.create_dataset(dname, data=data[key].pop('value'))
        attach_attributes(fid[dname], attrs=data[key])

    return fid


def aggregate_ancillary(ancillary_fnames, out_fname):
    """
    If the acquisition is part of a `tiled` scene such as Sentinel-2a,
    then we need to average the point measurements gathereed from
    all tiles.
    """
    # Initialise the output files
    if out_fname is None:
        fid1 = h5py.File('ancillary.h5', driver='core', backing_store=False)
    else:
        fid1 = h5py.File(out_fname, 'w')

    # initialise the mean result
    ozone = vapour = aerosol = elevation = 0.0

    n_tiles = len(ancillary_fnames)

    for fname in ancillary_fnames:
        with h5py.File(fname, 'r') as fid2:

            ozone += fid2['ozone'][()]
            vapour += fid2['water-vapour'][()]
            aerosol += fid2['aerosol'][()]
            elevation += fid2['elevation'][()]

    ozone /= n_tiles
    vapour /= n_tiles
    aerosol /= n_tiles
    elevation /= n_tiles

    description = ("The {} value is an average from all the {} values "
                   " retreived for each Granule.")
    attrs = {'data_source': 'granule_average'}

    dset = fid1.create_dataset('ozone', data=ozone)
    attrs['Description'] = description.format('Ozone')
    attach_attributes(dset, attrs)

    dset = fid1.create_dataset('water-vapour', data=vapour)
    attrs['Description'] = description.format('Water Vapour')
    attach_attributes(dset, attrs)

    dset = fid1.create_dataset('aerosol', data=aerosol)
    attrs['Description'] = description.format('Aerosol')
    attach_attributes(dset, attrs)

    dset = fid1.create_dataset('elevation', data=elevation)
    attrs['Description'] = description.format('Elevation')
    attach_attributes(dset, attrs)

    return fid1


def get_aerosol_data_v2(acquisition, aerosol_fname):
    """
    Extract the aerosol value for an acquisition.
    The version 2 retrieves the data from a HDF5 file, and provides
    more control over how the data is selected geo-metrically.
    Better control over timedeltas.
    """

    dt = acquisition.scene_center_datetime
    geobox = acquisition.gridded_geo_box()
    roi_poly = Polygon([geobox.ul_lonlat, geobox.ur_lonlat,
                        geobox.lr_lonlat, geobox.ll_lonlat])

    descr = ['AATSR_PIX', 'AATSR_CMP_YEAR_MONTH', 'AATSR_CMP_MONTH']
    names = ['ATSR_LF_%Y%m', 'aot_mean_%b_%Y_All_Aerosols',
             'aot_mean_%b_All_Aerosols']
    exts = ['/pix', '/cmp', '/cmp']
    pathnames = [ppjoin(ext, dt.strftime(n)) for ext, n in zip(exts, names)]

    fid = h5py.File(aerosol_fname, 'r')

    delta_tolerance = datetime.timedelta(days=0.5)

    data = None
    for pathname, description in zip(pathnames, descr):
        if pathname in fid:
            df = read_table(fid, pathname)
            aerosol_poly = wkt.loads(fid[pathname].attrs['extents'])

            if aerosol_poly.intersects(roi_poly):
                if description == 'AATSR_PIX':
                    abs_diff = (df['timestamp'] - dt).abs()
                    df = df[abs_diff < delta_tolerance]
                    df.reset_index(inplace=True, drop=True)

                if df.shape[0] == 0:
                    continue

                intersection = aerosol_poly.intersection(roi_poly)
                pts = GeoSeries([Point(x, y) for x, y in
                                 zip(df['lon'], df['lat'])])
                idx = pts.intersects(intersection)
                data = df[idx]['aerosol'].mean()

                if numpy.isfinite(data):
                    metadata = {'data_source': description,
                                'dataset_name': pathname,
                                'query_date': dt,
                                'extents': wkt.dumps(intersection)}

                    # ancillary metadata tracking
                    md = gaip.extract_ancillary_metadata(aerosol_fname)
                    for key in md:
                        metadata[key] = md[key]

                    fid.close()
                    return data, metadata

    # default aerosol value
    # assumes we are only processing Australia in which case it it should
    # be a coastal scene
    data = 0.06
    metadata = {'data_source': 'Default value used; Assumed a coastal scene'}

    fid.close()
    return data, metadata


def get_aerosol_data(acquisition, aerosol_path, aot_loader_path=None):
    """Extract the aerosol value for an acquisition.
    """

    dt = acquisition.scene_center_datetime
    geobox = acquisition.gridded_geo_box()
    ll_lon, ll_lat = geobox.ll_lonlat
    ur_lon, ur_lat = geobox.ur_lonlat

    descr = ['AATSR_PIX', 'AATSR_CMP_YEAR_MONTH', 'AATSR_CMP_MONTH']
    names = ['ATSR_LF_%Y%m.pix', 'aot_mean_%b_%Y_All_Aerosols.cmp',
             'aot_mean_%b_All_Aerosols.cmp']
    filenames = [pjoin(aerosol_path, dt.strftime(n)) for n in names]

    for filename, description in zip(filenames, descr):
        data = run_aot_loader(filename, dt, ll_lat, ll_lon, ur_lat,
                              ur_lon, aot_loader_path)
        if data:
            metadata = {'data_source': description,
                        'data_file': filename,
                        'Date used for querying': dt}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(filename)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    # default aerosol value
    # assumes we are only processing Australia in which case it it should
    # be a coastal scene
    data = 0.06
    metadata = {'data_source': 'Default value used; Assumed a coastal scene'}

    return data, metadata


def run_aot_loader(filename, dt, ll_lat, ll_lon, ur_lat, ur_lon,
                   aot_loader_path=None):
    """Load aerosol data for a specified `AATSR.
    <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_ data file.  This uses
    the executable ``aot_loader``.
    :param filename:
        The full path to the `AATSR
        <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_ file to load the
        data from.
    :type filename:
        :py:class:`str`
    :param dt:
        The date and time to extract the value for.
    :type dt:
        :py:class:`datetime.datetime`
    :param ll_lat:
        The latitude of the lower left corner of the region ('ll' for 'Lower
        Left').
    :type ll_lat:
        :py:class:`float`
    :param ll_lon:
        The longitude of the lower left corner of the region ('ll' for 'Lower
        Left').
    :type ll_lon:
        :py:class:`float`
    :param ur_lat:
        The latitude of the upper right corner of the region ('ur' for 'Upper
        Right').
    :type ur_lat:
        :py:class:`float`
    :param ur_lon:
        The longitude of the upper right corner of the region ('ur' for 'Upper
        Right').
    :type ur_lon:
        :py:class:`float`
    :param aot_loader_path:
        The directory where the executable ``aot_loader`` can be found.
    :type aot_loader_path:
        :py:class:`str`
    """
    filetype = splitext(filename)[1][1:]
    if not exists(filename):
        log.warning('Aerosol %s file (%s) not found', filetype, filename)
        return None

    if not aot_loader_path:
        aot_loader_path = abspath(pjoin(dirname(__file__), pardir, 'bin'))

    cmd = pjoin(aot_loader_path, 'aot_loader')
    if not exists(cmd):
        log.error('%s not found.', cmd)
    task = [cmd, '--' + filetype, filename, '--west', str(ll_lon),
            '--east', str(ur_lon), '--south', str(ll_lat),
            '--north', str(ur_lat), '--date', dt.strftime('%Y-%m-%d'),
            '--t', dt.strftime('%H:%M:%S'), '--print-values']
    task = ' '.join(task)
    print task
    result = subprocess.check_output(task, shell=True)

    m = re.search(r'AOT AATSR value:\s+(.+)$', result, re.MULTILINE)
    if m and m.group(1):
        return float(m.group(1).rstrip())

    log.warning('Aerosol file %s could not be parsed', filename)
    return None


def get_elevation_data(lonlat, dem_path):
    """
    Get elevation data for a scene.

    :param lon_lat:
        The latitude, longitude of the scene center.
    :type lon_lat:
        float (2-tuple)

    :dem_dir:
        The directory in which the DEM can be found.
    :type dem_dir:
        str
    """
    datafile = pjoin(dem_path, "DEM_one_deg.tif")
    data = gaip.get_pixel(datafile, lonlat) * 0.001  # scale to correct units

    metadata = {'data_source': 'Elevation',
                'data_file': datafile}

    # ancillary metadata tracking
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def get_ozone_data(ozone_path, lonlat, time):
    """
    Get ozone data for a scene. `lonlat` should be the (x,y) for the centre
    the scene.
    """
    filename = time.strftime('%b').lower() + '.tif'
    datafile = pjoin(ozone_path, filename)
    data = gaip.get_pixel(datafile, lonlat)

    metadata = {'data_source': 'Ozone',
                'data_file': datafile,
                'Date used for querying': time}

    # ancillary metadata tracking
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def get_water_vapour(acquisition, vapour_path, scale_factor=0.1):
    """
    Retrieve the water vapour value for an `acquisition` and the
    path for the water vapour ancillary data.
    """
    dt = acquisition.scene_center_datetime
    geobox = acquisition.gridded_geo_box()

    year = dt.strftime('%Y')
    filename = "pr_wtr.eatm.{year}.tif".format(year=year)
    datafile = os.path.join(vapour_path, filename)

    # calculate the water vapour band number based on the datetime

    doy = dt.timetuple().tm_yday
    hour = dt.timetuple().tm_hour
    band = (int(doy) - 1) * 4 + int((hour + 3) / 6)

    # Check for boundary condition: 1 Jan, 0-3 hours
    if band == 0 and doy == 1:
        band = 1

    # Get the number of bands
    with rasterio.open(datafile) as src:
        n_bands = src.count

    # Enable NBAR Near Real Time (NRT) processing
    if band > (n_bands + 1):
        rasterdoy = (((n_bands) - (int((hour + 3) / 6))) / 4) + 1
        if (doy - rasterdoy) < 7:
            band = (int(rasterdoy) - 1) * 4 + int((hour + 3) / 6)

    try:
        data = gaip.get_pixel(datafile, geobox.centre_lonlat, band=band)
    except IndexError:
        msg = "Invalid water vapour band number: {band}".format(band=band)
        raise IndexError(msg)

    data = data * scale_factor

    metadata = {'data_source': 'Water Vapour',
                'data_file': datafile,
                'Date used for querying': dt}

    # ancillary metadata tracking
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def ecwmf_elevation(datafile, lonlat):
    """
    Retrieve a pixel from the ECWMF invariant geo-potential
    dataset.
    Converts to Geo-Potential height in KM.
    """
    data = gaip.get_pixel(datafile, lonlat) / 9.80665 / 1000.0

    metadata = {'data_source': 'ECWMF Invariant Geo-Potential',
                'data_file': datafile}

    # ancillary metadata tracking
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def ecwmf_temperature_2metre(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF 2 metre Temperature
    collection.
    """
    product = 'temperature-2metre'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f))[0].split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= end:
            data = gaip.get_pixel(f, lonlat, time.day)

            metadata = {'data_source': 'ECWMF 2 metre Temperature',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        raise IOError('No ECWMF 2 metre Temperature ancillary data found.')


def ecwmf_dewpoint_temperature(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF 2 metre Dewpoint
    Temperature collection.
    """
    product = 'dewpoint-temperature'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f))[0].split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= end:
            data = gaip.get_pixel(f, lonlat, time.day)

            metadata = {'data_source': 'ECWMF 2 metre Dewpoint Temperature ',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        msg = 'No ECWMF 2 metre Dewpoint Temperature ancillary data found.'
        raise IOError(msg)


def ecwmf_surface_pressure(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Surface Pressure
    collection.
    Scales the result by 100 before returning.
    """
    product = 'surface-pressure'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f))[0].split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= end:
            data = gaip.get_pixel(f, lonlat, time.day) / 100.0

            metadata = {'data_source': 'ECWMF Surface Pressure',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        raise IOError('No Surface Pressure ancillary data found.')


def ecwmf_water_vapour(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Total Column Water Vapour
    collection.
    """
    product = 'water-vapour'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f)).split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= time:
            data = gaip.get_pixel(f, lonlat, time.day)

            metadata = {'data_source': 'ECWMF Total Column Water Vapour',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        msg = 'No ECWMF Total Column Water Vapour ancillary data found.'
        raise IOError(msg)


def ecwmf_temperature(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Temperature collection
    across 37 height pressure levels, for a given longitude,
    latitude and time.

    Reverses the order of elements
    (1000 -> 1 mb, rather than 1 -> 1000 mb) before returning.
    """
    product = 'temperature'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f))[0].split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= end:
            bands = get_4d_idx(time.day)
            data = gaip.get_pixel(f, lonlat, bands)[::-1]

            metadata = {'data_source': 'ECWMF Temperature',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_meatadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'Temperature', data)

            return df, metadata

    if data is None:
        raise IOError('No Temperature ancillary data found.')


def ecwmf_geo_potential(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Geo-Potential collection
    across 37 height pressure levels, for a given longitude,
    latitude and time.

    Converts to geo-potential height in KM, and reverses the order of
    the elements (1000 -> 1 mb, rather than 1 -> 1000 mb) before
    returning.
    """
    product = 'geo-potential'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f))[0].split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= end:
            bands = get_4d_idx(time.day)
            data = gaip.get_pixel(f, lonlat, bands)[::-1]
            scaled_data = data / 9.80665 / 1000.0

            metadata = {'data_source': 'ECWMF Geo-Potential',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # ancillary metadata tracking
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_meatadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'GeoPotential_Height', data)
            df.insert(1, 'GeoPotential', scaled_data)

            return df, md

    if data is None:
        raise IOError('No Geo-Potential ancillary data found.')


def ecwmf_relative_humidity(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Relative Humidity collection
    across 37 height pressure levels, for a given longitude,
    latitude and time.

    Reverses the order of elements
    (1000 -> 1 mb, rather than 1 -> 1000 mb) before returning.
    """
    product = 'relative-humidity'
    files = glob.glob(pjoin(input_path, '{}_*.grib'.format(product)))
    data = None
    for f in files:
        start, end = splitext(basename(f))[0].split('_')[1:]
        start = datetime.datetime.strptime(start, '%Y-%m-%d')
        end = datetime.datetime.strptime(end, '%Y-%m-%d')
        if start <= time <= end:
            bands = get_4d_idx(time.day)
            data = gaip.get_pixel(f, lonlat, bands)[::-1]

            metadata = {'data_source': 'ECWMF Relative Humidity',
                        'data_file': f,
                        'Date used for querying ECWMF': time}

            # file level metadata
            md = gaip.extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_meatadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'Relative_Humidity', data)

            return df, metadata

    if data is None:
        raise IOError('No Relative Humidity ancillary data found.')
