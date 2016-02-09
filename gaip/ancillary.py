"""Ancillary datasets."""

import logging
import subprocess
import re
import os
from os.path import join as pjoin, splitext, exists, abspath, dirname, pardir
from posixpath import join as ppjoin
import pwd
from datetime import datetime as dtime
from datetime import timedelta

import pandas
from geopandas import GeoSeries
import numpy
import rasterio
from shapely.geometry import Point
from shapely.geometry import Polygon
import gaip

log = logging.getLogger()


def extract_ancillary_metadata(fname):
    """
    Extracts the change (last metadata change), modified,
    accessed, and owner user id.

    :param fname:
        A string containing the full file pathname to a file
        on disk.

    :return:
        A `dictionary` with keys `change`, `modified`, `accessed`,
        and `user`.
    """
    res = {}
    fstat = os.stat(fname)
    res['change'] = dtime.utcfromtimestamp(fstat.st_ctime)
    res['modified'] = dtime.utcfromtimestamp(fstat.st_mtime)
    res['accessed'] = dtime.utcfromtimestamp(fstat.st_atime)
    res['user'] = pwd.getpwuid(fstat.st_uid).pw_gecos
    return res


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

    store = pandas.HDFStore(aerosol_fname, 'r')

    delta_tolerance = timedelta(days=0.5)

    value = None
    for pathname, description in zip(pathnames, descr):
        if pathname in store.keys():
            df = store[pathname]
            node = store.get_node(pathname)
            aerosol_poly = Polygon(node._v_attrs.extents)
            if aerosol_poly.intersects(roi_poly):
                if description == 'AATSR_PIX':
                    abs_diff = (df['timestamp'] - dt).abs()
                    df = df[abs_diff < delta_tolerance]
                intersection = aerosol_poly.intersection(roi_poly)
                pts = GeoSeries([Point(x, y) for x, y in
                                 zip(df['lon'], df['lat'])])
                idx = pts.intersects(intersection)
                value = df[idx]['aerosol'].mean()
                if numpy.isfinite(value):
                    res = {'data_source': description,
                           'data_file': pathname,
                           'value': value}

                    # ancillary metadata tracking
                    md = extract_ancillary_metadata(aerosol_fname)
                    for key in md:
                        res[key] = md[key]

                    store.close()
                    return res
    store.close()

    raise IOError('No aerosol ancillary data found.')


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
        value = run_aot_loader(filename, dt, ll_lat, ll_lon, ur_lat,
                               ur_lon, aot_loader_path)
        if value:
            res = {'data_source': description,
                   'data_file': filename,
                   'value': value}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(filename)
            for key in md:
                res[key] = md[key]

            return res

    raise IOError('No aerosol ancillary data found.')


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
            '--t', dt.strftime('%H:%M:%S')]
    task = ' '.join(task)
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
    value = gaip.get_pixel(datafile, lonlat) * 0.001  # scale to correct units

    res = {'data_source': 'Elevation',
           'data_file': datafile,
           'value': value}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        res[key] = md[key]

    return res


def get_ozone_data(ozone_path, lonlat, datetime):
    """
    Get ozone data for a scene. `lonlat` should be the (x,y) for the centre
    the scene.
    """
    filename = datetime.strftime('%b').lower() + '.tif'
    datafile = pjoin(ozone_path, filename)
    value = gaip.get_pixel(datafile, lonlat)

    res = {'data_source': 'Ozone',
           'data_file': datafile,
           'value': value}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        res[key] = md[key]

    return res


def get_solar_irrad(acquisitions, solar_path):
    """
    Extract solar irradiance values from the specified file. One for each band

    """
    acqs = [a for a in acquisitions if a.band_type == gaip.REF]
    bands = [a.band_num for a in acqs]

    with open(pjoin(solar_path, acqs[0].solar_irrad_file), 'r') as infile:
        header = infile.readline()
        if 'band solar irradiance' not in header:
            raise IOError('Cannot load solar irradiance file')

        irrads = {}
        for line in infile.readlines():
            band, value = line.strip().split()
            band, value = int(band), float(value)  # parse
            if band in bands:
                irrads[band] = value

        return irrads


def get_solar_dist(acquisition, sundist_path):
    """
    Extract Earth-Sun distance for this day of the year (varies during orbit).

    """
    doy = acquisition.scene_center_datetime.timetuple().tm_yday

    with open(sundist_path, 'r') as infile:
        for line in infile.readlines():
            index, dist = line.strip().split()
            index = int(index)
            if index == doy:
                res = {'data_source': 'Solar Distance',
                       'data_file': sundist_path,
                       'value': float(dist)}

                # ancillary metadata tracking
                md = extract_ancillary_metadata(sundist_path)
                for key in md:
                    res[key] = md[key]

                return res

    raise IOError('Cannot load Earth-Sun distance')


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
        value = gaip.get_pixel(datafile, geobox.centre_lonlat, band=band)
    except IndexError:
        msg = "Invalid water vapour band number: {band}".format(band=band)
        raise IndexError(msg)

    value = value * scale_factor

    water_vapour_data = {
        'data_source': 'Water Vapour',
        'data_file': datafile,
        'value': value
    }

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        water_vapour_data[key] = md[key]

    return water_vapour_data
