"""Ancillary datasets."""

import logging
import subprocess
import re
import os
from os.path import join as pjoin, splitext, exists, abspath, dirname, pardir
from posixpath import join as ppjoin
import datetime

import h5py
import pandas
from geopandas import GeoSeries
import numpy
import rasterio
from shapely.geometry import Point
from shapely.geometry import Polygon
import gaip

log = logging.getLogger()


def _collect_ancillary_data(acquisition, aerosol_path, water_vapour_path,
                            ozone_path, dem_path, brdf_path,
                            brdf_premodis_path, out_fname):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    def _write_scalar(fid, dataset_path, data, factor):
        scattering_names = {'iso': 'isometric',
                            'vol': 'volumetric',
                            'geo': 'geometric'}

        description = ("Bidirectional Reflectance Distribution Function "
                       "for the {} scattering fraction."
        data['Description'] = description.format(scattering_names[factor])

        value = data.pop('value')
        dset = fid.create_dataset(dataset_path, data=value)
        for key in data:
            if isinstance(data[key], datetime.datetime):
                data[key] = data[key].isoformat()

        attach_attrs(dset, data)
        fid.flush()

        return


    with h5py.File(out_fname, 'w') as fid:
        geobox = acquisition.gridded_geo_box()

        _write_scalar(fid, 'aerosol',
                      get_aerosol_data(acquisition, aerosol_path))
        _write_scalar(fid, 'water-vapour',
                      get_water_vapour(acquisition, water_vapour_path))
        _write_scalar(fid, 'ozone',
                      get_ozone_data(ozone_path, geobox.centre_lonlat,
                                     acquisition.scene_center_datetime))
        _write_scalar(fid, 'elevation',
                      get_elevation_data(geobox.centre_lonlat, dem_path))

        # brdf
        group = fid.create_group('brdf-image-datasets')
        data = gaip.get_brdf_data(acquisition, brdf_path, brdf_premodis_path,
                                  group)
        dname_format = "BRDF-Band-{band}-{factor}"
        for key in data:
            band, factor = key
            _write_scalar(fid, dname_format.format(band=band, factor=factor),
                          data[key], factor)

    return


def aggregate_ancillary(ancillary_fnames, out_fname):
    """
    If the acquisition is part of a `tiled` scene such as Sentinel-2a,
    then we need to average the point measurements gathereed from
    all tiles.
    """
    # initialise the mean result
    ozone = vapour = aerosol = elevation = 0.0

    n_tiles = len(ancillary_fnames)

    with h5py.File(out_fname, 'w') as fid1:
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
        attach_attrs(dset, attrs)

        dset = fid1.create_dataset('water-vapour', data=vapour)
        attrs['Description'] = description.format('Water Vapour')
        attach_attrs(dset, attrs)

        dset = fid1.create_dataset('aerosol', data=aerosol)
        attrs['Description'] = description.format('Aerosol')
        attach_attrs(dset, attrs)

        dset = fid1.create_dataset('elevation', data=elevation)
        attrs['Description'] = description.format('Elevation')
        attach_attrs(dset, attrs)

    return


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

    delta_tolerance = datetime.timedelta(days=0.5)

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
                    md = gaip.extract_ancillary_metadata(aerosol_fname)
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
            md = gaip.extract_ancillary_metadata(filename)
            for key in md:
                res[key] = md[key]

            return res

    # default aerosol value
    # assumes we are only processing Australia in which case it it should
    # be a coastal scene
    res = {'value': 0.06}

    return res


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
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        res[key] = md[key]

    return res


def get_ozone_data(ozone_path, lonlat, time):
    """
    Get ozone data for a scene. `lonlat` should be the (x,y) for the centre
    the scene.
    """
    filename = time.strftime('%b').lower() + '.tif'
    datafile = pjoin(ozone_path, filename)
    value = gaip.get_pixel(datafile, lonlat)

    res = {'data_source': 'Ozone',
           'data_file': datafile,
           'value': value}

    # ancillary metadata tracking
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        res[key] = md[key]

    return res


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
    md = gaip.extract_ancillary_metadata(datafile)
    for key in md:
        water_vapour_data[key] = md[key]

    return water_vapour_data
