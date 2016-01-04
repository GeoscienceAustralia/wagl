"""Ancillary datasets."""

import logging
import rasterio
import subprocess
import gaip
import re
import os
from datetime import timedelta

from os.path import join as pjoin
from posixpath import join as ppjoin
import pandas
from geopandas import GeoSeries
from shapely.geometry import Point
from shapely.geometry import Polygon
import numpy

log = logging.getLogger()


def get_aerosol_data(acquisition, aerosol_path, aot_loader_path=None):
    """Extract the aerosol value for an acquisition.
    """

    dt = acquisition.scene_center_datetime
    geobox = acquisition.gridded_geo_box()
    roi_poly = Polygon([geobox.ul, geobox.ur, geobox.lr, geobox.ll])

    descr = ['AATSR_PIX', 'AATSR_CMP_YEAR_MONTH', 'AATSR_CMP_MONTH']
    names = ['ATSR_LF_%Y%m', 'aot_mean_%b_%Y_All_Aerosols',
             'aot_mean_%b_All_Aerosols']
    exts = ['/pix', '/cmp', '/cmp']
    pathnames = [ppjoin(ext, dt.strftime(n)) for ext, n in zip(exts, names)]

    store = pandas.HDFStore(aerosol_path, 'r')

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
                    return {'data_source': description,
                            'data_file': pathname,
                            'value': value}

    store.close()


    # TODO write a shapefile of the intersected polygon
    # and optionally the points used for determining the mean value???
    # from fiona.crs import from_epsg
    # from shapely.geometry import mapping, Polygon, Point
    # schema = {'geometry': 'Polygon', 'properties': {'ancillary': 'str'}}
    # crs = from_epsg(4326)
    # with fiona.open('test-fiona-write.shp', 'w', 'ESRI Shapefile',
    #                 schema, crs=crs) as src:
    #     src.write({'geometry': mapping(poly),
    #                'properties': {'ancillary': 'aerosol'}})


    raise IOError('No aerosol ancillary data found.')


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
    return {'data_source': 'Elevation',
            'data_file': datafile,
            'value': value}


def get_ozone_data(ozone_path, lonlat, datetime):
    """
    Get ozone data for a scene. `lonlat` should be the (x,y) for the centre
    the scene.
    """
    filename = datetime.strftime('%b').lower() + '.tif'
    datafile = pjoin(ozone_path, filename)
    value = gaip.get_pixel(datafile, lonlat)
    return {'data_source': 'Ozone',
            'data_file': datafile,
            'value': value}


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
                return {'data_source': 'Solar Distance',
                        'data_file': sundist_path,
                        'distance': float(dist)}

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

    return water_vapour_data
