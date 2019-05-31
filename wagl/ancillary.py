#!/usr/bin/env python

"""
Ancillary dataset retrieval and storage
"""

from __future__ import absolute_import, print_function
from os.path import join as pjoin, basename, splitext
import datetime
import glob
from urllib.parse import urlparse

from posixpath import join as ppjoin
import numpy
import h5py
import pandas
from geopandas import GeoSeries
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely import wkt
from wagl.brdf import get_brdf_data
from wagl.data import get_pixel
from wagl.hdf5 import attach_attributes, write_scalar, write_dataframe
from wagl.hdf5 import read_h5_table, H5CompressionFilter
from wagl.hdf5 import attach_table_attributes
from wagl.metadata import extract_ancillary_metadata, read_metadata_tags, current_h5_metadata
from wagl.constants import DatasetName, POINT_FMT, GroupName, BandType
from wagl.constants import AerosolTier, WaterVapourTier, OzoneTier
from wagl.satellite_solar_angles import create_vertices


ECWMF_LEVELS = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200,
                225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
                775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]


class AncillaryError(Exception):

    """
    Specific error handle for ancillary retrieval
    """


def get_4d_idx(day):
    """
    A small utility function for indexing into a 4D dataset
    represented as a 3D dataset.
    [month, level, y, x], where level contains 37 levels, and day
    contains 28, 29, 30 or 31 days.
    """
    start = 1 + 37 * (day - 1)
    stop = start + 37
    return list(range(start, stop, 1))


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

    rh = 100 * ((112.0 - 0.1 * surf_t + dew_t) / (112.0 + 0.9 * surf_t)) ** 8

    return rh


def _collect_ancillary(container, satellite_solar_fname, nbar_paths,
                       sbt_path=None, invariant_fname=None, vertices=(3, 3),
                       out_fname=None, compression=H5CompressionFilter.LZF,
                       filter_opts=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_fname, 'r') as fid,\
        h5py.File(out_fname, 'w') as out_fid:

        sat_sol_grp = fid[GroupName.SAT_SOL_GROUP.value]
        collect_ancillary(container, sat_sol_grp, nbar_paths, sbt_path,
                          invariant_fname, vertices, out_fid, compression)


def collect_ancillary(container, satellite_solar_group, nbar_paths,
                      sbt_path=None, invariant_fname=None, vertices=(3, 3),
                      out_group=None, compression=H5CompressionFilter.LZF,
                      filter_opts=None):
    """
    Collects the ancillary required for NBAR and optionally SBT.
    This could be better handled if using the `opendatacube` project
    to handle ancillary retrieval, rather than directory passing,
    and filename grepping.

    :param container:
        An instance of an `AcquisitionsContainer` object.
        The container should consist of a single Granule or None,
        only. Use `AcquisitionsContainer.get_granule` method prior to
        calling this function.

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the solar zenith and
        solar azimuth datasets specified by the pathnames given by:

        * DatasetName.BOXLINE

    :param nbar_paths:
        A `dict` containing the ancillary pathnames required for
        retrieving the NBAR ancillary data. Required keys:

        * aerosol_data
        * water_vapour_data
        * ozone_path
        * dem_path
        * brdf_dict

    :param sbt_path:
        A `str` containing the base directory pointing to the
        ancillary products required for the SBT workflow.

    :param invariant_fname:
        A `str` containing the file path name to the invariant
        geopotential image file.

    :param vertices:
        An integer 2-tuple indicating the number of rows and columns
        of sample-locations ("coordinator") to produce.
        The vertex columns should be an odd number.
        Default is (3, 3).

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

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
    # Initialise the output files
    if out_group is None:
        fid = h5py.File('ancillary.h5', driver='core', backing_store=False)
    else:
        fid = out_group

    if filter_opts is None:
        filter_opts = {}

    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    group = fid.create_group(GroupName.ANCILLARY_GROUP.value)

    acquisition = container.get_highest_resolution()[0][0]

    boxline_dataset = satellite_solar_group[DatasetName.BOXLINE.value][:]
    coordinator = create_vertices(acquisition, boxline_dataset, vertices)
    lonlats = zip(coordinator['longitude'], coordinator['latitude'])

    desc = ("Contains the row and column array coordinates used for the "
            "atmospheric calculations.")
    attrs = {'description': desc, 'array_coordinate_offset': 0}
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    dset_name = DatasetName.COORDINATOR.value
    coord_dset = group.create_dataset(dset_name, data=coordinator, **kwargs)
    attach_table_attributes(coord_dset, title='Coordinator', attrs=attrs)

    if sbt_path:
        collect_sbt_ancillary(acquisition, lonlats, sbt_path, invariant_fname,
                              out_group=group, compression=compression,
                              filter_opts=filter_opts)

    collect_nbar_ancillary(container, out_group=group,
                           compression=compression, filter_opts=filter_opts,
                           **nbar_paths)

    if out_group is None:
        return fid


def collect_sbt_ancillary(acquisition, lonlats, ancillary_path,
                          invariant_fname=None, out_group=None,
                          compression=H5CompressionFilter.LZF,
                          filter_opts=None):
    """
    Collects the ancillary data required for surface brightness
    temperature.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param lonlats:
        A `list` of tuples containing (longitude, latitude) coordinates.

    :param ancillary_path:
        A `str` containing the directory pathname to the ECMWF
        ancillary data.

    :param invariant_fname:
        A `str` containing the file pathname to the invariant geopotential
        data.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

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
    # Initialise the output files
    if out_group is None:
        fid = h5py.File('sbt-ancillary.h5', driver='core', backing_store=False)
    else:
        fid = out_group

    fid.attrs['sbt-ancillary'] = True

    dt = acquisition.acquisition_datetime

    description = ('Combined Surface and Pressure Layer data retrieved from '
                   'the ECWMF catalogue.')
    attrs = {'description': description,
             'Date used for querying ECWMF': dt}

    for i, lonlat in enumerate(lonlats):
        pnt = POINT_FMT.format(p=i)
        # get data located at the surface
        dew = ecwmf_dewpoint_temperature(ancillary_path, lonlat, dt)
        t2m = ecwmf_temperature_2metre(ancillary_path, lonlat, dt)
        sfc_prs = ecwmf_surface_pressure(ancillary_path, lonlat, dt)
        sfc_hgt = ecwmf_elevation(invariant_fname, lonlat)
        sfc_rh = relative_humdity(t2m[0], dew[0])

        # output the scalar data along with the attrs
        dname = ppjoin(pnt, DatasetName.DEWPOINT_TEMPERATURE.value)
        write_scalar(dew[0], dname, fid, dew[1])

        dname = ppjoin(pnt, DatasetName.TEMPERATURE_2M.value)
        write_scalar(t2m[0], dname, fid, t2m[1])

        dname = ppjoin(pnt, DatasetName.SURFACE_PRESSURE.value)
        write_scalar(sfc_prs[0], dname, fid, sfc_prs[1])

        dname = ppjoin(pnt, DatasetName.SURFACE_GEOPOTENTIAL.value)
        write_scalar(sfc_hgt[0], dname, fid, sfc_hgt[1])

        dname = ppjoin(pnt, DatasetName.SURFACE_RELATIVE_HUMIDITY.value)
        attrs = {'description': 'Relative Humidity calculated at the surface'}
        write_scalar(sfc_rh, dname, fid, attrs)

        # get the data from each of the pressure levels (1 -> 1000 ISBL)
        gph = ecwmf_geo_potential(ancillary_path, lonlat, dt)
        tmp = ecwmf_temperature(ancillary_path, lonlat, dt)
        rh = ecwmf_relative_humidity(ancillary_path, lonlat, dt)

        dname = ppjoin(pnt, DatasetName.GEOPOTENTIAL.value)
        write_dataframe(gph[0], dname, fid, compression, attrs=gph[1],
                        filter_opts=filter_opts)

        dname = ppjoin(pnt, DatasetName.TEMPERATURE.value)
        write_dataframe(tmp[0], dname, fid, compression, attrs=tmp[1],
                        filter_opts=filter_opts)

        dname = ppjoin(pnt, DatasetName.RELATIVE_HUMIDITY.value)
        write_dataframe(rh[0], dname, fid, compression, attrs=rh[1],
                        filter_opts=filter_opts)

        # combine the surface and higher pressure layers into a single array
        cols = ['GeoPotential_Height', 'Pressure', 'Temperature',
                'Relative_Humidity']
        layers = pandas.DataFrame(columns=cols, index=range(rh[0].shape[0]),
                                  dtype='float64')

        layers['GeoPotential_Height'] = gph[0]['GeoPotential_Height'].values
        layers['Pressure'] = ECWMF_LEVELS[::-1]
        layers['Temperature'] = tmp[0]['Temperature'].values
        layers['Relative_Humidity'] = rh[0]['Relative_Humidity'].values

        # define the surface level
        df = pandas.DataFrame({'GeoPotential_Height': sfc_hgt[0],
                               'Pressure': sfc_prs[0],
                               'Temperature': kelvin_2_celcius(t2m[0]),
                               'Relative_Humidity': sfc_rh}, index=[0])

        # MODTRAN requires the height to be ascending
        # and the pressure to be descending
        wh = ((layers['GeoPotential_Height'] > sfc_hgt[0]) &
              (layers['Pressure'] < sfc_prs[0].round()))
        df = df.append(layers[wh])
        df.reset_index(drop=True, inplace=True)

        dname = ppjoin(pnt, DatasetName.ATMOSPHERIC_PROFILE.value)
        write_dataframe(df, dname, fid, compression, attrs=attrs,
                        filter_opts=filter_opts)

        fid[pnt].attrs['lonlat'] = lonlat

    if out_group is None:
        return fid


def collect_nbar_ancillary(container, aerosol_dict=None,
                           water_vapour_dict=None, ozone_path=None,
                           dem_path=None, brdf_dict=None,
                           out_group=None,
                           compression=H5CompressionFilter.LZF,
                           filter_opts=None):
    """
    Collects the ancillary information required to create NBAR.

    :param container:
        An instance of an `AcquisitionsContainer` object.

    :param aerosol_dict:
        A `dict` defined as either of the following:

        * {'user': <value>}
        * {'pathname': <value>}

    :param water_vapour_dict:
        A `dict` defined as either of the following:

        * {'user': <value>}
        * {'pathname': <value>}

    :param ozone_path:
        A `str` containing the full file pathname to the directory
        containing the ozone data.

    :param dem_path:
        A `str` containing the full file pathname to the directory
        containing the digital elevation model data.

    :param brdf_dict:
        A `dict` defined as either of the following:

        * {'user': {<band-alias>: {'alpha_1': <value>, 'alpha_2': <value>}, ...}}
        * {'brdf_path': <path-to-BRDF>, 'brdf_premodis_path': <path-to-average-BRDF>}

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

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

    :notes:
        The keywords compression and filter_opts aren't used as we no
        longer save the BRDF imagery. However, we may need to store
        tables in future, therefore they can remain until we know
        for sure they'll never be used.
    """
    # Initialise the output files
    if out_group is None:
        fid = h5py.File('nbar-ancillary.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    acquisition = container.get_highest_resolution()[0][0]
    dt = acquisition.acquisition_datetime
    geobox = acquisition.gridded_geo_box()

    aerosol = get_aerosol_data(acquisition, aerosol_dict)
    write_scalar(aerosol[0], DatasetName.AEROSOL.value, fid, aerosol[1])

    wv = get_water_vapour(acquisition, water_vapour_dict)
    write_scalar(wv[0], DatasetName.WATER_VAPOUR.value, fid, wv[1])

    ozone = get_ozone_data(ozone_path, geobox.centre_lonlat, dt)
    write_scalar(ozone[0], DatasetName.OZONE.value, fid, ozone[1])

    elev = get_elevation_data(geobox.centre_lonlat, dem_path)
    write_scalar(elev[0], DatasetName.ELEVATION.value, fid, elev[1])

    # brdf
    dname_format = DatasetName.BRDF_FMT.value
    for group in container.groups:
        for acq in container.get_acquisitions(group=group):
            if acq.band_type is not BandType.REFLECTIVE:
                continue
            data = get_brdf_data(acq, brdf_dict,
                                 compression)

            # output
            for param in data:
                dname = dname_format.format(parameter=param.value,
                                            band_name=acq.band_name)
                brdf_value = data[param].pop('value')
                write_scalar(brdf_value, dname, fid, data[param])

    if out_group is None:
        return fid


def _aggregate_ancillary(ancillary_fnames, write_access):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    # a horrible mechanism to ensure an unneeded logic contiunes; sigh....
    fnames = ancillary_fnames.copy()
    fnames.pop(fnames.index(write_access))

    # get file ids
    fids = [h5py.File(fname, 'r') for fname in fnames]
    fids.append(h5py.File(write_access, 'a'))
    aggregate_ancillary(fids)

    # close
    for fid in fids:
        fid.close()


def aggregate_ancillary(granule_groups):
    """
    If the acquisition is part of a `tiled` scene such as Sentinel-2a,
    then we need to average the point measurements gathered from
    all granules.
    """
    # initialise the mean result
    ozone = vapour = aerosol = elevation = 0.0

    # number of granules in the scene
    n_tiles = len(granule_groups)

    for granule in granule_groups:
        group = granule[GroupName.ANCILLARY_GROUP.value]

        ozone += group[DatasetName.OZONE.value][()]
        vapour += group[DatasetName.WATER_VAPOUR.value][()]
        aerosol += group[DatasetName.AEROSOL.value][()]
        elevation += group[DatasetName.ELEVATION.value][()]

    # average
    ozone /= n_tiles
    vapour /= n_tiles
    aerosol /= n_tiles
    elevation /= n_tiles

    description = ("The {} value is an average from all the {} values "
                   "retreived for each Granule.")
    attrs = {'data_source': 'granule_average'}

    # output each average value back into the same granule ancillary group
    group_name = ppjoin(GroupName.ANCILLARY_GROUP.value,
                        GroupName.ANCILLARY_AVG_GROUP.value)
    for granule in granule_groups:
        # for the multifile workflow, we only want to write to one granule
        try:
            group = granule.create_group(group_name)
        except ValueError:
            continue

        dset = group.create_dataset(DatasetName.OZONE.value, data=ozone)
        attrs['description'] = description.format(*(2*['Ozone']))
        attach_attributes(dset, attrs)

        dset = group.create_dataset(DatasetName.WATER_VAPOUR.value, data=vapour)
        attrs['description'] = description.format(*(2*['Water Vapour']))
        attach_attributes(dset, attrs)

        dset = group.create_dataset(DatasetName.AEROSOL.value, data=aerosol)
        attrs['description'] = description.format(*(2*['Aerosol']))
        attach_attributes(dset, attrs)

        dset = group.create_dataset(DatasetName.ELEVATION.value, data=elevation)
        attrs['description'] = description.format(*(2*['Elevation']))
        attach_attributes(dset, attrs)


def get_aerosol_data(acquisition, aerosol_dict):
    """
    Extract the aerosol value for an acquisition.
    The version 2 retrieves the data from a HDF5 file, and provides
    more control over how the data is selected geo-metrically.
    Better control over timedeltas.
    """

    dt = acquisition.acquisition_datetime
    geobox = acquisition.gridded_geo_box()
    roi_poly = Polygon([geobox.ul_lonlat, geobox.ur_lonlat,
                        geobox.lr_lonlat, geobox.ll_lonlat])

    descr = ['AATSR_PIX', 'AATSR_CMP_YEAR_MONTH', 'AATSR_CMP_MONTH']
    names = ['ATSR_LF_%Y%m', 'aot_mean_%b_%Y_All_Aerosols',
             'aot_mean_%b_All_Aerosols']
    exts = ['/pix', '/cmp', '/cmp']
    pathnames = [ppjoin(ext, dt.strftime(n)) for ext, n in zip(exts, names)]

    # temporary until we sort out a better default mechanism
    # how do we want to support default values, whilst still support provenance
    if 'user' in aerosol_dict:
        tier = AerosolTier.USER
        metadata = {
            'id': [None],
            'tier': tier
        }

        return aerosol_dict['user'], metadata

    aerosol_fname = aerosol_dict['pathname']

    fid = h5py.File(aerosol_fname, 'r')

    delta_tolerance = datetime.timedelta(days=0.5)

    data = None
    for pathname, description in zip(pathnames, descr):
        tier = AerosolTier['description']
        if pathname in fid:
            df = read_h5_table(fid, pathname)
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
                idx = pts.within(intersection)
                data = df[idx]['aerosol'].mean()

                if numpy.isfinite(data):
                    # ancillary metadata tracking
                    md = current_h5_metadata(fid, dataset_path=pathname)
                    metadata = {
                        'id': [md['id']],
                        'tier': tier
                    }

                    fid.close()
                    return data, metadata

    # default aerosol value
    data = 0.06
    metadata = {
        'id': [None],
        'tier': AerosolTier.FALLBACK_DEFAULT
    }

    fid.close()
    return data, metadata


def get_elevation_data(lonlat, pathname):
    """
    Get elevation data for a scene.

    :param lon_lat:
        The latitude, longitude of the scene center.
    :type lon_lat:
        float (2-tuple)

    :pathname:
        The pathname of the DEM with a ':' to seperate the
        dataset name.
    :type dem_dir:
        str
    """
    fname, dname = pathname.split(':')

    try:
        data, md_uuid = get_pixel(fname, dname, lonlat) * 0.001  # scale to correct units
    except ValueError:
        raise AncillaryError("No Elevation data")

    return data, md_uuid


def get_ozone_data(ozone_fname, lonlat, time):
    """
    Get ozone data for a scene. `lonlat` should be the (x,y) for the centre
    the scene.
    """
    dname = time.strftime('%b').lower()

    try:
        data, md_uuid = get_pixel(ozone_fname, dname, lonlat)
    except ValueError:
        raise AncillaryError("No Ozone data")

    metadata = {
        'id': [md_uuid],
        'tier': OzoneTier.DEFINITIVE
    }

    return data, metadata


def get_water_vapour(acquisition, water_vapour_dict, scale_factor=0.1,
                     tolerance=1):
    """
    Retrieve the water vapour value for an `acquisition` and the
    path for the water vapour ancillary data.
    """
    dt = acquisition.acquisition_datetime
    geobox = acquisition.gridded_geo_box()

    year = dt.strftime('%Y')
    hour = dt.timetuple().tm_hour
    filename = "pr_wtr.eatm.{year}.h5".format(year=year)

    if 'user' in water_vapour_dict:
        metadata = {
            'id': [None],
            'tier': WaterVapourTier.USER
        }
        return water_vapour_dict['user'], metadata

    water_vapour_path = water_vapour_dict['pathname']

    datafile = pjoin(water_vapour_path, filename)

    with h5py.File(datafile, 'r') as fid:
        index = read_h5_table(fid, 'INDEX')

    # set the tolerance in days to search back in time
    day_zero = dt - dt
    max_tolerance = day_zero - datetime.timedelta(days=tolerance)

    # only look for observations that have occured in the past
    time_delta = index.timestamp - dt
    result = time_delta[(time_delta < day_zero) & (time_delta > max_tolerance)]
    if result.shape[0] == 0:
        if 'fallback_dataset' not in water_vapour_dict:
            raise AncillaryError("No actual or fallback water vapour data.")

        tier = WaterVapourTier.FALLBACK_DATASET
        month = dt.strftime('%B-%d').upper()

        # closest previous observation
        # i.e. observations are at 0000, 0600, 1200, 1800
        # and an acquisition hour of 1700 will use the 1200 observation
        observations = numpy.array([0, 6, 12, 18])
        hr = observations[numpy.argmin(numpy.abs(hour - observations))]
        dataset_name = 'AVERAGE/{}/{:02d}00'.format(month, hr)
        datafile = water_vapour_dict['fallback_data']
    else:
        tier = WaterVapourTier.DEFINITIVE
        # get the index of the closest water vapour observation
        # which would be the maximum timedelta
        # as we're only dealing with negative timedelta's here
        idx = result.argmax()
        record = index.iloc[idx]
        dataset_name = record.band_name

    try:
        data, md_uuid = get_pixel(datafile, dataset_name, geobox.centre_lonlat)
    except ValueError:
        # h5py raises a ValueError not an IndexError for out of bounds
        raise AncillaryError("No Water Vapour data")

    # the metadata from the original file says (Kg/m^2)
    # so multiply by 0.1 to get (g/cm^2)
    data = data * scale_factor
    metadata = {
        'id': [md_uuid],
        'tier': tier
    }

    return data, metadata


# TODO; have swfo convert the files to HDF5
def ecwmf_elevation(datafile, lonlat):
    """
    Retrieve a pixel from the ECWMF invariant geo-potential
    dataset.
    Converts to Geo-Potential height in KM.
    2 metres is added to the result before returning.
    """
    try:
        data = get_pixel(datafile, lonlat) / 9.80665 / 1000.0 + 0.002
    except IndexError:
        raise AncillaryError("No Invariant Geo-Potential data")

    url = urlparse(datafile, scheme='file').geturl()

    metadata = {'data_source': 'ECWMF Invariant Geo-Potential',
                'url': url}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


# TODO; have swfo convert the files to HDF5
def ecwmf_temperature_2metre(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF 2 metre Temperature
    collection.
    """
    product = DatasetName.TEMPERATURE_2M.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat)

            metadata = {'data_source': 'ECWMF 2 metre Temperature',
                        'url': url,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        raise AncillaryError("No ECWMF 2 metre Temperature data")


# TODO; have swfo convert the files to HDF5
def ecwmf_dewpoint_temperature(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF 2 metre Dewpoint
    Temperature collection.
    """
    product = DatasetName.DEWPOINT_TEMPERATURE.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat)

            metadata = {'data_source': 'ECWMF 2 metre Dewpoint Temperature ',
                        'url': url,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        raise AncillaryError("No ECWMF 2 metre Dewpoint Temperature data")


# TODO; have swfo convert the files to HDF5
def ecwmf_surface_pressure(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Surface Pressure
    collection.
    Scales the result by 100 before returning.
    """
    product = DatasetName.SURFACE_PRESSURE.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat) / 100.0

            metadata = {'data_source': 'ECWMF Surface Pressure',
                        'url': url,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        raise AncillaryError("No ECWMF Surface Pressure data")


# TODO; have swfo convert the files to HDF5
def ecwmf_water_vapour(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Total Column Water Vapour
    collection.
    """
    product = DatasetName.WATER_VAPOUR.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat)

            metadata = {'data_source': 'ECWMF Total Column Water Vapour',
                        'url': url,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            return data, metadata

    if data is None:
        raise AncillaryError("No ECWMF Total Column Water Vapour data")


# TODO; have swfo convert the files to HDF5
def ecwmf_temperature(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Temperature collection
    across 37 height pressure levels, for a given longitude,
    latitude and time.

    Reverses the order of elements
    (1000 -> 1 mb, rather than 1 -> 1000 mb) before returning.
    """
    product = DatasetName.TEMPERATURE.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            bands = list(range(1, 38))
            data = get_pixel(f, lonlat, bands)[::-1]

            metadata = {'data_source': 'ECWMF Temperature',
                        'url': url,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_metadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'Temperature', data)

            return df, metadata

    if data is None:
        raise AncillaryError("No ECWMF Temperature profile data")


# TODO; have swfo convert the files to HDF5
def ecwmf_geo_potential(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Geo-Potential collection
    across 37 height pressure levels, for a given longitude,
    latitude and time.

    Converts to geo-potential height in KM, and reverses the order of
    the elements (1000 -> 1 mb, rather than 1 -> 1000 mb) before
    returning.
    """
    product = DatasetName.GEOPOTENTIAL.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            bands = list(range(1, 38))
            data = get_pixel(f, lonlat, bands)[::-1]
            scaled_data = data / 9.80665 / 1000.0

            metadata = {'data_source': 'ECWMF Geo-Potential',
                        'url': url,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_metadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'GeoPotential', data)
            df.insert(1, 'GeoPotential_Height', scaled_data)

            return df, md

    if data is None:
        raise AncillaryError("No ECWMF Geo-Potential profile data")


# TODO; have swfo convert the files to HDF5
def ecwmf_relative_humidity(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF Relative Humidity collection
    across 37 height pressure levels, for a given longitude,
    latitude and time.

    Reverses the order of elements
    (1000 -> 1 mb, rather than 1 -> 1000 mb) before returning.
    """
    product = DatasetName.RELATIVE_HUMIDITY.value.lower()
    search = pjoin(input_path, DatasetName.ECMWF_PATH_FMT.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        url = urlparse(f, scheme='file').geturl()
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            bands = list(range(1, 38))
            data = get_pixel(f, lonlat, bands)[::-1]

            metadata = {'data_source': 'ECWMF Relative Humidity',
                        'url': url,
                        'query_date': time}

            # file level metadata
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_metadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'Relative_Humidity', data)

            return df, metadata

    if data is None:
        raise AncillaryError("No ECWMF Relative Humidity profile data")
