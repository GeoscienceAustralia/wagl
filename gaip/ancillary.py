#!/usr/bin/env python

"""
Ancillary dataset retrieval and storage
"""

from __future__ import absolute_import, print_function
import logging
from os.path import join as pjoin, basename, splitext
from posixpath import join as ppjoin
import datetime
import glob
import numpy
import h5py
import pandas
from geopandas import GeoSeries
import rasterio
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely import wkt
from gaip.brdf import get_brdf_data
from gaip.data import get_pixel
from gaip.hdf5 import attach_attributes, write_scalar, write_dataframe
from gaip.hdf5 import read_h5_table, dataset_compression_kwargs
from gaip.hdf5 import attach_table_attributes
from gaip.metadata import extract_ancillary_metadata, read_meatadata_tags
from gaip.constants import DatasetName, POINT_FMT, GroupName, BandType
from gaip.satellite_solar_angles import create_vertices


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

    rh = 100 * ((112.0 - 0.1 * surf_t + dew_t) / (112.0 + 0.9 * surf_t)) **8

    return rh


def _collect_ancillary(container, satellite_solar_fname, nbar_paths,
                       sbt_path=None, invariant_fname=None, vertices=(3, 3),
                       out_fname=None, compression='lzf'):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_fname, 'r') as fid,\
        h5py.File(out_fname, 'w') as out_fid:

        sat_sol_grp = fid[GroupName.sat_sol_group.value]
        collect_ancillary(container, sat_sol_grp, nbar_paths, sbt_path,
                          invariant_fname, vertices, out_fid, compression)

    return


def collect_ancillary(container, satellite_solar_group, nbar_paths,
                      sbt_path=None, invariant_fname=None, vertices=(3, 3),
                      out_group=None, compression='lzf'):
    """
    Collects the ancillary required for NBAR and optionally SBT.
    This could be better handled if using the `opendatacube` project
    to handle ancillary retrieval, rather than directory passing,
    and filename grepping.

    :param container:
        An instance of an `AcquisitionsContainer` object.

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the solar zenith and
        solar azimuth datasets specified by the pathnames given by:

        * DatasetName.boxline

    :param nbar_paths:
        A `dict` containing the ancillary pathnames required for
        retrieving the NBAR ancillary data. Required keys:

        * aerosol_fname
        * water_vapour_path
        * ozone_path
        * dem_path
        * brdf_path
        * brdf_premodis_path

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
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # Initialise the output files
    if out_group is None:
        fid = h5py.File('ancillary.h5', driver='core', backing_store=False)
    else:
        fid = out_group

    group = fid.create_group(GroupName.ancillary_group.value)

    acquisition = container.get_acquisitions()[0]

    boxline_dataset = satellite_solar_group[DatasetName.boxline.value][:]
    coordinator = create_vertices(acquisition, boxline_dataset, vertices)
    lonlats = zip(coordinator['longitude'], coordinator['latitude'])

    desc = ("Contains the row and column array coordinates used for the "
            "atmospheric calculations.")
    attrs = {'Description': desc, 'array_coordinate_offset': 0}
    kwargs = dataset_compression_kwargs(compression=compression)
    dset_name = DatasetName.coordinator.value
    coord_dset = group.create_dataset(dset_name, data=coordinator, **kwargs)
    attach_table_attributes(coord_dset, title='Coordinator', attrs=attrs)


    if sbt_path:
        collect_sbt_ancillary(acquisition, lonlats, sbt_path, invariant_fname,
                              out_group=group, compression=compression)

    collect_nbar_ancillary(container, out_group=group, 
                           compression=compression, **nbar_paths)

    if out_group is None:
        return fid


def collect_sbt_ancillary(acquisition, lonlats, ancillary_path,
                          invariant_fname=None, out_group=None,
                          compression='lzf'):
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
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

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
    attrs = {'Description': description,
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
        dname = ppjoin(pnt, DatasetName.dewpoint_temperature.value)
        write_scalar(dew[0], dname, fid, dew[1])

        dname = ppjoin(pnt, DatasetName.temperature_2m.value)
        write_scalar(t2m[0], dname, fid, t2m[1])

        dname = ppjoin(pnt, DatasetName.surface_pressure.value)
        write_scalar(sfc_prs[0], dname, fid, sfc_prs[1])

        dname = ppjoin(pnt, DatasetName.surface_geopotential.value)
        write_scalar(sfc_hgt[0], dname, fid, sfc_hgt[1])

        dname = ppjoin(pnt, DatasetName.surface_relative_humidity.value)
        attrs = {'Description': 'Relative Humidity calculated at the surface'}
        write_scalar(sfc_rh, dname, fid, attrs)

        # get the data from each of the pressure levels (1 -> 1000 ISBL)
        gph = ecwmf_geo_potential(ancillary_path, lonlat, dt)
        tmp = ecwmf_temperature(ancillary_path, lonlat, dt)
        rh = ecwmf_relative_humidity(ancillary_path, lonlat, dt)

        dname = ppjoin(pnt, DatasetName.geopotential.value)
        write_dataframe(gph[0], dname, fid, compression, attrs=gph[1])

        dname = ppjoin(pnt, DatasetName.temperature.value)
        write_dataframe(tmp[0], dname, fid, compression, attrs=tmp[1])

        dname = ppjoin(pnt, DatasetName.relative_humidity.value)
        write_dataframe(rh[0], dname, fid, compression, attrs=rh[1])

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

        dname = ppjoin(pnt, DatasetName.atmospheric_profile.value)
        write_dataframe(df, dname, fid, compression, attrs=attrs)

        fid[pnt].attrs['lonlat'] = lonlat

    if out_group is None:
        return fid


def collect_nbar_ancillary(container, aerosol_fname=None,
                           water_vapour_path=None, ozone_path=None,
                           dem_path=None, brdf_path=None,
                           brdf_premodis_path=None, out_group=None,
                           compression='lzf'):
    """
    Collects the ancillary information required to create NBAR.

    :param container:
        An instance of an `AcquisitionsContainer` object.

    :param aerosol_fname:
        A `str` containing the full file pathname to the `HDF5` file
        containing the aerosol data.

    :param water_vapour_path:
        A `str` containing the full file pathname to the directory
        containing the water vapour data.

    :param ozone_path:
        A `str` containing the full file pathname to the directory
        containing the ozone data.

    :param dem_path:
        A `str` containing the full file pathname to the directory
        containing the digital elevation model data.

    :param brdf_path:
        A `str` containing the full file pathname to the directory
        containing the BRDF image mosaics.

    :param brdf_premodis_path:
        A `str` containing the full file pathname to the directory
        containing the premodis BRDF image mosaics.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

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
    """
    # Initialise the output files
    if out_group is None:
        fid = h5py.File('nbar-ancillary.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    acquisition = container.get_acquisitions()[0]
    dt = acquisition.acquisition_datetime
    geobox = acquisition.gridded_geo_box()

    aerosol = get_aerosol_data(acquisition, aerosol_fname)
    write_scalar(aerosol[0], DatasetName.aerosol.value, fid, aerosol[1])

    wv = get_water_vapour(acquisition, water_vapour_path)
    write_scalar(wv[0], DatasetName.water_vapour.value, fid, wv[1])

    ozone = get_ozone_data(ozone_path, geobox.centre_lonlat, dt)
    write_scalar(ozone[0], DatasetName.ozone.value, fid, ozone[1])

    elev = get_elevation_data(geobox.centre_lonlat, dem_path)
    write_scalar(elev[0], DatasetName.elevation.value, fid, elev[1])

    # brdf
    group = fid.create_group('brdf-image-datasets')
    dname_format = DatasetName.brdf_fmt.value
    for group in container.groups:
        for acq in container.get_acquisitions(group=group):
            if acq.band_type is not BandType.Reflective:
                continue
            data = get_brdf_data(acq, brdf_path, brdf_premodis_path, group,
                                 compression)

            # output
            for param in data:
                dname = dname_format.format(parameter=param.name,
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
        group = granule[GroupName.ancillary_group.value]

        ozone += group[DatasetName.ozone.value][()]
        vapour += group[DatasetName.water_vapour.value][()]
        aerosol += group[DatasetName.aerosol.value][()]
        elevation += group[DatasetName.elevation.value][()]

    # average
    ozone /= n_tiles
    vapour /= n_tiles
    aerosol /= n_tiles
    elevation /= n_tiles

    description = ("The {} value is an average from all the {} values "
                   "retreived for each Granule.")
    attrs = {'data_source': 'granule_average'}

    # output each average value back into the same granule ancillary group
    group_name = ppjoin(GroupName.ancillary_group.value,
                        GroupName.ancillary_avg_group.value)
    for granule in granule_groups:
        # for the multifile workflow, we only want to write to one granule
        try:
            group = granule.create_group(group_name)
        except ValueError:
            continue

        dset = group.create_dataset(DatasetName.ozone.value, data=ozone)
        attrs['Description'] = description.format(*(2*['Ozone']))
        attach_attributes(dset, attrs)

        dset = group.create_dataset(DatasetName.water_vapour.value, data=vapour)
        attrs['Description'] = description.format(*(2*['Water Vapour']))
        attach_attributes(dset, attrs)

        dset = group.create_dataset(DatasetName.aerosol.value, data=aerosol)
        attrs['Description'] = description.format(*(2*['Aerosol']))
        attach_attributes(dset, attrs)

        dset = group.create_dataset(DatasetName.elevation.value, data=elevation)
        attrs['Description'] = description.format(*(2*['Elevation']))
        attach_attributes(dset, attrs)


def get_aerosol_data(acquisition, aerosol_fname):
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

    fid = h5py.File(aerosol_fname, 'r')

    delta_tolerance = datetime.timedelta(days=0.5)

    data = None
    for pathname, description in zip(pathnames, descr):
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
                    metadata = {'data_source': description,
                                'dataset_pathname': pathname,
                                'query_date': dt,
                                'data_file': aerosol_fname,
                                'extents': wkt.dumps(intersection)}

                    # ancillary metadata tracking
                    md = extract_ancillary_metadata(aerosol_fname)
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
    data = get_pixel(datafile, lonlat) * 0.001  # scale to correct units

    metadata = {'data_source': 'Elevation',
                'data_file': datafile}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
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
    data = get_pixel(datafile, lonlat)

    metadata = {'data_source': 'Ozone',
                'data_file': datafile,
                'query_date': time}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def get_water_vapour(acquisition, vapour_path, scale_factor=0.1):
    """
    Retrieve the water vapour value for an `acquisition` and the
    path for the water vapour ancillary data.
    """
    dt = acquisition.acquisition_datetime
    geobox = acquisition.gridded_geo_box()

    year = dt.strftime('%Y')
    filename = "pr_wtr.eatm.{year}.tif".format(year=year)
    datafile = pjoin(vapour_path, filename)

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
        data = get_pixel(datafile, geobox.centre_lonlat, band=band)
    except IndexError:
        msg = "Invalid water vapour band number: {band}".format(band=band)
        raise IndexError(msg)

    data = data * scale_factor

    metadata = {'data_source': 'Water Vapour',
                'data_file': datafile,
                'query_date': dt}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def ecwmf_elevation(datafile, lonlat):
    """
    Retrieve a pixel from the ECWMF invariant geo-potential
    dataset.
    Converts to Geo-Potential height in KM.
    2 metres is added to the result before returning.
    """
    data = get_pixel(datafile, lonlat) / 9.80665 / 1000.0 + 0.002

    metadata = {'data_source': 'ECWMF Invariant Geo-Potential',
                'data_file': datafile}

    # ancillary metadata tracking
    md = extract_ancillary_metadata(datafile)
    for key in md:
        metadata[key] = md[key]

    return data, metadata


def ecwmf_temperature_2metre(input_path, lonlat, time):
    """
    Retrieve a pixel value from the ECWMF 2 metre Temperature
    collection.
    """
    product = 'temperature-2metre'
    files = glob.glob(pjoin(input_path, '{}_*.tif'.format(product)))
    product = DatasetName.temperature_2m.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat)

            metadata = {'data_source': 'ECWMF 2 metre Temperature',
                        'data_file': f,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
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
    product = DatasetName.dewpoint_temperature.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat)

            metadata = {'data_source': 'ECWMF 2 metre Dewpoint Temperature ',
                        'data_file': f,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
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
    product = DatasetName.surface_pressure.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat) / 100.0

            metadata = {'data_source': 'ECWMF Surface Pressure',
                        'data_file': f,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
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
    product = DatasetName.water_vapour.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            data = get_pixel(f, lonlat)

            metadata = {'data_source': 'ECWMF Total Column Water Vapour',
                        'data_file': f,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
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
    product = DatasetName.temperature.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            bands = list(range(1, 38))
            data = get_pixel(f, lonlat, bands)[::-1]

            metadata = {'data_source': 'ECWMF Temperature',
                        'data_file': f,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
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
    product = DatasetName.geopotential.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            bands = list(range(1, 38))
            data = get_pixel(f, lonlat, bands)[::-1]
            scaled_data = data / 9.80665 / 1000.0

            metadata = {'data_source': 'ECWMF Geo-Potential',
                        'data_file': f,
                        'query_date': time}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_meatadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'GeoPotential', data)
            df.insert(1, 'GeoPotential_Height', scaled_data)

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
    product = DatasetName.relative_humidity.value
    search = pjoin(input_path, DatasetName.ecmwf_path_fmt.value)
    files = glob.glob(search.format(product=product, year=time.year))
    data = None
    required_ymd = datetime.datetime(time.year, time.month, time.day)
    for f in files:
        ymd = splitext(basename(f))[0].split('_')[1]
        ancillary_ymd = datetime.datetime.strptime(ymd, '%Y-%m-%d')
        if ancillary_ymd == required_ymd:
            bands = list(range(1, 38))
            data = get_pixel(f, lonlat, bands)[::-1]

            metadata = {'data_source': 'ECWMF Relative Humidity',
                        'data_file': f,
                        'query_date': time}

            # file level metadata
            md = extract_ancillary_metadata(f)
            for key in md:
                metadata[key] = md[key]

            # internal file metadata (and reverse the ordering)
            df = read_meatadata_tags(f, bands).iloc[::-1]
            df.insert(0, 'Relative_Humidity', data)

            return df, metadata

    if data is None:
        raise IOError('No Relative Humidity ancillary data found.')
