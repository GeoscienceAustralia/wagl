"""
Contains the base implementations for the acquisition and AcquisitionsContainer objects
"""
from os.path import join as pjoin
from functools import total_ordering
from pkg_resources import resource_stream

import numpy
import rasterio

from ..geobox import GriddedGeoBox
from ..modtran import read_spectral_response
from ..tiling import generate_tiles
from ..constants import BandType

class AcquisitionsContainer(object):

    """
    A container for dealing with a hierarchial structure
    of acquisitions from different resolution groups, granules,
    but all part of the same geospatial area or scene.

    Note: Assuming that each granule contains the same resolution
    groups.
    """

    def __init__(self, label, granules):
        self._granules = granules
        self._label = label

    def __repr__(self):
        fmt = ("****Granule ID's****:\n{granules}\n"
               "\n****Resolution Groups****:\n{groups}")
        granules = "\n".join(self.granules)
        groups = "\n".join(self.groups)
        return fmt.format(granules=granules, groups=groups)

    @property
    def label(self):
        """
        Return the scene label.
        """
        return self._label

    @property
    def granules(self):
        """
        Lists the available granules within a scene.
        """
        return sorted(list(self._granules.keys()))

    @property
    def groups(self):
        """
        Lists the available groups within a scene.
        """
        return list(self._granules.get(self.granules[0]).keys())

    def get_acquisitions(self, group=None, granule=None, only_supported_bands=True):
        """
        Return a list of acquisitions for a given granule and group.

        :param group:
            A `str` defining the group layer from which to retrieve
            the acquisitions from. If `None` (default), return the
            acquisitions from the first group in the
            `AcquisitionsContainer.groups` list.

        :param granule:
            A `str` defining the granule layer from which to retrieve
            the acquisitions from. If `None` (default), return the
            acquisitions from the the first granule in the
            `AcquisitionsContainer.granule` list.

        :param only_supported_bands:
            boolean if set to True will return all bands that are
            defined in acquisition/sensors.json for the related platform.
            If set to False it will return all acquisitions

        :return:
            A `list` of `Acquisition` objects.
        """
        groups = self.get_granule(granule=granule)
        if group is None:
            acqs = groups[list(groups.keys())[0]]
        else:
            acqs = groups[group]

        if only_supported_bands:
            return list(
                filter(lambda acq: getattr(acq, 'supported_band', False) is True,
                acqs)
            )
        else:
            return acqs

    def get_granule(self, granule=None, container=False):
        """
        Return a granule containing groups of `Acquisition` objects.

        :param granule:
            A `str` defining the granule layer from which to retrieve
            groups of `Acquisition` objects. Default is `None`, which
            returns the the first granule in the
            `AcquisitionsContainer.granule` list.

        :param container:
            A boolean indicating whether to return the granule as an
            `AcquisitionsContainer` containing a single granule.
            Default is False.

        :return:
            A `dict` containing the groups of `Acquisition` objects
            for a given scene, unless container=True in which a new
            instance of an `AcquisitionsContainer` is returned.
        """
        if granule is None:
            granule = self.granules[0]

        if container:
            grps = {granule: self._granules[granule]}
            return AcquisitionsContainer(label=self.label, granules=grps)

        return self._granules[granule]

    def get_root(self, path='/', group=None, granule=None):
        """
        Get the root level file system path for a granule and/or group
        within the `AcquisitionsContainer` object.

        :param path:
            A `str` containing the root path on which to join the
            granule and/or group layers onto.

        :param group:
            A `str` containing the group layer to be joined onto
            `path`. If group is `None` (default), or `not in`
            `AcquisitionsContainer.groups` then no path join occurs.

        :param granule:
            A `str` containing the granule layer to be joined onto
            `path`. If granule is `None` (default), or `not in`
            `AcquisitionsContainer.granules` then no path join occurs.

        :return:
            A `str` representing the combined path for group and/or
            granule layers.
        """
        if (granule is None) or (granule not in self.granules):
            root = path
        else:
            root = pjoin(path, granule)

        if (group is not None) and (group in self.groups):
            root = pjoin(root, group)

        return root

    def get_highest_resolution(self, granule=None):
        """
        Retrieve a list of supported acquisitions from the highest
        resolution group for a given granule. The default is to
        return the first listed granule, and the first resolution
        group containing a supported band.
        Meaning that while there might exist a higher resolution
        group, the bands might not be supported, in which case
        the next highest resolution group will be searched.

        :return:
            A `tuple` of (`list`, group_name) where the `list`
            contains the acquisitions of supported bands from
            the highest resolution group, and group_name is the name
            of the group that the acquisitions came from.
        """
        for group in self.groups:
            acqs = self.get_acquisitions(group, granule)
            if acqs:
                break
        return acqs, group

    def get_mode_resolution(self, granule=None):
        """
        Retrieve the resolution group name that contains the mode
        resolution. 

        :return:
            A `tuple` of (`list`, group_name) where the `list`
            contains the acquisitions of supported bands from
            the mode resolution group, and group_name is the name
            of the group that the acquisitions came from.
        """
        dtype = numpy.dtype([('group', 'O'), ('frequency', 'int64')])
        data = numpy.zeros(len(self.groups), dtype=dtype)

        # insert data
        data['group'] = self.groups
        counts = [len(self.get_acquisitions(grp, None, False)) for grp in
                  self.groups]
        data['frequency'] = counts

        mode_grp = data['group'][data['frequency'].argmax()]
        acqs = self.get_acquisitions(mode_grp, granule)

        return acqs, mode_grp

    @property
    def supported_groups(self):
        """
        Return a list of resolution groups that have supported
        bands.
        """
        groups = []
        for group in self.groups:
            acqs = self.get_acquisitions(group=group)
            if acqs:
                groups.append(group)

        return groups

    def get_all_acquisitions(self, granule=None):
        """
        Retrieve all supported acquisitions from a given granule.
        This will ignore the fact that acquisitions could be from
        different resolution groups, and put them all in a single
        list.
        """
        acquisitions = []
        for group in self.groups:
            acqs =  self.get_acquisitions(group, granule)
            if acqs:
                acquisitions.extend(acqs)

        return acquisitions


@total_ordering
class Acquisition(object):

    """Acquisition metadata."""

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        self._pathname = pathname
        self._uri = uri
        self._acquisition_datetime = acquisition_datetime
        self._band_name = band_name
        self._band_id = band_id

        self._norad_id = None
        self._classification_type = None
        self._international_designator = None

        self._gps_file = False

        if metadata is not None:
            for key, value in metadata.items():
                if key == 'band_type':
                    value = BandType[value]
                setattr(self, key, value)

        self._open()

    def _open(self):
        """
        A private method for opening the dataset and
        retrieving the dimensional information.
        """
        with rasterio.open(self.uri) as ds:
            self._samples = ds.width
            self._lines = ds.height
            self._tile_size = ds.block_shapes[0]
            self._resolution = ds.res

    @property
    def pathname(self):
        """
        The pathname of the level1 dataset.
        """
        return self._pathname

    @property
    def uri(self):
        """
        The uri of the acquisition.
        """
        return self._uri

    @property
    def acquisition_datetime(self):
        """
        The acquisitions centre scantime.
        """
        return self._acquisition_datetime

    @property
    def band_name(self):
        """
        The band name, which goes by the format `BAND {}`.
        """
        return self._band_name

    @property
    def band_id(self):
        """
        The band id as given in the `sensors.json` file.
        """
        return self._band_id

    @property
    def norad_id(self):
        """The NORAD catalog id number."""
        return self._norad_id

    @property
    def classification_type(self):
        """The classification type; eg 'U' = unclassified."""
        return self._classification_type

    @property
    def international_designator(self):
        """The international designator."""
        return self._international_designator

    @property
    def samples(self):
        """The number of samples (aka. `width`)."""
        return self._samples

    @property
    def lines(self):
        """The number of lines (aka. `height`)."""
        return self._lines

    @property
    def tile_size(self):
        """
        The native tile size of the file on disk in
        (ysize, xsize) dimensions.
        """
        return self._tile_size

    @property
    def resolution(self):
        """
        The resolution of the file on disk reported as
        (y, x).
        """
        return self._resolution

    @property
    def gps_file(self):
        """
        Does the acquisition have an associated GPS file?
        """
        return self._gps_file

    def __eq__(self, other):
        return self.band_name == other.band_name

    def __lt__(self, other):
        return self.sortkey() < other.sortkey()

    def __repr__(self):
        return 'Acquisition(band_name=' + self.band_name + ')'

    def sortkey(self):
        """Representation used for sorting objects."""
        return self.band_name

    def data(self, out=None, window=None, masked=False):
        """
        Return `numpy.array` of the data for this acquisition.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        """
        with rasterio.open(self.uri) as ds:
            data = ds.read(1, out=out, window=window, masked=masked)

        return data

    def radiance_data(self, window=None, out_no_data=-999):
        """
        Return the data as radiance in watts/(m^2*micrometre).
        Override with a custom version for a specific sensor.
        """
        raise NotImplementedError

    def data_and_box(self, out=None, window=None, masked=False):
        """
        Return a tuple comprising the `numpy.array` of the data for this
        Acquisition and the `GriddedGeoBox` describing the spatial extent.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        for this acquisition.
        """
        with rasterio.open(self.uri) as ds:
            box = GriddedGeoBox.from_dataset(ds)
            if window is not None:
                rows = window[0][1] - window[0][0]
                cols = window[1][1] - window[1][0]
                prj = ds.crs.wkt
                res = ds.res

                # Get the new UL co-ordinates of the array
                ul_x, ul_y = ds.transform * (window[1][0], window[0][0])
                box = GriddedGeoBox(shape=(rows, cols), origin=(ul_x, ul_y),
                                    pixelsize=res, crs=prj)
            return (ds.read(1, out=out, window=window, masked=masked), box)

    def gridded_geo_box(self):
        """Return the `GriddedGeoBox` for this acquisition."""
        with rasterio.open(self.uri) as src:
            return GriddedGeoBox.from_dataset(src)

    def decimal_hour(self):
        """The time in decimal."""
        time = self.acquisition_datetime
        dec_hour = (time.hour + (time.minute + (time.second
                                                + time.microsecond / 1000000.0)
                                 / 60.0) / 60.0)
        return dec_hour

    def julian_day(self):
        """
        Return the Juilan Day of the acquisition_datetime.
        """
        return int(self.acquisition_datetime.strftime('%j'))

    @property
    def no_data(self):
        """
        Return the no_data value for this acquisition.
        Assumes that the acquisition is a single band file.
        """
        with rasterio.open(self.uri) as ds:
            nodata_list = ds.nodatavals
            return nodata_list[0]

    def spectral_response(self, as_list=False):
        """
        Reads the spectral response for the sensor.
        """
        fname = '../spectral_response/%s' % self.spectral_filter_file
        spectral_range = range(*self.spectral_range)
        with resource_stream(__name__, fname) as src:
            df = read_spectral_response(src, as_list, spectral_range)
        return df

    def close(self):
        """
        A simple additional utility for acquisitions that need
        to close any open files or set any specific properties to
        None. Used as a general cleanup.
        This utility might change over time as a better mechanism
        for handling various read methods is resolved.
        Override as needed.
        """
        pass

    def tiles(self):
        """
        Generate the tiling regime for this acquisition.
        """
        ysize, xsize = self.tile_size
        return generate_tiles(self.samples, self.lines, xsize, ysize)

    def close(self):
        """
        A method used by child classes to clear cached datasets
        """
        pass
