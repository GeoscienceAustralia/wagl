"""
Satellite TLE (two-line element) Loading
----------------------------------------
"""
from __future__ import absolute_import, print_function, unicode_literals
import datetime
import ephem
import re
import os

TLE_ENTRY_RE = (r'^([1])(\s+)([%(NUMBER)s%(CLASSIFICATION)s]+)(\s+)'
                r'([%(INTL_DESIGNATOR)s]+)(\s+)%(YYDDD)s(\.)(\d+)(\s+)'
                r'([\s\-]+)(\.)(\d+)(\s+)(\d+)([\-\+])(\d+)(\s+)(\d+)'
                r'([\-\+])(\d+)(\s+)(\d+)(\s+)(\d+)(\s)^([2])(.+)$')


def load_tle(acquisition, data_root, date_radius=45):
    """
    Loads satellite TLE (two-line element) for the given date and time.

    Arguments:
        centre_datetime: datetime.datetime instance
        data_root: root directory for TLE archive files
        date_radius: date radius for TLE search (days)

    Returns:
        ephem EarthSatellite instance

    """
    return (load_tle_from_archive(acquisition, data_root, date_radius) or
            load_tle_from_files(acquisition, data_root, date_radius))


def load_tle_from_archive(acquisition, data_root, day_radius=45):
    """Loads TLE (two-line element) for the satellite, date and time from
    an archive file.

    Arguments:
        centre_datetime: datetime.datetime instance
        data_root: directory containing TLE archive files
        date_radius: date radius for TLE search (days)

    Returns:
        ephem EarthSatellite instance

    """
    center_datetime = acquisition.scene_center_datetime

    offsets = sorted(range(-day_radius, day_radius),
                     cmp=lambda x, y: abs(x) - abs(y))
    tds = [datetime.timedelta(days=d) for d in offsets]
    yyddd_list = [(center_datetime + d).strftime('%02y%03j') for d in tds]

    name = acquisition.satellite_name.replace('-', '').upper()

    tle_archive_path = os.path.join(data_root, name,
                                    'TLE', '%s_ARCHIVE.txt' % acquisition.tag)

    text = ''
    try:
        with open(tle_archive_path, 'r') as fd:
            text = fd.read()
    except IOError:
        # no TLE archive file exists
        return None

    re_params = {
        'NUMBER': acquisition.number,
        'CLASSIFICATION': acquisition.classification,
        'INTL_DESIGNATOR': acquisition.intl_designator}

    for yyddd in yyddd_list:
        re_params['YYDDD'] = yyddd
        match = re.search(TLE_ENTRY_RE % re_params, text, re.MULTILINE)
        if match:
            # Reconstitute TLE entry from regex match groups.
            tle_text = (''.join(match.groups()[0:6])
                        + yyddd + ''.join(match.groups()[6:]))
            lines = tle_text.split('\n')
            return ephem.readtle(acquisition.spacecraft_id, lines[0], lines[1])

    return None


def load_tle_from_files(acquisition, data_root, day_range=45):
    """Load a TLE file for the specified datetime.

    Arguments:
        center_datetime: scene center datetime (datetime instance)
        data_root: ephemeris data root directory

    Returns:
        ephem EarthSatellite instance
    """

    name = acquisition.satellite_name.replace('-', '').upper()

    def open_tle(tle_path, center_datetime):
        """Open the TLE file and read."""
        with open(tle_path, 'r') as fd:
            tle_text = fd.readlines()
            if acquisition.tag == 'LS5':
                tle1, tle2 = tle_text[7:9]
            if acquisition.tag == 'LS7':
                tle1, tle2 = tle_text[1:3]
            return ephem.readtle(acquisition.satellite_name, tle1, tle2)

    center_datetime = acquisition.scene_center_datetime
    scene_doy = center_datetime.strftime('%j')  # Note format: '%03d'

    tle_dir = os.path.join(data_root,
                           name,
                           'TLE',
                           '%s_YEAR' % acquisition.tag,
                           '%4d' % center_datetime.year)

    tle_file = acquisition.tle_format % (center_datetime.year, scene_doy)
    tle_path = os.path.join(tle_dir, tle_file)

    if os.path.exists(tle_path):
        try:
            return open_tle(tle_path, center_datetime)
        except IOError:
            pass

    for d in range(1, day_range):
        ddelta = datetime.timedelta(days=d)
        for s in [-1, 1]:
            dt = center_datetime + (ddelta * s)
            tle_dir = os.path.join(data_root, name,
                                   'TLE',
                                   '%s_YEAR' % acquisition.tag,
                                   '%4d' % dt.year)
            tle_file = acquisition.tle_format % (dt.year, dt.strftime('%j'))
            tle_path = os.path.join(tle_dir, tle_file)
            if os.path.exists(tle_path):
                try:
                    return open_tle(tle_path, center_datetime)
                except IOError:
                    pass

    return None
