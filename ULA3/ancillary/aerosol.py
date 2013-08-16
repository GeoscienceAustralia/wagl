"""
Utilities for the extraction of aerosol data. These assume very specific file formats which are read by the
external program ``aot_loader`` (which must be located in the directory specified by the argument ``bin_dir``
of the functions contained herein).

More information can be found in TRIM document D2012-86155, which has been copied to
:download:`here <auxiliary/NBAR ANCILLARY DATA DOWNLOAD SCRIPTS Release Notes For 15 May.docx>` (which was
provided by `Paul Gardner <mailto:paul.gardner@ga.gov.au>`_) for your convenience.
"""
import logging, re, os
from datetime import datetime
from ULA3.utils import execute

logger = logging.getLogger('root.' + __name__)

SCENE_MONTH_NAMES = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')





def get_aerosol_value_for_region(
    dt,
    ll_lat, ll_lon,
    ur_lat, ur_lon,
    default_value,
    aerosol_dir,
    bin_dir,
    enable_aeronet=False,
    sc_lat=None, sc_lon=None,):
    """
    Extract the aerosol value for a region (typically a scene).

    :param dt:
        The date and time to extract the value for.
    :type dt:
        :py:class:`datetime.datetime`

    :param ll_lat:
        The latitude of the lower left corner of the region ('ll' for 'Lower Left').
    :type ll_lat:
        :py:class:`float`

    :param ll_lon:
        The longitude of the lower left corner of the region ('ll' for 'Lower Left').
    :type ll_lon:
        :py:class:`float`

    :param ur_lat:
        The latitude of the upper right corner of the region ('ur' for 'Upper Right').
    :type ur_lat:
        :py:class:`float`

    :param ur_lon:
        The longitude of the upper right corner of the region ('ur' for 'Upper Right').
    :type ur_lon:
        :py:class:`float`

    :param default_value:
        If no suitable value can be found from the AATSR or (possibly) AERONET data, then
        return this value instead.
    :type default_value:
        :py:class:`float`

    :param aerosol_dir:
        The directory containing the aerosol data files.
    :type aerosol_dir:
        :py:class:`str`

    :param bin_dir:
        The directory where the executable ``aot_loader`` can be found.
    :type bin_dir:
        :py:class:`str`

    :param enable_aeronet:
        Once upon a time in a version from long long ago,
        `AERONET <http://gcmd.nasa.gov/records/GCMD_AERONET_NASA.html>`_ aerosol data was also available (GA
        now only maintains `AATSR <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_ data for internal use).
        If this parameter is set to ``True``, then if a suitable (i.e. non-zero) value is not obtained from
        the `AATSR <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_ data, then the function will look for
        a value from AERONET data (you would have to ensure such data is available, which is not currently
        the case within GA).
    :type enable_aeronet:
        :py:class:`bool`

    :param sc_lat:
        The latitude of the center of the region ('sc' for 'Scene Center').
    :type sc_lat:
        :py:class:`float`

    :param sc_lon:
        The longitude of the center of the region ('sc' for 'Scene Center').
    :type sc_lon:
        :py:class:`float`

    While the natural place for the parameters sc_lat and sc_lon is just before the lower left and upper
    right locations, these are only required if ``enable_aeronet`` is ``True`` and hence, have been made
    optional.
    """

    # Parse the date string so we can use it to construct AATSR data paths.
    # self.nbar.centreDate should be of the form "2009-11-08".

    scene_month_name = SCENE_MONTH_NAMES[dt.month-1]

    # AATSR pix.
    data_desc = 'AATSR_PIX'
    atsr_file = os.path.join(aerosol_dir, 'ATSR_LF_%s%s.pix' % (dt.year, dt.month))
    value = get_aerosol_aatsr(
        atsr_file,
        dt,
        ll_lat, ll_lon,
        ur_lat, ur_lon,
        bin_dir)
    if value != 0.0:
        return {
            'data_source': data_desc,
            'data_file': atsr_file,
            'value': value}

    # AATSR year-month composite.
    data_desc = 'AATSR_CMP_YEAR_MONTH'
    atsr_file = os.path.join(
        aerosol_dir,
        'aot_mean_%s_%s_All_Aerosols.cmp' % (scene_month_name, dt.year))
    value = get_aerosol_aatsr(
        atsr_file, dt,
        ll_lat, ll_lon,
        ur_lat, ur_lon,
        bin_dir)
    if value != 0.0:
        return {
            'data_source': data_desc,
            'data_file': atsr_file,
            'value': value}

    # AATSR month composite.
    data_desc = 'AATSR_CMP_MONTH'
    atsr_file = os.path.join(
        aerosol_dir,
        'aot_mean_%s_All_Aerosols.cmp' % scene_month_name)
    value = get_aerosol_aatsr(
        atsr_file, dt,
        ll_lat, ll_lon,
        ur_lat, ur_lon,
        bin_dir)
    if value != 0.0:
        return {
            'data_source': data_desc,
            'data_file': atsr_file,
            'value': value}

    # Aeronet (CURRENTLY DISABLED)
    if enable_aeronet:
        assert sc_lat != None and sc_lon != None, "Both sc_lon and sc_lat must be specified (enable_aeronet is True and no suitable value found in AATSR data)"
        data_desc = 'AERONET'
        value = get_aerosol_aeronet(aerosol_dir, dt, sc_lat, sc_lon)
        if value != 0.0:
            return {
                'data_source': data_desc,
                'data_file': None,  # TODO aeronet data file path
                'value': value}

    # Default
    data_desc = 'DEFAULT'
    return {
        'data_source': data_desc,
        'data_file': None,
        'value': default_value}





def get_aerosol_aatsr(
    atsr_file,
    dt,
    ll_lat, ll_lon,
    ur_lat, ur_lon,
    bin_dir):
    """
    Load aerosol data for a specified `AATSR <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_ data file.
    This uses the executable ``aot_loader``.

    :param atsr_file:
        The full path to the `AATSR <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_ file to load the data
        from.
    :type atsr_file:
        :py:class:`str`

    :param dt:
        The date and time to extract the value for.
    :type dt:
        :py:class:`datetime.datetime`

    :param ll_lat:
        The latitude of the lower left corner of the region ('ll' for 'Lower Left').
    :type ll_lat:
        :py:class:`float`

    :param ll_lon:
        The longitude of the lower left corner of the region ('ll' for 'Lower Left').
    :type ll_lon:
        :py:class:`float`

    :param ur_lat:
        The latitude of the upper right corner of the region ('ur' for 'Upper Right').
    :type ur_lat:
        :py:class:`float`

    :param ur_lon:
        The longitude of the upper right corner of the region ('ur' for 'Upper Right').
    :type ur_lon:
        :py:class:`float`

    :param bin_dir:
        The directory where the executable ``aot_loader`` can be found.
    :type bin_dir:
        :py:class:`str`

    """
    ul_lon = ll_lon
    lr_lat = ll_lat
    lr_lon = ur_lon
    ul_lat = ur_lat

    s = re.search('\.(pix|cmp)', atsr_file)
    assert s, 'Invalid AATSR file name: ' + atsr_file
    filetype = s.group(1)

    if not os.path.exists(atsr_file):
        logger.debug('Aerosol %s file (%s) not found', filetype, atsr_file)
        return 0.0

    command = (os.path.join(bin_dir, 'aot_loader') +
        ' --' + filetype + ' ' + atsr_file +
        ' --west ' + str(min([ll_lon, ul_lon])) +
        ' --east ' + str(max([lr_lon, ur_lon])) +
        ' --south ' + str(min([ll_lat, lr_lat])) +
        ' --north ' + str(max([ul_lat, ur_lat])) +
        ' --date ' + dt.strftime('%Y-%m-%d') +
        ' --t ' + dt.strftime('%H:%M:%S'))

    result = execute(command)

    import pprint
    pprint.pprint(result) 

    def parse_aatsr_result_str(text, result_regex=r'AOT AATSR value:(\s+)(.+)$'):
        m = re.search(result_regex, text, re.MULTILINE)
        if m and m.group(2):
            return float(m.group(2).rstrip())
        return 0.0

    # TODO ...print something helpful when we get a non-zero return code...

    assert result['returncode'] == 0, 'get_aerosol_aatsr return code is non-zero'
    return parse_aatsr_result_str(result['stdout'])





def get_aerosol_from_aeronet_station(aeronet_file, dt):
    """
    Extract an aerosol visibility from a `AERONET <http://gcmd.nasa.gov/records/GCMD_AERONET_NASA.html>`_ file.

    :param aeronet_file:
        The file to extract the value from.
    :type aeronet_file:
        :py:class:`str`

    :param dt:
        The date and time to extract the value for.
    :type dt:
        :py:class:`datetime.datetime`

    """
    fpIn = open(aeronet_file)
    line = fpIn.readline()
    # print "Line:", line

    # Skip lines to the first one starting with " year"
    while line and line.find(" year") != 0:
        line = fpIn.readline()

    # Read through table of visibility values and find the closest one
    line = fpIn.readline()
    while line:
        temp = line.split()
        this_datetime = datetime()
        this_datetime.year = temp[0]
        this_datetime.month = temp[1]
        this_datetime.day = temp[2]
        this_datetime.hour = temp[3]
        this_datetime.minute = temp[4]
        this_datetime.second = temp[6]

        this_aero = float(temp[6])

        # print "thisDate",  thisDate
        if this_datetime > dt:
            # Found the first one AFTER
            break
        prev_datetime = this_datetime
        prev_aero = this_aero

        line = fpIn.readline()

    # logger.debug("Time interval: %s to %s" % (prevDate, thisDate))
    # Which is closer? This one or previous
    if abs(prev_datetime - dt) < abs(this_datetime - dt):
        logger.debug("Closest Time: %s {%f)" % (prev_datetime, prev_aero))
        return prev_aero
    else:
        logger.debug("Closest Time: %s (Aerosol: %f)" % (this_datetime, this_aero))
        return this_aero

    return this_aero





def get_closest_aeronet_station(stationFile, sc_lat, sc_lon):
    """
    Determine closest `AERONET <http://gcmd.nasa.gov/records/GCMD_AERONET_NASA.html>`_ ground station
    to a given location.

    :param sc_lat:
        The latitude of the center of the region ('sc' for 'Scene Center').
    :type sc_lat:
        :py:class:`float`

    :param sc_lon:
        The longitude of the center of the region ('sc' for 'Scene Center').
    :type sc_lon:
        :py:class:`float`

    """

    logger.debug("Reading Aerosol Stations file %s", stationFile)
    fpIn = open(stationFile)

    stationList = []

    line = fpIn.readline()
    while line:
        temp = line.split()
        if len(temp) < 3:
            break
        Long = float(temp[0])
        Lat = float(temp[1])
        Name = temp[2]

        # determine the stations distance
        distance = (sc_lat - Lat) * (sc_lat - Lat) + \
                   (sc_lon - Long) * (sc_lon - Long)

        stationList.append([distance, Name])
        line = fpIn.readline()

    if len(stationList) < 1:
        logger.error("Error reading station file %s" % stationFile )
        return ""

    # Need to determine the closest station by distance:
    # Sort them and return the first one
    stationList.sort()
    return stationList[0][1]





def get_aerosol_aeronet(aerosol_dir, dt, sc_lat, sc_lon):
    """
    Extract an `AERONET <http://gcmd.nasa.gov/records/GCMD_AERONET_NASA.html>`_ visibility value for a given
    date, time and location.

    :param dt:
        The date and time to extract the value for.
    :type dt:
        :py:class:`datetime.datetime`

    :param aerosol_dir:
        The directory containing the aerosol data files.
    :type aerosol_dir:
        :py:class:`str`

    :param sc_lat:
        The latitude of the center of the region ('sc' for 'Scene Center').
    :type sc_lat:
        :py:class:`float`

    :param sc_lon:
        The longitude of the center of the region ('sc' for 'Scene Center').
    :type sc_lon:
        :py:class:`float`

    """
    aeronet_stns_file = os.path.join(aerosol_dir, "stations.txt")
    aeronet_stn_name = get_closest_aeronet_station(aeronet_stns_file, sc_lat, sc_lon)
    aeronet_data_file = os.path.join(aerosol_dir, "visibility_%s.txt" % aeronet_stn_name)

    logger.debug("Closest Aerosol station: %s", aeronet_stn_name)
    logger.debug("Aerosol data file: %s", aeronet_data_file)

    if not os.path.exists(aeronet_data_file):
        print 'WARNING: aeronet data file (%s) not found' % aeronet_data_file
        return 0.0

    return get_aerosol_from_aeronet_station(aeronet_data_file, dt)
