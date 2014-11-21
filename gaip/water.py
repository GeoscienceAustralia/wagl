import rasterio
import gaip
import os


def get_water_vapour_data(water_vapour_path, lonlat, date_time,
                          scale_factor=0.1):
    """
    Retrieve the water vapour value given a longitude, latitude and
    date time.

    :param water_vapour_path:
        A string containing the full file path to a directory
        containing the water vapour ancillary data. Water vapour
        files should be of the form pr_wtr.eatm.{year}.tif with year
        indicating each years, i.e. 1999, 2000, 2001, worth of data.
        Each water vapour file contains a separate band for each
        day-of-year and hour.  Each band represents a continental
        scale 2D image.

    :param lonlat:
        The longitude and latitude are a 2 element floating point
        tuple (lon, lat) for the position of interest.

    :param date_time:
        A python datetime object containing the year, month, day,
        hour and second to be used in determining the correct band
        in the water vapour ancillary files.

    :param scale_factor:
        A scale factor to apply to the retrieve value. The water
        vapour ancillary files should be in kg/m^2, and the default
        action is to return the values in g/cm^2 for use within
        MODTRAN.
        Default is 0.1

    :return:
        A dictionary with 3 keys of the form:
        'data_source' -> Type of ancillary data.
        'data_file' -> File system path to the file of interest.
        'value' -> The retrieved data value.
    """

    year = date_time.strftime('%Y')
    filename = "pr_wtr.eatm.{year}.tif".format(year=year)
    datafile = os.path.join(water_vapour_path, filename)


    # calculate the water vapour band number based on the datetime

    doy = date_time.timetuple().tm_yday
    hour = date_time.timetuple().tm_hour
    band_number = (int(doy) - 1) * 4 + int((hour + 3) / 6)

    # Check for boundary condition: 1 Jan, 0-3 hours
    if band_number == 0 and doy == 1:
        band_number = 1

    # Get the number of bands
    with rasterio.open(datafile) as src:
        n_bands = src.count

    # Enable NBAR Near Real Time (NRT) processing
    if band_number > (n_bands + 1):
        rasterdoy = (((n_bands) - (int((hour + 3) / 6))) / 4) + 1
        if (doy - rasterdoy) < 7:
            band_idx = (int(rasterdoy) - 1) * 4 + int((hour + 3) / 6)

    try:
        value = gaip.get_pixel(datafile, lonlat, band=band_idx)
    except IndexError:
        msg = "Invalid water vapour band number: {band}".format(band=band_idx)
        raise IndexError(msg)

    value = value * scale_factor

    water_vapour_data = {
        'data_source': 'Water Vapour',
        'data_file': datafile,
        'value': value
    }

    return water_vapour_data
