import rasterio
import gaip
import os


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
        value = gaip.get_pixel(datafile, geobox.centre_lonlat, band=band_idx)
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
