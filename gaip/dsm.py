#!/usr/bin/env python

from gaip import Buffers
from gaip import filter_dsm

def get_dsm(acquisition, national_dsm, buffer, fname_subset, fname_smoothed):
    """
    Given an acquisition and a national Digitial Surface Model,
    extract a subset from the DSM based on the acquisition extents
    plus an x & y buffer. The subset is then smoothed with a 3x3
    gaussian filter.
    A square buffer is applied to the extents.
    """

    # Use the 1st acquisition to setup the geobox
    geobox = acquisition.gridded_geo_box()

    # Define Top, Bottom, Left, Right pixel buffers
    pixel_buf = Buffers(buffer)

    # Get the dimensions and geobox of the new image
    dem_cols = geobox.getShapeXY()[0] + pixel_buf.left + pixel_buf.right
    dem_rows = geobox.getShapeXY()[1] + pixel_buf.top + pixel_buf.bottom
    dem_shape = (dem_rows, dem_cols)
    dem_origin = geobox.convert_coordinates((0 - pixel_buf.left,
        0 - pixel_buf.top))
    dem_geobox = GriddedGeoBox(dem_shape, origin=dem_origin,
        pixelsize=geobox.pixelsize, crs=geobox.crs.ExportToWkt())

    # Retrive the DSM data
    dsm_data = reprojectFile2Array(national_dsm, dst_geobox=dem_geobox,
        resampling=RESAMPLING.bilinear)

    # Output the reprojected result
    write_img(dsm_data, fname_subset, geobox=dem_geobox)

    # Smooth the DSM
    dsm_data = filter_dsm(dsm_data)

    # Output the smoothed DSM
    write_img(dsm_data, fname_smoothed, geobox=dem_geobox)
