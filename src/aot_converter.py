#!/usr/bin/env python

from datetime import datetime as dt
import glob
from posixpath import join as pjoin
from posixpath import basename
from posixpath import splitext
import numpy
import pandas


def read_pix(filename):
    """
    The pix files are sparse 3D arrays.
    Will store as a Table and remove invalid data.
    """
    src = open(filename, 'rb')
    recs = numpy.fromfile(src, dtype='int32', count=3)
    xgrid = numpy.fromfile(src, dtype='float32', count=recs[0])
    ygrid = numpy.fromfile(src, dtype='float32', count=recs[1])
    idxlon = numpy.fromfile(src, dtype='int16', count=recs[2])
    idxlat = numpy.fromfile(src, dtype='int16', count=recs[2])
    date = numpy.fromfile(src, dtype='int16',
                          count=recs[2] * 3).reshape(3, recs[2])
    time = numpy.fromfile(src, dtype='int16',
                          count=recs[2] * 3).reshape(3, recs[2])
    lat = numpy.fromfile(src, dtype='float32', count=recs[2])
    lon = numpy.fromfile(src, dtype='float32', count=recs[2])
    aot = numpy.fromfile(src, dtype='float32', count=recs[2])
    src.close()

    obs_lon = xgrid[idxlon]
    obs_lat = ygrid[idxlat]
    timestamps = []
    for i in range(recs[2]):
        timestamps.append(dt(date[0, i], date[1, i], date[2, i],
                             time[0, i], time[1, i], time[2, i]))

    df = pandas.DataFrame({'timestamp': timestamps,
                           'lon': obs_lon,
                           'lat': obs_lat,
                           'aerosol': aot})

    # throw away bad data
    wh = (df['aerosol'] > 0.0) & (df['aerosol'] <= 1.0)
    df = df[wh]
    df.reset_index(inplace=True)

    # get the minimum bounding box
    ul = (df['lon'].min(), df['lat'].max())
    ur = (df['lon'].max(), df['lat'].max())
    lr = (df['lon'].max(), df['lat'].min())
    ll = (df['lon'].min(), df['lat'].min())
    extents = [ul, ur, lr, ll]

    return df, extents


def read_cmp(filename):
    """
    The cmp data is a 2D grid, but the pixel sizes are
    not constant (according to the lon and lat arrays.
    Will store as a Table and remove invalid data.
    """
    src = open(filename, 'rb')
    nx = numpy.fromfile(src, dtype='int32', count=1)[0]
    ny = numpy.fromfile(src, dtype='int32', count=1)[0]
    lon = numpy.fromfile(src, dtype='float32', count=nx)
    lat = numpy.fromfile(src, dtype='float32', count=ny)
    aot = numpy.fromfile(src, dtype='float32', count=nx*ny)

    obs_lon = []
    obs_lat = []
    for j in range(ny):
        for i in range(nx):
            obs_lon.append(lon[i])
            obs_lat.append(lat[j])

    df = pandas.DataFrame({'lon': obs_lon, 'lat': obs_lat, 'aerosol': aot})

    # throw away bad data
    wh = (df['aerosol'] > 0.0) & (df['aerosol'] <= 1.0)
    df = df[wh]
    df.reset_index(inplace=True)

    # get the minimum bounding box
    ul = (df['lon'].min(), df['lat'].max())
    ur = (df['lon'].max(), df['lat'].max())
    lr = (df['lon'].max(), df['lat'].min())
    ll = (df['lon'].min(), df['lat'].min())
    extents = [ul, ur, lr, ll]

    return df, extents


def main(aerosol_path, output_filename):
    # define a case switch
    func = {'pix': read_pix, 'cmp': read_cmp}

    # create the output file
    hstore = pandas.HDFStore(output_filename, 'w', complevel=4,
                             complib='bzip2')

    pattern = ['*.pix', '*.cmp']
    for p in pattern:
        search = pjoin(aerosol_path, p)
        files = glob.glob(search)
        for fname in files:
            pth, ext = splitext(fname)
            ext = ext.split(".")[-1]
            grp_name = basename(pth)
            out_path = pjoin(ext, grp_name)

            # read/write
            df, extents = func[ext](fname)
            hstore.append(out_path, df)
            hstore.flush()

            # assign extents attribute
            nd = hstore.get_node(out_path)
            nd._v_attrs.extents = extents
            hstore.flush()

    hstore.close()


if __name__ == '__main__':

    aerosol_path = '/g/data/v10/eoancillarydata/aerosol/AATSR/2.0/'
    out_fname = '/g/data/v10/eoancillarydata/aerosol/AATSR/2.0/aerosol.h5'
    main(aerosol_path, out_fname)
