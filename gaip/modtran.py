"""
MODTRAN drivers
---------------

"""
import os
from os.path import join as pjoin, exists, abspath, dirname
import subprocess

import numpy
import pandas
import rasterio
import gaip
from gaip import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO
from gaip import MIDLAT_SUMMER_TRANSMITTANCE, TROPICAL_TRANSMITTANCE

BIN_DIR = abspath(pjoin(dirname(__file__), '..', 'bin'))


def create_modtran_dirs(coords, albedos, modtran_root, modtran_exe_root,
                        workpath_format, input_format):
    """Create all modtran subdirectories. and input files."""

    if not exists(modtran_root):
        os.makedirs(modtran_root)

    data_dir = pjoin(modtran_exe_root, 'DATA')
    if not exists(data_dir):
        raise OSError('Cannot find MODTRAN')

    for coord in coords:
        for albedo in albedos:
            modtran_work = workpath_format.format(coord=coord, albedo=albedo)
            modtran_work = pjoin(modtran_root, modtran_work)
            mod5root_in = input_format.format(coord=coord, albedo=albedo)
            mod5root_in = pjoin(modtran_root, mod5root_in)

            if not exists(modtran_work):
                os.makedirs(modtran_work)

            with open(mod5root_in, 'w') as outfile:
                outfile.write(coord + '_alb_' + albedo + '\n')

            symlink_dir = pjoin(modtran_work, 'DATA')
            if exists(symlink_dir):
                os.unlink(symlink_dir)

            os.symlink(data_dir, symlink_dir)


def create_satellite_filter_file(acquisitions, satfilter_path, target):
    """Generate satellite filter input file."""
    refbands = [a for a in acquisitions if a.band_type == gaip.REF]
    filterfile = acquisitions[0].spectral_filter_file
    filterpath = os.path.join(satfilter_path, filterfile)

    with open(target, 'w') as outfile:
        outfile.write("%i\n" % len(refbands))
        outfile.write("%s\n" % filterpath)

    return target


# TODO: once validated, this can function can be deprecated
def write_modtran_input(acquisitions, modtran_input_file, ozone, vapour,
                        aerosol, elevation):
    """Generate modtran input file."""
    acq = acquisitions[0]
    filter_file = acq.spectral_filter_file
    cdate = acq.scene_centre_date
    altitude = acq.altitude / 1000.0  # in km
    dechour = acq.decimal_hour

    with open(modtran_input_file, 'w') as outfile:
        outfile.write("%f\n" % ozone)
        outfile.write("%f\n" % vapour)
        outfile.write("DATA/%s\n" % filter_file)
        outfile.write("-%f\n" % aerosol)
        outfile.write("%f\n" % elevation)
        outfile.write("Annotation, %s\n" % cdate.strftime('%Y-%m-%d'))
        outfile.write("%d\n" % altitude)
        outfile.write("%d\n" % int(cdate.strftime('%j')))
        outfile.write("%f\n" % dechour)


# TODO: once validated, this can function can be deprecated
# as we can write direct to the tp5 template
def write_modtran_inputs(acquisition, coordinator, view_fname, azi_fname,
                         lat_fname, lon_fname,  ozone, vapour, aerosol,
                         elevation, coords, albedos, out_fname_fmt):
    filter_file = acquisition.spectral_filter_file
    cdate = acquisition.scene_centre_date
    altitude = acquisition.altitude / 1000.0  # in km
    dechour = acquisition.decimal_hour
    coord = pandas.read_csv(coordinator, header=None, sep=r'\s+\s+',
                            engine='python', names=['row', 'col'])

    with rasterio.open(view_fname) as view_ds,\
        rasterio.open(azi_fname) as azi_ds,\
        rasterio.open(lat_fname) as lat_ds,\
        rasterio.open(lon_fname) as lon_ds:

        npoints = len(coords)
        view = numpy.zeros(npoints, dtype='float32')
        azi = numpy.zeros(npoints, dtype='float32')
        lat = numpy.zeros(npoints, dtype='float64')
        lon = numpy.zeros(npoints, dtype='float64')

        for i in range(1, npoints + 1):
            yidx = coord['row'][i]
            xidx = coord['col'][i]
            idx = ((yidx -1, yidx), (xidx -1, xidx))
            view[i-1] = view_ds.read(1, window=idx)[0, 0]
            azi[i-1] = azi_ds.read(1, window=idx)[0, 0]
            lat[i-1] = lat_ds.read(1, window=idx)[0, 0]
            lon[i-1] = lon_ds.read(1, window=idx)[0, 0]

    view_cor = 180 - view
    azi_cor = azi + 180
    rlon = 360 - lon
    
    # check if in western hemisphere
    wh = rlon >= 360
    rlon[wh] -= 360
    
    wh = (180 - view_cor) < 0.1
    view_cor[wh] = 180
    azi_cor[wh] = 0
    
    wh = azi_cor > 360
    azi_cor[wh] -= 360

    for i, p in enumerate(coords):
        for alb in albedos:
            out_fname = out_fname_fmt.format(coord=p, albedo=alb)
            with open(out_fname, 'w') as src:
                src.write("{:.8f}\n".format(float(alb)))
                src.write("{:.14f}\n".format(ozone))
                src.write("{:.14f}\n".format(vapour))
                src.write("DATA/{}\n".format(filter_file))
                src.write("-{:.14f}\n".format(aerosol))
                src.write("{:.14f}\n".format(elevation))
                src.write("Annotation, {}\n".format(cdate.strftime('%Y-%m-%d')))
                src.write("{:.14f}\n".format(altitude))
                src.write("{:f}\n".format(view_cor[i]))
                src.write("{:d}\n".format(int(cdate.strftime('%j'))))
                src.write("{:.14f}\n".format(lat[i]))
                src.write("{:.14f}\n".format(rlon[i]))
                src.write("{:.14f}\n".format(dechour))
                src.write("{:f}\n".format(azi_cor[i]))


def write_modis_brdf_files(acquisitions, fname_format, brdf_data,
                           solar_irrad_data, solar_dist_data):
    """Generate brdf input file."""
    ref_acqs = [a for a in acquisitions if a.band_type == gaip.REF]

    for acq in ref_acqs:
        band = acq.band_num
        modis_brdf_filename = fname_format.format(band_num=band)
        with open(modis_brdf_filename, 'w') as outfile:
            msg = "{iso} {vol} {geo}\n"
            msg = msg.format(iso=brdf_data[(band, 'iso')]['value'],
                             vol=brdf_data[(band, 'vol')]['value'],
                             geo=brdf_data[(band, 'geo')]['value'])
            outfile.write(msg)

            msg = "{bias} {gain} {irrad} {dist}\n"
            msg = msg.format(bias=acq.bias, gain=acq.gain,
                             irrad=solar_irrad_data[band],
                             dist=solar_dist_data)
            outfile.write(msg)


# TODO: remove
def generate_modtran_inputs(modtran_input, coordinator, sat_view_zenith,
                            sat_azimuth, lon_grid, lat_grid, coords, albedos,
                            fname_format, workdir):
    """Generate MODTRAN input files."""
    cmd = pjoin(BIN_DIR, 'generate_modtran_input')

    args = [cmd, modtran_input, coordinator, sat_view_zenith, sat_azimuth,
            lat_grid, lon_grid]

    targets = []
    for coord in coords:
        for albedo in albedos:
            target = fname_format.format(coord=coord, albedo=albedo)
            targets.append(pjoin(workdir, target))

    args.extend(targets)

    subprocess.check_call(args)

    return targets


def write_tp5(acquisition, coordinator, view_fname, azi_fname,
              lat_fname, lon_fname, ozone, vapour, aerosol, elevation,
              coords, albedos, out_fname_fmt):
    """Writes the tp5 files for the albedo (0, 1) and transmittance (t)."""
    geobox = acquisition.gridded_geo_box()
    filter_file = acquisition.spectral_filter_file
    cdate = acquisition.scene_centre_date
    doy = int(cdate.strftime('%j'))
    altitude = acquisition.altitude / 1000.0  # in km
    dechour = acquisition.decimal_hour
    coord = pandas.read_csv(coordinator, header=None, sep=r'\s+\s+',
                            engine='python', names=['row', 'col'])

    with rasterio.open(view_fname) as view_ds,\
        rasterio.open(azi_fname) as azi_ds,\
        rasterio.open(lat_fname) as lat_ds,\
        rasterio.open(lon_fname) as lon_ds:

        npoints = len(coords)
        view = numpy.zeros(npoints, dtype='float32')
        azi = numpy.zeros(npoints, dtype='float32')
        lat = numpy.zeros(npoints, dtype='float64')
        lon = numpy.zeros(npoints, dtype='float64')

        for i in range(1, npoints + 1):
            yidx = coord['row'][i]
            xidx = coord['col'][i]
            idx = ((yidx -1, yidx), (xidx -1, xidx))
            view[i-1] = view_ds.read(1, window=idx)[0, 0]
            azi[i-1] = azi_ds.read(1, window=idx)[0, 0]
            lat[i-1] = lat_ds.read(1, window=idx)[0, 0]
            lon[i-1] = lon_ds.read(1, window=idx)[0, 0]

    view_cor = 180 - view
    azi_cor = azi + 180
    rlon = 360 - lon
    
    # check if in western hemisphere
    wh = rlon >= 360
    rlon[wh] -= 360
    
    wh = (180 - view_cor) < 0.1
    view_cor[wh] = 180
    azi_cor[wh] = 0
    
    wh = azi_cor > 360
    azi_cor[wh] -= 360

    # get the modtran profiles to use based on the centre latitude 
    centre_lon, centre_lat = geobox.centre_lonlat
    if centre_lat < -23.0:
        albedo_profile = MIDLAT_SUMMER_ALBEDO
        trans_profile = MIDLAT_SUMMER_TRANSMITTANCE
    else:
        albedo_profile = TROPICAL_ALBEDO
        trans_profile = TROPICAL_TRANSMITTANCE

    # write the tp5 files required for input into MODTRAN
    for i, p in enumerate(coords):
        for alb in albedos:
            out_fname = out_fname_fmt.format(coord=p, albedo=alb)
            if alb == 't':
                data = trans_profile.format(albedo=0.0,
                                            water=vapour,
                                            ozone=ozone,
                                            filter_function=filter_file,
                                            visibility=-aerosol,
                                            elevation=elevation,
                                            sat_height=altitude,
                                            sat_view=view_cor[i],
                                            doy=doy,
                                            sat_view_offset=180.0-view_cor[i])
            else:
                data = albedo_profile.format(albedo=float(alb),
                                             water=vapour,
                                             ozone=ozone,
                                             filter_function=filter_file,
                                             visibility=-aerosol,
                                             elevation=elevation,
                                             sat_height=altitude,
                                             sat_view=view_cor[i],
                                             doy=doy,
                                             lat=lat[i],
                                             lon=rlon[i],
                                             time=dechour,
                                             sat_azimuth=azi_cor[i])
            with open(out_fname, 'w') as src:
                src.write(data)


def reformat_as_tp5(coords, albedos, profile, input_format, output_format,
                    workdir, cmd=pjoin(BIN_DIR, 'reformat_tp5_albedo')):
    """Reformat the MODTRAN input files in `tp5` format."""

    targets = []
    for coord in coords:
        for albedo in albedos:
            src = input_format.format(coord=coord, albedo=albedo)
            dst = output_format.format(coord=coord, albedo=albedo)
            targets.append(pjoin(workdir, dst))

            args = [cmd, pjoin(workdir, src), profile, pjoin(workdir, dst)]

            subprocess.check_call(args)

    return targets


def reformat_as_tp5_trans(coords, albedos, profile, input_format,
                          output_format, workdir):
    """Reformat the MODTRAN input files in `tp5` format in the trans case."""
    cmd = pjoin(BIN_DIR, 'reformat_tp5_transmittance')
    return reformat_as_tp5(coords, albedos, profile, input_format,
                           output_format, workdir, cmd)


def run_modtran(modtran_exe, workpath):
    """Run MODTRAN."""
    subprocess.check_call([modtran_exe], cwd=workpath)


def extract_flux(coords, albedos, input_format, output_format, satfilter):
    """Extract the flux data."""
    cmd = pjoin(BIN_DIR, 'read_flux_albedo')

    for coord in coords:
        for albedo in albedos:
            src = input_format.format(coord=coord, albedo=albedo)
            dst = output_format.format(coord=coord, albedo=albedo)
            args = [cmd, src, satfilter, dst]

            subprocess.check_call(args)


def extract_flux_trans(coords, input_format, output_format, satfilter):
    """Extract the flux data in the transmissive case."""
    cmd = pjoin(BIN_DIR, 'read_flux_transmittance')

    for coord in coords:
        src = input_format.format(coord=coord)
        dst = output_format.format(coord=coord)
        args = [cmd, src, satfilter, dst]

        subprocess.check_call(args)


def calc_coefficients(coords, chn_input_fmt, dir_input_fmt,
                      output_fmt, satfilter, cwd):
    """Calculate the coefficients from the MODTRAN output."""

    cmd = pjoin(BIN_DIR, 'calculate_coefficients')

    for coord in coords:
        args = [cmd, satfilter,
                pjoin(cwd, chn_input_fmt.format(coord=coord, albedo=0)),
                pjoin(cwd, chn_input_fmt.format(coord=coord, albedo=1)),
                pjoin(cwd, dir_input_fmt.format(coord=coord, albedo=0)),
                pjoin(cwd, dir_input_fmt.format(coord=coord, albedo=1)),
                pjoin(cwd, dir_input_fmt.format(coord=coord, albedo='t')),
                pjoin(cwd, output_fmt.format(coord=coord))]

        subprocess.check_call(args, cwd=cwd)


def reformat_atmo_params(acqs, coords, satfilter, factors, input_fmt,
                         output_fmt, workpath):
    """Reformat atmospheric parameters."""

    cmd = pjoin(BIN_DIR, 'reformat_modtran_output')

    bands = [str(a.band_num) for a in acqs]

    args = [cmd, satfilter]
    for coord in coords:
        args.append(input_fmt.format(coord=coord))

    for band in bands:
        for factor in factors:
            args.append(output_fmt.format(factor=factor, band=band))

    subprocess.check_call(args, cwd=workpath)


def bilinear_interpolate(acqs, factors, coordinator, boxline, centreline,
                         input_fmt, output_fmt, workpath):
    """Perform bilinear interpolation."""

    bands = [a.band_num for a in acqs]
    geobox = gaip.gridded_geo_box(acqs[0])
    cols, rows = geobox.get_shape_xy()

    # dataframes for the coords, scene centreline, boxline
    coords = pandas.read_csv(coordinator, header=None, sep=r'\s+\s+',
                             engine='python', skiprows=1, names=['row', 'col'])
    cent = pandas.read_csv(centreline, skiprows=2, header=None, sep=r'\s+\s+',
                           engine='python',
                           names=['line', 'centre', 'npoints', 'lat', 'lon'])
    box = pandas.read_csv(boxline, header=None, sep=r'\s+\s+', engine='python',
                          names=['line', 'cstart', 'cend'])

    coord = numpy.zeros((9, 2), dtype='int')
    coord[:, 0] = coords.row.values
    coord[:, 1] = coords.col.values
    centre = cent.centre.values
    start = box.cstart.values
    end = box.cend.values

    # Initialise the dict to store the locations of the bilinear outputs
    bilinear_outputs = {}

    for band in bands:
        for factor in factors:
            fname = output_fmt.format(factor=factor, band=band)
            fname = pjoin(workpath, fname)
            atmospheric_fname = input_fmt.format(factor=factor, band=band)
            bilinear_outputs[(band, factor)] = fname

            # atmospheric paramaters
            atmos = pandas.read_csv(atmospheric_fname, header=None,
                                    sep=r'\s+\s+', engine='python',
                                    names=['s1', 's2', 's3', 's4'])

            # get the individual atmospheric components
            s1 = atmos.s1.values
            s2 = atmos.s2.values
            s3 = atmos.s3.values
            s4 = atmos.s4.values

            res = numpy.zeros((rows, cols), dtype='float32')
            gaip.bilinear(cols, rows, coord, s1, s2, s3, s4, start, end,
                          centre, res.transpose())

            # Output the result to disk
            gaip.write_img(res, fname, geobox=geobox, nodata=-999)

    return bilinear_outputs
