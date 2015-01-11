"""
MODTRAN drivers
---------------

"""
import os
import gaip
import subprocess

from os.path import join as pjoin, exists, abspath, dirname

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


def write_modtran_input(acquisitions, modtran_input_file, ozone, vapour,
                        aerosol, elevation):
    """Generate modtran input file."""
    acq = acquisitions[0]
    geobox = acq.gridded_geo_box()
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
                             dist=solar_dist_data['distance'])
            outfile.write(msg)


def run_read_modtrancor_ortho(centreline, sat_view_zenith, coordinator,
                              boxline, cwd):
    """run read_modtrancor_ortho executable."""
    cmd = pjoin(BIN_DIR, 'read_modtrancor_ortho')

    args = [cmd, centreline, sat_view_zenith, coordinator, boxline]

    subprocess.check_call(args)


def generate_modtran_inputs(modtran_input, coordinator, sat_view_zenith,
                            sat_azimuth, lon_grid, lat_grid, coords, albedos,
                            fname_format, workdir):
    """Generate MODTRAN input files."""
    cmd = pjoin(BIN_DIR, 'input_modtran_ortho')

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


def reformat_as_tp5(coords, albedos, profile, input_format, output_format,
                    workdir, cmd = pjoin(BIN_DIR, 'refort_tp5_ga')):
    """Reformat the MODTRAN input files in `tp5` format."""

    targets = []
    for coord in coords:
        for albedo in albedos:
            src = input_format.format(coord=coord, albedo=albedo)
            dst = output_format.format(coord=coord, albedo=albedo)
            targets.append(pjoin(workdir, dst))

            args = [cmd, pjoin(workdir, src), profile, pjoin(workdir, dst)]

            print ' '.join(args)

            subprocess.check_call(args)

    return targets


def reformat_as_tp5_trans(coords, albedos, profile, input_format,
                          output_format, workdir):
    """Reformat the MODTRAN input files in `tp5` format in the trans case."""
    cmd = pjoin(BIN_DIR, 'refort_tp5_ga_trans')
    return reformat_as_tp5(coords, albedos, profile, input_format,
                           output_format, workdir, cmd)


def run_modtran(modtran_exe, workpath):
    """Run MODTRAN."""
    subprocess.check_call([modtran_exe], cwd=workpath)


def extract_flux(coords, albedos, input_format, output_format, satfilter):
    cmd = pjoin(BIN_DIR, 'read_flx_ga')

    for coord in coords:
        for albedo in albedos:
            src = input_format.format(coord=coord, albedo=albedo)
            dst = output_format.format(coord=coord, albedo=albedo)
            args = [cmd, src, satfilter, dst]

            subprocess.check_call(args)


def extract_flux_trans(coords, input_format, output_format, satfilter):
    """Extract the flux data in the transmissive case."""
    cmd = pjoin(BIN_DIR, 'read_flx_ga_trans')

    for coord in coords:
        src = input_format.format(coord=coord)
        dst = output_format.format(coord=coord)
        args = [cmd, src, satfilter, dst]

        subprocess.check_call(args)


def calc_coefficients(coords, chn_input_fmt, dir_input_fmt,
                      output_fmt, satfilter, cwd):
    """Calculate the coefficients from the MODTRAN output."""

    cmd = pjoin(BIN_DIR, 'coefficient')

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

    cmd = pjoin(BIN_DIR, 'read_modtran')

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

    cmd = pjoin(BIN_DIR, 'binear_ortho')

    bands = [a.band_num for a in acqs]

    # Initialise the dict to store the locations of the bilinear outputs
    bilinear_outputs = {}

    for band in bands:
        for factor in factors:
            fname = output_fmt.format(factor=factor, band=band)
            fname = pjoin(workpath, fname)
            bilinear_outputs[(band, factor)] = fname
            args = [cmd, coordinator,
                    input_fmt.format(factor=factor, band=band),
                    boxline, centreline,
                    fname]

            subprocess.check_call(args, cwd=workpath)

    return bilinear_outputs
