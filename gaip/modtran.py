"""
MODTRAN drivers
"""
import os
import subprocess

from os.path import join as pjoin, exists


class modtran_config(object):
    FULL_COORD_LIST = ['TL', 'TM', 'TR', 'ML', 'MM', 'MR', 'BL', 'BM', 'BR']
    FULL_ALBEDO_LIST = ['0', '1', 't']
    FULL_FACTOR_LIST = ['fv', 'fs', 'a', 'b', 's', 'dir', 'dif', 'ts']

FULL_COORD_LIST = ['TL', 'TM', 'TR', 'ML', 'MM', 'MR', 'BL', 'BM', 'BR']
FULL_ALBEDO_LIST = ['0', '1', 't']
FULL_FACTOR_LIST = ['fv', 'fs', 'a', 'b', 's', 'dir', 'dif', 'ts']


def create_modtran_dirs(workpath, modtran_root):
    """Create all modtran subdirectories."""
    paths = []
    for coord in FULL_COORD_LIST:
        for albedo in [albedo.lower() for albedo in FULL_ALBEDO_LIST]:
            modtran_work = pjoin(workpath, 'mod', coord, 'alb_' + albedo)
            mod5root_in = pjoin(modtran_work, 'mod5root.in')
            data_dir = pjoin(modtran_root, 'DATA')
            symlink_dir = pjoin(modtran_work, 'DATA')

            if not exists(modtran_workpath):
                os.makedirs(modtran_workpath)
                paths.append(modtran_workpath)

            if exists(symlink_dir):
                os.unlink(symlink_dir)
            os.symlink(data_dir, symlink_dir)
            paths.append(symlink_dir)

            with open(mod5root_in, 'w') as outfile:
                outfile.write(coord + '_alb_' + albedo + '\n')

    symlink_dir = pjoin(workpath, 'mod', 'DATA')
    if exists(symlink_dir):
        os.unlink(symlink_dir)
    os.symlink(data_dir, symlink_dir)
    paths.append(symlink_dir)

    return paths


def create_header_angle_file(max_view_angle):
    """Create header angle file.

    Arguments:
        ha_file: header angle file path
        max_view_angle: view angle cutoff (degrees)
    """
    # TODO Change of format
    # Fuqin uses a different format in her calculate angles stuff:
    # year month day hours
    # lines columns
    # centre_lat centre_lon
    # satellite_semi_mjr_radius orbital_inc angular_velocity
    # Also Fuqin says that it isn't actually used by MODTRAN, so maybe this
    # should be moved to a new area, potentially the calculate_angles file
    # The file itself is not used by the F2Py version (unlike the original
    # Fortran version), but can be used as an aid for validation purposes.
    # JS 20141205
    # ha_file, l1t_input_dataset, max_view_angle=None):
    
    # TODO remove hardwiring for LANDSAT orbital parameters

    max_view_angle = max_view_angle or CONFIG.VIEW_ANGLE_MAX

    # Centre of UL pixel (degrees)
    UL = l1t_input_dataset.lonlats['UL']

    # Resolution -- hardwired per original NBAR code
    # TODO calculate based on spatial res metadata
    resolution_deg = satellite.NOMINAL_PIXEL_DEGREES

    ha_text = '\n'.join([
        '%d %s' % (l1t_input_dataset.DOY, l1t_input_dataset.decimal_day),
        '%d %d' % l1t_input_dataset.shape,
        '%f' % resolution_deg,
        '%10.6f %10.6f' % (UL[1], UL[0]),
        '%f %f' % (l1t_input_dataset.lonlats['CENTRE'][
                   1], l1t_input_dataset.lonlats['CENTRE'][0]),
        '%f %f %f' % (satellite.SEMI_MAJOR_AXIS,
                      # Convert to degrees
                      math.degrees(satellite.INCLINATION),
                      satellite.OMEGA),  # Leave in radians
        '%f\n' % max_view_angle,
    ])

    try:
        fd = open(ha_file, 'w')
        fd.write(ha_text)
    finally:
        fd.close()


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
     ul_lon, ul_lat = geobox.ul_lonlat
     pixel_degrees = acq.nominal_pixel_degrees
     filter_file = acq.spectral_filter_file
     cdate = acq.scene_centre_date
     altitude = acq.altitude/1000.0  # in km
     dechour = acq.decimal_hour
     
     with open(modtran_input_file, 'w') as outfile:
        outfile.write("%f %f\n" % (ul_lat, ul_lon))
        outfile.write("%f\n" % pixel_degrees)
        outfile.write("%f\n" % ozone)
        outfile.write("%f\n" % vapour)
        outfile.write("DATA/%s\n" % filter_file)
        outfile.write("-%f\n" % aerosol)
        outfile.write("%f\n" % elevation)
        outfile.write("Annotation, %s\n" % cdate.strftime('%Y-%m-%d'))
        outfile.write("%d\n" % altitude)
        outfile.write("%d\n" % cdate.strftime('%j'))
        outfile.write("%f\n" % dechour)


def write_modis_brdf_files(acquisitions, prefix, brdf_data):
    """Generate brdf input file"""
    ref_acqs = [a for a in acquisitions if a.band_type == gaip.REF]

    for acq in ref_acqs:
        band = str(acq.band_num)
        modis_brdf_filename = prefix + band + ".txt"

        with open(modis_brdf_filename, 'w') as outfile:
            outfile.write("%f %f %f\n" % \
                    (brdf_data[(band, 'iso')]['value'],
                     brdf_data[(band, 'vol')]['value'],
                     brdf_data[(band, 'geo')]['value']))

            outfile.write(str(l1t_input_dataset.bias[band_number]) + " " +
                           str(l1t_input_dataset.gain[band_number]) + " " +
                           str(solar_irrad_data[band_number]['value']) + " " +
                           str(solar_dist_data['value']) + "\n")


def run_read_modtrancor_ortho():
    """run read_modtrancor_ortho executable."""
    command = (os.path.join(CONFIG.BIN_DIR, "read_modtrancor_ortho") + ' ' +
               os.path.join(CONFIG.work_path, 'CENTRELINE') + ' ' +
               os.path.join(CONFIG.work_path, 'SAT_V.bin') + ' ' +
               os.path.join(CONFIG.work_path, 'COORDINATOR') + ' ' +
               os.path.join(CONFIG.work_path, 'BOXLINE'))

    result = execute(command_string=command,
                     cwd=os.path.join(CONFIG.work_path, 'mod'))
    if result['returncode']:
        raise Exception('read_modtrancor_ortho failed')


def generate_modtran_inputs(modtran_input, coordinator, sat_view_zenith,
                            sat_azimuth, lon_grid, lat_grid, coords, albedos,
                            fname_format, workdir):
    """Generate MODTRAN input files."""
    cmd = pjoin(BIN_DIR, 'input_modtran_ortho_ula')

    args = [cmd, modtran_input, coordinator, sat_view_zenith, sat_azimuth]

    targets = []
    for coord in coordinators:
        for albedo in albedos:
            target = fname_format.format(coord=coord, albedo=albedo)
            targets.append(pjoin(workdir, target))

    args.extend(targets)
    args.append(lon_grid)
    args.append(lat_grid)

    subprocess.check_call(args)

    return targets

def reformat_as_tp5(coords, albedos, profile, input_format, output_format,
                    workdir):
    """Reformat the MODTRAN input files in `tp5` format."""
    cmd = pjoin(BIN_DIR, 'refort_tp5_ga')

    targets = []
    for coord in coordinators:
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
    cmd = pjoin(BIN_DIR, 'refort_tp5_ga_trans')

    targets = []
    for coord in coordinators:
        src = input_format.format(coord=coord)
        dst = output_format.format(coord=coord)
        targets.append(pjoin(workdir, dst))

        args = [cmd, pjoin(workdir, src), profile, pjoin(workdir, dst)]

        subprocess.check_call(args)

    return targets


def run_modtran(modtran_exe, workpath):
    """Run MODTRAN."""
    subprocess.check_call([modtran_exe], cwd=workpath)




def extract_flux(coords, albedos, input_format, output_format, satfilter):
    cmd = pjoin(BIN_DIR, 'read_flx_ga')

    for coord in coordinators:
        for albedo in albedos:
            src = input_format.format(coord=coord, albedo=albedo)
            dst = output_format.format(coord=coord, albedo=albedo)
            args = [cmd, src, satfilter, dst]

            subprocess.check_call(args)


def extract_flux_trans(coords, input_format, output_format, satfilter):
    cmd = pjoin(BIN_DIR, 'read_flx_ga_trans')

    for coord in coordinators:
        src = input_format.format(coord=coord)
        dst = output_format.format(coord=coord)
        args = [cmd, src, satfilter, dst]

        subprocess.check_call(args)



def calc_coefficients(coords, chn_input_fmt, dir_input_fmt,
                      output_fmt, satfilter, workpath):
    cmd = pjoin(BIN_DIR, 'coefficient')

    for coord in coords:
        args = [cmd, satfilter,
                pjoin(workpath, chn_input_fmt.format(coord=coord, albedo=0),
                pjoin(workpath, chn_input_fmt.format(coord=coord, albedo=1),
                pjoin(workpath, dir_input_fmt.format(coord=coord, albedo=0),
                pjoin(workpath, dir_input_fmt.format(coord=coord, albedo=1),
                pjoin(workpath, dir_input_fmt.format(coord=coord, albedo='t'),
                pjoin(workpath, output_fmt.format(coord=coord)]

        subprocess.check_call(args, cwd=workpath)

def reformat_atmo_params(acqs, coords, satfilter, factors, input_fmt,
                         output_fmt, workpath):

    cmd = pjoin(BIN_DIR, 'read_modtran')
    
    bands = [str(a.band_num) for a in acqs]

    args = [cmd, satfilter]
    for coord in coords:
        args.append(input_fmt.format(coord=coord))

    for band in bands:
        for factor in factors:
            args.append(output_fmt.format(factor=factor, band=band))

    subprocess.check_call(args, cwd=workpath)

def bilinear_interpolate(arqs, factors, coordinator, boxline, centreline,
                         input_fmt, output_fmt, workpath):
    
    cmd = pjoin(BIN_DIR, 'binear_ortho')

    bands = [str(a.band_num) for a in acqs]

    for coord in coords:
        args.append(input_fmt.format(coord=coord))

    for band in bands:
        for factor in factors:
            args = [cmd, coordinator,
                    input_fmt.format(coord=coord, band=band),
                    boxline,
                    centreline,
                    output_fmt.format(factor=factor, band=band)]

            subprocess.check_call(args, cwd=workpath)



def run_brdf_sim_bin(l1t_input_dataset, band_number, sublogger, work_path, bin_dir, DEBUG, L7_SLC_DATE):
    """
    Perform atmospheric correction (NBAR). This runs either ``brdf_sim_bin`` or ``brdf_sim_bin_slc`` depending
    on a few conditions (please see code).

    :param l1t_input_dataset:
        The L1 dataset to perform NBAR on

    :type l1t_input_dataset: :py:class:`ULA3.dataset.SceneDataset`.


    :param band_number:
        The band to run NBAR for.

    :param sublogger:
        A logger that information about the call, and possibly information about any errors which occur, will
        be written to. This should support the same interface as a 'standard' logger (i.e. as returned by
        :py:func:`logging.getlogger`).

    :param work_path:
        Base directory for input output directories.

    :param bin_dir:
        The directory where the programs ``coefficient`` is located.

    :param DEBUG:
        Should intermediate files be removed.

    :type DEBUG:
        logical

    :param L7_SLC_DATE:
        The date beyond which ``brdf_sim_bin_slc`` may be used.

    :todo:
        This function assumes the location of the MODTRAN parameter files on disk, rather than having them
        passed in as arguments (which is the approach taken in the :py:func:`ULA3.tc.run_brdfterrain`). This
        ties should be changed so that it can be useful outside the context of the image processor.

    """

    # Drop the factors that are only required for terrain correction.
    FULL_FACTOR_LIST = modtran_config.FULL_FACTOR_LIST[:-3]

    def write_band_file(band_number):
        """Function to save band data as a raw binary file and return filename
        """
        band_string = str(band_number)

        # TODO: make this more efficient for FST files - create symlinks
        # instead?
        band_filename = os.path.join(work_path, 'band' + band_string + '.dat')
        array = l1t_input_dataset.band_read_as_array(band_number)
        array.tofile(band_filename)
        return band_filename

    band_string = str(band_number)

    modtran_work_dir = os.path.join(work_path, 'mod')

    if (l1t_input_dataset.satellite.NAME == 'Landsat7' and
            l1t_input_dataset.sensor.startswith('ETM') and
            l1t_input_dataset.scene_centre_date > L7_SLC_DATE):
        brdf_sim_bin_exe = os.path.join(bin_dir, 'brdf_sim_bin_slc')
    elif (l1t_input_dataset.satellite.NAME == 'Landsat-8'):
        brdf_sim_bin_exe = os.path.join(bin_dir, 'brdf_sim_bin_LS8')
    else:
        brdf_sim_bin_exe = os.path.join(bin_dir, 'brdf_sim_bin')

    band_filename = write_band_file(band_number)
    command_string = (brdf_sim_bin_exe + ' ' +
                      os.path.join(work_path, 'brdf_modis_band' + band_string + '.txt') + ' ' +
                      os.path.join(work_path, 'STARTEND') + ' ' +
                      os.path.join(work_path, 'COORDINATOR') + ' ' +
                      band_filename + ' ' +
                      os.path.join(work_path, 'SOL_Z.bin') + ' ' +
                      os.path.join(work_path, 'SAT_V.bin') + ' ' +
                      os.path.join(work_path, 'REL_AZ.bin')
                      )

    for factor in FULL_FACTOR_LIST:
        command_string += ' ' + \
            os.path.join(
                modtran_work_dir, factor + '_band' + band_string + '.bin')

    command_string += ' ' + (os.path.join(work_path, 'ref_top_b' + band_string + '.bin') + ' ' +
                             os.path.join(work_path, 'ref_nobrdf_b' + band_string + '.bin') + ' ' +
                             os.path.join(
                                 work_path, 'ref_wbrdf_b' + band_string + '.bin')
                             )

    result = execute(command_string=command_string, cwd=work_path)

    if result['returncode']:
        raise Exception('%s failed for band %s' %
                        (brdf_sim_bin_exe, band_string))

    # Clean up band file if not required for debugging
    if not DEBUG:
        os.remove(band_filename)
