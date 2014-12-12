import os
import errno
import math
import numpy
from os.path import join as pjoin, exists
from osgeo import gdal
from ULA3.utils import execute


"""
Provides an interface to MODTRAN and some of the 'post-processing' steps that would be applied
to the outputs thereof when undertaking in the process of performing a BRDF correction.
"""

# TODO: not sure if the sublogger passed in makes sense - might have to
# change this.


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


def create_satellite_filter_file(acquisitions, work_path, satfilter_path):
    """Generate satellite filter input file."""
    sat_filter_file = pjoin(work_path, 'SATELLITEFILTER')
    refbands = [a for a in acquisitions if a.band_type == gaip.REF]
    filterfile = acquisitions[0].spectral_filter_file
    filterpath = os.path.join(satfilter_path, filterfile)

    with open(sat_filter_file, 'w') as outfile:
        outfile.write("%i\n" % len(refbands))
        outfile.write("%s\n" % filterpath)

    return sat_filter_file


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





def run_modtran(coordinator, albedo, sublogger, work_path, modtran_exe):
    """
    Run MODTRAN for a specified coordinator and albedo. This assumes that the directory structure and files have
    been setup prior to the call. (see :py:func:`prepare_modtran_input`). Note that if this is done as part of BRDF
    or NBAR correction, some of the outputs from :py:func:`prepare_modtran_input` will be required for other steps,
    even if they are not used here).

    This makes a system call (using :py:func:`ULA3.utils.execute`) to ``modtran_exe`` in the working
    directory ``os.path.join(work_path, 'mod', coordinator, 'alb_' + albedo)``.

    :param coordinator:
        A string denoting the location within the grid that MODTRAN should be executed for. Used
        to specify a directory (containing appropriate inputs) which the outputs will be written to.

    :param albedo:
        A string specifying the albedo to use. Used to specify a directory (containing appropriate inputs)
        which the outputs will be written to.

    :param sublogger:
        A logger that information about the call, and possibly information about any errors which occur, will
        be written to. This should support the same interface as a 'standard' logger (i.e. as returned by
        :py:func:`logging.getlogger`).

    :param work_path:
        Base directory for input output directories.

    :param modtran_exe:
        The path of the MODTRAN program.
    """
    modtran_work_dir = os.path.join(
        work_path, 'mod', coordinator, 'alb_' + albedo)

    result = execute(command_string=modtran_exe, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed for coordinator %s, albedo %s' %
                        (modtran_exe, coordinator, albedo))


def run_flux(coordinator, albedo, sublogger, work_path, bin_dir):
    """
    Extract information from MODTRAN output file \*.flx.

    This makes a system call to either ``read_flx_ga_trans`` or ``read_flx_ga`` (using
    :py:func:`ULA3.utils.execute`) in the working directory
    ``os.path.join(work_path, 'mod', coordinator, 'alb_' + albedo)``, passing the appropriate arguments.
    It assumes that :py:func:`prepare_modtran_input` and :py:func:`run_modtran` have been run in the appropriate
    location.

    :param coordinator:
        A string denoting the location within the grid that MODTRAN should be executed for. Used
        to specify a directory (containing appropriate inputs) which the outputs will be written to.

    :param albedo:
        A string specifying the albedo to use. Used to specify a directory (containing appropriate inputs)
        which the outputs will be written to.

    :param sublogger:
        A logger that information about the call, and possibly information about any errors which occur, will
        be written to. This should support the same interface as a 'standard' logger (i.e. as returned by
        :py:func:`logging.getlogger`).

    :param work_path:
        Base directory for input output directories.

    :param bin_dir:
        The directory where the programs ``read_flx_ga_trans`` and ``read_flx_ga`` are located.
    """
    runtrans_exe = 'read_flx_ga_trans' if albedo == 't' else 'read_flx_ga'
    modtran_work_dir = os.path.join(
        work_path, 'mod', coordinator, 'alb_' + albedo)

    command_string = (os.path.join(bin_dir, runtrans_exe) + ' ' +
                      os.path.join(modtran_work_dir, coordinator + '_alb_' + albedo + '.flx') + ' ' +
                      os.path.join(work_path, 'SATELLITEFILTER') + ' ' +
                      os.path.join(
                          work_path, 'mod', coordinator + '_alb_' + albedo + '.dir')
                      )

    result = execute(command_string=command_string, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed for coordinator %s, albedo %s' %
                        (runtrans_exe, coordinator, albedo))


def run_coefficient(coordinator, sublogger, work_path, bin_dir):
    """
    Compute the atmospheric parameters needed by BRDF and atmospheric correction couple model.

    This makes a system call to ``os.path.join(bin_dir, 'coefficient')`` (using :py:func:`ULA3.utils.execute`)
    in the working directory ``os.path.join(work_path, 'mod')``, passing the appropriate arguments.
    It assumes that :py:func:`run_flux` and its dependencies have been run in the appropriate location.

    :param coordinator:
        A string denoting the location within the grid that MODTRAN should be executed for. Used
        to specify a directory (containing appropriate inputs) which the outputs will be written to.

    :param sublogger:
        A logger that information about the call, and possibly information about any errors which occur, will
        be written to. This should support the same interface as a 'standard' logger (i.e. as returned by
        :py:func:`logging.getlogger`).

    :param work_path:
        Base directory for input output directories.

    :param bin_dir:
        The directory where the programs ``coefficient`` is located.
    """
    runcoefficient_exe = os.path.join(bin_dir, 'coefficient')
    modtran_work_dir = os.path.join(work_path, 'mod')

    command_string = (runcoefficient_exe + ' ' +
                      os.path.join(work_path, 'SATELLITEFILTER') + ' ' +
                      os.path.join(modtran_work_dir, coordinator, 'alb_0', coordinator + '_alb_0.chn') + ' ' +
                      os.path.join(modtran_work_dir, coordinator, 'alb_1', coordinator + '_alb_1.chn') + ' ' +
                      os.path.join(modtran_work_dir, coordinator + '_alb_0.dir') + ' ' +
                      os.path.join(modtran_work_dir, coordinator + '_alb_1.dir') + ' ' +
                      os.path.join(modtran_work_dir, coordinator + '_alb_t.dir') + ' ' +
                      os.path.join(modtran_work_dir, coordinator + '_alb.txt')
                      )

    result = execute(command_string=command_string, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed for coordinator %s' %
                        (runcoefficient_exe, coordinator))


def read_modtran(band_strings, sublogger, work_path, bin_dir):
    """
    Reformat the atmospheric parameters produced by MODTRAN for four boxes. These are needed to conduct bilinear analysis.

    This makes a system call to ``os.path.join(bin_dir, 'read_modtran')`` (using :py:func:`ULA3.utils.execute`)
    in the working directory ``os.path.join(work_path, 'mod')``, passing the appropriate arguments.
    It assumes that :py:func:`run_coefficient` has been run in the appropriate location.

    :param band_strings:
        An iterable over the list of bands for which ??? will be produced.

    :param sublogger:
        A logger that information about the call, and possibly information about any errors which occur, will
        be written to. This should support the same interface as a 'standard' logger (i.e. as returned by
        :py:func:`logging.getlogger`).

    :param work_path:
        Base directory for input output directories.

    :param bin_dir:
        The directory where the programs ``coefficient`` is located.
    """
    read_modtran_exe = os.path.join(bin_dir, 'read_modtran')
    modtran_work_dir = os.path.join(work_path, 'mod')

    command_string = (
        read_modtran_exe + ' ' + os.path.join(work_path, 'SATELLITEFILTER'))

    for coordinator in modtran_config.FULL_COORD_LIST:
        command_string += ' ' + \
            os.path.join(modtran_work_dir, coordinator + '_alb.txt')

    for band_string in band_strings:
        # NOTE: The last three factors in this list are only required for terrain correction. we might
        #       want to exclude these if we are only doing NBAR.
        # NOTE: The order of these strings is important as read_modtran.f depends on it... this is
        #       why they are hard coded here.
        # NOTE: This order is different to that in modtran_config.FULL_FACTOR_LIST, which is the order
        #       which brdf_sim_bin.f depends on.
        # TODO: It would be good to remove this dependency from the fortran
        # codes.
        for factor in ['fv', 'fs', 'b', 's', 'a', 'dir', 'dif', 'ts']:
            command_string += ' ' + \
                os.path.join(
                    modtran_work_dir, factor + '_out_b' + band_string + '.txt')

    result = execute(command_string=command_string, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed' % read_modtran_exe)


def run_bilinear_ortho(band_number, factor, sublogger, work_path, bin_dir):
    """
    Runs the (interpolation) program "binear_ortho" (which does contain a typo in its name).

    :todo: Someone whom understands what this program is used for should document it.
    """
    band_string = str(band_number)

    bilinear_ortho_exe = os.path.join(
        bin_dir, 'binear_ortho')  # Note typo in name
    modtran_work_dir = os.path.join(work_path, 'mod')

    output_filename = os.path.join(
        modtran_work_dir, factor + '_band' + band_string + '.bin')

    command_string = (bilinear_ortho_exe + ' ' +
                      os.path.join(work_path, 'COORDINATOR') + ' ' +
                      os.path.join(modtran_work_dir, factor + '_out_b' + band_string + '.txt') + ' ' +
                      os.path.join(work_path, 'BOXLINE') + ' ' +
                      os.path.join(work_path, 'CENTRELINE') + ' ' +
                      output_filename
                      )

    result = execute(command_string=command_string, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed' % bilinear_ortho_exe)

    return output_filename


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


def prepare_modtran_input(
        l1t_input_dataset,
        ha_file,
        modtran_input_file,
        brdf_data,
        solar_irrad_data,
        solar_dist_data,
        ozone_data,
        water_vapour_data,
        aerosol_data,
        elevation_data,
        max_view_angle,
        CONFIG):
    """
    Prepare all the input files for ModTran.

    :param l1t_input_dataset:
        The input dataset for which MODTRAN will run (or more precisely, for which NBAR or TC will be run). This is used
        to determine extents, times etc. which will be used when extracting ancilliary data.

    :type l1t_input_dataset: :py:class:`ULA3.dataset.SceneDataset`

    :param ha_file:
        The name of the header angle file to create for MODTRAN.

    :param modtran_input_file:
        The name of the input file to create for MODTRAN.

    :param brdf_data: ???

    :param solar_irrad_data: ???

    :param solar_dist_data: ???

    :param ozone_data: ???

    :param water_vapour_data: ???

    :param aerosol_data: ???

    :param elevation_data: ???

    :param max_view_angle: ???

    :param CONFIG:
        Contains a bunch of stuff that should really be specified as individual parameters. Please see the code
        to find out exactly what.

    :type CONFIG: :py:class:`ULA3.image_processor.ProcessorConfig`.
    """
    def run_read_modtrancor_ortho():
        """run read_modtrancor_ortho executable
        """
        command = (os.path.join(CONFIG.BIN_DIR, "read_modtrancor_ortho") + ' ' +
                   os.path.join(CONFIG.work_path, 'CENTRELINE') + ' ' +
                   os.path.join(CONFIG.work_path, 'SAT_V.bin') + ' ' +
                   os.path.join(CONFIG.work_path, 'COORDINATOR') + ' ' +
                   os.path.join(CONFIG.work_path, 'BOXLINE'))

        result = execute(command_string=command,
                         cwd=os.path.join(CONFIG.work_path, 'mod'))
        if result['returncode']:
            raise Exception('read_modtrancor_ortho failed')

    def run_input_modtran_ortho_ula():
        """ Run input_modtran_ortho_ula executable
        """
        command = (os.path.join(CONFIG.BIN_DIR, "input_modtran_ortho_ula") + ' ' +
                   os.path.join(CONFIG.work_path, 'MODTRANINPUT') + ' ' +
                   os.path.join(CONFIG.work_path, 'COORDINATOR') + ' ' +
                   os.path.join(CONFIG.work_path, 'SAT_V.bin') + ' ' +
                   os.path.join(CONFIG.work_path, 'SAT_AZ.bin'))

        print "\n\ncommand: %s\n\n" % command

        # Only run input_modtran_ortho_ula for albedos 0 & 1
        for coordinator in modtran_config.FULL_COORD_LIST:
            for albedo in [albedo.lower() for albedo in modtran_config.FULL_ALBEDO_LIST if albedo.lower() != 't']:
                command += ' ' + os.path.join(CONFIG.work_path, 'mod',  # coordinator, 'alb_' + albedo,
                                              coordinator + '_alb_' + albedo + '.txt')

        command += (' ' + os.path.join(CONFIG.work_path, 'LON.bin') + ' ' +
                    os.path.join(CONFIG.work_path, 'LAT.bin'))

        result = execute(command_string=command,
                         cwd=os.path.join(CONFIG.work_path, 'mod'))
        if result['returncode']:
            raise Exception('input_modtran_ortho_ula failed')

    def run_refort_tp5_ga():
        """ Run refort_tp5_ga / refort_tp5_ga_trans executable
        """
        # Set MODTRAN default profile
        modtranProfile = "tropical"
        if l1t_input_dataset.lonlats['CENTRE'][1] < -23.0:
            modtranProfile = "midlat_summer"

        for coordinator in modtran_config.FULL_COORD_LIST:
            for albedo in [albedo.lower() for albedo in modtran_config.FULL_ALBEDO_LIST]:
                if albedo == 't':
                    executable = 'refort_tp5_ga_trans'
                    text_file = coordinator + '_alb_0.txt'
                else:
                    executable = 'refort_tp5_ga'
                    text_file = coordinator + '_alb_' + albedo + '.txt'

                command = (os.path.join(CONFIG.BIN_DIR, executable) + ' ' +
                           os.path.join(CONFIG.work_path, 'mod', text_file) + ' ' +
                           os.path.join(CONFIG.ANCI_ROOT, 'lookup_tables', 'default_profile', modtranProfile + '.tp5') + ' ' +
                           os.path.join(CONFIG.work_path, 'mod', coordinator, 'alb_' + albedo,
                                        coordinator + '_alb_' + albedo + '.tp5'))

                result = execute(command_string=command,
                                 cwd=os.path.join(CONFIG.work_path, 'mod'))
                if result['returncode']:
                    raise Exception('%s failed' % executable)

    # TODO The STARTEND file is no longer used in the NBAR/TC code.
    # JS 20141205
    def create_startend_file():
        """Create startend file.
        """
        startend_path = os.path.join(CONFIG.work_path, 'STARTEND')

        bdata = l1t_input_dataset.band_read_as_array(
            l1t_input_dataset.root_band_number)
        startend_record_fmt = '%12d%12d%12d%12d%12d'

        record_list = []
        for irow in xrange(l1t_input_dataset.shape[0]):
            # Array of populated pixel indices across the row.
            j = numpy.where(bdata[irow, :] != 0)[0]
            if j.size == 0:
                # Empty line (no data)
                tup = ((irow + 1), 1, 1, 0, 0)
            else:
                # Populated line (all values 1-base)
                #  <row> <first> <middle> <last> <span>
                j0 = j.min()
                j1 = j.max()
                jspan = j1 - j0
                # Round to match mid values produced by angle_ortho_bin.
                jc = j0 + int(round(float(jspan) / 2))
                tup = ((irow + 1), (j0 + 1), (jc + 1), (j1 + 1), (jspan + 1))
            record_list.append(startend_record_fmt % tup)

        text = '\n'.join(record_list) + '\n'

        try:
            fd = open(startend_path, 'w')
            fd.write(text)
        finally:
            fd.close()

    satellite = l1t_input_dataset.satellite
    assert satellite, 'Unable to find Satellite object'

    create_header_angle_file(max_view_angle)
    #write_modtran_input()
    #write_modis_brdf_files()
    #create_startend_file()
    #run_read_modtrancor_ortho()
    #run_input_modtran_ortho_ula()
    run_refort_tp5_ga()
