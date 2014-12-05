import os, errno, math, numpy
from osgeo import gdal
from ULA3.utils import execute


"""
Provides an interface to MODTRAN and some of the 'post-processing' steps that would be applied
to the outputs thereof when undertaking in the process of performing a BRDF correction.
"""

# TODO: not sure if the sublogger passed in makes sense - might have to change this.
class modtran_config(object):
    FULL_COORD_LIST = ['TL','TM','TR','ML','MM','MR','BL','BM','BR']
    FULL_ALBEDO_LIST = ['0','1','t']
    FULL_FACTOR_LIST = ['fv', 'fs', 'a', 'b', 's', 'dir', 'dif', 'ts']





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
    modtran_work_dir = os.path.join(work_path, 'mod', coordinator, 'alb_' + albedo)

    result = execute(command_string=modtran_exe, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed for coordinator %s, albedo %s' % (modtran_exe, coordinator, albedo))





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
    modtran_work_dir = os.path.join(work_path, 'mod', coordinator, 'alb_' + albedo)

    command_string=(os.path.join(bin_dir, runtrans_exe) + ' ' +
            os.path.join(modtran_work_dir, coordinator + '_alb_' + albedo + '.flx') + ' ' +
            os.path.join(work_path, 'SATELLITEFILTER') + ' ' +
            os.path.join(work_path, 'mod', coordinator + '_alb_' + albedo + '.dir')
            )

    result = execute(command_string=command_string, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed for coordinator %s, albedo %s' % (runtrans_exe, coordinator, albedo))





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

    command_string=(runcoefficient_exe + ' ' +
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
        raise Exception('%s failed for coordinator %s' % (runcoefficient_exe, coordinator))





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

    command_string = (read_modtran_exe + ' ' + os.path.join(work_path, 'SATELLITEFILTER'))

    for coordinator in modtran_config.FULL_COORD_LIST:
        command_string += ' ' + os.path.join(modtran_work_dir, coordinator + '_alb.txt')

    for band_string in band_strings:
        # NOTE: The last three factors in this list are only required for terrain correction. we might
        #       want to exclude these if we are only doing NBAR.
        # NOTE: The order of these strings is important as read_modtran.f depends on it... this is
        #       why they are hard coded here.
        # NOTE: This order is different to that in modtran_config.FULL_FACTOR_LIST, which is the order
        #       which brdf_sim_bin.f depends on.
        # TODO: It would be good to remove this dependency from the fortran codes.
        for factor in ['fv', 'fs', 'b', 's', 'a', 'dir', 'dif', 'ts']:
            command_string += ' ' + os.path.join(modtran_work_dir, factor + '_out_b' + band_string + '.txt')

    result = execute(command_string=command_string, cwd=modtran_work_dir)

    if result['returncode']:
        raise Exception('%s failed' % read_modtran_exe)





def run_bilinear_ortho(band_number, factor, sublogger, work_path, bin_dir):
    """
    Runs the (interpolation) program "binear_ortho" (which does contain a typo in its name).

    :todo: Someone whom understands what this program is used for should document it.
    """
    band_string = str(band_number)

    bilinear_ortho_exe = os.path.join(bin_dir, 'binear_ortho') # Note typo in name
    modtran_work_dir = os.path.join(work_path, 'mod')

    output_filename = os.path.join(modtran_work_dir, factor + '_band' + band_string + '.bin')

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

        # TODO: make this more efficient for FST files - create symlinks instead?
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
        command_string += ' ' + os.path.join(modtran_work_dir, factor + '_band' + band_string + '.bin')

    command_string += ' ' + (os.path.join(work_path, 'ref_top_b' + band_string + '.bin') + ' ' +
        os.path.join(work_path, 'ref_nobrdf_b' + band_string + '.bin') + ' ' +
        os.path.join(work_path, 'ref_wbrdf_b' + band_string + '.bin')
        )

    result = execute(command_string=command_string, cwd=work_path)

    if result['returncode']:
        raise Exception('%s failed for band %s' % (brdf_sim_bin_exe, band_string))

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
    def create_modtran_dirs():
        # Create all modtran subdirectories
        for coordinator in modtran_config.FULL_COORD_LIST:
            for albedo in [albedo.lower() for albedo in modtran_config.FULL_ALBEDO_LIST]:
                modtran_work_dir = os.path.join(
                    CONFIG.work_path, 'mod', coordinator, 'alb_' + albedo)
                mod5root_in = os.path.join(modtran_work_dir, 'mod5root.in')
                data_dir = os.path.join(CONFIG.MODTRAN_ROOT, 'DATA')
                symlink_dir = os.path.join(modtran_work_dir, 'DATA')

                try:
                    os.makedirs(modtran_work_dir)
                except OSError, e:
                    if e.errno != errno.EEXIST:
                        raise

                try:
                    os.symlink(data_dir, symlink_dir)
                except OSError, e:
                    if e.errno != errno.EEXIST:
                        raise

                out_file = open(mod5root_in, 'w')
                out_file.write(coordinator + '_alb_' + albedo + '\n')
                out_file.close()

        symlink_dir = os.path.join(CONFIG.work_path, 'mod', 'DATA')
        try:
            os.symlink(data_dir, symlink_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise

    # TODO Re-Write so the correct data is written and the formatting
    # is more human readable. JS 20141205
    def create_header_angle_file(max_view_angle):#ha_file, l1t_input_dataset, max_view_angle=None):
        """Create header angle file.

        Arguments:
            ha_file: header angle file path
            max_view_angle: view angle cutoff (degrees)
        """
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
        #    '%f %f' % (l1t_input_dataset.scene_centre_lat, l1t_input_dataset.scene_centre_long),
            '%f %f' % (l1t_input_dataset.lonlats['CENTRE'][1],l1t_input_dataset.lonlats['CENTRE'][0]),
            '%f %f %f' % (satellite.SEMI_MAJOR_AXIS,
                          math.degrees(satellite.INCLINATION), # Convert to degrees
                          satellite.OMEGA),# Leave in radians
            '%f\n' % max_view_angle,
        ])

        try:
            fd = open(ha_file, 'w')
            fd.write(ha_text)
        finally:
            fd.close()


    def write_modtran_input():
        """Generate modtran input file"""

        out_file = open(modtran_input_file, 'w')

        out_file.write("%f %f\n" % (l1t_input_dataset.lonlats['UL'][1],l1t_input_dataset.lonlats['UL'][0]))
        out_file.write("%f\n" % satellite.NOMINAL_PIXEL_DEGREES)
        out_file.write("%s\n" % ozone_data['value'])

        # NOTE: The water vapour value for MODTRAN5 input should be the
        # scaled value, not the raw value from the GeoTIFF file. The processor
        # module should scale the value before storing it in nbar.vapourVal.
        out_file.write("%f\n" % water_vapour_data['value'])

        out_file.write("DATA/%s\n" % satellite.SPECTRAL_FILTER_FILE)

        # NOTE: Processor module scales raw value by -1 for MODTRAN 5.
        out_file.write("%f \n" % -aerosol_data['value'])

        # NOTE: The processor scales the raw DEM value by 0.001.
        out_file.write("%f\n" % elevation_data['value'])

        out_file.write("Annotation, %s\n" % l1t_input_dataset.scene_centre_date.strftime('%Y-%m-%d'))
        out_file.write("%d\n" % (satellite.ALTITUDE / 1000.0)) # Convert from metres to kilometres
        out_file.write("%d\n" % l1t_input_dataset.DOY)
        # Decimal hour
        out_file.write("%f\n" % l1t_input_dataset.decimal_hour)

        out_file.close()

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
                command += ' ' + os.path.join(CONFIG.work_path, 'mod', #coordinator, 'alb_' + albedo,
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

    def writeSatFilterFile():
        """Generate satellite filter input file"""
        satFilterFile = os.path.join(CONFIG.work_path, 'SATELLITEFILTER')

        out_file = open(satFilterFile, 'w')

        out_file.write("%s\n" % len(l1t_input_dataset.bands('REFLECTIVE')))
        satFilterPath = os.path.join(CONFIG.DIR_SatFilter, satellite.SPECTRAL_FILTER_FILE)
        out_file.write("%s\n" % satFilterPath)

        out_file.close()

    def write_modis_brdf_files():
        """Generate modtran input file"""
        for band_number in l1t_input_dataset.bands('REFLECTIVE'):
            band_string = str(band_number)

            modis_brdf_filename = os.path.join(CONFIG.work_path, "brdf_modis_band" + band_string + ".txt")

            out_file = open(modis_brdf_filename, 'w')

            out_file.write("%f %f %f\n" % (brdf_data[(band_number, 'iso')]['value'],
                                        brdf_data[(band_number, 'vol')]['value'],
                                        brdf_data[(band_number, 'geo')]['value']))
            out_file.write(str(l1t_input_dataset.bias[band_number]) + " " +
                        str(l1t_input_dataset.gain[band_number]) + " " +
                        str(solar_irrad_data[band_number]['value']) + " " +
                        str(solar_dist_data['value']) + "\n")
            out_file.close()

    # TODO REMOVE as no longer needed. JS 20141205
    def create_startend_file():
        """Create startend file.
        """
        startend_path = os.path.join(CONFIG.work_path, 'STARTEND')

        bdata = l1t_input_dataset.band_read_as_array(l1t_input_dataset.root_band_number)
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
                jc = j0 + int(round(float(jspan)/2))
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

    create_modtran_dirs()
    writeSatFilterFile()
    create_header_angle_file(max_view_angle)
    write_modtran_input()
    write_modis_brdf_files()
    create_startend_file()
    run_read_modtrancor_ortho()
    run_input_modtran_ortho_ula()
    run_refort_tp5_ga()

