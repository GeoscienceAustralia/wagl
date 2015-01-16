"""
NBAR Workflow
-------------

Workflow settings can be configured in `nbar.cfg` file.

"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals

import luigi
import gaip
import cPickle as pickle
import os

from os.path import join as pjoin, dirname, exists

CONFIG = luigi.configuration.get_config()
CONFIG.add_config_path(pjoin(dirname(__file__), 'nbar.cfg'))


def save(target, value):
    """Save `value` to `target` where `target` is a `luigi.Target` object. If
    the target filename ends with `pkl` then pickle the data. Otherwise, save
    as text."""
    with target.open('w') as outfile:
        if target.fn.endswith('pkl'):
            pickle.dump(value, outfile)
        else:
            print >>outfile, value


def load(target):
    """Load data from `target` where `target` is a `luigi.Target`."""
    if not target.fn.endswith('pkl'):
        raise IOError('Cannot load non-pickled object')
    with target.open('r') as infile:
        return pickle.load(infile)


def load_value(target):
    """Load the value from `target`."""
    if isinstance(target, str):
        target = luigi.LocalTarget(target)
    data = load(target)
    try:
        return data['value']
    except KeyError:
        return data


class GetElevationAncillaryData(luigi.Task):

    """Get ancillary elevation data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'dem_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        dem_path = CONFIG.get('ancillary', 'dem_path')
        value = gaip.get_elevation_data(geobox.centre_lonlat, dem_path)
        save(self.output(), value)


class GetOzoneAncillaryData(luigi.Task):

    """Get ancillary ozone data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'ozone_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        ozone_path = CONFIG.get('ancillary', 'ozone_path')
        centre = geobox.centre_lonlat
        dt = acqs[0].scene_center_datetime
        value = gaip.get_ozone_data(ozone_path, centre, dt)
        save(self.output(), value)


class GetSolarIrradianceAncillaryData(luigi.Task):

    """Get ancillary solar irradiance data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'irrad_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        solar_path = CONFIG.get('ancillary', 'solarirrad_path')
        value = gaip.get_solar_irrad(acqs, solar_path)
        save(self.output(), value)


class GetSolarDistanceAncillaryData(luigi.Task):

    """Get ancillary solar distance data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'sundist_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        sundist_path = CONFIG.get('ancillary', 'sundist_path')
        value = gaip.get_solar_dist(acqs[0], sundist_path)
        save(self.output(), value)


class GetWaterVapourAncillaryData(luigi.Task):

    """Get ancillary water vapour data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'vapour_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        vapour_path = CONFIG.get('ancillary', 'vapour_path')
        value = gaip.get_water_vapour(acqs[0], vapour_path)
        save(self.output(), value)


class GetAerosolAncillaryData(luigi.Task):

    """Get ancillary aerosol data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'aerosol_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        aerosol_path = CONFIG.get('ancillary', 'aerosol_path')
        value = gaip.get_aerosol_data(acqs[0], aerosol_path)
        save(self.output(), value)


class GetBrdfAncillaryData(luigi.Task):

    """Get ancillary BRDF data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'brdf_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        brdf_path = CONFIG.get('ancillary', 'brdf_path')
        brdf_premodis_path = CONFIG.get('ancillary', 'brdf_premodis_path')
        value = gaip.get_brdf_data(acqs[0], brdf_path, brdf_premodis_path,
                                   out_path)
        save(self.output(), value)


class GetAncillaryData(luigi.Task):

    """Get all ancillary data. This a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GetElevationAncillaryData(self.l1t_path, self.out_path),
                GetOzoneAncillaryData(self.l1t_path, self.out_path),
                GetSolarDistanceAncillaryData(self.l1t_path, self.out_path),
                GetSolarIrradianceAncillaryData(self.l1t_path, self.out_path),
                GetWaterVapourAncillaryData(self.l1t_path, self.out_path),
                GetAerosolAncillaryData(self.l1t_path, self.out_path),
                GetBrdfAncillaryData(self.l1t_path, self.out_path)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class CalculateLonGrid(luigi.Task):

    """Calculate the longitude grid."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'lon_grid_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output()
        gaip.create_lon_grid(acqs[0], target.fn)


class CalculateLatGrid(luigi.Task):

    """Calculate the latitude grid."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'lat_grid_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output()
        gaip.create_lat_grid(acqs[0], target.fn)


class CalculateLatLonGrids(luigi.Task):

    """Calculate the longitude and latitude grids. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]


class CalculateSatelliteAndSolarGrids(luigi.Task):

    """Calculate the satellite and solar grids."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        targets = [CONFIG.get('work', 'sat_view_target'),
                   CONFIG.get('work', 'sat_azimuth_target'),
                   CONFIG.get('work', 'solar_zenith_target'),
                   CONFIG.get('work', 'solar_azimuth_target'),
                   CONFIG.get('work', 'relative_azimuth_target'),
                   CONFIG.get('work', 'time_target'),
                   CONFIG.get('work', 'centreline_target'),
                   CONFIG.get('work', 'header_angle_target')]
        return [luigi.LocalTarget(pjoin(out_path, t)) for t in targets]

    def run(self):
        out_path = self.out_path
        targets = [CONFIG.get('work', 'sat_view_target'),
                   CONFIG.get('work', 'sat_azimuth_target'),
                   CONFIG.get('work', 'solar_zenith_target'),
                   CONFIG.get('work', 'solar_azimuth_target'),
                   CONFIG.get('work', 'relative_azimuth_target'),
                   CONFIG.get('work', 'time_target')]
        targets = [pjoin(out_path, t) for t in targets]
        centreline_target = pjoin(out_path,
                                  CONFIG.get('work', 'centreline_target'))
        header_angle_target = pjoin(out_path,
                                    CONFIG.get('work', 'header_angle_target'))
        lon_target = pjoin(out_path, CONFIG.get('work', 'lon_grid_target'))
        lat_target = pjoin(out_path, CONFIG.get('work', 'lat_grid_target'))

        acqs = gaip.acquisitions(self.l1t_path)

        geobox = acqs[0].gridded_geo_box()
        cols = acqs[0].samples

        (satellite_zenith, satellite_azimuth, solar_zenith, solar_azimuth,
         relative_azimuth, time, y_cent, x_cent, n_cent) = \
            gaip.calculate_angles(acqs[0], lon_target, lat_target,
                                  npoints=12, to_disk=targets)

        gaip.create_centreline_file(geobox, y_cent, x_cent, n_cent, cols,
                                    view_max=9.0, outfname=centreline_target)

        gaip.create_header_angle_file(acqs[0], view_max=9.0,
                                      outfname=header_angle_target)


class CalculateGridsTask(luigi.Task):

    """Calculate all the grids. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatLonGrids(self.l1t_path, self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path)]


class CreateModtranDirectories(luigi.Task):

    """Create the MODTRAN work directories and input driver files."""

    out_path = luigi.Parameter()

    def output(self):
        out_path = self.out_path
        input_format = CONFIG.get('modtran', 'input_format')
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        targets = []
        for coord in coords:
            for albedo in albedos:
                targets.append(input_format.format(coord=coord,
                                                   albedo=albedo))
        return [luigi.LocalTarget(pjoin(modtran_root, t)) for t in targets]

    def run(self):
        out_path = self.out_path
        modtran_exe_root = CONFIG.get('modtran', 'root')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('modtran', 'input_format')
        workpath_format = CONFIG.get('modtran', 'workpath_format')
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')

        gaip.create_modtran_dirs(coords, albedos, modtran_root,
                                 modtran_exe_root,
                                 workpath_format,
                                 input_format)


class CreateSatelliteFilterFile(luigi.Task):

    """Create the satellite filter file."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))
        return luigi.LocalTarget(target)

    def run(self):
        out_path = self.out_path
        acqs = gaip.acquisitions(self.l1t_path)
        satfilterpath = CONFIG.get('ancillary', 'satfilter_path')
        target = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))
        gaip.create_satellite_filter_file(acqs, satfilterpath,
                                          target)


class CreateModtranInputFile(luigi.Task):

    """Create the MODTRAN input file."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GetAncillaryData(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'modtran_input_target'))
        return luigi.LocalTarget(target)

    def run(self):
        out_path = self.out_path
        ozone_target = pjoin(out_path, CONFIG.get('work', 'ozone_target'))
        vapour_target = pjoin(out_path, CONFIG.get('work', 'vapour_target'))
        aerosol_target = pjoin(out_path, CONFIG.get('work', 'aerosol_target'))
        elevation_target = pjoin(out_path, CONFIG.get('work', 'dem_target'))
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output().fn
        ozone = load_value(ozone_target)
        vapour = load_value(vapour_target)
        aerosol = load_value(aerosol_target)
        elevation = load_value(elevation_target)
        gaip.write_modtran_input(acqs, target, ozone, vapour, aerosol,
                                 elevation)


class CreateModisBrdfFiles(luigi.Task):

    """Create the Modis BRDF files."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GetAncillaryData(l1t_path, self.out_path)]

    def output(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        modis_brdf_format = pjoin(out_path,
            CONFIG.get('brdf', 'modis_brdf_format'))

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        targets = []
        for acq in acqs:
            band = acq.band_num
            if band not in bands_to_process:
                continue
            modis_brdf_filename = modis_brdf_format.format(band_num=band)
            target = pjoin(out_path, modis_brdf_filename)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        outdir = self.out_path
        modis_brdf_format = pjoin(outdir,
            CONFIG.get('brdf', 'modis_brdf_format'))
        brdf_target = pjoin(outdir, CONFIG.get('work', 'brdf_target'))
        brdf_data = load_value(brdf_target)
        irrad_target = pjoin(outdir, CONFIG.get('work', 'irrad_target'))
        solar_irrad_data = load_value(irrad_target)
        solar_dist_target = pjoin(outdir, CONFIG.get('work', 'sundist_target'))
        solar_dist_data = load_value(solar_dist_target)
        gaip.write_modis_brdf_files(acqs, modis_brdf_format, brdf_data,
                                    solar_irrad_data, solar_dist_data)


class RunModtranCorOrtho(luigi.Task):

    """Run `modtran_cor_ortho` binary."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        targets = [CONFIG.get('work', 'coordinator_target'),
                   CONFIG.get('work', 'boxline_target')]
        return [luigi.LocalTarget(pjoin(out_path, t)) for t in targets]

    def run(self):
        out_path = self.out_path
        # sources
        centreline_target = pjoin(out_path,
                                  CONFIG.get('work', 'centreline_target'))
        sat_view_zenith_target = pjoin(out_path,
                                       CONFIG.get('work', 'sat_view_target'))
        # targets
        coordinator_target = pjoin(out_path,
                                   CONFIG.get('work', 'coordinator_target'))
        boxline_target = pjoin(out_path,
                               CONFIG.get('work', 'boxline_target'))
        cwd = pjoin(out_path, CONFIG.get('work', 'read_modtrancor_ortho_cwd'))

        gaip.run_read_modtrancor_ortho(centreline_target,
                                       sat_view_zenith_target,
                                       coordinator_target,
                                       boxline_target,
                                       cwd)


class GenerateModtranInputFiles(luigi.Task):

    """Generate the MODTRAN input files by running the Fortran binary
    `input_modtran_ortho`."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [RunModtranCorOrtho(self.l1t_path, self.out_path),
                CreateModtranDirectories(self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                CreateModtranInputFile(self.l1t_path, self.out_path),
                CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('input_modtran', 'coords').split(',')
        albedos = CONFIG.get('input_modtran', 'albedos').split(',')
        output_format = CONFIG.get('input_modtran', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'input_modtran_cwd'))
        output_format = pjoin(workdir, output_format)

        targets = []
        for coord in coords:
            for albedo in albedos:
                targets.append(output_format.format(coord=coord,
                                                    albedo=albedo))
        return [luigi.LocalTarget(t) for t in targets]

    def run(self):
        out_path = self.out_path
        # sources
        modtran_input_target = pjoin(out_path,
                                     CONFIG.get('work',
                                                'modtran_input_target'))
        coordinator_target = pjoin(out_path,
                                   CONFIG.get('work', 'coordinator_target'))
        sat_view_zenith_target = pjoin(out_path,
                                       CONFIG.get('work', 'sat_view_target'))
        sat_azimuth_target = pjoin(out_path,
                                   CONFIG.get('work', 'sat_azimuth_target'))
        lon_grid_target = pjoin(out_path,
                                CONFIG.get('work', 'lon_grid_target'))
        lat_grid_target = pjoin(out_path,
                                CONFIG.get('work', 'lat_grid_target'))

        coords = CONFIG.get('input_modtran', 'coords').split(',')
        albedos = CONFIG.get('input_modtran', 'albedos').split(',')
        fname_format = CONFIG.get('input_modtran', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'input_modtran_cwd'))

        gaip.generate_modtran_inputs(modtran_input_target,
                                     coordinator_target,
                                     sat_view_zenith_target,
                                     sat_azimuth_target,
                                     lon_grid_target,
                                     lat_grid_target,
                                     coords,
                                     albedos,
                                     fname_format,
                                     workdir)


class ReformatAsTp5(luigi.Task):

    """Reformat the MODTRAN input files in `tp5` format. This runs the
    Fortran binary `refort_tp5_ga` multiple times."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GenerateModtranInputFiles(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('reformat_tp5', 'coords').split(',')
        albedos = CONFIG.get('reformat_tp5', 'albedos').split(',')
        output_format = CONFIG.get('reformat_tp5', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'reformat_tp5_cwd'))
        targets = []
        for coord in coords:
            for albedo in albedos:
                targets.append(output_format.format(coord=coord,
                                                    albedo=albedo))
        return [luigi.LocalTarget(pjoin(workdir, t)) for t in targets]

    def run(self):
        out_path = self.out_path
        modtran_profile_path = CONFIG.get('ancillary', 'modtran_profile_path')
        profile_format = CONFIG.get('modtran', 'profile_format')
        input_format = CONFIG.get('reformat_tp5', 'input_format')
        output_format = CONFIG.get('reformat_tp5', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'reformat_tp5_cwd'))
        coords = CONFIG.get('reformat_tp5', 'coords').split(',')
        albedos = CONFIG.get('reformat_tp5', 'albedos').split(',')

        # determine modtran profile
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        centre_lon, centre_lat = geobox.centre_lonlat
        profile = 'tropical'
        if centre_lat < -23.0:
            profile = 'midlat_summer'

        profile = pjoin(modtran_profile_path,
                        profile_format.format(profile=profile))

        gaip.reformat_as_tp5(coords, albedos, profile,
                             input_format, output_format,
                             workdir)


class ReformatAsTp5Trans(luigi.Task):

    """Reformat the MODTRAN input files in `tp5` format in the transmissive
    case. This runs the Fortran binary `refort_tp5_ga_trans` multiple
    times."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GenerateModtranInputFiles(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('reformat_tp5_trans', 'coords').split(',')
        albedos = CONFIG.get('reformat_tp5_trans', 'albedos').split(',')
        workdir = pjoin(out_path,
                        CONFIG.get('work', 'reformat_tp5_trans_cwd'))
        output_format = CONFIG.get('reformat_tp5_trans', 'output_format')
        output_format = pjoin(workdir, output_format)
        targets = []
        for coord in coords:
            for albedo in albedos:
                target = output_format.format(coord=coord, albedo=albedo)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        modtran_profile_path = CONFIG.get('ancillary', 'modtran_profile_path')
        profile_format = CONFIG.get('modtran', 'profile_format')
        input_format = CONFIG.get('reformat_tp5_trans', 'input_format')
        output_format = CONFIG.get('reformat_tp5_trans', 'output_format')
        workdir = pjoin(out_path,
                        CONFIG.get('work', 'reformat_tp5_trans_cwd'))
        coords = CONFIG.get('reformat_tp5_trans', 'coords').split(',')
        albedos = CONFIG.get('reformat_tp5_trans', 'albedos').split(',')

        # determine modtran profile
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        centre_lon, centre_lat = geobox.centre_lonlat
        profile = 'tropical'
        if centre_lat < -23.0:
            profile = 'midlat_summer'

        profile = pjoin(modtran_profile_path,
                        profile_format.format(profile=profile))

        gaip.reformat_as_tp5_trans(coords, albedos, profile,
                                   input_format, output_format,
                                   workdir)


class PrepareModtranInput(luigi.Task):

    """Prepare MODTRAN inputs. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CreateModtranDirectories(self.out_path),
                CreateSatelliteFilterFile(self.l1t_path, self.out_path),
                GenerateModtranInputFiles(self.l1t_path, self.out_path),
                ReformatAsTp5(self.l1t_path, self.out_path),
                ReformatAsTp5Trans(self.l1t_path, self.out_path)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class RunModtranCase(luigi.Task):

    """Run MODTRAN for a specific `coord` and `albedo`. This task is
    parameterised this way to allow parallel instances of MODTRAN to run."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()
    coord = luigi.Parameter()
    albedo = luigi.Parameter()

    def requires(self):
        return [PrepareModtranInput(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        flux_format = CONFIG.get('modtran', 'flx_output_format')
        flux_format = pjoin(modtran_root, flux_format)
        coef_format = CONFIG.get('modtran', 'chn_output_format')
        coef_format = pjoin(modtran_root, coef_format)
        flx_target = flux_format.format(coord=self.coord, albedo=self.albedo)
        chn_target = coef_format.format(coord=self.coord, albedo=self.albedo)
        return [luigi.LocalTarget(flx_target),
                luigi.LocalTarget(chn_target)]

    def run(self):
        out_path = self.out_path
        modtran_exe = CONFIG.get('modtran', 'exe')
        workpath_format = CONFIG.get('modtran', 'workpath_format')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        workpath = workpath_format.format(coord=self.coord, albedo=self.albedo)
        gaip.run_modtran(modtran_exe, pjoin(modtran_root, workpath))


class RunModtran(luigi.Task):

    """Run MODTRAN for all coords and albedos. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        reqs = [PrepareModtranInput(self.l1t_path, self.out_path)]
        for coord in coords:
            for albedo in albedos:
                reqs.append(RunModtranCase(self.l1t_path, self.out_path, 
                                           coord, albedo))
        return reqs

    def complete(self):
        return all([t.complete() for t in self.requires()])


class ExtractFlux(luigi.Task):

    """Extract the flux data from the MODTRAN outputs."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux', 'coords').split(',')
        albedos = CONFIG.get('extract_flux', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('extract_flux', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        targets = []
        for coord in coords:
            for albedo in albedos:
                target = output_format.format(coord=coord, albedo=albedo)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux', 'coords').split(',')
        albedos = CONFIG.get('extract_flux', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('extract_flux', 'input_format')
        input_format = pjoin(modtran_root, input_format)
        output_format = CONFIG.get('extract_flux', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))

        gaip.extract_flux(coords, albedos, input_format, output_format,
                          satfilter)


class ExtractFluxTrans(luigi.Task):

    """Extract the flux data from the MODTRAN output in the transmissive
    case."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux_trans', 'coords').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('extract_flux_trans', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        targets = []
        for coord in coords:
            target = output_format.format(coord=coord)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux_trans', 'coords').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('extract_flux_trans', 'input_format')
        input_format = pjoin(modtran_root, input_format)
        output_format = CONFIG.get('extract_flux_trans', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))

        gaip.extract_flux_trans(coords, input_format, output_format,
                                satfilter)


class CalculateCoefficients(luigi.Task):

    """Calculate the atmospheric parameters needed by BRDF and atmospheric
    correction model."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [ExtractFlux(self.l1t_path, self.out_path),
                ExtractFluxTrans(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('coefficients', 'coords').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('coefficients', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        targets = []
        for coord in coords:
            target = output_format.format(coord=coord)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('coefficients', 'coords').split(',')
        chn_input_format = CONFIG.get('coefficients', 'chn_input_format')
        dir_input_format = CONFIG.get('coefficients', 'dir_input_format')
        output_format = CONFIG.get('coefficients', 'output_format')
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))

        gaip.calc_coefficients(coords, chn_input_format, dir_input_format,
                               output_format, satfilter, workpath)


class ReformatAtmosphericParameters(luigi.Task):

    """Reformat the atmospheric parameters produced by MODTRAN for four boxes.
    These are needed to conduct bilinear interpolation. This runs the binary
    `read_modtran`. """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateCoefficients(self.l1t_path, self.out_path),
                CreateSatelliteFilterFile(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        factors = CONFIG.get('read_modtran', 'factors').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('read_modtran', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        bands = [a.band_num for a in acqs]
        targets = []
        for factor in factors:
            for band in bands:
                if band not in bands_to_process:
                    # Skip
                    continue
                target = output_format.format(factor=factor, band=band)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('read_modtran', 'coords').split(',')
        factors = CONFIG.get('read_modtran', 'factors').split(',')
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('read_modtran', 'input_format')
        input_format = pjoin(workpath, input_format)
        output_format = CONFIG.get('read_modtran', 'output_format')
        output_format = pjoin(workpath, output_format)
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))

        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # Initialise the list to contain the acquisitions we wish to process
        acqs_to_process = []
        for acq in acqs:
            band_number = acq.band_num
            if band_number in bands_to_process:
                acqs_to_process.append(acq)

        gaip.reformat_atmo_params(acqs_to_process, coords, satfilter, factors,
                                  input_format, output_format, workpath)


class BilinearInterpolation(luigi.Task):

    """Perform the bilinear interpolation."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [ReformatAtmosphericParameters(self.l1t_path, self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        factors = CONFIG.get('bilinear', 'factors').split(',')
        output_format = CONFIG.get('bilinear', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        bands = [a.band_num for a in acqs]
        targets = []
        target = pjoin(out_path,
                       CONFIG.get('work', 'bilinear_outputs_target'))
        targets.append(luigi.LocalTarget(target))
        for factor in factors:
            for band in bands:
                if band not in bands_to_process:
                    # Skip
                    continue
                target = output_format.format(factor=factor, band=band)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        factors = CONFIG.get('bilinear', 'factors').split(',')
        coordinator = pjoin(out_path,
                            CONFIG.get('work', 'coordinator_target'))
        boxline = pjoin(out_path,
                        CONFIG.get('work', 'boxline_target'))
        centreline = pjoin(out_path,
                           CONFIG.get('work', 'centreline_target'))
        input_format = CONFIG.get('bilinear', 'input_format')
        output_format = CONFIG.get('bilinear', 'output_format')
        workpath = pjoin(out_path,
                         CONFIG.get('work', 'modtran_root'))

        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # Initialise the list to contain the acquisitions we wish to process
        acqs_to_process = []
        for acq in acqs:
            band_number = acq.band_num
            if band_number in bands_to_process:
                acqs_to_process.append(acq)

        bilinear_fnames = gaip.bilinear_interpolate(acqs_to_process, factors,
                                                    coordinator, boxline,
                                                    centreline, input_format,
                                                    output_format, workpath)

        save(self.output()[0], bilinear_fnames)


class CreateTCRflDirs(luigi.Task):
    """
    Setup the directories to contain the Intermediate files
    produced for terrain corection.
    """

    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        rfl_path = pjoin(out_path, CONFIG.get('work', 'rfl_output_dir'))

        targets = [luigi.LocalTarget(tc_path), luigi.LocalTarget(rfl_path)]
        return targets

    def run(self):
        out_path = self.out_path
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        rfl_path = pjoin(out_path, CONFIG.get('work', 'rfl_output_dir'))
        if not exists(tc_path):
            os.makedirs(tc_path)
        if not exists(rfl_path):
            os.makedirs(rfl_path)

class DEMExctraction(luigi.Task):

    """
    Extract the DEM covering the acquisition extents plus an
    arbitrary buffer. The subset is then smoothed with a gaussian
    filter.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CreateTCRflDirs(self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        subset_target = pjoin(work_path,
                              CONFIG.get('extract_dsm', 'dsm_subset'))
        smoothed_target = pjoin(work_path,
                                CONFIG.get('extract_dsm', 'dsm_smooth_subset'))
        targets = [luigi.LocalTarget(subset_target),
                   luigi.LocalTarget(smoothed_target)]
        return targets

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        national_dsm = CONFIG.get('ancillary', 'dem_tc')
        subset_target = CONFIG.get('extract_dsm', 'dsm_subset')
        smoothed_target = CONFIG.get('extract_dsm', 'dsm_smooth_subset')
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        dsm_subset_fname = pjoin(work_path, subset_target)
        dsm_subset_smooth_fname = pjoin(work_path, smoothed_target)

        gaip.get_dsm(acqs[0], national_dsm, buffer, dsm_subset_fname,
                     dsm_subset_smooth_fname)


class SlopeAndSelfShadow(luigi.Task):

    """
    Compute the slope, aspect, incident, azimuth incident, exiting,
    azimuth exiting, relative slope angles, and a self shadow mask.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # These could've been under the work config, but i thought it might be
        # an ok idea to separate these targets into their own section
        slope_target = pjoin(work_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(work_path,
                              CONFIG.get('self_shadow', 'aspect_target'))
        incident_target = pjoin(work_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        azi_incident_target = pjoin(work_path,
                                    CONFIG.get('self_shadow',
                                               'azimuth_incident_target'))
        exiting_target = pjoin(work_path,
                               CONFIG.get('self_shadow',
                                          'exiting_target'))
        azi_exiting_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_exiting_target'))
        relative_slope_target = pjoin(work_path,
                                      CONFIG.get('self_shadow',
                                                 'relative_slope_target'))
        self_shadow_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'self_shadow_target'))

        targets = [luigi.LocalTarget(self_shadow_target),
                   luigi.LocalTarget(slope_target),
                   luigi.LocalTarget(aspect_target),
                   luigi.LocalTarget(incident_target),
                   luigi.LocalTarget(exiting_target),
                   luigi.LocalTarget(azi_incident_target),
                   luigi.LocalTarget(azi_exiting_target),
                   luigi.LocalTarget(relative_slope_target)]
        return targets

    def run(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        acqs = gaip.acquisitions(self.l1t_path)

        # Input targets
        satellite_view_fname = pjoin(out_path,
                                     CONFIG.get('work', 'sat_view_target'))
        satellite_azimuth_fname = pjoin(out_path,
                                        CONFIG.get('work',
                                                   'sat_azimuth_target'))
        solar_zenith_fname = pjoin(out_path,
                                   CONFIG.get('work', 'solar_zenith_target'))
        solar_azimuth_fname = pjoin(out_path,
                                    CONFIG.get('work',
                                               'solar_azimuth_target'))
        smoothed_dsm_fname = pjoin(work_path, CONFIG.get('extract_dsm',
                                                         'dsm_smooth_subset'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))

        # Output targets
        self_shadow_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'self_shadow_target'))
        slope_target = pjoin(work_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(work_path,
                              CONFIG.get('self_shadow', 'aspect_target'))
        incident_target = pjoin(work_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        exiting_target = pjoin(work_path,
                               CONFIG.get('self_shadow', 'exiting_target'))
        azi_incident_target = pjoin(work_path,
                                    CONFIG.get('self_shadow',
                                               'azimuth_incident_target'))
        azi_exiting_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_exiting_target'))
        relative_slope_target = pjoin(work_path,
                                      CONFIG.get('self_shadow',
                                                 'relative_slope_target'))
        header_slope_target = pjoin(work_path,
                                    CONFIG.get('work', 'header_slope_target'))

        out_targets = [self_shadow_target, slope_target, aspect_target,
                       incident_target, exiting_target, azi_incident_target,
                       azi_exiting_target, relative_slope_target]

        gaip.calculate_self_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   solar_zenith_fname, solar_azimuth_fname,
                                   satellite_view_fname,
                                   satellite_azimuth_fname, out_targets,
                                   header_slope_target)


class CalculateCastShadow(luigi.Task):

    """Calculate cast shadow masks. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateCastShadowSun(self.l1t_path, self.out_path),
                CalculateCastShadowSatellite(self.l1t_path, self.out_path)]

    def complete(self):
        return all([t.complete() for t in self.requires()])

class CalculateCastShadowSun(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path)]

    def output(self):
        out_path = pjoin(self.out_path,
                         CONFIG.get('work', 'tc_intermediates'))
        sun_target = pjoin(out_path,
                           CONFIG.get('cast_shadow', 'sun_direction_target'))

        target = luigi.LocalTarget(sun_target)

        return target

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        tc_work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        smoothed_dsm_fname = pjoin(tc_work_path,
                                   CONFIG.get('extract_dsm',
                                              'dsm_smooth_subset'))
        solar_zenith_target = pjoin(out_path,
                                    CONFIG.get('work', 'solar_zenith_target'))
        solar_azimuth_target = pjoin(out_path,
                                     CONFIG.get('work',
                                                'solar_azimuth_target'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        window_height = int(CONFIG.get('terrain_correction',
                                       'shadow_sub_matrix_height'))
        window_width = int(CONFIG.get('terrain_correction',
                                      'shadow_sub_matrix_width'))

        # Output targets
        sun_target = pjoin(tc_work_path,
                           CONFIG.get('cast_shadow', 'sun_direction_target'))

        gaip.calculate_cast_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   window_height, window_width,
                                   solar_zenith_target, solar_azimuth_target,
                                   sun_target)


class CalculateCastShadowSatellite(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path)]

    def output(self):
        out_path = pjoin(self.out_path,
                         CONFIG.get('work', 'tc_intermediates'))
        satellite_target = pjoin(out_path,
                                 CONFIG.get('cast_shadow',
                                            'satellite_direction_target'))
        target = luigi.LocalTarget(satellite_target)

        return target

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        tc_work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        smoothed_dsm_fname = pjoin(tc_work_path,
                                   CONFIG.get('extract_dsm',
                                              'dsm_smooth_subset'))
        satellite_view_target = pjoin(out_path,
                                      CONFIG.get('work', 'sat_view_target'))
        satellite_azimuth_target = pjoin(out_path,
                                         CONFIG.get('work',
                                                    'sat_azimuth_target'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        window_height = int(CONFIG.get('terrain_correction',
                                       'shadow_sub_matrix_height'))
        window_width = int(CONFIG.get('terrain_correction',
                                      'shadow_sub_matrix_width'))

        # Output targets
        satellite_target = pjoin(tc_work_path,
                                 CONFIG.get('cast_shadow',
                                            'satellite_direction_target'))

        gaip.calculate_cast_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   window_height, window_width,
                                   satellite_view_target,
                                   satellite_azimuth_target, satellite_target)


class TerrainCorrection(luigi.Task):

    """Perform the terrain correction."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [BilinearInterpolation(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path),
                SlopeAndSelfShadow(self.l1t_path, self.out_path),
                CalculateCastShadow(self.l1t_path, self.out_path),
                CreateModisBrdfFiles(self.l1t_path, self.out_path)]

    def output(self):
        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # Get the reflectance levels and base output format
        rfl_levels = CONFIG.get('terrain_correction', 'rfl_levels').split(',')
        output_format = CONFIG.get('terrain_correction', 'output_format')

        # Output directory
        out_path = pjoin(self.out_path,
                         CONFIG.get('work', 'rfl_output_dir'))

        # Create the targets
        targets = []
        for level in rfl_levels:
            for band in bands_to_process:
                target = pjoin(out_path,
                               output_format.format(level=level, band=band))
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path

        # Get the necessary config params
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        outdir = pjoin(out_path, CONFIG.get('work', 'rfl_output_dir'))
        bilinear_target = pjoin(out_path, 
                                CONFIG.get('work', 'bilinear_outputs_target'))
        bilinear_target = load_value(bilinear_target)
        rori = float(CONFIG.get('terrain_correction', 'rori'))
        modis_brdf_format = pjoin(out_path,
                                  CONFIG.get('brdf', 'modis_brdf_format'))
        new_modis_brdf_format = pjoin(tc_path,
                                      CONFIG.get('brdf',
                                                 'new_modis_brdf_format'))

        # Get the reflectance levels and base output format
        rfl_levels = CONFIG.get('terrain_correction', 'rfl_levels').split(',')
        output_format = CONFIG.get('terrain_correction', 'output_format')

        # Input targets (images)
        self_shadow_target = pjoin(tc_path,
                                   CONFIG.get('self_shadow',
                                              'self_shadow_target'))
        slope_target = pjoin(tc_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(tc_path,
                              CONFIG.get('self_shadow', 'aspect_target'))
        incident_target = pjoin(tc_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        exiting_target = pjoin(tc_path,
                               CONFIG.get('self_shadow', 'exiting_target'))
        relative_slope_target = pjoin(tc_path,
                                      CONFIG.get('self_shadow',
                                                 'relative_slope_target'))
        sun_target = pjoin(tc_path,
                           CONFIG.get('cast_shadow', 'sun_direction_target'))
        satellite_target = pjoin(tc_path,
                                 CONFIG.get('cast_shadow',
                                            'satellite_direction_target'))
        solar_zenith_target = pjoin(out_path,
                                    CONFIG.get('work', 'solar_zenith_target'))
        solar_azimuth_target = pjoin(out_path,
                                     CONFIG.get('work',
                                                'solar_azimuth_target'))
        satellite_view_target = pjoin(out_path,
                                      CONFIG.get('work', 'sat_view_target'))
        relative_angle_target = pjoin(out_path,
                                      CONFIG.get('work',
                                                 'relative_azimuth_target'))

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # Initialise the list to contain the acquisitions we wish to process
        acqs_to_process = []
        for acq in acqs:
            band_number = acq.band_num
            if band_number in bands_to_process:
                acqs_to_process.append(acq)

        # Output targets
        # Create a dict of filenames per reflectance level per band
        rfl_lvl_fnames = {}
        for level in rfl_levels:
            for band in bands_to_process:
                outfname = output_format.format(level=level, band=band)
                rfl_lvl_fnames[(band, level)] = pjoin(outdir, outfname)

        gaip.calculate_reflectance(acqs_to_process, bilinear_target, rori,
                                   self_shadow_target, sun_target,
                                   satellite_target, solar_zenith_target,
                                   solar_azimuth_target, satellite_view_target,
                                   relative_angle_target, slope_target,
                                   aspect_target, incident_target,
                                   exiting_target, relative_slope_target,
                                   rfl_lvl_fnames, modis_brdf_format,
                                   new_modis_brdf_format)


if __name__ == '__main__':  # FIXME
    l1t_path = '/g/data1/v10/NBAR_validation_reference/Nov2013/L1T_Input/LS7_90-84_2000-09-13/UTM/LS7_ETM_OTH_P51_GALPGS01-002_090_084_20000913'
    out_path = '/g/data1/v10/testing_ground/jps547/test_gaip/test_build'
    luigi.build([TerrainCorrection(l1t_path, out_path)],
                local_scheduler=True)
