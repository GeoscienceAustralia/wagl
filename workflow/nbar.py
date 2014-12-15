import luigi
import gaip
import os
import cPickle as pickle

from os.path import join as pjoin, dirname, exists

CONFIG = luigi.configuration.get_config()
CONFIG.add_config_path(pjoin(dirname(__file__), 'nbar.cfg'))


def save(target, value):
    with target.open('w') as outfile:
        if target.fn.endswith('pkl'):
            pickle.dump(value, outfile)
        else:
            print >>outfile, value


def load(target):
    if not target.fn.endswith('pkl'):
        raise IOError('Cannot load non-pickled object')
    with target.open('r') as infile:
        return pickle.load(infile)


def load_value(target)
    data = load(target)
    try:
        return data['value']
    except KeyError:
        return data


class GetElevationAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'dem_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        dem_path = CONFIG.get('ancillary', 'dem_path')
        value = gaip.get_elevation_data(geobox.centre_lonlat, dem_path)
        save(self.output(), value)


class GetOzoneAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'ozone_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        ozone_path = CONFIG.get('ancillary', 'ozone_path')
        centre = geobox.centre_lonlat
        dt = acqs[0].scene_center_datetime
        value = gaip.get_ozone_data(ozone_path, centre, dt)
        save(self.output(), value)


class GetSolarIrradianceAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'irrad_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        solar_path = CONFIG.get('ancillary', 'solarirrad_path')
        value = gaip.get_solar_irrad(acqs, solar_path)
        save(self.output(), value)


class GetSolarDistanceAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'sundist_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        sundist_path = CONFIG.get('ancillary', 'sundist_path')
        value = gaip.get_solar_dist(acqs[0], sundist_path)
        save(self.output(), value)


class GetWaterVapourAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'vapour_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        vapour_path = CONFIG.get('ancillary', 'vapour_path')
        value = gaip.get_water_vapour(acqs[0], vapour_path)
        save(self.output(), value)


class GetAerosolAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'aerosol_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        aerosol_path = CONFIG.get('ancillary', 'aerosol_path')
        value = gaip.get_aerosol_data(acqs[0], aerosol_path)
        save(self.output(), value)


class GetBrdfAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'brdf_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        brdf_path = CONFIG.get('ancillary', 'brdf_path')
        brdf_premodis_path = CONFIG.get('ancillary', 'brdf_premodis_path')
        work_path = CONFIG.get('work', 'path')
        value = gaip.get_brdf_data(acqs[0], brdf_path, brdf_premodis_path,
                                   work_path)
        save(self.output(), value)


class GetAncillaryData(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [GetElevationAncillaryDataTask(self.l1t_path),
                GetOzoneAncillaryDataTask(self.l1t_path),
                GetSolarDistanceAncillaryDataTask(self.l1t_path),
                GetWaterVapourAncillaryDataTask(self.l1t_path),
                GetAerosolAncillaryDataTask(self.l1t_path),
                GetBrdfAncillaryDataTask(self.l1t_path)]

    def complete(self):
        return all([d.complete() for d in self.requires()])


class CalculateLonGrid(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'lon_grid_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output()
        gaip.create_lon_grid(acqs[0], target.fn)


class CalculateLatGrid(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'lat_grid_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output()
        gaip.create_lat_grid(acqs[0], target.fn)


class CalculateLatLonGrids(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatGrid(self.l1t_path),
                CalculateLonGrid(self.l1t_path)]

    def complete(self):
        #FIXME
        return all([d.complete() for d in self.requires()])


class CalculateSatelliteAndSolarGrids(luigi.Task):

    l1t_path = luigi.Parameter()
    targets = None

    def requires(self):
        return [CalculateLatGrid(self.l1t_path),
                CalculateLonGrid(self.l1t_path)]

    def run(self):
        self.targets = [
            CONFIG.get('work', 'sat_view_zenith_target'),
            CONFIG.get('work', 'sat_azimuth_target'),
            CONFIG.get('work', 'solar_zenith_target'),
            CONFIG.get('work', 'solar_azimuth_target'),
            CONFIG.get('work', 'relative_azimuth_target'),
            CONFIG.get('work', 'time_target')]
        work_path = CONFIG.get('work', 'path')
        lon_target = CONFIG.get('work', 'lon_target')
        lat_target = CONFIG.get('work', 'lat_target')

        acqs = gaip.acquisitions(self.l1t_path)

        geobox = acqs[0].gridded_geo_box()
        cols = acqs[0].samples
        rows = acqs[0].lines

        (satellite_zenith, satellite_azimuth, solar_zenith, solar_azimuth,
         relative_azimuth, time, Y_cent, X_cent, N_cent) = \
            gaip.calculate_angles(acqs[0], lon_target, lat_target,
                                  npoints=12, to_disk=self.targets)

        gaip.create_centreline_file(geobox, Y_cent, X_cent, N_cent, cols,
                                    view_max=9.0, outdir=work_path)

    def complete(self):
        #FIXME
        if self.targets:
            return all([exists(t) for t in self.targets])
        else:
            return False


class CalculateGridsTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatLonGrids(self.l1t_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path)]

    def complete(self):
        #FIXME
        return all([d.complete() for d in self.requires()])


class CreateModtranDirectories(luigi.Task):

    created = None

    def run(self):
        workpath = CONFIG.get('work', 'path')
        modtran_root = CONFIG.get('modtran', 'root')
        modtran_exe = CONFIG.get('modtran', 'exe')
        self.created = gaip.create_modtran_dirs(workpath, modtran_root)

    def complete(self):
        #FIXME
        if self.created:
            return all([exists(p) for p in self.created])
        else:
            return False


class CreateSatelliteFilterFile(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        target = CONFIG.get('work', 'sat_filter_target')
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        satfilterpath = CONFIG.get('ancillary', 'satfilter_path')
        target = CONFIG.get('work', 'sat_filter_target')
        self.created = gaip.create_satellite_filter_file(acqs, satfilterpath,
                                                         target)


class WriteModtranInputFile(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [GetAncillaryData(self.l1t_path)]

    def output(self):
        target = CONFIG.get('work', 'modtran_input_target')
        return luigi.LocalTarget(target)

    def run(self):
        ozone_target = CONFIG.get('work', 'ozone_target')
        vapour_target = CONFIG.get('work', 'vapour_target')
        aerosol_target = CONFIG.get('work', 'aerosol_target')
        elevation_target = CONFIG.get('work', 'elevation_target')
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output().fn
        ozone = load_value(ozone_target)
        vapour = load_value(vapour_target)
        aerosol = load_value(aerosol_target)
        elevation = load_value(elevation_target)
        gaip.write_modtran_input(acqs, target, ozone, vapour, aerosol,
                                 elevation)


class WriteModisBrdfFiles(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires():
        return []

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        modis_brdf_prefix = CONFIG.get('work', 'modis_brdf_prefix')
        brdf_target = CONFIG.get('work', 'brdf_target')
        brdf_data = load(brdf_target)
        # TODO
        gaip.write_modis_brdf_files(acqs, prefix, brdf_data)


class RunModtranCorOrtho(luigi.Task):

    l1t_path = luigi.Parameter()
    created = None

    def requires():
        return [CalculateSatelliteAndSolarGrids(self.l1t_path)]

    def run(self):
        # sources
        centreline_target = CONFIG.get('work', 'centreline_target')
        sat_view_zenith_target = CONFIG.get('work', 'sat_view_zenith_target')
        # targets
        coordinator_target = CONFIG.get('work', 'coordinator_target')
        boxline_target = CONFIG.get('work', 'boxline_target')

        cwd = CONFIG.get('work', 'read_modtrancor_ortho_cwd')
        gaip.run_read_modtrancor_ortho([centreline_target,
                                        sat_view_zenith_target,
                                        coordinator_target,
                                        boxline_target],
                                       cwd=cwd)

        self.created = [coordinator_target, boxline_target]

    def complete(self):
        if self.created:
            return all([exists(p) for p in self.created])
        else:
            return False


class GenerateModtranInputFiles(luigi.task):

    """Generate the MODTRAN input files by running the Fortran binary
    `input_modtran_ortho_ula`."""

    l1t_path = luigi.Parameter()
    created = None

    def requires():
        return [RunModtranCorOrtho(self.l1t_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path),
                WriteModtranInputFile(self.l1t_path),
                CalculateLatGrid(self.l1t_path),
                CalculateLonGrid(self.l1t_path)]

    def run(self):
        # sources
        modtran_input_target = CONFIG.get('work', 'modtran_input_target')
        coordinator_target = CONFIG.get('work', 'coordinator_target')
        sat_view_zenith_target = CONFIG.get('work', 'sat_view_zenith_target')
        sat_azimuth_target = CONFIG.get('work', 'sat_azimuth_target')
        lon_grid_target = CONFIG.get('work', 'lon_grid_target')
        lat_grid_target = CONFIG.get('work', 'lat_grid_target')

        section = 'input_modtran'
        coords = CONFIG.get(section, 'coords').split(',')
        albedos = CONFIG.get(section, 'albedos').split(',')
        fname_format = CONFIG.get(section, 'filename_format')
        workdir = CONFIG.get('work', 'input_modtran_cwd')

        self.created = gaip.generate_modtran_inputs(modtran_input_target,
                                                    coordinator_target,
                                                    sat_view_zenith_target,
                                                    sat_azimuth_target,
                                                    lon_grid_target,
                                                    lat_grid_target,
                                                    coords,
                                                    albedos,
                                                    fname_format,
                                                    workdir)

    def complete(self):
        #FIXME
        if self.created:
            return all([exists(p) for p in self.created])
        else:
            return False


class ReformatAsTp5(luigi.task):

    """Reformat the MODTRAN input files in `tp5` format. This runs the
    Fortran binary `refort_tp5_ga` multiple times."""

    l1t_path = luigi.Parameter()
    created = None

    def requires():
        return [GenerateModtranInputFiles(self.l1t_path)]

    def run(self):
        modtran_profile_path = CONFIG.get('ancillary', 'modtran_profile_path')
        profile_format = CONFIG.get('modtran', 'profile_format')
        input_format = CONFIG.get('reformat_tp5', 'input_format')
        output_format = CONFIG.get('reformat_tp5', 'output_format')
        workdir = CONFIG.get('work', 'reformat_tp5_cwd')
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

        self.created = gaip.reformat_as_tp5(coords, albedos, profile,
                                            input_format, output_format,
                                            workdir)

    def complete(self):
        if self.created:
            return all([exists(p) for p in self.created])
        else:
            return False


class ReformatAsTp5Trans(luigi.task):

    """Reformat the MODTRAN input files in `tp5` format in the transmissive
    case. This runs the Fortran binary `refort_tp5_ga_trans` multiple
    times."""

    l1t_path = luigi.Parameter()

    def requires():
        return [GenerateModtranInputFiles(self.l1t_path)]

    def output(self):
        coords = CONFIG.get('reformat_tp5_trans', 'coords').split(',')
        albedos = CONFIG.get('reformat_tp5_trans', 'albedos').split(',')
        output_format = CONFIG.get('reformat_tp5_trans', 'output_format')
        targets = []
        for coord in coords:
            for albedo in albedos:
                target = output_format.format(coord=coord, albedo=albedo)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        modtran_profile_path = CONFIG.get('ancillary', 'modtran_profile_path')
        profile_format = CONFIG.get('modtran', 'profile_format')
        input_format = CONFIG.get('reformat_tp5_trans', 'input_format')
        output_format = CONFIG.get('reformat_tp5_trans', 'output_format')
        workdir = CONFIG.get('work', 'reformat_tp5_trans_cwd')
        coords = CONFIG.get('reformat_tp5_trans', 'coords').split(',')

        # determine modtran profile
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        centre_lon, centre_lat = geobox.centre_lonlat
        profile = 'tropical'
        if centre_lat < -23.0:
            profile = 'midlat_summer'

        profile = pjoin(modtran_profile_path,
                        profile_format.format(profile=profile))

        self.created = gaip.reformat_as_tp5_trans(coords, profile,
                                                  input_format, output_format,
                                                  workdir)


class PrepareModtranInput(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [CreateModtranDirectories(),
                CreateSatelliteFilterFile(self.l1t_path),
                GenerateModtranInputFiles(self.l1t_path),
                ReformatAsTp5(self.l1t_path),
                ReformatAsTp5Trans(self.l1t_path)]


class RunModtranCase(luigi.Task):

    l1t_path = luigi.Parameter()
    coord = luigi.Parameter()
    albedo = luigi.Parameter()

    def requires(self):
        return [PrepareModtranInput(self.l1t_path)]

    def output(self):
        flx_format = CONFIG.get('modtran', 'flx_output_format')
        chn_format = CONFIG.get('modtran', 'chn_output_format')
        flx_target = flux_format.format(coord=self.coord, albedo=self.albedo)
        chn_target = coef_format.format(coord=self.coord, albedo=self.albedo)
        return [luigi.LocalTarget(flx_target),
                luigi.LocalTarget(chn_target)]

    def run(self):
        modtran_exe = CONFIG.get('modtran', 'exe')
        workpath_format = CONFIG.get('modtran', 'workpath_format')
        workpath = workdir_format.format(coord=self.coord, albedo=self.albedo)
        gaip.run_modtran(modtran_exe, workpath)


class RunModtran(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        reqs = [PrepareModtranInput(self.l1t_path)]
        for coord in coords:
            for albedo in albedos:
                reqs.append(RunModtranCase(self.l1t_path, coord, albedo))
        return reqs


class ExtractFlux(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.l1t_path)]

    def output(self):
        coords = CONFIG.get('extract_flux', 'coords').split(',')
        albedos = CONFIG.get('extract_flux', 'albedos').split(',')
        output_format = CONFIG.get('extract_flux', 'output_format')
        targets = []
        for coord in coords:
            for albedo in albedos:
                target = output_format.format(coord=coord, albedo=albedo)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        coords = CONFIG.get('extract_flux', 'coords').split(',')
        albedos = CONFIG.get('extract_flux', 'albedos').split(',')
        input_format = CONFIG.get('extract_flux', 'input_format')
        output_format = CONFIG.get('extract_flux', 'output_format')
        satfilter = CONFIG.get('work', 'satfilter_target')

        gaip.extract_flux(coords, albedos, input_format, output_format,
                          satfilter)


class ExtractFluxTrans(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.l1t_path)]

    def output(self):
        coords = CONFIG.get('extract_flux_trans', 'coords').split(',')
        output_format = CONFIG.get('extract_flux_trans', 'output_format')
        targets = []
        for coord in coords:
            target = output_format.format(coord=coord)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        coords = CONFIG.get('extract_flux_trans', 'coords').split(',')
        input_format = CONFIG.get('extract_flux_trans', 'input_format')
        output_format = CONFIG.get('extract_flux_trans', 'output_format')
        satfilter = CONFIG.get('work', 'satfilter_target')

        gaip.extract_flux_trans(coords, input_format, output_format,
                                satfilter)


class CalculateCoefficients(luigi.Task):

    """Calculate the atmospheric parameters needed by BRDF and atmospheric
    correction model."""

    l1t_path = luigi.Parameter()

    def requires(self):
        return [ExtractFlux(self.l1t_path), ExtractFluxTrans(self.l1t_path)]

    def output(self):
        coords = CONFIG.get('coefficients', 'coords').split(',')
        output_format = CONFIG.get('coefficients', 'output_format')
        targets = []
        for coord in coords:
            target = output_format.format(coord=coord)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        coords = CONFIG.get('coefficients', 'coords').split(',')
        chn_input_format = CONFIG.get('coefficients', 'chn_input_format')
        dir_input_format = CONFIG.get('coefficients', 'dir_input_format')
        output_format = CONFIG.get('coefficients', 'output_format')
        satfilter = CONFIG.get('work', 'satfilter_target')
        workpath = CONFIG.get('work', 'modtran_root')

        gaip.calc_coefficients(coords, chn_input_format, dir_input_format,
                               output_format, satfilter, workpath)



class ReformatAtmosphericParameters(luigi.Task):

    """Reformat the atmospheric parameters produced by MODTRAN for four boxes.
    These are needed to conduct bilinear interpolation. This runs the binary
    `read_modtran`. """

    l1t_path = luigi.Parameter()

    def requires(self):
        return [CalculateCoefficients(self.l1t_path)]

    def output(self):
        factors = CONFIG.get('read_modtran', 'factors').split(',')
        output_format = CONFIG.get('read_modtran', 'output_format')
        acqs = gaip.acquisitions(self.l1t_path)
        bands = [str(a.band_num) for a in acqs]
        targets = []
        for factor in factors:
            for band in bands:
                target = output_format.format(factor=factor, band=band)
                targets.append(luigi.LocalTarget(target))
        return targets


    def run(self):
        coords = CONFIG.get('read_modtran', 'coords').split(',')
        factors = CONFIG.get('read_modtran', 'factors').split(',')
        input_format = CONFIG.get('read_modtran', 'input_format')
        output_format = CONFIG.get('read_modtran', 'output_format')
        workpath = CONFIG.get('work', 'modtran_root')

        acqs = gaip.acquisitions(self.l1t_path)

        gaip.reformat_atmo_params(acqs, coords, satfilter, factors,
                                  input_format, output_format, workpath)


class BilinearInterpolation(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return [ReformatAtmosphericParameters(self.l1t_path)]

    def output(self):
        factors = CONFIG.get('bilinear', 'factors').split(',')
        output_format = CONFIG.get('bilinear', 'output_format')
        acqs = gaip.acquisitions(self.l1t_path)
        bands = [str(a.band_num) for a in acqs]
        targets = []
        for factor in factors:
            for band in bands:
                target = output_format.format(factor=factor, band=band)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        factors = CONFIG.get('bilinear', 'factors').split(',')
        coordinator = CONFIG.get('work', 'coordinator_target')
        boxline = CONFIG.get('work', 'boxline_target')
        centreline = CONFIG.get('work', 'centreline_target')
        input_format = CONFIG.get('bilinear', 'input_format')
        output_format = CONFIG.get('bilinear', 'output_format')
        workpath = CONFIG.get('work', 'modtran_root')

        acqs = gaip.acquisitions(self.l1t_path)

        gaip.bilinear_interpolate(acqs, factors, coordinator, boxline,
                                  centreline, input_format, output_format,
                                  workpath)


class RadiativeTransferPrepTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class BrdfTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class WriteTiffFilesTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class RadiativeTransferNbarTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class RadiativeTransferPostprocessingTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


if __name__ == '__main__':
    l1t_path = '../gaip/tests/data/L1T/LS7_90-81_2009-04-15/UTM/LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415'
    luigi.build([GetElevationAncillaryDataTask(l1t_path)],
                local_scheduler=True)
