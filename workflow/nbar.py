import luigi
import logging

config = luigi.configuration.get_config()
config.add_config_path('workflow.cfg')

HOSTNAME = 'raijin'

log = logging.getLogger('root.' + __name__)


class Parameter(luigi.Parameter):

    def __init__(self, default=_no_value, is_list=False, is_boolean=False,
                 is_global=False, significant=True, description=None,
                 config_path=None, default_from_config=None, name=None):
        path = dict(section=HOSTNAME, name=name)
        super(Parameter, self).__init__(default=default,
                                        is_list=is_list,
                                        is_boolean=is_boolean,
                                        is_global=is_global,
                                        significant=significant,
                                        description=description,
                                        config_path=path,
                                        default_from_config=None)


class CheckNBAROutputTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass

    def is_complete(self):
        return True


class GetElevationAncillaryDataTask(luigi.Task):

    l1t_path = Parameter(name='l1t_path')

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        return luigi.LocalTarget('data/thresh_%s_%s.txt' % (cell, self.value))

    def run(self):
        from ULA3.image_processor import ProcessorConfig
        from ULA3.ancillary.elevation import get_elevation_data

        CONFIG = ProcessorConfig()
        l1t_input_dataset = DATA.get_item(
            CONFIG.input['l1t']['path'], SceneDataset)
        elevation_data = get_elevation_data(
            l1t_input_dataset.lonlats['CENTRE'], CONFIG.DIR_DEM)
        pass


class GetOzoneAncillaryDataTask(luigi.Task):

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        pass

    def run(self):
        pass


class GetSolarIrradienceAncillaryDataTask(luigi.Task):

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        pass

    def run(self):
        pass


class GetSolarDistanceAncillaryDataTask(luigi.Task):

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        pass

    def run(self):
        pass


class GetWaterVapourAncillaryDataTask(luigi.Task):

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        pass

    def run(self):
        pass


class GetAerosolAncillaryDataTask(luigi.Task):

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        pass

    def run(self):
        pass


class GetBrdfAncillaryDataTask(luigi.Task):

    def requires(self):
        return [CheckNBAROutputTask()]

    def output(self):
        pass

    def run(self):
        pass


class GetAncillaryData(luigi.Task):

    def requires(self):
        return [GetElevationAncillaryDataTask(),
                GetOzoneAncillaryDataTask(),
                GetSolarIrradienceAncillaryDataTask(),
                GetSolarDistanceAncillaryDataTask(),
                GetWaterVapourAncillaryDataTask(),
                GetAerosolAncillaryDataTask(),
                GetBrdfAncillaryDataTask()]

    def output(self):
        pass

    def run(self):
        pass


class CalculateLatLonGrids(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class CalculateSatelliteGrids(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class CalculateSolarGrids(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class CalculateGridsTask(luigi.Task):

    def requires(self):
        return [CalculateLatLonGrids(),
                CalculateSatelliteGrids(),
                CalculateSolarGrids()]

    def output(self):
        pass

    def run(self):
        pass


class PrepareModtranInputTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class RunModtranTask(luigi.Task):

    albedo = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class RunFluxTask(luigi.Task):

    albedo = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class RunCoefficientTask(luigi.Task):

    coef = luigi.Parameter()
    # TL, TM, TR, ML, MM, MR, BL, BM, BR

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class ReadModtranTask(luigi.Task):

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


class BilinearOrthoTask(luigi.Task):

    param = luigi.Parameter()  # fv, fs, b, s, a

    def requires(self):
        return []

    def output(self):
        pass

    def run(self):
        pass


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
