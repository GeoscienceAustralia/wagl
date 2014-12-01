import luigi
import gaip
import os

from os.path import join as pjoin, dirname

CONFIG = luigi.configuration.get_config()
CONFIG.add_config_path(pjoin(dirname(__file__), 'nbar.cfg'))

WORK_PATH = CONFIG.get('work', 'path')
DEM_PATH = CONFIG.get('ancillary', 'dem_path')


def save(target, value):
    with target().open('w') as outfile:
        print >> outfile, value


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

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(pjoin(WORK_PATH, 'elevation.txt'))

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        centre = (149.8356150617209, -34.610552280145946)  # FIXME
        value = gaip.get_elevation_data(centre, DEM_PATH)
        save(self.output, value)


class GetOzoneAncillaryDataTask(luigi.Task):

    l1t_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(pjoin(WORK_PATH, 'ozone.txt'))

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        centre = (149.8356150617209, -34.610552280145946)  # FIXME
        value = gaip.get_ozone_data(OZONE_PATH, centre, dt)
        save(self.output, value)


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


if __name__ == '__main__':
    l1t_path = '../gaip/tests/data/L1T/LS7_90-81_2009-04-15/UTM/LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415'
    luigi.build([GetElevationAncillaryDataTask(l1t_path)],
                local_scheduler=True)
