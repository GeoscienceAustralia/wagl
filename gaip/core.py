"""
Core code.
"""
import xml.dom.minidom as xml

DATASETS_CONFIG = 'datasets.xml'


class Acquisition(object):

    def __init__(self, wave_length, sensor, satellite, time):
        self.wave_length = wave_length
        self.sensor = sensor
        self.satellite = satellite
        self.time = time


def get_acquisition(path, acquisition):
    """
    Get acquisition data.
    """
    pass


class Filetype(object):

    def __init__(self, ext, suffix, product):
        self.ext = ext
        self.suffix = suffix
        self.product = product

    def __repr__(self):
        return str((self.ext, self.suffix, self.product))


def acquisitions(path):
    """
    Return a list of available acquisitions.
    """
    dom = xml.parse(DATASETS_CONFIG)
    classes = dom.documentElement.getElementsByTagName('CLASS')[0]
    datas = classes.getElementsByTagName('DATA')[0]
    datadirs = datas.getElementsByTagName('DATADIR')
    paths = [n.getAttribute('PATH') for n in datadirs]

    filetypes = []
    for filetype in datas.getElementsByTagName('FILETYPE'):
        ext = filetype.getAttribute('EXTENSION')
        suffix = filetype.getAttribute('ROOTSUFFIX')
        product = filetype.getAttribute('PRODUCTFORMAT')
        filetypes.append(Filetype(ext, suffix, product))

    return filetypes


if __name__ == '__main__':
    import pprint
    pprint.pprint(acquisitions(''))

