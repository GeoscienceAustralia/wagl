"""
Core code.
"""
import os
import re
import copy
import datetime

from functools import total_ordering
from os.path import isdir, join as pjoin


@total_ordering
class Acquisition(object):

    """Acquisition information.
    """

    def __init__(self, metadata):
        for k, v in metadata.iteritems():
            self.__dict__.update(v)

    def __eq__(self, other):
        return self.band_name == other.band_name

    def __lt__(self, other):
        return self.band_name < other.band_name

    def __repr__(self):
        return 'Acquisition(band_name=' + self.band_name + ')'


def parse_type(s):
    """Parse the string `s` and return a native python object."""

    strptime = datetime.datetime.strptime

    def yesno(s):
        if len(s) == 1:
            if s == 'Y':
                return True
            if s == 'N':
                return False
        raise ValueError

    def none(s):
        if len(s) == 4 and s == 'NONE':
            return None
        raise ValueError

    parsers = [int,
               float,
               lambda x: strptime(x, '%Y-%m-%dT%H:%M:%SZ'),
               lambda x: strptime(x, '%Y-%m-%d').date(),
               lambda x: strptime(x[0:15], '%H:%M:%S.%f').time(),
               lambda x: yesno(x.strip('"')),
               lambda x: none(x.strip('"')),
               lambda x: str(x.strip('"'))]

    for parser in parsers:
        try:
            return parser(s)
        except ValueError:
            pass
    raise ValueError


def load_mtl(filename, root='L1_METADATA_FILE', pairs=r'(\w+)\s=\s(.*)'):
    """Parse an MTL file and return dict-of-dict's containing the metadata"""

    def parse(lines, tree, level=0):
        while lines:
            line = lines.pop(0)
            match = re.findall(pairs, line)
            if match:
                key, value = match[0]
                if key == 'GROUP':
                    tree[value] = {}
                    parse(lines, tree[value], level + 1)
                elif key == 'END_GROUP':
                    break
                else:
                    tree[key.lower()] = parse_type(value)

    tree = {}
    with open(filename, 'r') as fo:
        parse(fo.readlines(), tree)

    return tree[root]


def acquisitions(filename):
    """Obtain a list of Acquisition objects from an MTL file."""

    data = load_mtl(filename)
    bands = [k[:-10] for k in data['PRODUCT_METADATA'].keys() if
             k.startswith("band") and k.endswith('file_name')]

    # We now create an acquisition object for each band and make the
    # parameters names nice.

    acquisitions = []
    for band in bands:
        # create a new copy
        new = copy.deepcopy(data)

        # set band name and number
        new['PRODUCT_METADATA']['band_name'] = band
        new['PRODUCT_METADATA']['band_num'] = int(band.replace('band', ''))

        # replace filename
        product = new['PRODUCT_METADATA']
        filename = product[band + '_file_name']
        for k in product.keys():
            if k.startswith('band') and k.endswith('file_name'):
                product.pop(k, None)
            product['file_name'] = filename

        def replace_minmax(section_name, tag):
            section = new[section_name]
            minimum = section[tag + 'min_' + band]
            maximum = section[tag + 'max_' + band]
            for k in section.keys():
                if k.startswith(tag + 'min_') or k.startswith(tag + 'max_'):
                    section.pop(k, None)
            section[tag + 'min'] = minimum
            section[tag + 'max'] = maximum

        replace_minmax('MIN_MAX_RADIANCE', 'l')
        replace_minmax('MIN_MAX_PIXEL_VALUE', 'qcal')

        def replace_param(section_name):
            section = new[section_name]
            for k in section.keys():
                if k.endswith(band):
                    newkey = k.replace('_' + band, '')
                    section[newkey] = section[k]
                elif k.startswith(band):
                    newkey = k.replace(band + '_', '')
                    section[newkey] = section[k]
                else:
                    section.pop(k, None)

        replace_param('PRODUCT_PARAMETERS')
        replace_param('CORRECTIONS_APPLIED')

        acquisitions.append(Acquisition(new))

    return sorted(acquisitions)
