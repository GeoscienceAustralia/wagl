"""
Utilities to handle MTL files.
"""

import datetime
import re

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
    """Parse an MTL file and return dict-of-dict's containing the metadata."""

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

