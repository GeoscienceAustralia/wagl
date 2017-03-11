"""
Setup gaip
"""

from __future__ import absolute_import

import versioneer
from numpy.distutils.core import setup


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
    )

    config.add_subpackage('gaip')
    return config


setup(
    name='gaip',
    configuration=configuration,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url='https://github.com/GeoscienceAustralia/gaip',
    maintainer='Josh Sixsmith',
    maintainer_email='joshua.sixsmith@ga.gov.au',
    scripts=['scripts/test_calculate_angles.py',
             'scripts/test_dsm.py',
             'scripts/test_exiting_angles.py',
             'scripts/test_incident_angles.py',
             'scripts/test_relative_slope.py',
             'scripts/test_shadow_masks.py',
             'scripts/test_slope_aspect.py']
)
