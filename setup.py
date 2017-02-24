"""
Setup gaip
"""

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
    # config.get_version('gaip/version.py') # sets config.version
    return config


setup(
    name='gaip',
    configuration=configuration,
)
