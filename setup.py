"""
Setup gaip
"""
from __future__ import absolute_import

import setuptools

import versioneer
from numpy.distutils.core import setup


tests_require = [
    'pytest',
]

install_requires = [
    'luigi>=2.4.0',
    'numpy>=1.8',
    'scipy>=0.14',
    'numexpr>=2.4.6',
    'ephem>=3.7.5.3',
    'pyproj>1.9.5',
    'scikit-image>=0.8.2',
    'GDAL>=1.9.2',
    'rasterio>0.9', # Hack to get the alpha release
    'fiona>=1.7.0',
    'shapely>=1.5.13',
    'h5py>=2.5.0',
    'tables>=3.4.2',
    'pandas>=0.17.1',
    'geopandas>=0.1.1',
    'pyyaml>=3.11',
    'nested_lookup>=0.1.3',
    'python-dateutil>=2.6.1',
    'structlog>=16.1.0',
    'idl-functions>=0.5.2', # custom package
]

dependency_links = [
    'git+git://github.com/sixy6e/idl-functions.git@master#egg=idl-functions-0.5.2',
]

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
    license='CC0 1.0 Universal',
    author='The gaip authors',
    maintainer='gaip developers',
    scripts=['utils/test_satellite_solar_angles',
             'utils/test_dsm',
             'utils/test_exiting_angles',
             'utils/test_incident_angles',
             'utils/test_relative_slope',
             'utils/test_terrain_shadow_masks',
             'utils/test_slope_aspect',
             'utils/aot_converter',
             'utils/gaip_convert',
             'utils/gaip_ls',
             'utils/gaip_residuals',
             'utils/gaip_pbs'],
    setup_requires=['pytest-runner'],
    tests_require=tests_require,
    install_requires=install_requires,
    dependency_links=dependency_links,
)
