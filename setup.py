"""
Setup wagl
"""
from __future__ import absolute_import

import setuptools

from numpy.distutils.core import setup


tests_require = [
    "pytest",
]

install_requires = [
    "luigi>=2.7.3",
    "numpy>=1.8",
    "scipy>=0.14",
    "numexpr>=2.4.6",
    "ephem>=3.7.5.3",
    "pyproj>1.9.5",
    "scikit-image>=0.8.2",
    "GDAL>=1.9.2",
    "rasterio>1,!=1.0.3.post1,!=1.0.3",  # issue with /vsizip/ reader
    "fiona>=1.7.0",
    "shapely>=1.5.13",
    "h5py>=2.5.0",
    "tables>=3.4.2",
    "pandas>=0.17.1",
    "geopandas>=0.1.1",
    "pyyaml>=3.11",
    "nested_lookup>=0.1.3",
    "python-dateutil>=2.6.1",
    "structlog>=16.1.0",
    "idl-functions>=0.5.2",  # custom package
    "attrs>=17.4.0",
    "importlib-metadata;python_version<'3.8'",
]

dependency_links = [
    "git+git://github.com/sixy6e/idl-functions.git@master#egg=idl-functions-0.5.2",
]


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
    )

    config.add_subpackage("wagl")
    return config


setup(
    name="wagl",
    configuration=configuration,
    use_scm_version=True,
    url="https://github.com/GeoscienceAustralia/wagl",
    license="CC0 1.0 Universal",
    author="The wagl authors",
    author_email="earth.observation@ga.gov.au",
    maintainer="wagl developers",
    packages=setuptools.find_packages(exclude=("tests", "env_tests")),
    scripts=[
        "utils/test_satellite_solar_angles",
        "utils/test_dsm",
        "utils/test_exiting_angles",
        "utils/test_incident_angles",
        "utils/test_relative_slope",
        "utils/test_terrain_shadow_masks",
        "utils/test_slope_aspect",
        "utils/aot_converter",
        "utils/wagl_convert",
        "utils/wagl_ls",
        "utils/wagl_residuals",
        "utils/wagl_pbs",
    ],
    setup_requires=["pytest-runner", "setuptools_scm"],
    tests_require=tests_require,
    install_requires=install_requires,
    dependency_links=dependency_links,
)
