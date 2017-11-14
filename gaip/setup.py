"""
Setup
-----

This compiles all the Fortran extensions.
"""


from __future__ import absolute_import

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('gaip', parent_package, top_path)
    config.add_subpackage('scripts')
    config.add_data_files('sensors.json')
    config.add_data_dir('spectral_response')
    config.add_extension(
        '__cast_shadow_mask',
        [
            'f90_sources/sys_variables.f90',
            'f90_sources/cast_shadow_main.f90',
            'f90_sources/terrain_border_margins.f90',
            'f90_sources/cast_shadow_mask.f90',
            'f90_sources/terrain_occlusion.f90',
            'f90_sources/geo2metres_pixel_size.f90',
        ]
    ),
    config.add_extension(
        '__exiting_angle',
        [
            'f90_sources/sys_variables.f90',
            'f90_sources/exiting_angle.f90',
            'f90_sources/earth_rotation.f90',
        ]
    ),
    config.add_extension(
        '__incident_angle',
        [
            'f90_sources/incident_angle.f90',
            'f90_sources/earth_rotation.f90',
        ]
    ),
    config.add_extension(
        '__slope_aspect',
        [
            'f90_sources/sys_variables.f90',
            'f90_sources/slope_aspect.f90',
            'f90_sources/geo2metres_pixel_size.f90',
        ]
    ),
    config.add_extension(
        '__surface_reflectance',
        [
            'f90_sources/surface_reflectance.f90',
            'f90_sources/white_sky.f90',
            'f90_sources/black_sky.f90',
            'f90_sources/brdf_shape.f90',
        ]
    ),
    config.add_extension(
        '__satellite_model',
        [
            'f90_sources/sys_variables.f90',
            'f90_sources/geo2metres_pixel_size.f90',
            'f90_sources/satellite_model.f90',
        ]
    ),
    config.add_extension(
        '__track_time_info',
        [
            'f90_sources/sys_variables.f90',
            'f90_sources/geod2geo.f90',
            'f90_sources/q_cal.f90',
            'f90_sources/geo2metres_pixel_size.f90',
            'f90_sources/satellite_track.f90',
            'f90_sources/track_time_info.f90',
        ]
    ),
    config.add_extension(
        '__sat_sol_angles',
        [
            'f90_sources/sys_variables.f90',
            'f90_sources/solar_angle.f90',
            'f90_sources/geod2geo.f90',
            'f90_sources/q_cal.f90',
            'f90_sources/compute_angles.f90',
            'f90_sources/satellite_solar_angles_main.f90',
        ]
    ),
    config.add_extension(
        '__bilinear_interpolation',
        [
            'f90_sources/bilinear_interpolation.f90',
        ]
    )

    return config
