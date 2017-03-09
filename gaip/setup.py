"""
Setup
-----

This compiles all the Fortran extensions.
"""


from __future__ import absolute_import
def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('gaip', parent_package, top_path)
    config.add_data_files('sensors.json')
    config.add_data_dir('spectral_response')
    config.add_extension(
        '__cast_shadow_mask',
        [
            'sys_variables.f90',
            'cast_shadow_main.f90',
            'terrain_border_margins.f90',
            'cast_shadow_mask.f90',
            'terrain_occlusion.f90',
            'geo2metres_pixel_size.f90',
        ]
    ),
    config.add_extension(
        '__exiting_angle',
        [
            'sys_variables.f90',
            'exiting_angle.f90',
            'earth_rotation.f90',
        ]
    ),
    config.add_extension(
        '__incident_angle',
        [
            'incident_angle.f90',
            'earth_rotation.f90',
        ]
    ),
    config.add_extension(
        '__slope_aspect',
        [
            'sys_variables.f90',
            'slope_aspect.f90',
            'geo2metres_pixel_size.f90',
        ]
    ),
    config.add_extension(
        '__surface_reflectance',
        [
            'surface_reflectance.f90',
            'white_sky.f90',
            'black_sky.f90',
            'brdf_shape.f90',
        ]
    ),
    config.add_extension(
        '__satellite_model',
        [
            'sys_variables.f90',
            'geo2metres_pixel_size.f90',
            'satellite_model.f90',
        ]
    ),
    config.add_extension(
        '__track_time_info',
        [
            'sys_variables.f90',
            'geod2geo.f90',
            'q_cal.f90',
            'geo2metres_pixel_size.f90',
            'satellite_track.f90',
            'track_time_info.f90',
        ]
    ),
    config.add_extension(
        '__sat_sol_angles',
        [
            'sys_variables.f90',
            'solar_angle.f90',
            'geod2geo.f90',
            'q_cal.f90',
            'compute_angles.f90',
            'satellite_solar_angles_main.f90',
        ]
    ),
    config.add_extension(
        '__interpolation',
        [
            'bilinear_interpolation.f90',
        ]
    )
    config.add_extension(
        'unmiximage',
        [
            'unmiximage.f90',
            'constants_NSWC.f90',
            'nnls.f90',
            'unmiximage.pyf',
        ]
    )

    return config
