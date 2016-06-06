"""
Setup
-----

This compiles all the Fortran extensions.
"""

from numpy.distutils.core import Extension, setup

setup(name='gaip',
      ext_modules=[Extension(name='_cast_shadow_mask',
                             sources=['cast_shadow_main.f90',
                                      'terrain_border_margins.f90',
                                      'cast_shadow_mask.f90',
                                      'terrain_occlusion.f90',
                                      'geo2metres_pixel_size.f90']),
                   Extension(name='_exiting_angle',
                             sources=['exiting_angle.f90',
                                      'earth_rotation.f90']),
                   Extension(name='_incident_angle',
                             sources=['incident_angle.f90',
                                      'earth_rotation.f90']),
                   Extension(name='_slope_aspect',
                             sources=['slope_aspect.f90',
                                      'geo2metres_pixel_size.f90']),
                   Extension(name='_surface_reflectance',
                             sources=['surface_reflectance.f90',
                                      'reflectance_mod.f90']),
                   Extension(name='_satellite_model',
                             sources=['geo2metres_pixel_size.f90',
                                      'satellite_model.f90']),
                   Extension(name='_track_time_info',
                             sources=['geod2geo.f90',
                                      'q_cal.f90',
                                      'geo2metres_pixel_size.f90',
                                      'satellite_track.f90',
                                      'track_time_info.f90']),
                   Extension(name='_sat_sol_angles',
                             sources=['solar_angle.f90',
                                      'geod2geo.f90',
                                      'q_cal.f90',
                                      'compute_angles.f90',
                                      'satellite_solar_angles_main.f90']),
                   Extension(name='_interpolation',
                             sources=['bilinear_interpolation.f90'])])
