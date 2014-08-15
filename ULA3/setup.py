from numpy.distutils.core import Extension, setup

setup(name='ULA3',
      ext_modules=[
          Extension(name='filter',
                    sources=['read_array_int32.f90',
                             'read_array_int16.f90',
                             'read_array_int8.f90',
                             'read_array_float32.f90',
                             'filter_float.f90']),
           Extension(name='_shade_main_landsat_pixel',
                    sources=['shade_main_landsat_pixel.f90',
                             'set_borderf.f90',
                             'get_proj_shadows.f90',
                             'proj_terrain.f90',
                             'pixelsize.f90']),
          Extension(name='_slope_pixelsize_newpole',
                    sources=['slope_pixelsize_newpole.f90',
                             'cal_pole.f90',
                             'pixelsize.f90']),
          Extension(name='_brdf_terrain_newdiff_all',
                    sources=['terrain_correction.f90',
                             'white_sky.f90',
                             'black_sky.f90',
                             'rl_brdf.f90'])
      ])
