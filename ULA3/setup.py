from numpy.distutils.core import Extension, setup

setup(name='ULA3',
      ext_modules=[
          Extension(name='filter',
                    sources=['filter.f90']),
           Extension(name='_shade_main_landsat_pixel',
                    sources=['shade_main_landsat_pixel.f90']),
          Extension(name='_slope_pixelsize_newpole',
                    sources=['slope_pixelsize_newpole.f90']),
          Extension(name='_brdf_terrain_newdiff_all',
                    sources=['brdf_terrain_newdiff_all.f90'])
      ])
