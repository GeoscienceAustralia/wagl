all: filter.so _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so .fc

.fc:
	$(MAKE) --directory=fc --makefile=Makefile all

filter.so _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so:
	python setup.py build_ext -i

clean: rootclean fcclean

rootclean:
	-rm *.so

fcclean:
	$(MAKE) --directory=fc clean
