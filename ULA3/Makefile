all: filter.so _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so .fc sys_variables.o sys_variables.mod angle_all.so set_satmod.so set_times.so

.fc:
	$(MAKE) --directory=fc --makefile=Makefile all

filter.so _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so sys_variables.o sys_variables.mod angle_all.so set_satmod.so set_times.so:
	gfortran -c sys_variables.f90
	python setup.py build_ext -i

clean: rootclean fcclean

rootclean:
	-rm *.so
	-rm *.o
	-rm *.mod

fcclean:
	$(MAKE) --directory=fc clean
