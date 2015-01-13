FC := gfortran

all: filter.so _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so .fc sys_variables.o sys_variables.mod angle_all.so set_satmod.so set_times.so

filter.so _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so sys_variables.o sys_variables.mod angle_all.so set_satmod.so set_times.so:
	$$(module load gcc/4.9.0)
	$(FC) -c sys_variables.f90
	#python setup.py config_fc --f90flags="-assume byterecl" build_ext -i
	python setup.py build_ext -i

clean:
	-rm *.so
	-rm *.o
	-rm *.mod
