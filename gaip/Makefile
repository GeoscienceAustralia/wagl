FC := gfortran
F2PY:= f2py

all: _shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so sys_variables.o sys_variables.mod angle_all.so set_satmod.so set_times.so unmiximage.pyf unmiximage.so

_shade_main_landsat_pixel.so _slope_pixelsize_newpole.so _brdf_terrain_newdiff_all.so sys_variables.o sys_variables.mod angle_all.so set_satmod.so set_times.so:
	$$(module load gcc/4.9.0)
	$(FC) -c sys_variables.f90
	python setup.py build_ext -i

unmiximage.pyf unmiximage.so:
	$(F2PY) -h unmiximage.pyf -m unmiximage unmiximage.f90
	$(F2PY) -c -m unmiximage unmiximage.pyf constants_NSWC.f90 nnls.f90 unmiximage.f90

clean:
	-rm *.so
	-rm *.o
	-rm *.mod
	-rm *.pyf
