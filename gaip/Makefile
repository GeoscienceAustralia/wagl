FC := gfortran
F2PY:= f2py

all: _cast_shadow_mask.so _slope_self_shadow.so _surface_reflectance.so sys_variables.o sys_variables.mod _sat_sol_angles.so _satellite_model.so _track_time_info.so unmiximage.pyf unmiximage.so _interpolation.so

_cast_shadow_mask.so _slope_self_shadow.so _surface_reflectance.so sys_variables.o sys_variables.mod _sat_sol_angles.so _satellite_model.so _track_time_info.so:
	$$(module load gcc/5.2.0)
	$(FC) -c sys_variables.f90
	python setup.py build_ext -i

unmiximage.pyf unmiximage.so:
	$(F2PY) -h unmiximage.pyf -m unmiximage unmiximage.f90
	$(F2PY) -c -m unmiximage unmiximage.pyf constants_NSWC.f90 nnls.f90 unmiximage.f90

clean:
	@-rm -f *.so
	@-rm -f *.o
	@-rm -f *.mod
	@-rm -f *.pyf
	@-rm -fr build
