SHELL := /bin/bash

F90_OPTS = --f90flags="-assume byterecl"

all: .fc filter shade slope terrain

.fc:
	$(MAKE) --directory=fc --makefile=Makefile all

filter:
	f2py -m filter $(F90_OPTS) -c filter.f90

shade:
	source ~/.profile; \
	module unload intel-fc; \
	module unload intel-cc; \
	module load gcc/4.4.4; \
	f2py -m _shade_main_landsat_pixel --fcompiler=gnu95 --opt="-O0" -c shade_main_landsat_pixel.f90

slope:
	source ~/.profile; \
	module unload intel-fc; \
	module unload intel-cc; \
	module load gcc/4.4.4; \
	f2py -m _slope_pixelsize_newpole --fcompiler=gnu95 --opt="-O0" -c slope_pixelsize_newpole.f90

terrain:
	source ~/.profile; \
	module unload intel-fc; \
	module unload intel-cc; \
	module load gcc/4.4.4; \
	f2py -m _brdf_terrain_newdiff_all --fcompiler=gnu95 --opt="-O0" -c brdf_terrain_newdiff_all.f90





clean: rootclean fcclean

rootclean:
	-rm *.so

fcclean:
	$(MAKE) --directory=fc clean
