# make file for compiling the fortran code unmiximage.f90
# in such a way that it interfaces to python.

# The recipe is:
# 1. Create the target unmiximage.pyf (make unmiximage.pyf)
# 2. Edit the target file $(UNMIX_IMAGE).pyf so that
#    all variables except image, endmembermatrix, innull, outnull
#    and fractionsimage are hidden (use intent(hide)).
#    image (or pixelrefl), endmembermatrix, innull and outnull are input variables (use intent(in))
#    fractionsimage is an output variable (use intent(out)).
# 3. Create the binary file. (make unmiximage.so)
# 4. make install (simply moves the binary to a subdirectory)
# 5. Test by running unmiximage_test.py. The output will be printed
# to the screen and must be the same as:
# [ 0.44815772  0.25296948  0.21086879  0.06933006]
# [ 0.58656233  0.          0.2432272   0.05688818]
# [ 0.66786787  0.06846847  0.          0.11460207]
# [ 0.          0.99988303  0.107227    0.12129509]
# [ 0.          0.8780015   0.23776422  0.09579347]
# [-10. -10. -10. -10.]
# [ 0.44815772  0.25296948  0.21086879  0.06933006]
# 5. Run a second test to make sure that the numbers match a
# set of reference pixels - unmix_image_compare_with_reference_fractions.py


F2PY = f2py
MODULEDIR = .
UNMIX_IMAGE = unmiximage
UNMIX_IMAGE_PYF = $(UNMIX_IMAGE).pyf
UNMIX_IMAGE_SO = $(UNMIX_IMAGE).so

default:	all

all:	$(UNMIX_IMAGE_PYF) $(UNMIX_IMAGE_SO) install

$(UNMIX_IMAGE_PYF): unmiximage.f90
	echo $(F2PY) -h $@ -m $(UNMIX_IMAGE) $<
	$(F2PY) -h $@ -m $(UNMIX_IMAGE) $<

$(UNMIX_IMAGE_SO): $(UNMIX_IMAGE_PYF) constants_NSWC.f90 nnls.f90 unmiximage.f90
	echo $(F2PY) $(F2PYFLAGS) -c -m $(UNMIX_IMAGE) $^
	$(F2PY) $(F2PYFLAGS) -c -m $(UNMIX_IMAGE) $^
	#mv $(UNMIX_IMAGE_SO) $(MODULEDIR)

install:
	# RSC_FRACCOV_PYTHONBINDIR is defined in the fractionalcover environment module.
	# cp $(MODULEDIR)/$(UNMIX_IMAGE_SO) $(RSC_FRACCOV_PYTHONBINDIR)

clean:
	-rm $(UNMIX_IMAGE_PYF) $(MODULEDIR)/$(UNMIX_IMAGE_SO)

