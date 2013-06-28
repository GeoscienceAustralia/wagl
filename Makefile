SHELL := /bin/bash

# Required modules (NCI)
#   - hdf4/4.2.6_2012   ==> for hdf_extractor (DEPRECATED)
#   - intel-fc/11.1.073 ==> for FORTRAN sources (using ifort)
#   - intel-mkl/10.3.0  ==> for FORTRAN sources (using ifort)

include Makefile.platform

all: binaries .doc


#
# Executables
#

binaries: cbin fbin ula3

cbin:
	$(MAKE) --directory=src install

fbin:
	$(MAKE) --directory=fortran --makefile=$(FMAKEFILE) install

ula3:
	$(MAKE) --directory=ULA3

gen-doc:
	docs/build-sphinx.sh

.doc: gen-doc
	$(MAKE) --directory=docs html

#
# Cleanup
#

clean: cclean fclean binclean pyclean ula3clean .docsclean

cclean:
	$(MAKE) --directory=src clean

fclean:
	$(MAKE) --directory=fortran --makefile=$(FMAKEFILE) clean

ula3clean:
	$(MAKE) --directory=ULA3 clean

binclean:
	-rm -f ./bin/*

pyclean:
	-find . -name "*.pyc" -delete

.docsclean:
	$(MAKE) --directory=docs clean

