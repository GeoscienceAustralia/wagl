SHELL := /bin/bash

FC := gfortran

FCFLAGS := -Wall -O3 -fbounds-check

EXES := binear_ortho \
        brdf_sim_bin \
        brdf_sim_bin_slc \
        brdf_sim_bin_LS8 \
        input_modtran_ortho_ula \
        coefficient \
        read_flx_ga \
        read_modtrancor_ortho \
        input_modtran_ortho \
        refort_tp5_ga \
        refort_tp5_ga_trans \
        read_modtran \
        read_flx_ga_trans \

BINDIR         := ../bin
INSTALLED_EXES := $(addprefix $(BINDIR)/,$(EXES))

all: $(EXES)

install: $(EXES) $(INSTALLED_EXES)

$(EXES): %: %.f
	$(FC) $(FCFLAGS) $< -o $@

$(INSTALLED_EXES): $(BINDIR)/%: %
	cp -p $< $@

clean:
	-rm $(EXES)

