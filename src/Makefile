SHELL := /bin/bash
CC := gcc
FC := gfortran
FCFLAGS := -Wall -O3 -fbounds-check
CFLAGS := -O3

EXES := binear_ortho \
        brdf_sim_bin \
        brdf_sim_bin_slc \
        brdf_sim_bin_LS8 \
        coefficient \
        read_flx_ga \
        read_modtrancor_ortho \
        input_modtran_ortho \
        refort_tp5_ga \
        refort_tp5_ga_trans \
        read_modtran \
        read_flx_ga_trans \
	aot_loader

OBJS := $(addsuffix .o,$(EXES))

BINDIR         := ../bin
INSTALLED_EXES := $(addprefix $(BINDIR)/,$(EXES))

all: $(EXES)

install: $(EXES) $(INSTALLED_EXES)

$(INSTALLED_EXES): $(BINDIR)/%: %
	cp -p $< $@

%: %.c
	$$(module load gcc)
	$(CC) $(CFLAGS) $< -o $@

%: %.f
	$$(module load gcc)
	$(FC) $(FCFLAGS) $< -o $@

clean:
	-rm -f $(OBJS) $(EXES)
