SHELL := /bin/bash
CC := gcc
FC := gfortran
FCFLAGS := -Wall -O3 -fbounds-check
CFLAGS := -O3

EXES := binear_ortho \
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
	$$(module load gcc/4.9.0)
	$(CC) $(CFLAGS) $< -o $@

%: %.f
	$$(module load gcc/4.9.0)
	$(FC) $(FCFLAGS) $< -o $@

clean:
	-rm -f $(OBJS) $(EXES)
