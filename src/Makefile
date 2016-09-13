SHELL := /bin/bash
CC := gcc
FC := gfortran
FCFLAGS := -Wall -O3 -fbounds-check
CFLAGS := -O3

EXES := aot_loader

OBJS := $(addsuffix .o,$(EXES))

BINDIR         := ../bin
INSTALLED_EXES := $(addprefix $(BINDIR)/,$(EXES))

all: $(EXES)

install: $(EXES) $(INSTALLED_EXES)

$(INSTALLED_EXES): $(BINDIR)/%: %
	@cp -p $< $@

%: %.c
	$$(module load gcc/5.2.0)
	$(CC) $(CFLAGS) $< -o $@

%: %.f
	$$(module load gcc/5.2.0)
	$(FC) $(FCFLAGS) $< -o $@

clean:
	@rm -f $(OBJS) $(EXES)
