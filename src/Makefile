SHELL := /bin/bash

include ../Makefile.platform


CC      := gcc
CFLAGS  := -O3

# DEPRECATED (hdf_extractor)
#LDFLAGS := -ldf -lmfhdf -ltiff -ljpeg -lz -lm

# The PLATFORM variable is required to accommodate minor differences
# in HDF4 configuration at NCI. For NEO/VCLOUD use PLATFORM=DCC.
# TODO: refactor conditional compilation mechanism in source code

### OLD (without modules) ###
#CFLAGS   := $(CFLAGS) -DPLATFORM_$(PLATFORM) -I$(HDF_ROOT)/include
#LDFLAGS  := $(LDFLAGS) -L$(HDF_ROOT)/lib

### NEW (with modules) ###
CFLAGS := $(CFLAGS) -DPLATFORM_$(PLATFORM)

EXES := aot_loader
OBJS := $(addsuffix .o,$(EXES))

BINDIR         := ../bin
INSTALLED_EXES := $(addprefix $(BINDIR)/,$(EXES))

all: $(EXES)
install: $(EXES) $(INSTALLED_EXES)

$(INSTALLED_EXES): $(BINDIR)/%: %
	cp -p $< $@

$(OBJS): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(EXES): %: %.o
	$(CC) $(LDFLAGS) $< -o $@

clean:
	-rm -f $(OBJS) $(EXES)




# HDF test program

#DEVOPT    := /data/nbar/dev/optbuild/__rel__
DEVOPT    := /data/nbar/opt
DEVOPTINC := $(DEVOPT)/include
DEVOPTLIB := $(DEVOPT)/lib

hdfx: hdfx.o
	$(CC) $< -o $@ -L$(DEVOPTLIB) -lmfhdf -ldf

hdfx.o: hdfx.c
	$(CC) -c $< -o $@ -O2 -Wall -I$(DEVOPTINC)

hclean:
	-rm hdf_test hdf_test_l hdf_test.o hdf_test_l.o


