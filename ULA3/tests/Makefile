SHELL := /bin/bash

all: test_histogram

test_histogram:
	f2py -m _test_IDL_Histogram -c test_IDL_Histogram.f90

clean:
	-rm *.so
