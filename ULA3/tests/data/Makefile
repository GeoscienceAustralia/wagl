FC=gfortran

all: float32.bin int32.bin int16.bin int8.bin

genbindata: genbindata.f90
	$(FC) $< -o $@

float32.bin int32.bin int16.bin int8.bin: genbindata
	./genbindata
