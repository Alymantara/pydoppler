# PyDoppler compiles this solver automatically (``make -B dopp``) and caches
# the executable by source and compiler hashes. Arrays are allocated at run
# time, so one executable serves every data and map size.

FC = gfortran
FFLAGS = -O2
DEBUG_FFLAGS = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

SOURCES = dop_core.f90 dop.f90

.PHONY: debug clean

dopp: $(SOURCES)
	$(FC) $(FFLAGS) -o dopp $(SOURCES)

debug:
	$(FC) $(DEBUG_FFLAGS) -o dopp $(SOURCES)

clean:
	rm -f dopp *.mod *.o
