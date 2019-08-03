# IBM aix
# with error traps and bounds checks:
# compile=cp -f iclock.f clock.f ; xlf -D -C -g -q flttrap
# optimized:
# compile=cp -f iclock.f clock.f ; xlf -O
# optimized, with extended memory allocation of 2048 MB:
#compile=cp -f iclock.f clock.f ; xlf -bmaxdata:2048000000 -bmaxstack:2048000000 -O

# DEC alpha with debug and bound checks
# compile=cp -f cclock.f clock.f ; f90 -C -g -ladebug
# optimized:
# compile=cp -f cclock.f clock.f ; g77  -O

# HP HP/UX
# compile=cp -f hclock.f clock.f ; f77 -O +U77
# with error traps and bounds checks:
# compile=cp -f hclock.f clock.f ; f77 +FPVZO -g -C -O +U77

# SUN sparc
#compile=cp -f cclock.f clock.f ; gfortran -O
#compile=cp -f cclock.f clock.f ; g77 -O -fno-globals -Wno-globals
#-fno-globals

compile=cp -f cclock.f clock.f ; gfortran -O -w

dop.out: dopin dop.in dopp
	dopp
dopp: dop.f emap.par
	$(compile) -o dopp dop.f clock.f

lobe.out: lobe lobe.in
	lobe
lobe: lobe.f
	$(compile) -o lobe lobe.f
