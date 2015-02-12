objects=fas_com_mod.o fs2-xtal.o fs3p-gps.o fs3p-matl.o fs3p-swp.o fs3p-bkc.o fs3p-io_nowin.o fs3p-rcm.o fs3pre.o
FC=ifort
FFLAGS=-i-static -O3
%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

fs3pre:$(objects)
	ifort $(FFLAGS) -o fs3pre $(objects)
fs3pre.o:fs3pre.f90
fs3p-swp.o:fs3p-swp.f90
fs3p-bkc.o:fs3p-bkc.f90
fs3p-gps.o:fs3p-gps.f90
fs3p-matl.o:fs3p-matl.f90
fs3p-io_nowin.o:fs3p-io_nowin.f90
fs3p-rcm.o:fs3p-rcm.f90
fs2-xtal.o:fs2-xtal.f90
fas_com_mod.o:fas_com_mod.f90

.PHONY:clean
clean:
	rm fs3pre $(objects) fas_com.mod  

install:
	cp fs3pre ~/bin/
