objects=fas-com-mod.o Bconds.o Cpstrs.o Elvals.o \
Faslt_io_nowin.o Fasolt3.o Shapef.o Sig_eps.o \
Stif_sol.o crys-cons.o crys_math.o ela-plas.o

#ifeq(machine,'i486') 
#FC=ifort -i-static -g -fpe0 -traceback
#FC=ifort -i-static -O3 -p -g
FC=ifort -i-static -O3
#endif
#LIBS=-L/opt/intel/fce/9.1.036/lib/libimf.so

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

fasolt3:$(objects)
	$(FC) -o fasolt3 $(objects)

fasolt3.o:fasolt.f90

Bconds.o:Bconds.f90

Cpstrs.o:Cpstrs.f90

Elvals.o:Elvals.f90

Shapef.o:Shapef.f90

Sig_eps.o:Sig_eps.f90

Stif_sol.o:Stif_sol.f90

crys-cons.o:crys-cons.f90

crys_math.o:crys_math.f90

ela_plas.o:ela_plas.f90

Faslt-io_nowin.o:Faslt-io_nowin.f90

fas_com_mod.o:fas_com_mod.f90

.PHONY:clean
clean:
	rm fasolt3 $(objects) fas_com.mod

install:
	cp fasolt3 /Users/joao/bin/
