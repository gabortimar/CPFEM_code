objects= fas-con-mod.o data_in.o data_out.o fspl_orir.o csvmkr.o
	
FC=ifort
FFLAGS= -i-static -O3

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

csvmkr:$(objects)
	ifort $(FFLAGS) -o csvmkr $(objects)

fas-con-mod.o:fas-con-mod.f90

fspl_orir.o:fspl_orir.f90

csvmkr.o:csvmkr.f90

data_in.o:data_in.f90

data_out.o:data_out.f90

.PHONY: clean

clean:
	rm csvmkr $(objects) fas_com.mod 
install:
	cp csvmkr /home/fonseca/bin
