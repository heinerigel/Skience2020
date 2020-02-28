#gfortran
F90=gfortran
Opt=-O3 -ffree-line-length-none
#intel fortran
#F90=ifort
#Opt=-O3
#
all: sem1D
#
OBJ=  nrtype.o  nrutil.o module_spline.o  speclib.o sem_module.o sem1D.o

sem1D: $(OBJ)
	$(F90) $(Opt) $(OBJ)   -o sem1D
	mv sem1D ../bin
clean:
	\rm $(OBJ) *.mod

dfour1.o: dfour1.f
	$(F90) $(Opt) -c dfour1.f
gaussj.o: gaussj.f
	$(F90) $(Opt) -c gaussj.f
%.o: %.f90
	$(F90) $(Opt) -c $*.f90
