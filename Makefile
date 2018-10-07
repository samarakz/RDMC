LIBDIR = -L$(LAPACKPATH)
LDLIBS = -llapack -lrefblas
FCOMPILER = gfortran
#F90FLAGS = -g -fbacktrace -fcheck=all
F90FLAGS = -g -fbacktrace

%.o: %.F90       # if any file with .o looking at, also consider .F90 files with the same name
	$(FCOMPILER) $(F90FLAGS) -c $<

%.o: %.f90
	$(FCOMPILER) $(F90FLAGS) -c $<

seedsig:	seedsig_selfenergy.o seedsig_harmonics.o seedsig.o
	$(FCOMPILER) $^ $(LIBDIR) $(LDLIBS) -o $(@)

seedsig.o: seedsig_selfenergy.o seedsig_harmonics.o # dependence of each routine to other routines
seedsig_selfenergy.o:
seedsig_harmonics.o:

clean:
	rm -f *.o seedsig
