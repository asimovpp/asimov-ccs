#
#        FILE: Makefile
# DESCRIPTION: Makefile for PETSc ex3 prototype of ASiMoV-CCS.
#

FC = gfortran
FFLAGS =

EXE = ex3

SRC = $(EXE).f90
OBJ = $(SRC:%.f90=%.o)

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^

clean:
	rm -f $(EXE) *.o
