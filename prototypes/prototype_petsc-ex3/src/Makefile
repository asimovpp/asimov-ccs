#
#        FILE: Makefile
# DESCRIPTION: Makefile for PETSc ex3 prototype of ASiMoV-CCS.
#

FC = gfortran
FFLAGS = -cpp

EXE = ex3

SRC = $(EXE).f90
OBJ = $(SRC:%.f90=%.o)

INC = -I$(PETSC_DIR)/include
LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LIB)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^ $(INC)

clean:
	rm -f $(EXE) *.o
