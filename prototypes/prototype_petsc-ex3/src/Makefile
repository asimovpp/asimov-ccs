#
#        FILE: Makefile
# DESCRIPTION: Makefile for PETSc ex3 prototype of ASiMoV-CCS.
#

FC = gfortran
FFLAGS = -cpp -std="f2018"
MODFLAGS = -fsyntax-only

EXE = ex3

MOD_SRC = $(wildcard *_mod.f90)
SRC =	$(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
MOD = $(MOD_SRC:%_mod.f90=%.mod)

INC = -I$(PETSC_DIR)/include
LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(INC) $(LIB)

# $(OBJ): $(MOD)

%.mod: %_mod.f90
	$(FC) $(FFLAGS) $(MODFLAGS) $< $(INC)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^ $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
