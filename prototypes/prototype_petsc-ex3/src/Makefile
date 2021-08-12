#
#        FILE: Makefile
# DESCRIPTION: Makefile for PETSc ex3 prototype of ASiMoV-CCS.
#

FC = gfortran
FFLAGS = -cpp -std="f2018"
MODFLAGS = -fsyntax-only

EXE = ex3

MOD_SRC = accs_types_mod.f90 accs_petsctypes_mod.f90 accsvec_mod.f90 accs_utils_mod.f90
SRC =	$(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
MOD = $(MOD_SRC:%_mod.f90=%.mod)

INC = -I$(PETSC_DIR)/include
LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(MOD) $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(INC) $(LIB)

%.mod: %_mod.f90
	$(FC) $(FFLAGS) $(MODFLAGS) $< $(INC)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^ $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
