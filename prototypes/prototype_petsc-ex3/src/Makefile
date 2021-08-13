#
#        FILE: Makefile
# DESCRIPTION: Makefile for PETSc ex3 prototype of ASiMoV-CCS.
#

.SUFFIXES:

FC = gfortran
FFLAGS = -cpp -std="f2018"
MODFLAGS = -fsyntax-only
DEPFLAGS = -M -MD -MF

EXE = ex3

SRC = $(wildcard *.f90)
OBJ = $(SRC:%.f90=%.o)
DEP = $(SRC:%.f90=%.d)

FFLAGS += -DACCS_PETSC
INC = -I$(PETSC_DIR)/include
LIB = -L$(PETSC_DIR)/lib -lpetsc

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(INC) $(LIB)

%.mod %.smod %._mod.o: %_mod.f90 %_mod.d
	$(FC) $(FFLAGS) -c $< $(INC)

%.mod %._mod.o: %_mod.f90 %_mod.d
	$(FC) $(FFLAGS) -c $< $(INC)

%.o: %.f90 %.d
	$(FC) $(FFLAGS) -o $@ -c $< $(INC)

%.d: %.f90
	-@$(FC) $(FFLAGS) $(DEPFLAGS) $@ $< $(INC)

clean:
	rm -f $(EXE) *.o *.mod *.smod
clean-full: clean
	rm -f *.d

$(DEP):

-include $(DEP)
